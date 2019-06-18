import time
from fractions import Fraction
import itertools
from functools import lru_cache
from collections import Counter, namedtuple
from heapq import heappop, heappush

from fast_fractions import FastFraction
from base_map import BaseMap, gen_constraint_cube_maps
import curve_sat_adapter 


@lru_cache(maxsize=2**20)
def get_int_cube_with_cache(dim, N, cubes):
    # Все x умножаем на N*l
    x = [0] * dim
    Npower = 1
    for cube in reversed(cubes):
        for j in range(dim):
            x[j] += cube[j] * Npower
        Npower *= N
    return x

@lru_cache(maxsize=2**20)
def get_int_time_with_cache(dim, N, cnums):
    G = N**dim
    # все t на G^l
    # t = c0/G + c1/G**2 + ... = (c_{l-1} + c_{l-2}*G + ..) / G^l
    t = 0
    Gpower = 1
    for cnum in reversed(cnums):
        t += cnum * Gpower
        Gpower *= G
    return t


class PartialFractalCurve:
    # base_maps - list as in FractalCurve, may contain None
    # repr_maps - список представителей (coset representatives) base_map-ов, которые сохраняют вход-выход
    # symmetries - симметрии кривой
    def __init__(self, dim, div, proto, base_maps, repr_maps, symmetries):
        self.dim = dim
        self.div = div
        self.proto = tuple(proto)
        self.base_maps = tuple(base_maps)
        self.repr_maps = tuple(repr_maps)
        self.symmetries = tuple(symmetries)
        self.genus = self.div ** self.dim

    # создать кривую с другим прототипом/base_maps/whatever
    def changed(self, proto=None, base_maps=None, repr_maps=None, symmetries=None):
        return type(self)(
            dim=self.dim,
            div=self.div,
            proto=proto if proto is not None else self.proto,
            base_maps=base_maps if base_maps is not None else self.base_maps,
            repr_maps=repr_maps if repr_maps is not None else self.repr_maps,
            symmetries=symmetries if symmetries is not None else self.symmetries,
        )

    def get_fraction(self, cnum):
        """Get fraction as a curve."""
        return self.apply_base_map(self.base_maps[cnum])

    def reverse(self):
        """Reverse time in a curve."""

        kwargs = {}
        if hasattr(self, 'repr_maps'):
            kwargs['repr_maps'] = reversed(self.repr_maps)

        # симметрии не меняются!

        return self.changed(
            # прототип проходится в обратном порядке
            proto = reversed(self.proto),

            # базовые преобразования проходятся в обратном порядке
            # сами по себе они не меняются:
            #   - если обращения времени не было, то его и не будет
            #   - изометрия куба не меняется, т.к. время не играет роли
            base_maps = reversed(self.base_maps),

            **kwargs,
        )

    # кандидаты в self.base_maps[cnum]
    def gen_allowed_maps(self, cnum):
        if self.base_maps[cnum] is not None:
            # базовое преобразование уже определено!
            yield self.base_maps[cnum]
            return

        repr_map = self.repr_maps[cnum]
        for symm in self.symmetries:
            yield repr_map * symm

    def get_piece_position(self, cnum):
        return CurvePiecePosition(
            dim=self.dim,
            div=self.div,
            cnums=[cnum],
            cubes=[self.proto[cnum]],
        )

    def specify(self, cnum, base_map):
        if base_map not in self.gen_allowed_maps(cnum):
            raise Exception("Can't specify curve")

        new_base_maps = list(self.base_maps)
        new_base_maps[cnum] = base_map
        return self.changed(base_maps=new_base_maps)

    def apply_base_map(self, base_map):
        """Apply base map to a fractal curve, return new curve."""
        if base_map.dim != self.dim:
            raise Exception("Incompatible base map!")

        # можно разложить базовое преобразование в произведение (коммутирующих) 
        # преобразований: обращение времени с тождественной изометрией  +  изометрии куба
        if base_map.time_rev:
            return self.reverse().apply_base_map(base_map.cube_map())

        # применяем изометрию куба

        # прототип подвергается изометрии
        proto = [base_map.apply_cube(self.div, cube) for cube in self.proto]

        # базовые преобразования сопрягаются: действительно, чтобы получить
        # из преобразованной кривой её фракцию, можно сделать так:
        # - сначала вернуться к исходной кривой (inv)
        # - применить преобразование исходной кривой для перехода к фракции (bm)
        # - перейти к преобразованной кривой (base_map)
        inv = base_map.inverse()
        conj_cache = {}
        def conjugate(bm):
            if bm not in conj_cache:
                conj_cache[bm] = base_map * bm * inv
            return conj_cache[bm]

        new_maps = [conjugate(bm) if bm is not None else None for bm in self.base_maps]
        kwargs = {}
        if hasattr(self, 'symmetries'):
            kwargs['symmetries'] = [conjugate(bm) for bm in self.symmetries]
        if hasattr(self, 'repr_maps'):
            kwargs['repr_maps'] = [conjugate(bm) for bm in self.repr_maps]

        return self.changed(proto=proto, base_maps=new_maps, **kwargs)

    def bm_info(self):
        return {cnum: bm for cnum, bm in enumerate(self.base_maps) if bm is not None}

    def is_specialization(self, tmpl):
        return all(self.base_maps[cnum] == bm for cnum, bm in tmpl.bm_info().items())

    #
    # Про отношение
    #

    def init_pairs_tree(self):
        G = self.genus
        for cnum1 in range(G):
            for cnum2 in range(cnum1 + 2, G):
                pos1 = self.get_piece_position(cnum1)
                pos2 = self.get_piece_position(cnum2)
                yield CurvePieceBalancedPair(self, None, pos1, pos2)

        for junc in self.get_junctions_info().keys():
            last_cnum1 = 0 if junc.time_rev else G - 1
            first_cnum2 = G - 1 if junc.base_map.time_rev else 0
            for cnum1 in range(G):
                for cnum2 in range(G):
                    if (cnum1, cnum2) == (last_cnum1, first_cnum2):
                        continue
                    pos1 = self.get_piece_position(cnum1)
                    pos2 = self.get_piece_position(cnum2)
                    yield CurvePieceBalancedPair(self, junc, pos1, pos2)

    # upper_bound - если кривая лучше - нам подходит
    def estimate_ratio(self, ratio_func, lower_bound, upper_bound, max_iter=10**9, log_pack=100, sat_pack=100, find_model=False, verbose=False):
        adapter = curve_sat_adapter.CurveSATAdapter(dim=self.dim)
        adapter.init_curve(self)

        pairs_tree = PairsTree(ratio_func)
        pairs_tree.set_good_threshold(upper_bound)
        pairs_tree.set_bad_threshold(lower_bound)

        for pair in self.init_pairs_tree():
            pairs_tree.add_pair(pair)

        it = 0
        while it < max_iter:
            it += 1
            pairs_tree.divide()
            while pairs_tree.bad_pairs:
                bad_pair = pairs_tree.bad_pairs.pop()
                adapter.add_forbid_clause(bad_pair.junc, bad_pair.curve)

            if not pairs_tree.data:
                break

            if it % log_pack == 0 and verbose:
                worst_node = pairs_tree.data[0]
                worst_pair = worst_node.pair
                print({
                    'iter': it,
                    'pairs': len(pairs_tree.data),
                    'pqstats:': pairs_tree.stats,
                    'up': float(worst_node.up),
                    'depth': (worst_pair.pos1.depth, worst_pair.pos2.depth),
                })

            if it % sat_pack == 0:
                if verbose:
                    print('iter:', it, 'adapter stats:', adapter.stats())
                if not adapter.solve():
                    print('no SAT model')
                    return False

        if it == max_iter:
            print('used all iterations...')
            return False

        if not adapter.solve():
            print('no SAT model')
            return False

        print('SAT model exists!')
        if not find_model:
            return True

        # это если попросят модель
        model = adapter.get_model()
        return {
            "model": model,
            "curve": adapter.get_curve_from_model(self, model),
            "pairs_tree": pairs_tree,
        }


    #
    # Стыки.
    #

    # словарь {junc: curve_list} кривых, приводящих к данному стыку
    def get_junctions_info(self):
        # строим конфигурации:
        # это набор (i, curve), где в curve заданы bm[0], bm[-1], bm[i], bm[i+1]
        configs = []
        G = self.genus
        for bm_first in self.gen_allowed_maps(0):
            for bm_last in self.gen_allowed_maps(G - 1):
                for cnum in range(G - 1):
                    for bm_i in self.gen_allowed_maps(cnum):
                        if cnum == 0 and bm_i != bm_first:
                            continue
                        for bm_ii in self.gen_allowed_maps(cnum + 1):
                            if cnum + 1 == G - 1 and bm_ii != bm_last:
                                continue
                            curve = self.specify(0, bm_first)\
                                .specify(G - 1, bm_last)\
                                .specify(cnum, bm_i)\
                                .specify(cnum + 1, bm_ii)

                            configs.append((cnum, curve))

        # конфигурации, приводящие к данному стыку
        junc_curves = {}

        for cnum, curve in configs:
            base_junc = curve.get_base_junction(cnum)
            for junc in curve.gen_junctions_from_base([base_junc]):
                junc_curves.setdefault(junc, []).append(curve)

        return junc_curves

    def get_base_junction(self, cnum):
        delta = [c2j - c1j for c1j, c2j in zip(self.proto[cnum], self.proto[cnum + 1])]
        return Junction(delta, self.base_maps[cnum], self.base_maps[cnum + 1])

    # возвращает стыки вместе с производными
    def gen_junctions_from_base(self, base_juncs):
        for junc in base_juncs:
            yield junc
        seen = set(base_juncs)
        to_derive = list(base_juncs)
        while to_derive:
            junc = to_derive.pop()
            dj = self.get_derived_junction(junc)
            if dj not in seen:
                yield dj
                seen.add(dj)
                to_derive.append(dj)

    def get_derived_junction(self, junc):
        base_map = junc.base_map
        cnum1 = 0 if junc.time_rev else -1
        cnum2 = -1 if base_map.time_rev else 0

        if self.base_maps[cnum1] is None or self.base_maps[cnum2] is None:
            raise Exception("Can't get derivative: base_map not defined")

        cube1 = self.proto[cnum1]
        bm1 = self.base_maps[cnum1]
        if junc.time_rev:
            bm1 = bm1.reverse_time()

        cube2 = base_map.apply_cube(self.div, self.proto[cnum2])
        bm2 = base_map * self.base_maps[cnum2]

        der_delta = tuple(c2j + dxj * self.div - c1j for c1j, c2j, dxj in zip(cube1, cube2, junc.delta_x))

        return Junction(der_delta, bm1, bm2)

# стык двух фракций кривой
# всегда приводим к стандартному виду:
# первая фракция стандартной пространственной ориентации, но, возможно, с обращением времени (self.time_rev)
# вторая фракция в ориентации self.base_map
# куб второй фракции получается из куба первого сдвигом на delta_x \in {0,1,-1}^d
class Junction:
    def __init__(self, delta, bm1, bm2):
        bm1_cube_inv = bm1.cube_map().inverse()
        self.base_map = bm1_cube_inv * bm2
        self.delta_x = bm1_cube_inv.apply_vec(delta)
        self.time_rev = bm1.time_rev

    def _data(self):
        return (self.delta_x, self.time_rev, self.base_map)

    def __eq__(self, other):
        return self._data() == other._data()

    def __hash__(self):
        return hash(self._data())

    def __repr__(self):
        return '1: ' + ('t->1-t' if self.time_rev else '') + ', 2: ' + str(self.base_map) + ', --> ' + str(self.delta_x)


# вся инфа о положении фракции в кривой
# cnums - последовательность номеров кубов
# cubes - последовательность кубов
# каждый куб - для той кривой, которая в соотв. фракции
class CurvePiecePosition:
    def __init__(self, dim, div, cnums, cubes):
        self.dim = dim
        self.div = div
        self.cnums = cnums
        self.cubes = cubes
        self.depth = len(self.cnums)
        self.sub_div = div**self.depth
        self.sub_genus = self.sub_div**dim

    def specify(self, cnum, cube):
        return type(self)(
            dim=self.dim,
            div=self.div,
            cnums=self.cnums + [cnum],
            cubes=self.cubes + [cube],
        )

    # естественные целочисленные координаты - время и нулевой угол куба
    # int time: t - время, умноженное на G^l, кусок - [t, t+1], abs time: [t/G^l, (t+1)/G^l]
    # int cube: cj <= xj <= cj+1, abs cube: cj/N^l <= xj <= (cj+1)/N^l
    # int curve: [0, N**l]^d, abs curve: [0,1]^d
    # l = depth
    # возвращает l,x,t
    def get_int_coords(self):
        return (
            self.depth,
            get_int_cube_with_cache(self.dim, self.div, tuple(self.cubes)),
            get_int_time_with_cache(self.dim, self.div, tuple(self.cnums)),
        )

class CurvePiece:
    # фракция кривой = кривая + позиция фракции
    def __init__(self, curve, pos):
        self.curve = curve
        self.pos = pos

    # делим кривую дальше всеми способами
    def divide(self):
        curve = self.curve
        dim, G = curve.dim, curve.genus

        # определим ориентацию предпоследнего кусочка
        # в последнем кубе ориентация не задана!
        prev_map = BaseMap.id_map(dim)
        for cnum in self.pos.cnums[:-1]:
            cnum = prev_map.apply_cnum(G, cnum)
            prev_map = prev_map * curve.base_maps[cnum]  # именно в таком порядке!

        prev_curve = curve.apply_base_map(prev_map)

        active_cnum = self.pos.cnums[-1]  # кубик, где всё происходит
        spec_cnum = prev_map.apply_cnum(G, active_cnum)

        # делим
        for bm in prev_curve.gen_allowed_maps(active_cnum):

            # сопряжение, как в apply:
            # для prev_curve.base_maps[active_cnum] = bm =>  orig_curve.base_maps[spec_cnum] = ...
            new_map = prev_map.inverse() * bm * prev_map
            specified_curve = curve.specify(spec_cnum, new_map)

            last_curve = prev_curve.apply_base_map(bm)
            for cnum, cube in enumerate(last_curve.proto):
                new_pos = self.pos.specify(cnum, cube)
                new_piece = CurvePiece(specified_curve, new_pos)
                yield new_piece


class CurvePieceBalancedPair:
    # сбалансированная пара фракций кривых
    # считаем, что t_2 > t_1
    def __init__(self, curve, junc, pos1, pos2):
        self.curve = curve
        self.junc = junc
        self.pos1 = pos1
        self.pos2 = pos2
        self.piece1 = CurvePiece(self.curve, self.pos1)
        self.piece2 = CurvePiece(self.curve, self.pos2)

    def divide(self):
        # при делении кусочка у него уточняется кривая, поэтому берём кривую из него
        # решаем, кого делить
        if self.pos1.depth > self.pos2.depth:
            for subpiece in self.piece2.divide():
                yield type(self)(subpiece.curve, self.junc, self.pos1, subpiece.pos)
        else:
            for subpiece in self.piece1.divide():
                yield type(self)(subpiece.curve, self.junc, subpiece.pos, self.pos2)

    # пространственные и временные расстояния
    def int_dist(self):
        if not hasattr(self, '_int_dist'):
            self._int_dist = self.get_int_dist()
        return self._int_dist

    def get_int_dist(self):
        dim = self.curve.dim
        N = self.curve.div
        G = self.curve.genus

        pos1, pos2 = self.pos1, self.pos2

        l1, x1, t1 = pos1.get_int_coords()
        l2, x2, t2 = pos2.get_int_coords()

        junc = self.junc
        if junc is None:
            junc_dt = 0
            junc_dx = (0,) * dim
        else:
            junc_dt = 1
            junc_dx = self.junc.delta_x

            # junc: time_rev
            if junc.time_rev:
                t1 = pos1.sub_genus - 1 - t1

            # junc: сначала поворот
            x2 = junc.base_map.apply_cube(pos2.sub_div, x2)
            t2 = junc.base_map.apply_cnum(pos2.sub_genus, t2)

        # приведение к единому масштабу
        if l1 == l2:
            mx2 = mt2 = 1
        elif l1 == l2 + 1:
            x2 = [x2j * N for x2j in x2]
            t2 *= G
            mx2 = N
            mt2 = G
        else:
            raise Exception("Bad coordinates!")

        mx = pos1.sub_div
        mt = pos1.sub_genus

        # мы привели целые координаты к следующим:
        # cube1: x1j <= xj <= x1j + 1         -- кубик внутри кривой [0, mx]^d
        # cube2: x2j <= xj <= x2j + mx2  -- кубик внутри кривой [0, mx2 * mx]^d, + сдвиг на junc_dx * mx
        #
        # time1: t1 <= t <= t1 + 1
        # time2: t2 <= t <= t2 + mt2, после сдвига: t2 + junc_dt * mt <= t <= t2 + mt2 + junc_dt * mt

        max_dx = [None] * dim
        for j in range(dim):
            x1j = x1[j]

            # junc: потом сдвиг
            x2j = x2[j] + junc_dx[j] * mx

            dxj = max(abs(x1j - x2j + 1), abs(x1j - x2j - mx2))
            max_dx[j] = dxj

        max_dt = t2 + junc_dt * mt + mt2 - t1  # max(t_2 - t_1)
        min_dt = t2 + junc_dt * mt - (t1 + 1)  # min(t_2 - t_1)

        return {'max_dx': max_dx, 'min_dt': min_dt, 'max_dt': max_dt}

    def upper_bound(self, ratio_func):
        dist = self.int_dist()
        return FastFraction(*ratio_func(self.curve.dim, dist['max_dx'], dist['min_dt']))

    def lower_bound(self, ratio_func):
        dist = self.int_dist()
        return FastFraction(*ratio_func(self.curve.dim, dist['max_dx'], dist['max_dt']))

# пары фракций с приоритетом
class PairsTree:

    RichPair = namedtuple('RichPair', [
        'priority',  # первое поле, по нему идёт сравнение; для быстрой сортировки по upper_bound
        'pair',  # CurvePieceBalancedPair
        'up',  # upper_bound
        'lo',  # lower_bound
    ])

    def __init__(self, ratio_func):
        self.data = []
        self.stats = Counter()
        self.bad_pairs = []
        self.ratio_func = ratio_func
        self.good_threshold = None  # если отношение <= порога, пара хорошая
        self.bad_threshold = None  # если отношение > порога, пара плохая
        self._inc = 0
        self.LOG = []

    def set_good_threshold(self, thr):
        self.good_threshold = thr

    def set_bad_threshold(self, thr):
        self.bad_threshold = thr

    def add_pair(self, pair):
        up = pair.upper_bound(self.ratio_func)
        gthr = self.good_threshold
        if gthr is not None and float(up) < gthr:  # TOOD нет ли проблемы с округлением
            self.stats['good'] += 1
            return

        lo = pair.lower_bound(self.ratio_func)
        bthr = self.bad_threshold
        if bthr is not None and float(lo) > bthr:
            self.bad_pairs.append(pair)
            self.stats['bad'] += 1
            return

        self._inc += 1  # чтобы сравнение не проваливалось к парам
        node = PairsTree.RichPair(priority=(-float(up), self._inc), pair=pair, up=up, lo=lo)
        heappush(self.data, node)

    def divide(self):
        if not self.data:
            return

        worst_node = heappop(self.data)
        for pair in worst_node.pair.divide():
            self.add_pair(pair)
