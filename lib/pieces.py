from functools import lru_cache
from collections import Counter, namedtuple
from heapq import heappop, heappush
from fractions import Fraction

from fast_fractions import FastFraction
from base_maps import BaseMap, gen_constraint_cube_maps


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

    def get_last_map(self):
        curve = self.curve
        dim, G = curve.dim, curve.genus
        last_map = BaseMap.id_map(dim)
        for cnum in self.pos.cnums:
            cnum = last_map.apply_cnum(G, cnum)
            last_map = last_map * curve.base_maps[cnum]  # именно в таком порядке!
        return last_map

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

    # brkline - ломаная, последовательность пар (x,t); x \in [0,1]^d, t \in [0,1]
    def get_bounds(self, ratio_func, brkline=None):
        dim = self.curve.dim
        N = self.curve.div
        G = self.curve.genus

        pos1, pos2 = self.pos1, self.pos2

        l1, x1, t1 = pos1.get_int_coords()
        l2, x2, t2 = pos2.get_int_coords()

        use_brkline = (brkline is not None)
        if use_brkline:
            last1 = self.piece1.get_last_map()
            brkline1 = [(last1.apply_x(x), last1.apply_t(t)) for x, t in brkline]
            last2 = self.piece2.get_last_map()
            brkline2 = [(last2.apply_x(x), last2.apply_t(t)) for x, t in brkline]

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
                if use_brkline:
                    brkline1 = [(x, 1 - t) for x, t in brkline1]

            # junc: сначала поворот
            base_map = junc.base_map
            x2 = base_map.apply_cube(pos2.sub_div, x2)
            t2 = base_map.apply_cnum(pos2.sub_genus, t2)

            if use_brkline:
                brkline2 = [(base_map.apply_x(x), base_map.apply_t(t)) for x, t in brkline2]

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

        # junc: потом сдвиг
        t2 += junc_dt * mt
        x2 = [x2j + junc_dxj * mx for x2j, junc_dxj in zip(x2, junc_dx)]

        max_dx = [None] * dim
        for j in range(dim):
            x1j = x1[j]
            x2j = x2[j]
            dxj = max(abs(x1j - x2j + 1), abs(x1j - x2j - mx2))
            max_dx[j] = dxj

        max_dt = t2 + mt2 - t1  # max(t_2 - t_1)
        min_dt = t2 - (t1 + 1)  # min(t_2 - t_1)

        lo = FastFraction(*ratio_func(dim, max_dx, max_dt))
        up = FastFraction(*ratio_func(dim, max_dx, min_dt))


        argmax = None
        if use_brkline:
            for x1rel, t1rel in brkline1:
                t1_final = t1 + t1rel
                x1_final = [x1j + x1relj for x1j, x1relj in zip(x1, x1rel)]
                for x2rel, t2rel in brkline2:
                    t2_final = t2 + t2rel * mt2
                    x2_final = [x2j + x2relj * mx2 for x2j, x2relj in zip(x2, x2rel)]

                    dx = [x1j - x2j for x1j, x2j in zip(x1_final, x2_final)]
                    dt = t2_final - t1_final

                    lo_final = Fraction(*ratio_func(dim, dx, dt))
                    if lo_final > 80:
                        print(pos1.cnums, pos1.cubes, pos2.cnums, pos2.cubes)
                        print(junc)
                        print(x1, t1, x2, t2)
                        print(x1rel, t1rel, x2rel, t2rel)

                    lo_final = FastFraction(lo_final.numerator, lo_final.denominator)
                    if lo_final > lo:
                        x1_real = [Fraction(x1j, mx) for x1j in x1_final]
                        x2_real = [Fraction(x2j, mx) for x2j in x2_final]
                        t1_real = Fraction(t1_final, mt)
                        t2_real = Fraction(t2_final, mt)
                        argmax = {'x1': x1_real, 't1': t1_real, 'x2': x2_real, 't2': t2_real, 'junc': junc}
                        lo = lo_final

        return lo, up, argmax


# пары фракций с приоритетом
class PairsTree:

    RichPair = namedtuple('RichPair', [
        'priority',  # первое поле, по нему идёт сравнение; для быстрой сортировки по upper_bound
        'pair',  # CurvePieceBalancedPair
        'up',  # upper_bound
        'lo',  # lower_bound
        'argmax',
    ])

    def __init__(self, ratio_func, brkline=None):
        self.data = []
        self.stats = Counter()
        self.bad_pairs = []
        self.ratio_func = ratio_func
        self.good_threshold = None  # если отношение <= порога, пара хорошая
        self.bad_threshold = None  # если отношение > порога, пара плохая
        self.brkline = brkline
        self._inc = 0

    def set_good_threshold(self, thr):
        assert isinstance(thr, FastFraction)
        self.good_threshold = thr

    def set_bad_threshold(self, thr):
        assert isinstance(thr, FastFraction)
        self.bad_threshold = thr

    def add_pair(self, pair):
        lo, up, argmax = pair.get_bounds(self.ratio_func, brkline=self.brkline)
        gthr = self.good_threshold
        if gthr is not None and up < gthr:
            self.stats['good'] += 1
            return

        bthr = self.bad_threshold
        if bthr is not None and lo > bthr:
            self.bad_pairs.append(pair)
            self.stats['bad'] += 1
            return

        self._inc += 1  # чтобы сравнение не проваливалось к парам
        node = PairsTree.RichPair(priority=(-up, self._inc), pair=pair, lo=lo, up=up, argmax=argmax)
        heappush(self.data, node)

    def divide(self):
        if not self.data:
            return

        worst_node = heappop(self.data)
        for pair in worst_node.pair.divide():
            self.add_pair(pair)
