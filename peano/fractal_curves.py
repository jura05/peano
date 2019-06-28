# coding: utf-8

import itertools
from fractions import Fraction

from .utils import get_lcm
from .base_maps import BaseMap, gen_constraint_cube_maps
from .fast_fractions import FastFraction
from . import partial_fractal_curves
from . import pieces


class FractalCurve(partial_fractal_curves.PartialFractalCurve):
    """Class representing fractal peano curve in [0,1]^d.
    Params:
        div         positive integer, number of divisions of each side of the cube (characteristic of a curve)
        dim         dimension
        proto       curve prototype - sequence of "cubes" with integer coordinates (x_0,..,x_{d-1}), 0<=x_j<div
        base_maps   sequence of base maps (instances of BaseMap class)

    Immutable, hashable.
    For examples see get_hilbert_curve()
    """

    # В коде приняты следующие обозначения:
    # cube -- куб из прототипа
    # cnum -- номер куба в прототипе

    def __init__(self, proto, base_maps, dim = None, div = None):
        self.proto = tuple(tuple(cube) for cube in proto)
        self.base_maps = tuple(base_maps)
        self.dim = dim if dim is not None else self.get_dim()
        self.div = div if div is not None else self.get_div()
        self.genus = self.div ** self.dim

    def _data(self):
        return self.dim, self.div, self.proto, self.base_maps

    def __eq__(self, other):
        return self._data() == other._data()

    def __hash__(self):
        return hash(self._data())

    def get_div(self):
        return max(self.proto)[0]+1

    def get_dim(self):
        return len(self.proto[0])

    def changed(self, proto=None, base_maps=None):
        return type(self)(
            dim=self.dim,
            div=self.div,
            proto=(proto if proto is not None else self.proto),
            base_maps=(base_maps if base_maps is not None else self.base_maps),
        )
    
    #
    # Работа с базовыми преобразованиями и фракциями
    #

    def get_subdivision(self, k=1):
        """Get k-th subdivision of a curve."""
        N = self.div
        current_curve = self
        for _ in range(k):
            new_proto = []
            new_base_maps = []
            for cube, base_map in zip(current_curve.proto, current_curve.base_maps):
                proto = self.proto
                base_maps = self.base_maps

                if base_map.time_rev:
                    # в этой фракции прототип и базовые преобразования идут в обратном порядке
                    proto = reversed(proto)
                    base_maps = reversed(base_maps)

                for c in proto:
                    nc = base_map.apply_cube(N, c)
                    new_cube = [cj*N + ncj for cj, ncj in zip(cube, nc)]
                    new_proto.append(new_cube)

                # базовые преобразования для подраздедения:
                # пусть (cube, base_map) соответствуют i-й фракции
                # в ней мы взяли j-ю подфракцию (bm)
                # Какое преобразование переводит кривую в j-ю фракцию внутри i-й?
                # - сначала к исходной кривой мы применим bm, чтобы перевести её в j-ю фракцию,
                # - потом ко всей этой картинке применяем base_map, чтобы перевести всё в i-ю фракцию (base_map)
                # можно сделать наоборот:
                # - сначала кривую переводим в i-ю фракцию (base_map)
                # - применяем внутри i-й фракции преобразования для перехода в j-ю
                #   но там оно будет сопряженное: base_map * bm * base_map^{-1}, см. apply_base_map
                for bm in base_maps:
                    new_base_maps.append(base_map * bm)

            current_curve = type(self)(
                dim = self.dim,
                div = N*current_curve.div,
                proto = new_proto,
                base_maps = new_base_maps,
            )

        return current_curve


    #
    # Точки входа-выхода, время выхода на грань
    #

    def get_entrance(self):
        """Entrance of a curve, i.e. point f(0)."""
        start, period = self._get_cubes(0)
        return self._get_cube_limit(start, period)

    def get_exit(self):
        """Exit of a curve, i.e. point f(1)."""
        start, period = self._get_cubes(self.genus-1)
        return self._get_cube_limit(start, period)

    def _get_cubes(self, cnum):
        # находим последовательность вложенных кубов, которая получится, если в каждой фракции брать куб с номером cnum
        # возвращает пару (непериодическая часть, периодическая часть)
        cur_map = BaseMap.id_map(self.dim)  # current curve = cur_map * self
        cubes = []
        index = {}

        while True:
            cur_curve = self.apply_base_map(cur_map)
            cube = cur_curve.proto[cnum]

            cubes.append(cube)
            index[cur_map] = len(cubes)-1

            cur_map = cur_curve.base_maps[cnum] * cur_map
            if cur_map in index:
                idx = index[cur_map]
                return cubes[0:idx], cubes[idx:]

    def _get_cube_limit(self, start, period):
        # дана последовательность кубов, периодическая с некоторого момента, ищем предельную точку
        p = [0] * self.dim
        for j in range(self.dim):
            start_j = [x[j] for x in start]
            period_j = [x[j] for x in period]
            p[j] = self._get_periodic_sum(start_j, period_j, self.div)
        return tuple(p)

    def get_vertex_moments(self):
        # строим список всех вершин
        vertices = [[]]
        for j in range(self.dim):
            new_vertices = []
            for xj in [0,1]:
                for edge in vertices:
                    new_edge = edge + [xj]
                    new_vertices.append(new_edge)
            vertices = new_vertices
        moment = {}
        for v in vertices:
            moment[tuple(v)] = self.get_edge_touch(v)
        return moment

    def get_vertex_brkline(self):
        # строим ломаную (BRoKen LINE) из вершин фракций и моментов в порядке прохождения
        vm = self.get_vertex_moments()
        
        # ломаная для самой кривой
        brkline = [(v,vm[v]) for v in sorted(vm.keys(), key=lambda x:vm[x])]
        
        # локаная для кривой с обращенным временем
        brkline_rev = [(v,1-t) for v,t in reversed(brkline)]  # моменты при прохождении обратной кривой
        result = []
        for cnum, cube, base_map in zip(range(self.genus), self.proto, self.base_maps):
            if base_map.time_rev:
                curr_brkline = brkline_rev
            else:
                curr_brkline = brkline
            for v,t in curr_brkline:
                # есть вершина и момент
                # сначала поворачиваем, потом переносим во фракцию
                bv = base_map.apply_x(v)
                real_bv = [Fraction(cube[j] + bv[j], self.div) for j in range(self.dim)]
                real_t = Fraction(cnum + t, self.genus)
                result.append((real_bv, real_t))

        # удаляем дублирующиеся точки перехода
        new_result = []
        for i, r in enumerate(result):
            if i >= 1 and r == result[i-1]:
                # точка перехода
                continue
            new_result.append(r)
        assert len(new_result) == self.genus * (2**self.dim-1) + 1
        return result

    def get_edge_touch(self, edge):
        """Find moment of first edge touch.
        Edge is a tuple of {0,1,None} defining a set
        {(x_0,...,x_{d-1}): x_i==0 if e[i]==0, x_i==1 if e[i]==1, or arbitrary x[i] if e[i] is None.
        E.g., tuples (0,0,0) or (0,1,1) define vertices
        """
        curve = self
        cur_map = BaseMap.id_map(self.dim)  # curve = cur_map * self
        cnums = []
        index = {}
        while True:
            cnum = curve._get_edge_cnum(edge)
            bm = curve.base_maps[cnum]
            cnums.append(cnum)

            index[cur_map] = len(cnums)-1

            cur_map = bm * cur_map
            if cur_map in index:
                period_start = index[cur_map]
                break
            curve = curve.apply_base_map(bm)

        return self._get_time_limit(cnums[0:period_start], cnums[period_start:])

    def _get_edge_cnum(self, edge):
        # какой куб из прототипа первым касается грани
        N = self.div
        for cnum, cube in enumerate(self.proto):
            # проверяем, что куб касается грани
            touch = True
            for x, e in zip(cube, edge):
                if e is None:
                    continue
                elif e == 1 and x != (N-1):
                    touch = False
                    break
                elif e == 0 and x != 0:
                    touch = False
                    break
            if touch:
                return cnum

    def _get_time_limit(self, start, period):
        # задана начальная и периодическая последовательность номеров кубов, считаем время
        return self._get_periodic_sum(start, period, self.genus)

    # считаем сумму:
    #   s_0/d + s_1/d^2 + ... + s_{k-1}/d^k (непериодическая часть = start) +
    #       + s_k/d^{k+1} + ... + s_{k+m-1}/d^{k+m} + ... (периодическая часть = period)
    @staticmethod
    def _get_periodic_sum(start, period, d):
        n0 = 0
        d_power = 1
        for x in reversed(start):
            n0 += x * d_power
            d_power *= d
        t0 = Fraction(n0, d_power)

        np = 0
        dp_power = 1
        for x in reversed(period):
            np += x * dp_power
            dp_power *= d
        tp = Fraction(np, d_power*(dp_power - 1))

        return t0 + tp

    def check(self):
        """Check consistency of curve params."""
        d = self.dim
        n = self.div

        # dummy checks
        assert d > 0
        assert n > 0
        assert len(self.proto) == self.genus, 'bad proto length'

        for cube in self.proto:
            for j in range(d):
                assert 0 <= cube[j] < n, 'bad cube coordinates'
        assert len(set(self.proto)) == len(self.proto), 'non-unique cubes'

        curve_entrance = self.get_entrance()
        curve_exit = self.get_exit()

        # проверяем соответствие входов-выходов
        gates = [(None, curve_entrance)]  # пары (вход,выход)

        for cube, bm in zip(self.proto, self.base_maps):
            entrance_rel = bm.apply_x(curve_entrance)
            entrance = tuple(Fraction(c + e, n) for c, e in zip(cube, entrance_rel))

            exit_rel = bm.apply_x(curve_exit)
            exit = tuple(Fraction(c + e, n) for c, e in zip(cube, exit_rel))
            if bm.time_rev:
                gates.append((exit,entrance))
            else:
                gates.append((entrance,exit))
        gates.append((curve_exit, None))

        for i in range(len(gates)-1):
            if gates[i][1] != gates[i+1][0]:
                msg = 'exit does not correspond to entrance at ' + str(i)
                raise Exception(msg)

    def forget(self, allow_time_rev=False):
        entr = self.get_entrance()
        exit = self.get_exit()

        symmetries = []
        for bm in gen_constraint_cube_maps(self.dim, {entr: entr, exit: exit}):
            symmetries.append(bm)
            
        if allow_time_rev:
            for bm in gen_constraint_cube_maps(self.dim, {entr: exit, exit: entr}):
                symmetries.append(bm.reverse_time())

        return partial_fractal_curves.PartialFractalCurve(
            dim=self.dim,
            div=self.div,
            proto=self.proto,
            base_maps=[None for j in range(self.genus)],  # забыли BaseMap-ы!
            repr_maps=self.base_maps,  # base_map-ы стали лишь представителями
            symmetries=symmetries,
        )

    #
    # Стыки.
    #

    def gen_junctions(self):
        base_junctions = set(self.get_base_junction(cnum) for cnum in range(self.genus - 1))
        yield from self.gen_junctions_from_base(base_junctions)

    def get_junctions(self):
        return list(self.gen_junctions())

    #
    # Показатели гладкости кривой
    #

    def estimate_ratio_simple(self, ratio_func, subdivision=0):
        curve = self.get_subdivision(subdivision) if subdivision > 0 else self

        gates = []
        entr = [int(xj) for xj in self.get_entrance()]
        exit = [int(xj) for xj in self.get_exit()]
        for cnum, base_map in enumerate(curve.base_maps):
            piece_entr = base_map.apply_x(entr)
            piece_exit = base_map.apply_x(exit)
            if base_map.time_rev:
                piece_entr, piece_exit = piece_exit, piece_entr
            gates.append((piece_entr, piece_exit))
            
        max_r = None
        data = zip(range(curve.genus), curve.proto, gates)
        it = 0
        tot = curve.genus * (curve.genus - 1) // 2
        for pair in itertools.combinations(data, 2):
            it += 1
            if it % 100000 == 0:
                print('iter: {} of {} ({:.2f} %)'.format(it + 1, tot, 100 * (it + 1) / tot))
            cnum1, cube1, gate1 = pair[0]
            cnum2, cube2, gate2 = pair[1]

            t1 = cnum1  # only entrance
            x1 = [cube1j + entr1j for cube1j, entr1j in zip(cube1, gate1[0])]

            t2 = cnum2  # only entrance
            x2 = [cube2j + entr2j for cube2j, entr2j in zip(cube2, gate2[0])]

            dx = [x1j - x2j for x1j, x2j in zip(x1, x2)]

            r = FastFraction(*ratio_func(self.dim, dx, t2 - t1))
            if max_r is None or r > max_r:
                print('max_r:', max_r)
                max_r = r

        return max_r

    def init_pairs_tree(self):
        juncs = set(self.gen_junctions())
        for pair in super().init_pairs_tree():
            if pair.junc is None or pair.junc in juncs:
                yield pair

    def estimate_ratio(self, ratio_func, rel_tol_inv=100, upper_bound=None, max_iter=None, use_vertex_brkline=False, verbose=False):
        curr_up = None
        curr_lo = FastFraction(0, 1)

        pairs_tree_par = {}
        if use_vertex_brkline:
            vertex_brkline = [(x, t) for x, t in self.get_vertex_moments().items()]
            pairs_tree_par['brkline'] = IntegerBrokenLine(self.dim, vertex_brkline)
        pairs_tree = pieces.PairsTree(ratio_func, **pairs_tree_par)

        for pair in self.init_pairs_tree():
            pairs_tree.add_pair(pair)

        argmax = None
        for node in pairs_tree.data:
            if node.lo > curr_lo:
                curr_lo = node.lo
                argmax = node.argmax

        curr_up = pairs_tree.data[0].up
        pairs_tree.set_good_threshold(curr_lo)
        if verbose:
            print('start bounds: ', curr_lo, curr_up)

        tolerance = FastFraction(rel_tol_inv + 1, rel_tol_inv)
        iter_no = 0
        while curr_up > curr_lo * tolerance:
            iter_no += 1
            if not pairs_tree.data:
                break
            if max_iter is not None and iter_no > max_iter:
                # оставляем оценки, полученные на текущий момент
                break
            
            pairs_tree.divide()

            for node in pairs_tree.data:
                if node.lo > curr_lo:
                    curr_lo = node.lo
                    argmax = node.argmax
                    if verbose:
                        print('new lower bound: ', curr_lo, curr_up)
            pairs_tree.set_good_threshold(curr_lo)

            new_up = pairs_tree.data[0].up
            if new_up < curr_up:
                if verbose:
                    print('new upper bound: ', curr_lo, new_up)
                curr_up = new_up

        res = {'up': curr_up, 'lo': curr_lo}
        if argmax is not None:
            res['argmax'] = argmax

        return res


class IntegerBrokenLine:
    def __init__(self, dim, brkline):
        denoms = set()
        for x, t in brkline:
            if isinstance(t, Fraction):
                denoms.add(t.denominator)
            for xj in x:
                if isinstance(xj, Fraction):
                    denoms.add(xj.denominator)
        lcm = get_lcm(denoms)
        mx = lcm
        mt = lcm**dim
        self.points = [([int(xj * mx) for xj in x], int(t * mt)) for x, t in brkline]
        self.mx = mx
        self.mt = mt
