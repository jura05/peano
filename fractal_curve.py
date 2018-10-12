# coding: utf-8

from heapq import heappop, heappush
from fractions import Fraction
from math import gcd
from collections import Counter

from base_map import BaseMap, PieceMap

class FractalCurve:
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

    def __init__(self, dim, div, proto, base_maps):
        self.dim = dim
        self.div = div
        self.proto = tuple(tuple(cube) for cube in proto)
        self.base_maps = tuple(base_maps)

    def _data(self):
        return self.dim, self.div, self.proto, self.base_maps

    def __eq__(self, other):
        return self._data() == other._data()

    def __hash__(self):
        return hash(self._data())

    def genus(self):
        """Fractal genus of the curve."""
        return self.div ** self.dim

    #
    # Работа с базовыми преобразованиями и фракциями
    #

    def reverse(self):
        """Reverse time in a curve."""
        return type(self)(
            dim = self.dim,
            div = self.div,

            # прототип проходится в обратном порядке
            proto = reversed(self.proto),

            # базовые преобразования проходятся в обратном порядке
            # сами по себе они не меняются:
            #   - если обращения времени не было, то его и не будет
            #   - изометрия куба не меняется, т.к. время не играет роли
            base_maps = reversed(self.base_maps),
        )

    def apply_base_map(self, base_map):
        """Apply base map to a fractal curve, return new curve."""

        # можно разложить базовое преобразование в произведение (коммутирующих) 
        # преобразований: обращение времени с тождественной изометрией  +  изометрии куба
        if base_map.time_rev:
            curve = self.reverse()
            cube_map = base_map.cube_map()
        else:
            curve = self
            cube_map = base_map

        # применяем изометрию куба
        inv = cube_map.inverse()
        return type(self)(
            dim = self.dim,
            div = self.div,

            # прототип подвергается изометрии
            proto = [cube_map.apply_cube(curve.div, cube) for cube in curve.proto],

            # базовые преобразования сопрягаются: действительно, чтобы получить
            # из преобразованной кривой её фракцию, можно сделать так:
            # - сначала вернуться к исходной кривой (inv)
            # - применить преобразование исходной кривой для перехода к фракции (bm)
            # - перейти к преобразованной кривой (cube_map)
            base_maps = [cube_map * bm * inv for bm in curve.base_maps],
        )

    def get_fraction(self, cnum):
        """Get fraction as a curve."""
        return self.apply_base_map(self.base_maps[cnum])

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
        start, period = self._get_cubes(self.genus()-1)
        return self._get_cube_limit(start, period)

    def _get_cubes(self, cnum):
        # находим последовательность вложенных кубов, которая получится, если в каждой фракции брать куб с номером cnum
        # возвращает пару (непериодическая часть, периодическая часть)
        cur_map = BaseMap(dim=self.dim)  # current curve = cur_map * self
        cubes = []
        index = {}

        while True:
            if cur_map.time_rev:
                cur_cnum = -cnum
            else:
                cur_cnum = cnum
            cube = cur_map.apply_cube(self.div, self.proto[cur_cnum])
            cubes.append(cube)
            index[cur_map] = len(cubes)-1

            # сначала переходим из исходной кривой во фракцию, потом всё отображаем в текущую кривую
            # можно было бы хранить cur_curve и писать cur_map = cur_curve.base_maps[cnum] * cur_map
            cur_map = cur_map * self.base_maps[cur_cnum]

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
        for cnum, cube, base_map in zip(range(self.genus()), self.proto, self.base_maps):
            if base_map.time_rev:
                curr_brkline = brkline_rev
            else:
                curr_brkline = brkline
            for v,t in curr_brkline:
                # есть вершина и момент
                # сначала поворачиваем, потом переносим во фракцию
                bv = base_map.apply_x(v)
                real_bv = [Fraction(cube[j] + bv[j], self.div) for j in range(self.dim)]
                real_t = Fraction(cnum + t, self.genus())
                result.append((real_bv, real_t))

        # удаляем дублирующиеся точки перехода
        new_result = []
        for i, r in enumerate(result):
            if i >= 1 and r == result[i-1]:
                # точка перехода
                continue
            new_result.append(r)
        assert len(new_result) == self.genus() * (2**self.dim-1) + 1
        return result

    def get_edge_touch(self, edge):
        """Find moment of first edge touch.
        Edge is a tuple of {0,1,None} defining a set
        {(x_0,...,x_{d-1}): x_i==0 if e[i]==0, x_i==1 if e[i]==1, or arbitrary x[i] if e[i] is None.
        E.g., tuples (0,0,0) or (0,1,1) define vertices
        """
        curve = self
        cur_map = BaseMap(dim=self.dim)  # curve = cur_map * self
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
        return self._get_periodic_sum(start, period, self.genus())

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
        assert len(self.proto) == self.genus(), 'bad proto length'

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
            assert gates[i][1] == gates[i+1][0], 'exit does not correspond to entrance'

    #
    # Стыки. Реализовано для кривых без обращения времени
    #

    def get_junctions(self):
        """Junction is a pair (delta, base_map). Get all junctions of a curve."""
        if any(bm.time_rev for bm in self.base_maps):
            raise Exception("get_junctions not implemented for time reverse!")

        junctions = set()

        for i in range(self.genus()-1):
            cube = self.proto[i]
            next_cube = self.proto[i+1]
            delta = tuple(nc-c for nc, c in zip(next_cube, cube))
            junctions.add(self._get_std_junction(delta, self.base_maps[i], self.base_maps[i+1]))

        to_derive = list(junctions)
        while to_derive:
            junction = to_derive.pop()
            dj = self.get_derived_junction(junction)
            if dj not in junctions:
                junctions.add(dj)
                to_derive.append(dj)

        return junctions

    def get_derived_junction(self, junction):
        delta, base_map = junction
        cube1 = self.proto[-1]
        cube2 = base_map.apply_cube(self.div, self.proto[0])
        der_delta = tuple(delta[k]*self.div + cube2[k] - cube1[k] for k in range(self.dim))
        return self._get_std_junction(der_delta, self.base_maps[-1], base_map * self.base_maps[0])

    # поворачиваем, чтобы обеспечить тождественное преобразование на первой фракции
    @staticmethod
    def _get_std_junction(delta, bm1, bm2):
        bm1_inv = bm1.inverse()
        return bm1_inv.apply_vec(delta), bm1_inv * bm2

    #
    # Показатели гладкости кривой
    #

    class CurveBalancedPair:
        # работаем с парой фракций кривой
        # первую приводим к стандартной ориентации (как для стыка)
        # одна из фракций может быть "раздута" (не более чем на одно подразделение), другая - стандартная:
        # compare=0:  первая фракция [0,1]^d, вторая: delta_x + [0,1]^d
        # compare=-1: первая фракция [0,1]^d, вторая: delta_x + [0,N]^d
        # compare=1:  первая фракция [0,N]^d, вторая: delta_x + [0,1]^d <--- сейчас не используется
        def __init__(self, delta_x, delta_t, base_map, compare=0):
            self.delta_x = delta_x
            self.delta_t = delta_t  # время от начала первой фракции до начала второй (за 1 времени заметаем [0,1]^d)
            self.base_map = base_map  # как повёрнута вторая фракция
            self.compare = compare

    def divide_pair(self, pair):
        """Divide CurveBalancedPair.
        Will process orig_map attribute: it maps original pair to new one
        """

        d = self.dim
        N = self.div
        G = self.genus()

        delta_x = pair.delta_x
        delta_t = pair.delta_t
        base_map = pair.base_map
        has_orig_map = hasattr(pair, 'orig_map')
        orig_maps = []
        pairs = []
        if pair.compare == 0:
            # делим стандартную фракцию, нужно раздуть всё и стандартизовать!
            for cnum, cube, bm in zip(range(G), self.proto, self.base_maps):
                scaled_delta_x = tuple(delta_x[j] * N - cube[j] for j in range(d))
                scaled_delta_t = delta_t * G - cnum

                # применим inv_map, стандартизуя первую фракцию
                inv_map = bm.inverse()
                new_base_map = inv_map * base_map
                new_delta_x  = inv_map.apply_cube_start(cube_start=scaled_delta_x, cube_length=N)
                new_pair = self.CurveBalancedPair(
                    delta_x=new_delta_x,
                    delta_t=scaled_delta_t,
                    base_map=new_base_map,
                    compare=-1,
                )
                if has_orig_map:
                    orig_map = PieceMap(
                        base_map=inv_map,
                        shift=tuple(-cube[j] for j in range(d)),
                        scale=N,
                        time_scale=G,
                        time_shift=-cnum,
                    )
                    orig_maps.append(orig_map)
                pairs.append(new_pair)
        elif pair.compare == -1:
            # делим вторую, раздутую фракцию, стандартизовать не нужно!
            for cnum, cube, bm in zip(range(G), self.proto, self.base_maps):
                new_cube = base_map.apply_cube(N, cube)
                new_base_map = base_map * bm
                new_pair = self.CurveBalancedPair(
                    delta_x=tuple(delta_x[j] + new_cube[j] for j in range(d)),
                    delta_t=delta_t + cnum,
                    base_map=new_base_map,
                    compare=0,
                )
                if has_orig_map:
                    orig_maps.append(PieceMap.id_map(d))
                pairs.append(new_pair)
        elif pair.compare == 1:
            # делим первую, раздутую фракцию, нужно только повернуть
            for cnum, cube, bm in zip(range(G), self.proto, self.base_maps):
                new_delta_x = tuple(delta_x[j] - cube[j] for j in range(d))
                new_delta_t = delta_t - cnum

                # применим inv_map, стандартизуя первую фракцию
                inv_map = bm.inverse()
                new_base_map = inv_map * base_map
                new_delta_x  = inv_map.apply_cube_start(cube_start=new_delta_x, cube_length=1)
                new_pair = self.CurveBalancedPair(
                    delta_x=new_delta_x,
                    delta_t=new_delta_t,
                    base_map=new_base_map,
                    compare=0,
                )
                if has_orig_map:
                    orig_map = PieceMap(
                        base_map=inv_map,
                        shift=tuple(-cube[j] for j in range(d)),
                        scale=1,
                        time_scale=1,
                        time_shift=-cnum,
                    )
                    orig_maps.append(orig_map)
                pairs.append(new_pair)

        if has_orig_map:
            for new_pair, new_orig_map in zip(pairs, orig_maps):
                new_pair.orig_map = new_orig_map * pair.orig_map

        return pairs

    def estimate_ratio(self, ratio_func, junctions=None, upper_bound=None, max_iter=None, rel_tol=None, find_argmax=False, verbose=0):
        """Estimate ratio of a curve for given junctions.
        Params:
            ratio_func      generic function for ratio; args: d, delta_v, delta_t
                                we assume it is d-homogeneous
            junctions       list of junctions (None for all junctions, [None] for self reduced ratio)
            upper_bound     apriori upper bound for ratio; we stop if higher ratio is found
            max_iter        subj
            rel_tol         relative error tolerance
            verbose         print every iteration
            find_argmax     subj
        Returns dict with keys:
            lower_bound     bounds for ratio
            upper_bound
            argmax          dict with keys {v1,t1,v2,t2,junc}, describing ratio argmax (if options find_argmax is set)
            stats           counter with some statistics
        """

        # Алгоритм работы:
        # сравниваем пары фракций:
        # 1) оцениваем отношение снизу (через вершины, например)
        # если отношение больше текущей верхней оценки - выходим
        # 2) оцениваем отношение сверху (по прототипу)
        # - если меньше текущей нижней оценки - останавливаемся, не подразбиваем
        # - иначе, подразбиваем обе фракции

        if any(bm.time_rev for bm in self.base_maps):
            raise Exception("Not implemented for curves with time reversal!")
        if max_iter is None and rel_tol is None:
            raise Exception("Define max_iter or rel_tol!")

        if junctions is None:
            junctions = self.get_junctions()
            junctions.add(None)

        d = self.dim
        N = self.div
        G = self.genus()

        # тут можно использовать другие ломаные
        self_brkline = self.get_vertex_brkline()

        # приводим ломаную к целым числам; множители lcm_x, lcm_t применим лишь в самом конце вычислений
        denoms = []
        for v, t in self_brkline:
            denoms.append(t.denominator)
            denoms += [x.denominator for x in v]
        lcm_x = get_lcm(denoms)
        lcm_t = lcm_x**d
        int_brkline = []
        for v, t in self_brkline:
            new_t = int(t * lcm_t)
            new_v = tuple(int(x * lcm_x) for x in v)
            int_brkline.append((new_v, new_t))

        argmax = None
        stats = Counter()
        curr_lower_bound = 0
        curr_upper_bound = None

        rich_pairs = []  # пары кривых с доп. информацией
        entry_count = 0
        for junc in junctions:
            if junc is None:
                delta_x = (0,) * d
                delta_t = 0
                junc_base_map = BaseMap(dim=self.dim)
            else:
                delta_x, junc_base_map = junc
                delta_t = 1

            # считаем максимально возможное редуцированное отношение
            max_dv = tuple(1 + abs(delta_x[j]) for j in range(d))
            min_dt = Fraction(1, G)
            dummy_upper_bound = ratio_func(d, max_dv, min_dt)

            start_pair = self.CurveBalancedPair(
                delta_x=delta_x,
                delta_t=delta_t,
                base_map=junc_base_map,
                compare=0,
            )
            if find_argmax:
                start_pair.orig_map = PieceMap.id_map(d)

            # дробим каждую из частей (дробление только одной не уменьшает кол-во итераций)
            junc_pairs = []
            for pair in self.divide_pair(start_pair):
                for sub_pair in self.divide_pair(pair):
                    if sub_pair.delta_t <= 1:
                        continue
                    junc_pairs.append(sub_pair)

            for pair in junc_pairs:
                rich_pair = {'pair': pair, 'upper_bound': dummy_upper_bound, 'max_subdivision': 1}
                if find_argmax:
                    rich_pair['junc'] = junc
                entry_count += 1
                priority = -dummy_upper_bound
                heappush(rich_pairs, (priority, entry_count, rich_pair))

        seen_pairs = set()
        while rich_pairs:
            stats['iter'] += 1
            if max_iter is not None and stats['iter'] >= max_iter:
                break

            rich_pair = heappop(rich_pairs)[-1]
            pair = rich_pair['pair']
            delta_x = pair.delta_x
            delta_t = pair.delta_t
            base_map = pair.base_map

            # учитываем, что одна из фракций "раздута"
            if pair.compare == 0:
                mx1, mx2 = 1, 1
            elif pair.compare == 1:
                mx1, mx2 = N, 1
            else:
                mx1, mx2 = 1, N
            mt1, mt2 = mx1**d, mx2**d

            pair_position = (delta_x, delta_t, base_map, pair.compare)
            if pair_position in seen_pairs:
                # пара уже встречалась, новых оценок нам не даст
                stats['seen_pair'] += 1
                continue
            else:
                seen_pairs.add(pair_position)

            stats['max_subdivision'] = max(stats['max_subdivision'], rich_pair['max_subdivision'])

            if verbose:
                print('iter: {}; max_subdivision: {}; ratio: {} <= X <= {}'.format(
                    stats['iter'],
                    stats['max_subdivision'],
                    float(curr_lower_bound),
                    (float(curr_upper_bound) if curr_upper_bound is not None else '?'),
                ))

            # могла обновиться нижняя оценка, и пара больше не актуальна!
            if rich_pair['upper_bound'] <= curr_lower_bound:
                stats['stop_early'] += 1
                continue

            #
            # Обновим верхнюю границу (curr_upper_bound)
            #

            # нужно найти максимальное расстояние между точками двух кубов
            # точки первого куба: 0 <= x[j] <= mx1, второго: delta_x[j] <= x[j] <= delta_x[j] + mx2
            max_dv = tuple(max(abs(mx1 - dx), abs(dx + mx2)) for dx in delta_x)
            min_dt = delta_t - mt1
            up_ratio = ratio_func(d, max_dv, min_dt)

            # здесь мы используем, что rich_pairs это heap по priority = (-upper_bound)
            curr_upper_bound = max(up_ratio, rich_pairs[0][-1]['upper_bound']) if rich_pairs else up_ratio

            if up_ratio <= curr_lower_bound:
                # нам эта пара больше не интересна!
                stats['stop_divide'] += 1
                continue

            #
            # Обновим нижнюю оценку (curr_lower_bound)
            #

            # поскольку приводили к целым числам, нужно всё умножить на lcm_x (соотв., lcm_t)
            if mx1 == 1:
                brkline1 = int_brkline
            else:
                brkline1 = []
                for v, t in int_brkline:
                    scaled_v = tuple(v[j] * mx1 for j in range(d))
                    scaled_t = t * mt1
                    brkline1.append((scaled_v, scaled_t))

            brkline2 = []
            for v, t in int_brkline:
                # повернуть + масштабировать
                rot_v = base_map.apply_x2(v, lcm_x)
                if mx2 == 1:
                    scaled_v = rot_v
                    scaled_t = t
                else:
                    scaled_v = tuple(rot_v[j] * mx2 for j in range(d))
                    scaled_t = t * mt2
                brkline2.append((scaled_v, scaled_t))

            # так как в ломаных уже "сидят" lcm-ы, умножаем только дельты
            scaled_delta_x = tuple(dx * lcm_x for dx in delta_x)
            scaled_delta_t = delta_t * lcm_t
            for v2, t2 in brkline2:
                v2_shifted = tuple(v2[j] + scaled_delta_x[j] for j in range(d))
                t2_shifted = t2 + scaled_delta_t
                for v1, t1 in brkline1:
                    dv = tuple(v2_shifted[j] - v1[j] for j in range(d))
                    dt = t2_shifted - t1
                    ratio = ratio_func(d, dv, dt)
                    if ratio > curr_lower_bound:
                        curr_lower_bound = ratio
                        if find_argmax:
                            inv = pair.orig_map.inverse()
                            # нужно вернуться обратно
                            # сначала уберём lcm-ы
                            v1_frac = tuple(Fraction(v1[j], lcm_x) for j in range(d))
                            v2_frac = tuple(Fraction(v2_shifted[j], lcm_x) for j in range(d))

                            t1_frac = Fraction(t1, lcm_t)
                            t2_frac = Fraction(t2_shifted, lcm_t)

                            # теперь обратное преобразование
                            (orig_v1, orig_t1) = inv.apply(v1_frac, t1_frac)
                            (orig_v2, orig_t2) = inv.apply(v2_frac, t2_frac)

                            # v2, t2 считаем относительно второй фракции
                            junc = rich_pair['junc']
                            if junc is not None:
                                orig_v2 = tuple(orig_v2[j] - junc[0][j] for j in range(d))
                                orig_t2 -= 1
                            argmax = {
                                'v1': orig_v1, 't1': orig_t1,
                                'v2': orig_v2, 't2': orig_t2,
                                'junc': junc,
                            }

            if upper_bound is not None and curr_lower_bound > upper_bound:
                break

            if rel_tol is not None and curr_upper_bound <= (1 + rel_tol) * curr_lower_bound:
                break

            for sub_pair in self.divide_pair(pair):
                max_subdivision = rich_pair['max_subdivision']
                if pair.compare == 0:
                    max_subdivision += 1  # если фракции были равны, номер подразбиения увеличился
                new_rich_pair = {'pair': sub_pair, 'upper_bound': up_ratio, 'max_subdivision': max_subdivision}
                if 'junc' in rich_pair:
                    new_rich_pair['junc'] = rich_pair['junc']
                entry_count += 1
                heappush(rich_pairs, (-up_ratio, entry_count, new_rich_pair))

        return {
            'lower_bound': curr_lower_bound,
            'upper_bound': curr_upper_bound,
            'argmax': argmax,
            'stats': stats,
        }

def get_lcm(iterable):
    """Least common multiple of integer sequence."""
    lcm = 1
    for x in iterable:
        lcm = (lcm * x) // gcd(lcm, x)
    return lcm
