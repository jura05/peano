from fractions import Fraction

from .fast_fractions import FastFraction
from .base_maps import BaseMap, Spec
from .utils import get_lcm
from . import fuzzy_poly_curves
from . import pieces


class PolyCurve(fuzzy_poly_curves.FuzzyPolyCurve):

    #
    # Точки входа-выхода
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
        cur_spec = Spec(base_map=BaseMap.id_map(self.dim), pnum=self.pnum)  # current curve = cur_spec * self
        cubes = []
        index = {}

        while True:
            cur_curve = cur_spec * self
            cube = cur_curve.proto[cnum]

            cubes.append(cube)
            index[cur_spec] = len(cubes)-1

            cur_spec = cur_curve.specs[cnum] * cur_spec

            if cur_spec in index:
                idx = index[cur_spec]
                return cubes[0:idx], cubes[idx:]

    def _get_cube_limit(self, start, period):
        # дана последовательность кубов, периодическая с некоторого момента, ищем предельную точку
        p = [0] * self.dim
        for j in range(self.dim):
            start_j = [x[j] for x in start]
            period_j = [x[j] for x in period]
            p[j] = self._get_periodic_sum(start_j, period_j, self.div)
        return tuple(p)

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

    #
    # Стыки
    #
    # Для полностью заданной кривой можно определить все стыки

    def gen_base_junctions(self):
        seen = set()
        for pnum in range(self.pattern_count):
            for cnum in range(self.genus - 1):
                junc = self.get_base_junction(cnum=cnum, pnum=pnum)
                if junc not in seen:
                    yield junc
                    seen.add(junc)

    def gen_junctions(self):
        yield from self.gen_junctions_from_base(list(self.gen_base_junctions()))

    def get_junctions(self):
        return list(self.gen_junctions())

    #

    # кандидаты в self.base_maps[cnum]
    def gen_allowed_specs(self, pnum, cnum):
        yield self.patterns[pnum].specs[cnum]

    def init_pairs_tree(self):
        G = self.genus
        for junc in self.gen_auto_junctions():
            for cnum1 in range(G):
                for cnum2 in range(cnum1 + 2, G):
                    pos1 = self.get_piece_position(junc.spec1.pnum, cnum1)
                    pos2 = self.get_piece_position(junc.spec2.pnum, cnum2)
                    yield pieces.CurvePieceBalancedPair(self, junc, pos1, pos2)

        for junc in self.gen_junctions():
            last_cnum1 = 0 if junc.spec1.base_map.time_rev else G - 1
            first_cnum2 = G - 1 if junc.spec2.base_map.time_rev else 0
            for cnum1 in range(G):
                for cnum2 in range(G):
                    if (cnum1, cnum2) == (last_cnum1, first_cnum2):
                        continue
                    pos1 = self.get_piece_position(junc.spec1.pnum, cnum1)
                    pos2 = self.get_piece_position(junc.spec2.pnum, cnum2)
                    yield pieces.CurvePieceBalancedPair(self, junc, pos1, pos2)

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
