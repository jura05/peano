from fractions import Fraction

from . import fuzzy_poly_curves
from .base_maps import BaseMap, Spec


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
