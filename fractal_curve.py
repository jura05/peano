# coding: utf-8

from base_map import BaseMap
from fractions import Fraction

class FractalCurve:
    """Class representing fractal peano curve in [0,1]^d.
    Params:
        div         positive integer, number of divisions of each side of the cube (characteristic of a curve)
        dim         dimension
        proto       curve prototype - sequence of "cubes" with integer coordinates (x_0,..,x_{d-1}), 0<=x_j<div
          or
        chain_code  'ikjKJ'
        base_maps   sequence of base maps (instances of BaseMap class)

    For examples see get_hilbert_curve()
    """
    def __init__(self, dim, div, base_maps, proto=None, chain_code=None):
        self.div = div
        self.dim = dim
        if chain_code is not None:
            proto = _chain2proto(dim, chain_code)

        self.proto = tuple(proto)
        self.base_maps = tuple(base_maps)

    def get_entrance(self):
        return self._get_limit(0)

    def get_exit(self):
        return self._get_limit(self.genus()-1)

    def reverse(self):
        return type(self)(
            dim = self.dim,
            div = self.div,
            proto = reversed(self.proto),
            base_maps = reversed(self.base_maps),
        )

    def apply_base_map(self, base_map):
        """Apply base map to a fractal curve, return new curve."""
        if base_map.time_rev:
            curve = self.reverse()
            cube_map = base_map.cube_map()
        else:
            curve = self
            cube_map = base_map

        inv = cube_map.inverse()
        return type(self)(
            dim = self.dim,
            div = self.div,
            proto = [cube_map.apply_cube(curve.div, cube) for cube in curve.proto],
            base_maps = [cube_map * bm * inv for bm in curve.base_maps],
        )

    def get_subcurve(self, i):
        """Get fraction as a curve."""
        return self.apply_base_map(self.base_maps[i])

    def get_subdivision(self):
        """Get divided curve (with squared genus)."""
        N = self.div
        new_proto = []
        new_base_maps = []
        for cube, base_map in zip(self.proto, self.base_maps):
            if not base_map.time_rev:
                for c in self.proto:
                    nc = base_map.apply_cube(N, c)
                    new_cube = [cj*N + ncj for cj, ncj in zip(cube, nc)]
                    new_proto.append(new_cube)

                for bm in self.base_maps:
                    new_base_maps.append(base_map * bm)
            else:
                for c in reversed(self.proto):
                    nc = base_map.apply_cube(N, c)
                    new_cube = [cj*N + ncj for cj, ncj in zip(cube, nc)]
                    new_proto.append(new_cube)

                for bm in reversed(self.base_maps):
                    new_base_maps.append(base_map * bm)

        return type(self)(
            dim = self.dim,
            div = N**2,
            proto = new_proto,
            base_maps = new_base_maps,
        )

    def _get_edge_cnum(self, edge):
        """Edge is a tuple of {0,1,None} defining a set
        {(x_0,...,x_{d-1}): x_i==0 if e[i]==0, x_i==1 if e[i]==1, or arbitrary x[i] if e[i] is None."""

        N = self.div
        def check_touch(cube, edge, N):
            for x, e in zip(cube, edge):
                if e is None:
                    continue
                if e and x != (N-1):
                    return False
                if (not e) and x != 0:
                    return False
            return True

        for i, cube in enumerate(self.proto):
            if check_touch(cube, edge, N):
                return i


    def get_edge_touch(self, edge):
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

        return self._get_time_limit(cnums, period_start)


    def _get_limit(self, cnum):
        div = self.div
        genus = self.genus()

        cur_map = BaseMap(dim=self.dim)  # current curve = cur_map * self
        seen = set()
        cube_period = []

        # составляем номера под-кубов (в изначальной ориентации)
        # каждый следующий куб получается умножением слева на base_map
        while True:
            if cur_map.time_rev:
                cur_cnum = genus-1-cnum
            else:
                cur_cnum = cnum
            cube = cur_map.apply_cube(div, self.proto[cur_cnum])
            cube_period.append(cube)
            cur_map = cur_map * self.base_maps[cur_cnum]
            if cur_map in seen:
                break
            seen.add(cur_map)

        return self._get_periodic_limit(cube_period)

    def _get_time_limit(self, cnums, period_start):
        g = self.genus()
        t0 = Fraction(0, 1)

        cnum_base = cnums[0:period_start]
        for i, cnum in enumerate(cnum_base):
            t0 += Fraction(cnum, g**(i+1))

        cnum_period = cnums[period_start:]
        m = len(cnum_period)
        tp = 0
        for i, cnum in enumerate(cnum_period):
            tp += cnum * g**(m-1-i) 

        return t0 + Fraction(tp, g**period_start*(g**m - 1))

    # суммируем геометрическую прогрессию (левые нижние углы кубов) со знаменателем 1/div
    # ТУТ ТОЖЕ НЕ ПЕРИОДИЧЕСКИ МОЖЕТ БЫТЬ
    def _get_periodic_limit(self, cube_period):
        div = self.div
        m = len(cube_period)
        pt = [0] * self.dim
        for i, cube in enumerate(cube_period):
            for j, c in enumerate(cube):
                pt[j] += c * div**(m-1-i) 

        return tuple(Fraction(x, div**m - 1) for x in pt)

    def genus(self):
        """Fractal genus of the curve."""
        return self.div ** self.dim

    def check(self):
        """Check consistency of curve params."""
        n = self.div
        d = self.dim

        # dummy checks
        assert n > 0
        assert len(self.proto) == self.genus(), 'bad proto length'

        for cube in self.proto:
            for j in range(d):
                assert 0 <= cube[j] < n, 'bad cube coordinates'
        sqset = set(tuple(cube) for cube in self.proto)
        assert len(sqset) == len(self.proto), 'non-unique cubes'

        entrance = self.get_entrance()
        exit = self.get_exit()

        # проверяем соответствие входов-выходов
        gates = []  # пары (вход,выход), реальные координаты в кубе, умноженные на div
        entrance_n = tuple(n*entrance[j] for j in range(d))  # начальный вход
        gates.append((None, entrance_n))
        for cube, bm in zip(self.proto, self.base_maps):
            entrance_pos = bm.apply_x(entrance)
            entrance_n = tuple(cube[j] + entrance_pos[j] for j in range(d))
            exit_pos = bm.apply_x(exit)
            exit_n = tuple(cube[j] + exit_pos[j] for j in range(d))
            if bm.time_rev:
                gates.append((exit_n,entrance_n))
            else:
                gates.append((entrance_n,exit_n))
        exit_n = tuple(n*exit[j] for j in range(d))
        gates.append((exit_n, None))

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
        return (bm1_inv.apply_vec(delta), bm1_inv * bm2)


# some utility functions
def _chain2proto(dim, chain_code, start=None):
    """Convert chain code like 'ijK' to curve prototype."""
    assert dim <= 6
    letters = 'ijklmn'
    l2v = {}
    for k in range(dim):
        l = letters[k]
        v = [0] * dim
        v[k] = 1
        l2v[l] = v
        l2v[l.upper()] = [-x for x in v]

    if start is None:
        start = (0,) * dim

    cube = start
    proto = [cube]
    for l in chain_code:
        diff = l2v[l]
        cube = [c + d for c, d in zip(cube, diff)]
        proto.append(cube)

    return proto
