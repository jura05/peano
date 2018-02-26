# coding: utf-8

from fractions import Fraction

class FractalCurve:
    """Class representing fractal peano curve in [0,1]^d.
    Params:
        div         positive integer, number of divisions of each side of the cube (characteristic of a curve)
        dim         dimension (default: 2)
        proto       curve prototype - sequence of "cubes" with integer coordinates (x_0,..,x_{d-1}), 0<=x_j<div
        base_maps   sequence of base maps (instances of BaseMap class)
        entrance    entrance point (default: None, i.e. determine by proto & base_maps)
        exit        exit point (default: None, see entrance)

    For examples see get_hilbert_curve()
    """
    def __init__(self, div, proto, base_maps, entrance=None, exit=None, dim=2):
        self.div = div
        self.dim = dim
        self.proto = proto
        self.base_maps = base_maps
        self.entrance = entrance if entrance is not None else self.get_entrance()
        self.exit = exit if exit is not None else self.get_exit()

    def get_entrance(self):
        return self._get_limit(self.div, self.proto[0], self.base_maps[0])

    def get_exit(self):
        return self._get_limit(self.div, self.proto[-1], self.base_maps[-1])

    def apply_base_map(self, base_map):
        """Apply base map to a fractal curve, return new curve (TODO: use time_rev!)."""
        inv = base_map.inverse()
        return type(self)(
            div = self.div,
            dim = self.dim,
            proto = [base_map.apply_cube(self.div, cube) for cube in self.proto],
            base_maps = [base_map * bm * inv for bm in self.base_maps],
        )

    def get_subcurve(self, i):
        """Get fraction as a curve."""
        return self.apply_base_map(self.base_maps[i])

    def get_subdivision(self):
        """Get divided curve (with squared genius)."""
        N = self.div
        new_proto = []
        new_base_maps = []
        for cube, base_map in zip(self.proto, self.base_maps):
            for c in self.proto:
                nc = base_map.apply_cube(N, c)
                new_cube = [cj*N + ncj for cj, ncj in zip(cube, nc)]
                new_proto.append(new_cube)

            for bm in self.base_maps:
                new_base_maps.append(bm * base_map)

        return type(self)(
            dim = self.dim,
            div = N**2,
            proto = new_proto,
            base_maps = new_base_maps,
        )

    @staticmethod
    def _get_limit(div, cube, base_map):
        cube_period = [cube]
        cur_map = base_map
        # составляем номера под-кубов (в изначальной ориентации)
        # каждый следующий куб получается умножением слева на base_map
        while True:
            new_cube = cur_map.apply_cube(div, cube)
            cur_map = base_map * cur_map
            if tuple(new_cube) == tuple(cube):
                break
            cube_period.append(new_cube)

        # суммируем геометрическую прогрессию (левые нижние углы кубов) со знаменателем div**(-m)
        m = len(cube_period)
        x = 0
        y = 0
        for i, c in enumerate(cube_period):
            x += c[0] * div**(m-1-i) 
            y += c[1] * div**(m-1-i)

        lx = Fraction(x, div**m - 1)
        ly = Fraction(y, div**m - 1)
        return (lx,ly)

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

        for j in range(d):
            assert 0 <= self.entrance[j] < n , 'bad entrance'
            assert 0 <= self.exit[j] < n, 'bad exit'

        # проверяем соответствие входов-выходов
        gates = []  # пары (вход,выход), реальные координаты в кубе, умноженные на div
        entrance_n = tuple(n*self.entrance[j] for j in range(d))  # начальный вход
        gates.append((None, entrance_n))
        for cube, bm in zip(self.proto, self.base_maps):
            entrance_pos = bm.apply_x(self.entrance)
            entrance_n = tuple(cube[j] + entrance_pos[j] for j in range(d))
            exit_pos = bm.apply_x(self.exit)
            exit_n = tuple(cube[j] + exit_pos[j] for j in range(d))
            if bm.time_rev:
                gates.append((exit_n,entrance_n))
            else:
                gates.append((entrance_n,exit_n))
        exit_n = tuple(n*self.exit[j] for j in range(d))
        gates.append((exit_n, None))

        for i in range(len(gates)-1):
            assert gates[i][1] == gates[i+1][0], 'exit does not correspond to entrance'

    def get_junctions(self):
        """Junction is a pair (delta, base_map). Get all junctions of a curve."""
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


    # здесь пока не учитываем time_rev !!!
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
