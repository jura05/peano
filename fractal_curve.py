# coding: utf-8

from fractions import Fraction

class BaseMap:
    """Base map: isometry of cube and (possibly) time reversal.
    Params:
        perm, flip: params, which define isometry of cube
                      perm = [k_0,...,k_{d-1}]
                      flip = [b_0,...,b_{d-1}]
                      define the map
                      (x_0,...,x_{d-1}) -> (f(b_0,x_{k_0}),...,f(b_{d-1},x_{k_{d-1}})),
                      where f(b,x)=x if b is False, and 1-x otherwise
        time_rev    time reversal (boolean), default: False
    """

    @classmethod
    def id_map(cls, dim):
        perm = list(range(dim))
        flip = [False]*dim
        time_rev = False
        return cls(perm, flip, time_rev)

    def __init__(self, perm, flip, time_rev=False):
        assert len(perm) == len(flip)
        self.dim = len(perm)
        self.perm = perm
        self.flip = [bool(b) for b in flip]
        self.time_rev = bool(time_rev)

    def __mul__(self, other):
        """Composition of base maps."""
        assert self.dim == other.dim
        perm = []
        flip = []
        for i in range(self.dim):
            k = other.perm[self.perm[i]]
            b1 = self.flip[i]
            b2 = other.flip[self.perm[i]]
            perm.append(k)
            flip.append(b1 ^ b2)
        time_rev = self.time_rev ^ other.time_rev
        return self.__class__(perm, flip, time_rev)

    def apply_x(self, x):
        """Apply isometry to a point x of [0,1]^d."""
        return [1-x[ki] if bi else x[ki] for ki, bi in zip(self.perm, self.flip)]

    def apply_cube(self, div, cube):
        """Apply isometry to a sub-cube."""
        return [div-cube[ki]-1 if bi else cube[ki] for ki, bi in zip(self.perm, self.flip)]


class FractalCurve:
    """Class representing fractal peano curve in [0,1]^d.
    Params:
        div         positive integer, number of divisions of each side of the cube (characteristic of a curve)
        dim         dimension (default: 2)
        proto       curve prototype - sequence of "cubes" with coordinates (x_0,..,x_{d-1}), 0<=x_j<div
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
