from fractions import Fraction
from fast_fractions import FastFraction
import functools

from base_map import BaseMap, constraint_base_maps

def ratio_l2(d, dv, dt):
    return FastFraction(sum(x**2 for x in dv)**d, dt**2)

class PartialFractalCurve:
    # base_maps - list as in FractalCurve, may contain None
    # gates - (relative entrance, relative exit) pairs  <=>  brkline
    def __init__(self, dim, div, proto, base_maps, gates):
        self.dim = dim
        self.div = div
        self.proto = proto
        self.base_maps = base_maps
        self.gates = gates

    def genus(self):
        return self.div ** self.dim

    def get_entrance(self):
        cube = self.proto[0]
        rel_entrance = self.gates[0][0]
        return tuple(Fraction(cube[j] + rel_entrance[j], self.div) for j in range(self.dim))

    def get_exit(self):
        cube = self.proto[-1]
        rel_exit = self.gates[-1][1]
        return tuple(Fraction(cube[j] + rel_exit[j], self.div) for j in range(self.dim))

    def get_allowed_maps(self):
        entrance = self.get_entrance()
        exit = self.get_exit()
        return constraint_base_maps(self.dim, {entrance: entrance, exit: exit})

    def apply_base_map(self, base_map):
        """Apply base map to a fractal curve, return new curve."""

        if base_map.time_rev:
            raise Exception("Not implemented")

        cube_map = base_map

        # применяем изометрию куба
        inv = cube_map.inverse()
        # прототип подвергается изометрии
        proto = [cube_map.apply_cube(self.div, cube) for cube in self.proto]

        # базовые преобразования сопрягаются: действительно, чтобы получить
        # из преобразованной кривой её фракцию, можно сделать так:
        # - сначала вернуться к исходной кривой (inv)
        # - применить преобразование исходной кривой для перехода к фракции (bm)
        # - перейти к преобразованной кривой (cube_map)
        new_maps = []
        for bm in self.base_maps:
            if bm is None:
                new_bm = None
            else:
                new_bm = cube_map * bm * inv
            new_maps.append(new_bm)

        gates = []
        for entrance, exit in self.gates:
            new_entrance = cube_map.apply_x(entrance)
            new_exit = cube_map.apply_x(exit)
            gates.append((new_entrance, new_exit))

        return type(self)(
            dim = self.dim,
            div = self.div,
            proto = proto,
            base_maps = new_maps,
            gates = gates,
        )


class CurvePiece:
    # фракция кривой, задаётся следующим образом:
    # cnums - последовательность номеров в прототипе
    # base_maps - последовательность базовых преобразований (для перехода)
    def __init__(self, curve, cnums, cubes, base_maps):
        self.curve = curve
        self.cnums = cnums
        self.cubes = cubes
        self.base_maps = base_maps

    # делим кривую дальше всеми способами
    def divide(self):
        # определим ориентацию последней фракции
        prev_map = BaseMap(dim=self.curve.dim)
        for bm in self.base_maps:
            prev_map = bm * prev_map

        # определим допустимые basemap-ы
        prev_curve = self.curve.apply_base_map(prev_map)
        entrance = prev_curve.gates[self.cnums[-1]][0]
        exit = prev_curve.gates[self.cnums[-1]][1]
        allowed_maps = constraint_base_maps(
            self.curve.dim,
            {prev_curve.get_entrance(): entrance, prev_curve.get_exit(): exit}
        )

        # делим
        new_pieces = []
        for bm in allowed_maps:
            last_map = bm * prev_map
            last_curve = self.curve.apply_base_map(last_map)
            for cnum, cube in enumerate(last_curve.proto):
                new_piece = type(self)(
                    curve = self.curve,
                    cnums = self.cnums + [cnum],
                    cubes = self.cubes + [cube],
                    base_maps = self.base_maps + [bm],
                )
                new_pieces.append(new_piece)
        return new_pieces

    def coords(self):
        if not hasattr(self, "_coords"):
            self._coords = self._get_coords()
        return self._coords

    def _get_coords(self):
        curve = self.curve
        G = curve.genus()
        N = curve.div
        dim = curve.dim

        #
        # Все x умножаем на N*l, все t на G^l
        #

        # t = c0/G + c1/G**2 + ... = (c_{l-1} + c_{l-2}*G + ..) / G^l
        t = 0
        Gpower = 1
        for cnum in reversed(self.cnums):
            t += cnum * Gpower
            Gpower *= G

        x = [0 for j in range(dim)]
        Npower = 1
        for cube in reversed(self.cubes):
            for j in range(dim):
                x[j] += cube[j] * Npower
            Npower *= N

        return {'x': x, 't': t, 'depth': len(self.cnums)}


class CurvePieceBalancedPair:
    # сбалансированная пара фракций кривых
    # считаем, что t_2 > t_1
    def __init__(self, piece1, piece2):
        self.piece1 = piece1
        self.piece2 = piece2
        self.curve = self.piece1.curve  # == self.piece2.curve

    def divide(self):
        # решаем, кого делить
        if len(self.piece1.cnums) > len(self.piece2.cnums):
            # делить второй кусок
            new_pairs = []
            for subpiece in self.piece2.divide():
                new_pairs.append(type(self)(self.piece1, subpiece))
        else:
            new_pairs = []
            for subpiece in self.piece1.divide():
                new_pairs.append(type(self)(subpiece, self.piece2))

        return  new_pairs

    def upper_bound(self, verbose=False):
        dim = self.curve.dim
        N = self.curve.div
        G = self.curve.genus()

        coords1 = self.piece1.coords().copy()
        coords2 = self.piece2.coords().copy()

        if coords1['depth'] == coords2['depth']:
            coords1['mx'] = coords1['mt'] = coords2['mx'] = coords2['mt'] = 1
            # глубина одинакова, делать вообще ничего не нужно
        elif coords1['depth'] == (coords2['depth'] + 1):
            coords1['mx'] = coords1['mt'] = 1
            coords2['x'] = tuple(coords2['x'][j] * N for j in range(dim))
            coords2['t'] = coords2['t'] * G
            coords2['mx'] = N
            coords2['mt'] = G
        else:
            raise Exception("Bad coordinates!")

        if verbose:
            print(coords1, coords2)
        dim = self.curve.dim
        dx = []
        for j in range(dim):
            x1j = coords1['x'][j]
            x2j = coords2['x'][j]
            dxj = max(abs(x1j - x2j + coords1['mx']), abs(x1j - x2j - coords2['mx']))
            dx.append(dxj)

        dt = coords2['t'] - (coords1['t'] + coords1['mt'])

#        print('=== PIECE 1 ===')
#        print(self.piece1.cnums, self.piece1.cubes, coords1)
#        print('=== PIECE 2 ===')
#        print(self.piece2.cnums, self.piece2.cubes, coords2)
#        print(dx,dt)
#        print('')

        return ratio_l2(dim, dx, dt)

    # тут можно попробовать использовать gates
    def lower_bound(self):
        dim = self.curve.dim
        N = self.curve.div
        G = self.curve.genus()

        coords1 = self.piece1.coords().copy()
        coords2 = self.piece2.coords().copy()
        dx = []

        if coords1['depth'] == coords2['depth']:
            coords1['mx'] = coords1['mt'] = coords2['mx'] = coords2['mt'] = 1
            # глубина одинакова, делать вообще ничего не нужно
        elif coords1['depth'] == (coords2['depth'] + 1):
            coords1['mx'] = coords1['mt'] = 1

            coords2['x'] = tuple(coords2['x'][j] * N for j in range(dim))
            coords2['t'] = coords2['t'] * G
            coords2['mx'] = N
            coords2['mt'] = G
        else:
            raise Exception("Bad coordinates!")


        for j in range(dim):
            x1j = coords1['x'][j]
            x2j = coords2['x'][j]
            a1j = x1j - x2j + coords1['mx']
            a2j = x1j - x2j - coords2['mx']
            if (a1j >= 0 and a2j >= 0) or (a1j <= 0 and a2j <= 0):
                # одного знака, ноль не задели
                dxj = min(abs(a1j), abs(a2j))
            else:
                dxj = 0
            dx.append(dxj)

        dt = coords2['t'] + coords2['mt'] - coords1['t']

        return ratio_l2(dim, dx, dt)

def forget(curve):
    curve_entrance = curve.get_entrance()
    curve_exit = curve.get_exit()
    gates = []
    for bm in curve.base_maps:
        new_entrance = bm.apply_x(curve_entrance)
        new_exit = bm.apply_x(curve_exit)
        gates.append((new_entrance, new_exit))

    return PartialFractalCurve(
        dim=curve.dim,
        div=curve.div,
        proto=curve.proto,
        base_maps=[None for j in range(curve.genus())],  # забыли BaseMap-ы!
        gates=gates,
    )
