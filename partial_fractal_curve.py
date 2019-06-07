from fractions import Fraction
from fast_fractions import FastFraction

from base_map import BaseMap, constraint_base_maps


# time_rev not supported!
class PartialFractalCurve:
    # base_maps - list as in FractalCurve, may contain None
    # gates - (relative entrance, relative exit) pairs  <=>  brkline
    def __init__(self, dim, div, proto, base_maps, gates):
        self.dim = dim
        self.div = div
        self.proto = tuple(proto)

        if any(bm is not None and bm.time_rev for bm in base_maps):
            raise Exception("time_rev not implemented!")

        self.base_maps = tuple(base_maps)
        self.gates = tuple(gates)

    def genus(self):
        return self.div ** self.dim

    def get_entrance(self):
        cube = self.proto[0]
        rel_entrance = self.gates[0][0]
        return tuple(Fraction(cube[j] + rel_entrance[j], self.div) for j in range(self.dim))

    def get_fraction(self, cnum):
        return self.apply_base_map(self.base_maps[cnum])

    def get_exit(self):
        cube = self.proto[-1]
        rel_exit = self.gates[-1][1]
        return tuple(Fraction(cube[j] + rel_exit[j], self.div) for j in range(self.dim))

    # кандидаты в self.base_maps[cnum]
    def get_allowed_maps(self, cnum):
        if self.base_maps[cnum] is not None:
            # базовое преобразование уже определено!
            return [self.base_maps[cnum]]

        curve_entrance = self.get_entrance()
        curve_exit = self.get_exit()
        (piece_entrance, piece_exit) = self.gates[cnum]
        return constraint_base_maps(
            self.dim,
            {curve_entrance: piece_entrance, curve_exit: piece_exit},
        )

    def get_piece_position(self, cnum):
        return CurvePiecePosition(
            dim=self.dim,
            div=self.div,
            cnums=[cnum],
            cubes=[self.proto[cnum]],
        )

    def specify(self, cnum, base_map):
        if self.base_maps[cnum] is not None:
            return self

        new_base_maps = list(self.base_maps)
        new_base_maps[cnum] = base_map

        return type(self)(
            dim=self.dim,
            div=self.div,
            proto=self.proto,
            base_maps=new_base_maps,
            gates=self.gates,
        )

    def apply_base_map(self, base_map):
        """Apply base map to a fractal curve, return new curve."""

        if base_map.time_rev:
            raise Exception("Not implemented")

        cube_map = base_map

        # применяем изометрию куба

        # прототип подвергается изометрии
        proto = [cube_map.apply_cube(self.div, cube) for cube in self.proto]

        # базовые преобразования сопрягаются (см. fractal_curve.py)
        inv = cube_map.inverse()
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
            dim=self.dim,
            div=self.div,
            proto=proto,
            base_maps=new_maps,
            gates=gates,
        )

# вся инфа о положении фракции в кривой
class CurvePiecePosition:
    def __init__(self, dim, div, cnums, cubes):
        self.dim = dim
        self.div = div
        self.cnums = cnums
        self.cubes = cubes

    def depth(self):
        return len(self.cnums)

    def specify(self, cnum, cube):
        return type(self)(
            dim=self.dim,
            div=self.div,
            cnums=self.cnums + [cnum],
            cubes=self.cubes + [cube],
        )

    # естественные целочисленные координаты - время и нулевой угол куба
    # int time: [t, t+1], abs time: [t/G^l, (t+1)/G^l]
    # int cube: cj <= xj <= cj+1, abs cube: cj/N^l <= xj <= (cj+1)/N^l
    # l = depth
    def int_coords(self):
        if not hasattr(self, "_int_coords"):
            self._int_coords = self.get_int_coords()
        return self._int_coords

    def get_int_coords(self):
        dim = self.dim
        N = self.div
        G = N**dim

        # Все x умножаем на N*l, все t на G^l
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

        return {'x': x, 't': t, 'l': len(self.cnums)}

class CurvePiece:
    # фракция кривой = кривая + позиция фракции
    def __init__(self, curve, pos):
        self.curve = curve
        self.pos = pos

    def prev_map(self):
        if not hasattr(self, '_prev_map'):
            self._prev_map = self.get_prev_map()
        return self._prev_map

    # TODO: comment
    # в последнем кубе ориентация не задана!
    def get_prev_map(self):
        prev_map = BaseMap(dim=self.curve.dim)
        for cnum in self.pos.cnums[:-1]:
            prev_map = prev_map * self.curve.base_maps[cnum]  # именно в таком порядке!
        return prev_map

    def get_last_gates(self):
        prev_map = self.prev_map()
        (entrance, exit) = self.curve.gates[self.pos.cnums[-1]]
        new_entrance = prev_map.apply_x(entrance)
        new_exit = prev_map.apply_x(exit)
        return (new_entrance, new_exit)

    # делим кривую дальше всеми способами
    def divide(self):
        # определим ориентацию предпоследнего кусочка

        prev_map = self.prev_map()
        prev_curve = self.curve.apply_base_map(prev_map)

        active_cnum = self.pos.cnums[-1]  # кубик, где всё происходит
        allowed_maps = prev_curve.get_allowed_maps(active_cnum)

        # делим
        new_pieces = []
        for bm in allowed_maps:
            last_curve = prev_curve.apply_base_map(bm)

            # сопряжение, как в apply:
            # для prev_map.base_maps[active_cnum] = bm =>  orig_curve.base_maps[active_cnum] = ...
            new_map = prev_map.inverse() * bm * prev_map
            specified_curve = self.curve.specify(active_cnum, new_map)

            for cnum, cube in enumerate(last_curve.proto):
                new_pos = self.pos.specify(cnum, cube)
                new_piece = CurvePiece(specified_curve, new_pos)
                new_pieces.append(new_piece)

        return new_pieces


class CurvePieceBalancedPair:
    # сбалансированная пара фракций кривых
    # считаем, что t_2 > t_1
    def __init__(self, curve, pos1, pos2):
        self.curve = curve
        self.pos1 = pos1
        self.pos2 = pos2
        self.piece1 = CurvePiece(self.curve, self.pos1)
        self.piece2 = CurvePiece(self.curve, self.pos2)

    def divide(self):
        # при делении кусочка у него уточняется кривая, поэтому берём кривую из него
        # решаем, кого делить
        if self.pos1.depth() > self.pos2.depth():
            new_pairs = []
            for subpiece in self.piece2.divide():
                new_pairs.append(type(self)(subpiece.curve, self.pos1, subpiece.pos))
        else:
            new_pairs = []
            for subpiece in self.piece1.divide():
                new_pairs.append(type(self)(subpiece.curve, subpiece.pos, self.pos2))

        return new_pairs

    def upper_bound(self, ratio_func, verbose=False):
        dim = self.curve.dim
        N = self.curve.div
        G = self.curve.genus()

        coords1 = self.pos1.int_coords().copy()
        coords2 = self.pos2.int_coords().copy()

        if coords1['l'] == coords2['l']:
            coords1['mx'] = coords1['mt'] = coords2['mx'] = coords2['mt'] = 1
            # глубина одинакова, делать вообще ничего не нужно
        elif coords1['l'] == (coords2['l'] + 1):
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
        return FastFraction(*ratio_func(dim, dx, dt))

    # с учётом ворот!
    def lower_bound(self, ratio_func):
        dim = self.curve.dim
        N = self.curve.div
        G = self.curve.genus()

        coords1 = self.pos1.int_coords().copy()
        coords2 = self.pos2.int_coords().copy()
        dx = []

        # УЧЁТ ВОРОТ
        piece1_exit = self.piece1.get_last_gates()[1]
        piece2_entrance = self.piece2.get_last_gates()[0]
        coords1['x'] = [cj + gj for cj, gj in zip(coords1['x'], piece1_exit)]
        coords2['x'] = [cj + gj for cj, gj in zip(coords2['x'], piece2_entrance)]

        if coords1['l'] == coords2['l']:
            coords1['mx'] = coords1['mt'] = coords2['mx'] = coords2['mt'] = 1
            # глубина одинакова, делать вообще ничего не нужно
        elif coords1['l'] == (coords2['l'] + 1):
            coords1['mx'] = coords1['mt'] = 1

            coords2['x'] = tuple(coords2['x'][j] * N for j in range(dim))
            coords2['t'] = coords2['t'] * G
            coords2['mx'] = N
            coords2['mt'] = G
        else:
            raise Exception("Bad coordinates!")

        for j in range(dim):
            # всё просто, т.к. это реальные координаты точек!
            x1j = coords1['x'][j]
            x2j = coords2['x'][j]
            dxj = abs(x1j - x2j)
            dx.append(dxj)

        # т.к. ворота, то разница во времени минимальна!
        dt = coords2['t'] - (coords1['t'] + coords1['mt'])

        return FastFraction(*ratio_func(dim, dx, dt))


def forget(curve):
    curve_entrance = curve.get_entrance()
    curve_exit = curve.get_exit()

    # не создаём дробные ворота!
    def fix(x):
        return int(x) if x == int(x) else x

    curve_entrance = tuple(fix(ce) for ce in curve_entrance)
    curve_exit = tuple(fix(ce) for ce in curve_exit)

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
