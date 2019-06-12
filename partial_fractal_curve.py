from fractions import Fraction
import itertools
from functools import lru_cache

from fast_fractions import FastFraction
from base_map import BaseMap, constraint_base_maps, list_base_maps


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


class PartialFractalCurve:
    # base_maps - list as in FractalCurve, may contain None
    # gates - (relative entrance, relative exit) pairs  <=>  brkline
    def __init__(self, dim, div, proto, base_maps, gates):
        self.dim = dim
        self.div = div
        self.proto = tuple(proto)
        self.base_maps = tuple(base_maps)
        self.gates = tuple(gates)
        self.genus = self.div ** self.dim

    def get_entrance(self):
        cube = self.proto[0]
        rel_entrance = self.gates[0][0]
        return tuple(Fraction(cube[j] + rel_entrance[j], self.div) for j in range(self.dim))

    def get_fraction(self, cnum):
        """Get fraction as a curve."""
        return self.apply_base_map(self.base_maps[cnum])

    def get_exit(self):
        cube = self.proto[-1]
        rel_exit = self.gates[-1][1]
        return tuple(Fraction(cube[j] + rel_exit[j], self.div) for j in range(self.dim))

    def reverse(self):
        """Reverse time in a curve."""

        kwargs = {}
        if hasattr(self, 'gates'):
            # ворота идут в обратном порядке, вход и выход меняются местами
            gates = reversed(reversed(g) for g in gates)
            kwargs['gates'] = gates

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

            **kwargs,
        )

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
            if self.base_maps[cnum] == base_map:
                return self  # nothing to do
            else:
                raise Exception("Can't specify curve")

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

        # можно разложить базовое преобразование в произведение (коммутирующих) 
        # преобразований: обращение времени с тождественной изометрией  +  изометрии куба
        if base_map.time_rev:
            curve = self.reverse()
            cube_map = base_map.cube_map()
        else:
            curve = self
            cube_map = base_map

        # применяем изометрию куба

        # прототип подвергается изометрии
        proto = [cube_map.apply_cube(self.div, cube) for cube in self.proto]

        # базовые преобразования сопрягаются: действительно, чтобы получить
        # из преобразованной кривой её фракцию, можно сделать так:
        # - сначала вернуться к исходной кривой (inv)
        # - применить преобразование исходной кривой для перехода к фракции (bm)
        # - перейти к преобразованной кривой (cube_map)
        inv = cube_map.inverse()
        new_maps = []
        for bm in self.base_maps:
            if bm is None:
                new_bm = None
            else:
                new_bm = cube_map * bm * inv
            new_maps.append(new_bm)

        kwargs = {}
        if hasattr(self, 'gates'):
            gates = []
            for entrance, exit in self.gates:
                new_entrance = cube_map.apply_x(entrance)
                new_exit = cube_map.apply_x(exit)
                gates.append((new_entrance, new_exit))
            kwargs['gates'] = gates

        return type(self)(
            dim=self.dim,
            div=self.div,
            proto=proto,
            base_maps=new_maps,
            **kwargs,
        )

    def bm_info(self):
        return {cnum: bm for cnum, bm in enumerate(self.base_maps) if bm is not None}

    def is_specialization(self, tmpl):
        return all(self.base_maps[cnum] == bm for cnum, bm in tmpl.bm_info().items())

    #
    # Стыки.
    #

    # словарь {junc: curve_list} кривых, приводящих к данному стыку
    def get_junctions_info(self):
        # строим конфигурации:
        # это набор (i, curve), где в curve заданы bm[0], bm[-1], bm[i], bm[i+1]
        configs = []
        G = self.genus
        for bm_first in self.get_allowed_maps(0):
            for bm_last in self.get_allowed_maps(G - 1):
                for i in range(G - 1):
                    for bm_i in self.get_allowed_maps(i):
                        if i == 0 and bm_i != bm_first:
                            continue
                        for bm_ii in self.get_allowed_maps(i + 1):
                            if i == G - 2 and bm_ii != bm_last:
                                continue
                            curve = self.specify(0, bm_first)\
                                .specify(G-1, bm_last)\
                                .specify(i, bm_i)\
                                .specify(i+1, bm_ii)

                            configs.append((i, curve))

        # конфигурации, приводящие к данному стыку
        junc_curves = {}

        for i, curve in configs:
            seen_junc = set()
            cube = curve.proto[i]
            next_cube = curve.proto[i+1]
            delta = tuple(nc-c for nc, c in zip(next_cube, cube))
            base_junc = curve._get_std_junction(delta, curve.base_maps[i], curve.base_maps[i+1])
            seen_junc.add(base_junc)

            to_derive = [base_junc]
            while to_derive:
                junc = to_derive.pop()
                dj = curve.get_derived_junction(junc)
                if dj not in seen_junc:
                    seen_junc.add(dj)
                    to_derive.append(dj)

            for junc in seen_junc:
                if junc not in junc_curves:
                    junc_curves[junc] = []
                junc_curves[junc].append(curve)

        return junc_curves

    def get_derived_junction(self, junction):
        if self.base_maps[0] is None or self.base_maps[-1] is None:
            raise Exception("Can't get derivative")
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


# вся инфа о положении фракции в кривой
class CurvePiecePosition:
    def __init__(self, dim, div, cnums, cubes):
        self.dim = dim
        self.div = div
        self.cnums = cnums
        self.cubes = cubes
        self.depth = len(self.cnums)

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
    # int curve: [0, N**l]^d, abs curve: [0,1]^d
    # l = depth
    # возвращает [l,x,t]
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
        new_pairs = []
        if self.pos1.depth > self.pos2.depth:
            for subpiece in self.piece2.divide():
                new_pairs.append(type(self)(subpiece.curve, self.junc, self.pos1, subpiece.pos))
        else:
            for subpiece in self.piece1.divide():
                new_pairs.append(type(self)(subpiece.curve, self.junc, subpiece.pos, self.pos2))

        return new_pairs

    # пространственные и временные расстояния
    def int_dist(self):
        if not hasattr(self, '_int_dist'):
            self._int_dist = self.get_int_dist()
        return self._int_dist

    def get_int_dist(self):
        dim = self.curve.dim
        N = self.curve.div
        G = self.curve.genus

        l1, x1, t1 = self.pos1.get_int_coords()
        l2, x2, t2 = self.pos2.get_int_coords()

        if self.junc is None:
            junc_dt = 0
            junc_dx = (0,) * dim
        else:
            junc_dt = 1
            junc_dx, junc_bm = self.junc
            # junc: сначала поворот
            x2 = junc_bm.apply_cube(N**l2, x2)

        if l1 == l2:
            mx2 = mt2 = 1
        elif l1 == l2 + 1:
            x2 = [x2j * N for x2j in x2]
            t2 *= G
            mx2 = N
            mt2 = G
        else:
            raise Exception("Bad coordinates!")

        mx = N**l1
        mt = mx**dim

        # мы привели целые координаты к следующим:
        # cube1: x1j <= xj <= x1j + 1         -- кубик внутри кривой [0, mx]^d
        # cube2: x2j <= xj <= x2j + mx2  -- кубик внутри кривой [0, mx2 * mx]^d

        max_dx = [None] * dim
        for j in range(dim):
            x1j = x1[j]

            # junc: потом сдвиг
            x2j = x2[j] + junc_dx[j] * mx

            dxj = max(abs(x1j - x2j + 1), abs(x1j - x2j - mx2))
            max_dx[j] = dxj

        min_dt = junc_dt * mt + t2 - (t1 + 1)
        max_dt = junc_dt * mt + (t2 + mt2) - t1
        return {'max_dx': max_dx, 'min_dt': min_dt, 'max_dt': max_dt}

    def upper_bound(self, ratio_func):
        dist = self.int_dist()
        return FastFraction(*ratio_func(self.curve.dim, dist['max_dx'], dist['min_dt']))

    def lower_bound(self, ratio_func):
        dist = self.int_dist()
        return FastFraction(*ratio_func(self.curve.dim, dist['max_dx'], dist['max_dt']))


