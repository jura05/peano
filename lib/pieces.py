from functools import lru_cache
from collections import Counter, namedtuple
from heapq import heappop, heappush

from fast_fractions import FastFraction
from base_map import BaseMap, gen_constraint_cube_maps


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


# вся инфа о положении фракции в кривой
# cnums - последовательность номеров кубов
# cubes - последовательность кубов
# каждый куб - для той кривой, которая в соотв. фракции
class CurvePiecePosition:
    def __init__(self, dim, div, cnums, cubes):
        self.dim = dim
        self.div = div
        self.cnums = cnums
        self.cubes = cubes
        self.depth = len(self.cnums)
        self.sub_div = div**self.depth
        self.sub_genus = self.sub_div**dim

    def specify(self, cnum, cube):
        return type(self)(
            dim=self.dim,
            div=self.div,
            cnums=self.cnums + [cnum],
            cubes=self.cubes + [cube],
        )

    # естественные целочисленные координаты - время и нулевой угол куба
    # int time: t - время, умноженное на G^l, кусок - [t, t+1], abs time: [t/G^l, (t+1)/G^l]
    # int cube: cj <= xj <= cj+1, abs cube: cj/N^l <= xj <= (cj+1)/N^l
    # int curve: [0, N**l]^d, abs curve: [0,1]^d
    # l = depth
    # возвращает l,x,t
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

    # делим кривую дальше всеми способами
    def divide(self):
        curve = self.curve
        dim, G = curve.dim, curve.genus

        # определим ориентацию предпоследнего кусочка
        # в последнем кубе ориентация не задана!
        prev_map = BaseMap.id_map(dim)
        for cnum in self.pos.cnums[:-1]:
            cnum = prev_map.apply_cnum(G, cnum)
            prev_map = prev_map * curve.base_maps[cnum]  # именно в таком порядке!

        prev_curve = curve.apply_base_map(prev_map)

        active_cnum = self.pos.cnums[-1]  # кубик, где всё происходит
        spec_cnum = prev_map.apply_cnum(G, active_cnum)

        # делим
        for bm in prev_curve.gen_allowed_maps(active_cnum):

            # сопряжение, как в apply:
            # для prev_curve.base_maps[active_cnum] = bm =>  orig_curve.base_maps[spec_cnum] = ...
            new_map = prev_map.inverse() * bm * prev_map
            specified_curve = curve.specify(spec_cnum, new_map)

            last_curve = prev_curve.apply_base_map(bm)
            for cnum, cube in enumerate(last_curve.proto):
                new_pos = self.pos.specify(cnum, cube)
                new_piece = CurvePiece(specified_curve, new_pos)
                yield new_piece


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
        if self.pos1.depth > self.pos2.depth:
            for subpiece in self.piece2.divide():
                yield type(self)(subpiece.curve, self.junc, self.pos1, subpiece.pos)
        else:
            for subpiece in self.piece1.divide():
                yield type(self)(subpiece.curve, self.junc, subpiece.pos, self.pos2)

    # пространственные и временные расстояния
    def int_dist(self):
        if not hasattr(self, '_int_dist'):
            self._int_dist = self.get_int_dist()
        return self._int_dist

    def get_int_dist(self):
        dim = self.curve.dim
        N = self.curve.div
        G = self.curve.genus

        pos1, pos2 = self.pos1, self.pos2

        l1, x1, t1 = pos1.get_int_coords()
        l2, x2, t2 = pos2.get_int_coords()

        junc = self.junc
        if junc is None:
            junc_dt = 0
            junc_dx = (0,) * dim
        else:
            junc_dt = 1
            junc_dx = self.junc.delta_x

            # junc: time_rev
            if junc.time_rev:
                t1 = pos1.sub_genus - 1 - t1

            # junc: сначала поворот
            x2 = junc.base_map.apply_cube(pos2.sub_div, x2)
            t2 = junc.base_map.apply_cnum(pos2.sub_genus, t2)

        # приведение к единому масштабу
        if l1 == l2:
            mx2 = mt2 = 1
        elif l1 == l2 + 1:
            x2 = [x2j * N for x2j in x2]
            t2 *= G
            mx2 = N
            mt2 = G
        else:
            raise Exception("Bad coordinates!")

        mx = pos1.sub_div
        mt = pos1.sub_genus

        # мы привели целые координаты к следующим:
        # cube1: x1j <= xj <= x1j + 1         -- кубик внутри кривой [0, mx]^d
        # cube2: x2j <= xj <= x2j + mx2  -- кубик внутри кривой [0, mx2 * mx]^d, + сдвиг на junc_dx * mx
        #
        # time1: t1 <= t <= t1 + 1
        # time2: t2 <= t <= t2 + mt2, после сдвига: t2 + junc_dt * mt <= t <= t2 + mt2 + junc_dt * mt

        max_dx = [None] * dim
        for j in range(dim):
            x1j = x1[j]

            # junc: потом сдвиг
            x2j = x2[j] + junc_dx[j] * mx

            dxj = max(abs(x1j - x2j + 1), abs(x1j - x2j - mx2))
            max_dx[j] = dxj

        max_dt = t2 + junc_dt * mt + mt2 - t1  # max(t_2 - t_1)
        min_dt = t2 + junc_dt * mt - (t1 + 1)  # min(t_2 - t_1)

        return {'max_dx': max_dx, 'min_dt': min_dt, 'max_dt': max_dt}

    def upper_bound(self, ratio_func):
        dist = self.int_dist()
        return FastFraction(*ratio_func(self.curve.dim, dist['max_dx'], dist['min_dt']))

    def lower_bound(self, ratio_func):
        dist = self.int_dist()
        return FastFraction(*ratio_func(self.curve.dim, dist['max_dx'], dist['max_dt']))

# пары фракций с приоритетом
class PairsTree:

    RichPair = namedtuple('RichPair', [
        'priority',  # первое поле, по нему идёт сравнение; для быстрой сортировки по upper_bound
        'pair',  # CurvePieceBalancedPair
        'up',  # upper_bound
        'lo',  # lower_bound
    ])

    def __init__(self, ratio_func):
        self.data = []
        self.stats = Counter()
        self.bad_pairs = []
        self.ratio_func = ratio_func
        self.good_threshold = None  # если отношение <= порога, пара хорошая
        self.bad_threshold = None  # если отношение > порога, пара плохая
        self._inc = 0
        self.LOG = []

    def set_good_threshold(self, thr):
        self.good_threshold = thr

    def set_bad_threshold(self, thr):
        self.bad_threshold = thr

    def add_pair(self, pair):
        up = pair.upper_bound(self.ratio_func)
        gthr = self.good_threshold
        if gthr is not None and float(up) < gthr:  # TOOD нет ли проблемы с округлением
            self.stats['good'] += 1
            return

        lo = pair.lower_bound(self.ratio_func)
        bthr = self.bad_threshold
        if bthr is not None and float(lo) > bthr:
            self.bad_pairs.append(pair)
            self.stats['bad'] += 1
            return

        self._inc += 1  # чтобы сравнение не проваливалось к парам
        node = PairsTree.RichPair(priority=(-float(up), self._inc), pair=pair, up=up, lo=lo)
        heappush(self.data, node)

    def divide(self):
        if not self.data:
            return

        worst_node = heappop(self.data)
        for pair in worst_node.pair.divide():
            self.add_pair(pair)
