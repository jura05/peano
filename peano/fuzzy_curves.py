import time
from fractions import Fraction
import itertools

from .fast_fractions import FastFraction
from . import sat_adapters
from . import pieces
from . import curves
from . import fuzzy_poly_curves
from .common import Junction, Spec, Pattern


class FuzzyCurve(fuzzy_poly_curves.FuzzyPolyCurve):
    # base_maps - list as in Curve, may contain None
    # repr_maps - список представителей (coset representatives) base_map-ов, которые сохраняют вход-выход
    # symmetries - симметрии кривой
    def __init__(self, dim, div, proto, base_maps, repr_maps, symmetries):
        self.dim = dim
        self.div = div
        self.proto = tuple(proto)
        self.base_maps = tuple(base_maps)
        self.repr_maps = tuple(repr_maps)
        self.symmetries = tuple(symmetries)
        self.genus = self.div ** self.dim

        # to make parent methods work
        self.patterns = [Pattern(proto=proto, specs=[Spec(base_map=bm) for bm in base_maps])]
        self.pattern_count = 1

    # создать кривую с другим прототипом/base_maps/whatever
    def changed(self, proto=None, base_maps=None, repr_maps=None, symmetries=None):
        return type(self)(
            dim=self.dim,
            div=self.div,
            proto=proto if proto is not None else self.proto,
            base_maps=base_maps if base_maps is not None else self.base_maps,
            repr_maps=repr_maps if repr_maps is not None else self.repr_maps,
            symmetries=symmetries if symmetries is not None else self.symmetries,
        )

    def get_fraction(self, cnum):
        """Get fraction as a curve."""
        return self.apply_base_map(self.base_maps[cnum])

    def reverse(self):
        """Reverse time in a curve."""

        kwargs = {}
        if hasattr(self, 'repr_maps'):
            kwargs['repr_maps'] = reversed(self.repr_maps)

        # симметрии не меняются!

        return self.changed(
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
    def gen_allowed_maps(self, cnum):
        if self.base_maps[cnum] is not None:
            # базовое преобразование уже определено!
            yield self.base_maps[cnum]
            return

        repr_map = self.repr_maps[cnum]
        for symm in self.symmetries:
            yield repr_map * symm

    def gen_possible_curves(self):
        bm_variants = [self.gen_allowed_maps(cnum) for cnum in range(self.genus)]
        for base_maps in itertools.product(*bm_variants):
            yield curves.Curve(
                dim=self.dim,
                div=self.div,
                proto=self.proto,
                base_maps=base_maps,
            )

    def get_piece_position(self, cnum):
        return pieces.CurvePiecePosition(
            dim=self.dim,
            div=self.div,
            cnums=[cnum],
            cubes=[self.proto[cnum]],
        )

    def specify(self, cnum, base_map):
        if base_map not in self.gen_allowed_maps(cnum):
            raise Exception("Can't specify curve")

        new_base_maps = list(self.base_maps)
        new_base_maps[cnum] = base_map
        return self.changed(base_maps=new_base_maps)

    def apply_base_map(self, base_map):
        """Apply base map to a fractal curve, return new curve."""
        if base_map.dim != self.dim:
            raise Exception("Incompatible base map!")

        # можно разложить базовое преобразование в произведение (коммутирующих) 
        # преобразований: обращение времени с тождественной изометрией  +  изометрии куба
        if base_map.time_rev:
            return self.reverse().apply_base_map(base_map.cube_map())

        # применяем изометрию куба

        # прототип подвергается изометрии
        proto = [base_map.apply_cube(self.div, cube) for cube in self.proto]

        # базовые преобразования сопрягаются: действительно, чтобы получить
        # из преобразованной кривой её фракцию, можно сделать так:
        # - сначала вернуться к исходной кривой (inv)
        # - применить преобразование исходной кривой для перехода к фракции (bm)
        # - перейти к преобразованной кривой (base_map)
        inv = base_map.inverse()
        conj_cache = {}
        def conjugate(bm):
            if bm not in conj_cache:
                conj_cache[bm] = base_map * bm * inv
            return conj_cache[bm]

        new_maps = [conjugate(bm) if bm is not None else None for bm in self.base_maps]
        kwargs = {}
        if hasattr(self, 'symmetries'):
            kwargs['symmetries'] = [conjugate(bm) for bm in self.symmetries]
        if hasattr(self, 'repr_maps'):
            kwargs['repr_maps'] = [conjugate(bm) for bm in self.repr_maps]

        return self.changed(proto=proto, base_maps=new_maps, **kwargs)

    def bm_info(self):
        return {cnum: bm for cnum, bm in enumerate(self.base_maps) if bm is not None}

    def is_specialization(self, tmpl):
        return all(self.base_maps[cnum] == bm for cnum, bm in tmpl.bm_info().items())

    #
    # Про отношение
    #

    def init_pairs_tree(self):
        auto_junc = Junction.get_auto_junc(dim=self.dim)
        G = self.genus
        for cnum1 in range(G):
            for cnum2 in range(cnum1 + 2, G):
                pos1 = self.get_piece_position(cnum1)
                pos2 = self.get_piece_position(cnum2)
                yield pieces.CurvePieceBalancedPair(self, auto_junc, pos1, pos2)

        for junc in self.get_junctions_info().keys():
            last_cnum1 = 0 if junc.time_rev else G - 1
            first_cnum2 = G - 1 if junc.base_map.time_rev else 0
            for cnum1 in range(G):
                for cnum2 in range(G):
                    if (cnum1, cnum2) == (last_cnum1, first_cnum2):
                        continue
                    pos1 = self.get_piece_position(cnum1)
                    pos2 = self.get_piece_position(cnum2)
                    yield pieces.CurvePieceBalancedPair(self, junc, pos1, pos2)

    # если отношение кривой больше upper_bound, дальше не идём
    def estimate_ratio(self, ratio_func, rel_tol_inv, upper_bound=None, max_iter=None, sat_pack=100, find_model=False, verbose=False):
        curr_lo = FastFraction(0, 1)

        # в качестве начально верхней оценки возьмём оценку для первой полной кривой
        curve0 = next(self.gen_possible_curves())
        curr_up = curve0.estimate_ratio(ratio_func, rel_tol_inv=rel_tol_inv)['up']

        # инвариант: лучшая кривая в данном классе лежит в [curr_lo, curr_up]

        # делаем аналог bisect для поиска хорошей кривой
        tolerance = FastFraction(rel_tol_inv + 1, rel_tol_inv)
        while curr_up > curr_lo * tolerance:
            new_lo = FastFraction(2, 3) * curr_lo + FastFraction(1, 3) * curr_up
            new_up = FastFraction(1, 3) * curr_lo + FastFraction(2, 3) * curr_up
            print('best in ', float(curr_lo), float(curr_up), 'seek with thresholds:', float(new_lo), float(new_up))
            has_good = self.test_ratio(
                ratio_func,
                lower_bound=new_lo,
                upper_bound=new_up,
                max_iter=max_iter,
                sat_pack=sat_pack,
                find_model=False,
                verbose=verbose,
            )
            if has_good:
                curr_up = new_up
            else:
                curr_lo = new_lo

            if upper_bound is not None and curr_lo > upper_bound:
                return

        print('ratio in:', curr_lo, curr_up)
        data = self.test_ratio(
            ratio_func,
            lower_bound=curr_lo,
            upper_bound=curr_up,
            sat_pack=sat_pack,
            find_model=find_model,
            verbose=verbose,
        )
        return data


    # если возвращает True (или кривую), то есть кривая с отношением <= upper_bound
    # если возвращает False, то нет кривой с отношением < lower_bound
    def test_ratio(self, ratio_func, lower_bound, upper_bound, max_iter=None, sat_pack=100, find_model=False, verbose=False):
        adapter = sat_adapters.CurveSATAdapter(dim=self.dim)
        adapter.init_curve(self)

        pairs_tree = pieces.PairsTree(ratio_func)
        pairs_tree.set_good_threshold(upper_bound)
        pairs_tree.set_bad_threshold(lower_bound)

        for pair in self.init_pairs_tree():
            pairs_tree.add_pair(pair)

        iter_no = 0
        while True:
            iter_no += 1
            if max_iter is not None and iter_no > max_iter:
                print('used all iterations...')
                return None

            pairs_tree.divide()
            while pairs_tree.bad_pairs:
                bad_pair = pairs_tree.bad_pairs.pop()
                adapter.add_forbid_clause(bad_pair.junc, bad_pair.curve)

            if not pairs_tree.data:
                break

            if iter_no % 1000 == 0 and verbose:
                worst_node = pairs_tree.data[0]
                worst_pair = worst_node.pair
                print({
                    'iter': iter_no,
                    'pairs': len(pairs_tree.data),
                    'pair_tree_stats:': pairs_tree.stats,
                    'up': float(worst_node.up),
                    'depth': (worst_pair.pos1.depth, worst_pair.pos2.depth),
                })

            if iter_no % sat_pack == 0:
                if verbose:
                    print('iter:', iter_no, 'adapter stats:', adapter.stats())
                if not adapter.solve():
                    print('no SAT model')
                    return False

        if not adapter.solve():
            print('no SAT model')
            return False

        print('SAT model exists!')
        if not find_model:
            return True

        # это если попросят модель
        model = adapter.get_model()
        return {
            "model": model,
            "curve": adapter.get_curve_from_model(self, model),
            "pairs_tree": pairs_tree,
        }


    #
    # Стыки.
    #

    # словарь {junc: curve_list} кривых, приводящих к данному стыку
    def get_junctions_info(self):
        # строим конфигурации:
        # это набор (i, curve), где в curve заданы bm[0], bm[-1], bm[i], bm[i+1]
        configs = []
        G = self.genus
        for bm_first in self.gen_allowed_maps(0):
            for bm_last in self.gen_allowed_maps(G - 1):
                for cnum in range(G - 1):
                    for bm_i in self.gen_allowed_maps(cnum):
                        if cnum == 0 and bm_i != bm_first:
                            continue
                        for bm_ii in self.gen_allowed_maps(cnum + 1):
                            if cnum + 1 == G - 1 and bm_ii != bm_last:
                                continue
                            curve = self.specify(0, bm_first)\
                                .specify(G - 1, bm_last)\
                                .specify(cnum, bm_i)\
                                .specify(cnum + 1, bm_ii)

                            configs.append((cnum, curve))

        # конфигурации, приводящие к данному стыку
        junc_curves = {}

        for cnum, curve in configs:
            base_junc = curve.get_base_junction(cnum=cnum, pnum=0)
            for junc in curve.gen_junctions_from_base([base_junc]):
                junc_curves.setdefault(junc, []).append(curve)

        return junc_curves
