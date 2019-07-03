from collections import namedtuple
import itertools
from fractions import Fraction

from .common import Junction
from .base_maps import BaseMap, Spec, gen_constraint_cube_maps
from .fast_fractions import FastFraction
from . import pieces
from . import curves
from . import sat_adapters


# patterns - список пар (proto, specs)
# pnum - выделенный шаблон, задаёт реальную кривую [0,1]->[0,1]^d
class FuzzyCurve:
    Pattern = namedtuple('Pattern', ['proto', 'specs'])

    def __init__(self, dim, div, patterns, pnum=0):
        self.dim = dim
        self.div = div

        self.patterns = []
        for proto, specs in patterns:
            proto = (tuple(cube) for cube in proto)
            specs = (sp if isinstance(sp, Spec) else Spec(sp) if sp is not None else None for sp in specs)
            pattern = FuzzyCurve.Pattern(proto=tuple(proto), specs=tuple(specs))
            self.patterns.append(pattern)

        self.pattern_count = len(self.patterns)
        self.pnum = pnum

        self.proto = self.patterns[pnum].proto
        self.specs = self.patterns[pnum].specs

        self.genus = div**dim

    def get_fraction(self, cnum):
        return self.specs[cnum] * self

    # создать кривую с другими данным
    def changed(self, patterns=None, pnum=None, **kwargs):
        return type(self)(
            dim=self.dim,
            div=self.div,
            patterns=patterns if patterns is not None else self.patterns,
            pnum=pnum if pnum is not None else self.pnum,
            **kwargs,
        )

    def gen_possible_curves(self):
        # предполагается, что все gen_allowed_specs совместимы!
        sp_variants = []
        G = self.genus
        # делаем большой список генераторов вариантов, по pnum * genus + cnum; TODO - сделать по-красивше
        for pnum in range(self.pattern_count):
            sp_variants += [self.gen_allowed_specs(pnum, cnum) for cnum in range(G)]

        for all_specs in itertools.product(*sp_variants):
            patterns = []
            for pnum, pattern in enumerate(self.patterns):
                specs = all_specs[pnum * G:(pnum + 1) * G]
                patterns.append((pattern.proto, specs))

            yield curves.Curve(dim=self.dim, div=self.div, patterns=patterns, pnum=self.pnum)

    def reverse(self):
        """Reverse time in a curve."""
        new_patterns = []
        for pattern in self.patterns:
            # прототип проходится в обратном порядке
            new_proto = reversed(pattern.proto)

            # базовые преобразования проходятся в обратном порядке
            # сами по себе они не меняются:
            #   - если обращения времени не было, то его и не будет
            #   - изометрия куба не меняется, т.к. время не играет роли
            new_specs = reversed(pattern.specs)
            new_patterns.append((new_proto, new_specs))

        return self.changed(patterns=new_patterns)

    def apply_cube_map(self, base_map):
        new_patterns = []
        for pattern in self.patterns:
            # прототип подвергается изометрии
            new_proto = [base_map.apply_cube(self.div, cube) for cube in pattern.proto]

            # спеки преобразования сопрягаются: действительно, чтобы получить
            # из преобразованной кривой её фракцию, можно сделать так:
            # - сначала вернуться к исходной кривой (inv)
            # - применить преобразование исходной кривой для перехода к фракции (bm)
            # - перейти к преобразованной кривой (base_map)
            inv = base_map.inverse()
            conj_cache = {}
            def conjugate(sp):
                if sp not in conj_cache:
                    conj_cache[sp] = base_map * sp * inv
                return conj_cache[sp]

            new_specs = [conjugate(spec) if spec is not None else None for spec in pattern.specs]
            new_patterns.append((new_proto, new_specs))

        return self.changed(patterns=new_patterns)

    def __rmul__(self, other):
        """Apply base map/spec to a fractal curve, return new curve."""
        if isinstance(other, Spec):
            new_curve = other.base_map * self
            return new_curve.changed(pnum=other.pnum)
        elif not isinstance(other, BaseMap):
            return NotImplemented

        base_map = other

        if base_map.dim != self.dim:
            raise Exception("Incompatible base map!")

        # можно разложить базовое преобразование в произведение (коммутирующих) 
        # преобразований: обращение времени с тождественной изометрией  +  изометрии куба
        if base_map.time_rev:
            return base_map.cube_map() * self.reverse()

        return self.apply_cube_map(base_map.cube_map())

    def gen_allowed_specs(self, pnum, cnum):
        raise NotImplementedError("Define in child class")

    def sp_info(self):
        curve_info = []
        for pnum, pattern in enumerate(self.patterns):
            for cnum, spec in enumerate(pattern.specs):
                if spec is not None:
                    curve_info.append((pnum, cnum, spec))
        return curve_info

    def is_specialization(self, tmpl):
        for pnum, cnum, sp in tmpl.sp_info():
            if self.patterns[pnum].specs[cnum] != sp:
                return False
        return True

    def get_piece_position(self, pnum, cnum):
        return pieces.CurvePiecePosition(
            dim=self.dim,
            div=self.div,
            cnums=[cnum],
            cubes=[self.patterns[pnum].proto[cnum]],
        )

    def specify(self, pnum, cnum, spec):
        if spec not in self.gen_allowed_specs(pnum, cnum):
            raise Exception("Can't specify curve")

        pattern = self.patterns[pnum]
        new_specs = list(pattern.specs)
        new_specs[cnum] = spec
        new_pattern = (pattern.proto, new_specs)

        new_patterns = list(self.patterns)
        new_patterns[pnum] = new_pattern

        return self.changed(patterns=new_patterns)


    #
    # Стыки - общие методы
    #
    # Все стыки = автостыки + базовые стыки + производные базовых стыков
    #

    # автостыки - они учитываются отдельно!
    def gen_auto_junctions(self):
        for pnum in range(self.pattern_count):
            yield Junction.get_auto_junc(dim=self.dim, pnum=pnum)

    # возвращает стыки вместе с производными
    def gen_junctions_from_base(self, base_juncs):
        for junc in base_juncs:
            yield junc
        seen = set(base_juncs)
        to_derive = list(base_juncs)
        while to_derive:
            junc = to_derive.pop()
            dj = self.get_derived_junction(junc)
            if dj not in seen:
                yield dj
                seen.add(dj)
                to_derive.append(dj)

    def get_base_junction(self, pnum, cnum):
        pattern = self.patterns[pnum]
        delta_x = [c2j - c1j for c1j, c2j in zip(pattern.proto[cnum], pattern.proto[cnum + 1])]
        return Junction.get_junc(
            pattern.specs[cnum],
            pattern.specs[cnum + 1],
            delta_x,
        )

    # производные стыки - не для автостыков!
    def get_derived_junction(self, junc):
        if junc.delta_t != 1:
            raise Exception("Derivative is defined for dt=1 junctions!")

        spec1 = junc.spec1
        spec2 = junc.spec2

        p1 = self.patterns[spec1.pnum]
        p2 = self.patterns[spec2.pnum]

        cnum1 = 0 if spec1.base_map.time_rev else -1
        cnum2 = -1 if spec2.base_map.time_rev else 0

        if p1.specs[cnum1] is None or p2.specs[cnum2] is None:
            raise Exception("Can't get derivative: spec not defined")

        cube1 = spec1.base_map.apply_cube(self.div, p1.proto[cnum1])  # сейчас не нужно - нормализуем cube_map1 -> id
        cube2 = spec2.base_map.apply_cube(self.div, p2.proto[cnum2])
        der_delta = tuple(c2j + dxj * self.div - c1j for c1j, c2j, dxj in zip(cube1, cube2, junc.delta_x))

        return Junction.get_junc(
            spec1.base_map * p1.specs[cnum1],
            spec2.base_map * p2.specs[cnum2],
            der_delta,
        )

    def init_pairs_tree(self):
        auto_junc = Junction.get_auto_junc(dim=self.dim)
        G = self.genus
        for junc in self.gen_auto_junctions():
            for cnum1 in range(G):
                for cnum2 in range(cnum1 + 2, G):
                    pos1 = self.get_piece_position(junc.spec1.pnum, cnum1)
                    pos2 = self.get_piece_position(junc.spec2.pnum, cnum2)
                    yield pieces.CurvePieceBalancedPair(self, junc, pos1, pos2)

        for junc in self.get_junctions_info().keys():
            last_cnum1 = 0 if junc.spec1.base_map.time_rev else G - 1
            first_cnum2 = G - 1 if junc.spec2.base_map.time_rev else 0
            for cnum1 in range(G):
                for cnum2 in range(G):
                    if (cnum1, cnum2) == (last_cnum1, first_cnum2):
                        continue
                    pos1 = self.get_piece_position(pnum=junc.spec1.pnum, cnum=cnum1)
                    pos2 = self.get_piece_position(pnum=junc.spec2.pnum, cnum=cnum2)
                    yield pieces.CurvePieceBalancedPair(self, junc, pos1, pos2)

    # если отношение кривой больше upper_bound, дальше не идём
    def estimate_ratio(self, ratio_func, rel_tol_inv=1000, upper_bound=None, max_iter=None, sat_pack=100, find_model=False, verbose=False):
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
        # это набор (pnum, cnum, curve), где в curve заданы:
        # - все начальные и конечные спеки (для нахождения производных стыков)
        # - спека в (pnum,cnum), (pnum,cnum+1)
        configs = []
        G = self.genus
        P = self.pattern_count
        variants = []
        for pnum in range(P):
            variants.append(self.gen_allowed_specs(pnum=pnum, cnum=0))
            variants.append(self.gen_allowed_specs(pnum=pnum, cnum=G - 1))

        for variant in itertools.product(*variants):
            variant = list(variant)
            gate_specs = {}
            for pnum in range(P):
                gate_specs[(pnum, 0)] = variant.pop(0)
                gate_specs[(pnum, G - 1)] = variant.pop(0)

            for pnum in range(P):
                for cnum in range(G - 1):
                    pc1 = (pnum, cnum)
                    for sp1 in self.gen_allowed_specs(*pc1):
                        if pc1 in gate_specs and gate_specs[pc1] != sp1:
                            continue
                        pc2 = (pnum, cnum + 1)
                        for sp2 in self.gen_allowed_specs(*pc2):
                            if pc2 in gate_specs and gate_specs[pc2] != sp2:
                                continue
                            def_specs = gate_specs.copy()
                            def_specs[pc1] = sp1
                            def_specs[pc2] = sp2
                            curve = self
                            for pc, spec in def_specs.items():
                                p, c = pc
                                curve = curve.specify(pnum=p, cnum=c, spec=spec)

                            configs.append((pnum, cnum, curve))

        # конфигурации, приводящие к данному стыку
        junc_curves = {}

        for pnum, cnum, curve in configs:
            base_junc = curve.get_base_junction(pnum=pnum, cnum=cnum)
            for junc in curve.gen_junctions_from_base([base_junc]):
                junc_curves.setdefault(junc, []).append(curve)

        return junc_curves


# patterns_symm - список пар (для каждого паттерна):
#   repr_specs - список представителей (coset representatives) spec-ов
#   symmetries - симметрии кривой (пока только base_map-ы, pnum-ы не меняем!)
class SymmFuzzyCurve(FuzzyCurve):
    def __init__(self, *args, **kwargs):
        patterns_symm = kwargs.pop('patterns_symm')
        super().__init__(*args, **kwargs)
        self.patterns_symm = tuple(patterns_symm)

    def changed(self, *args, **kwargs):
        if 'patterns_symm' not in kwargs:
            kwargs['patterns_symm'] = self.patterns_symm
        return super().changed(*args, **kwargs)

    def reverse(self):
        new_patterns_symm = []
        for repr_specs, symmetries in self.patterns_symm:
            new_patterns_symm.append((tuple(reversed(repr_specs)), symmetries))
        return super().reverse().changed(patterns_symm=new_patterns_symm)

    def apply_cube_map(self, cube_map):
        curve = super().apply_cube_map(cube_map)
        inv = cube_map.inverse()
        conj_cache = {}
        def conjugate(sp):
            if sp not in conj_cache:
                conj_cache[sp] = cube_map * sp * inv
            return conj_cache[sp]

        new_patterns_symm = []
        for repr_specs, symmetries in curve.patterns_symm:
            new_repr_specs = [conjugate(sp) for sp in repr_specs]
            new_symmetries = [conjugate(bm) for bm in symmetries]
            new_patterns_symm.append((tuple(new_repr_specs), new_symmetries))

        return curve.changed(patterns_symm=new_patterns_symm)

    def gen_allowed_specs(self, pnum, cnum):
        pattern = self.patterns[pnum]
        if pattern.specs[cnum] is not None:
            # базовое преобразование уже определено!
            yield pattern.specs[cnum]
            return

        repr_specs, symmetries = self.patterns_symm[pnum]
        repr_spec = self.patterns_symm[pnum][0][cnum]
        # симметрии надо брать для другого шаблона!
        symmetries = self.patterns_symm[repr_spec.pnum][1]
        for symm in symmetries:
            yield repr_spec * symm  # в этом порядке?

    @classmethod
    def init_from_brkline(cls, dim, div, brkline, allow_time_rev):
        proto = [brk[0] for brk in brkline]

        cube_first, entr_first, _ = brkline[0]
        entr = tuple(Fraction(cj + ej, div) for cj, ej in zip(cube_first, entr_first))

        cube_last, _, exit_last = brkline[-1]
        exit = tuple(Fraction(cj + ej, div) for cj, ej in zip(cube_last, exit_last))

        symmetries = []
        for bm in gen_constraint_cube_maps(dim, {entr: entr, exit: exit}):
            symmetries.append(bm)
        if allow_time_rev:
            for bm in gen_constraint_cube_maps(dim, {entr: exit, exit: entr}):
                symmetries.append(bm.reverse_time())

        specs = [None] * len(proto)
        repr_specs = [None] * len(proto)

        for cnum, brk in enumerate(brkline):
            cube, rel_entr, rel_exit = brk
            repr_specs[cnum] = Spec(next(gen_constraint_cube_maps(dim, {entr: rel_entr, exit: rel_exit})))

        return cls(
            dim=dim, div=div,
            patterns=[(proto, specs)],
            patterns_symm=[(repr_specs, symmetries)],
        )
