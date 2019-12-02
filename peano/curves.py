from collections import namedtuple
import itertools

from .fast_fractions import FastFraction
from .base_maps import BaseMap, Spec, gen_constraint_cube_maps
from .utils import get_periodic_sum


class Proto:
    """Curve prototype -- sequence of cubes."""
    def __init__(self, dim, div, cubes):
        self.dim = dim
        self.div = div
        self._cubes = tuple(tuple(cube) for cube in cubes)

    def __iter__(self):
        return iter(self._cubes)

    def __getitem__(self, idx):
        return self._cubes[idx]

    def __len__(self):
        return len(self._cubes)

    def __eq__(self, other):
        return self._cubes == other._cubes

    def __rmul__(self, base_map):
        if not isinstance(base_map, BaseMap):
            return NotImplemented
        cubes = [base_map.apply_cube(self.div, cube) for cube in self._cubes]
        if base_map.time_rev:
            cubes = reversed(cubes)
        return type(self)(self.dim, self.div, cubes)


Pattern = namedtuple('Pattern', ['proto', 'specs'])



class FuzzyCurve:
    """
    Polyfractal peano curve, not fully specified.

    Object attributes:
    .proto -- prototype for selected pattern
    .specs -- specs for selected pattern
    .genus -- subj
    """

    def __init__(self, dim, div, patterns, pnum=0):
        """
        Create FuzzyCurve instance.

        Params:
        dim -- dimension d of image [0,1]^d of the curve
        div -- number of divisions for each of the coordinates, so genus = G = div**dim
        patterns -- list of patterns
                    each pattern is a tuple (proto, specs),
                    proto -- prototype, list of cubes (with integer coords) of length G
                    specs -- list of specs (Spec or None) of length G
        pnum -- selected pattern, to define actual curve f:[0,1]->[0,1]^d

        Note that we assume that each pattern has the same div.
        """
        self.dim = dim
        self.div = div

        self.patterns = []
        for proto, specs in patterns:
            proto = proto if isinstance(proto, Proto) else Proto(dim, div, proto)
            specs = tuple(sp if isinstance(sp, Spec) else Spec(sp) if sp is not None else None for sp in specs)
            pattern = Pattern(proto=proto, specs=specs)
            self.patterns.append(pattern)

        self.pattern_count = len(self.patterns)
        self.pnum = pnum

        self.proto = self.patterns[pnum].proto
        self.specs = self.patterns[pnum].specs

        self.genus = div**dim

    def get_fraction(self, cnum):
        """First-order fraction of a curve."""
        return self.specs[cnum] * self

    def changed(self, patterns=None, pnum=None, **kwargs):
        """Create curve with changed parameters."""
        return type(self)(
            dim=self.dim,
            div=self.div,
            patterns=patterns if patterns is not None else self.patterns,
            pnum=pnum if pnum is not None else self.pnum,
            **kwargs,
        )

    def reversed(self):
        """Reverse time in a curve."""
        # reversal of the curve implies reversal of all patterns
        # sequence of specs, and proto, are reversed
        # each base_map does not change:
        #   - pnums are the same
        #   - if there was not time_rev, there should not be after reversal
        #   - cube map does not change, because geometry does not change
        new_patterns = [(reversed(pattern.proto), reversed(pattern.specs)) for pattern in self.patterns]
        return self.changed(patterns=new_patterns)

    def apply_cube_map(self, base_map):
        """Apply cube isometry to a curve."""

        if base_map.dim != self.dim:
            raise Exception("Incompatible base map!")
        if base_map.time_rev:
            raise Exception("Do not use this method with time_rev!")

        new_patterns = []
        for pattern in self.patterns:
            # isometry for the prototype
            new_proto = base_map * pattern.proto

            # specs are conjugated: to get from mapped curve its fraction, we can:
            # - map curve to original (base_map^{-1})
            # - then map curve to the fraction (spec)
            # - then map the whole thing to mapped curve (base_map)
            new_specs = [spec.conjugate_by(base_map) if spec is not None else None for spec in pattern.specs]
            new_patterns.append((new_proto, new_specs))

        return self.changed(patterns=new_patterns)

    def __rmul__(self, other):
        """Apply base map or spec to a fractal curve, return new curve."""
        if isinstance(other, Spec):
            base_map = other.base_map
            pnum = other.pnum
        elif isinstance(other, BaseMap):
            base_map = other
            pnum = self.pnum
        else:
            return NotImplemented

        # base_map is the composition of commutating maps: cube_map and time_map
        curve = self
        if base_map.time_rev:
            curve = curve.reversed()
        curve = curve.apply_cube_map(base_map.cube_map())

        if curve.pnum != pnum:
            curve = curve.changed(pnum=pnum)

        return curve

    def compose_spec(self, spec, cnum):
        """
        Returns spec X, such that: (spec * self).specs[cnum] * (spec * self) = X * self.

        Does not require actual curve multiplication (it is slow).
        Method allows to get orientations of deep fractions of a curve.
        """

        active_pattern = self.patterns[spec.pnum]
        active_cnum = spec.base_map.apply_cnum(self.genus, cnum)
        last_spec = active_pattern.specs[active_cnum]
        return spec.base_map * last_spec * ~spec.base_map * spec

    def gen_allowed_specs(self, pnum, cnum):
        raise NotImplementedError("Define in child class")

    def gen_possible_curves(self):
        """
        Generate all curves, compatible with self.

        We use gen_allowed_specs (defined in child classes) and suppose
        that specs for different fractions are independent.
        This is very important condition, which is provided by continuity.
        """

        sp_variant_generators = []
        G = self.genus

        # large list of generators of specs for each (pnum, cnum); flatten by pnum * G + cnum
        for pnum in range(self.pattern_count):
            for cnum in range(G):
                sp_variant_generators.append(self.gen_allowed_specs(pnum, cnum))

        for all_specs in itertools.product(*sp_variant_generators):
            patterns = []
            for pnum, pattern in enumerate(self.patterns):
                specs = all_specs[(G * pnum):(G * (pnum + 1))]
                patterns.append((pattern.proto, specs))

            yield Curve(dim=self.dim, div=self.div, patterns=patterns, pnum=self.pnum)

    def sp_info(self):
        """List of triples (pnum, cnum, spec) of defined specs."""
        curve_info = []
        for pnum, pattern in enumerate(self.patterns):
            for cnum, spec in enumerate(pattern.specs):
                if spec is not None:
                    curve_info.append((pnum, cnum, spec))
        return curve_info

    def is_specialization(self, tmpl):
        """Check is self has all of defined specs of the given curve, and they are the same."""
        for pnum, cnum, sp in tmpl.sp_info():
            if self.patterns[pnum].specs[cnum] != sp:
                return False
        return True

    def specify(self, pnum, cnum, spec):
        """
        Check that we can set specs to spec at pnum, cnum, and return specified curve if so.

        This is the main method while dividing pairs_tree in estimators,
        so the efficiency is important here!
        """
        if spec not in self.gen_allowed_specs(pnum, cnum):
            raise Exception("Can't specify curve")

        pattern = self.patterns[pnum]
        if pattern.specs[cnum] is not None:
            return self  # optimization

        new_specs = list(pattern.specs)
        new_specs[cnum] = spec
        new_pattern = (pattern.proto, new_specs)

        new_patterns = list(self.patterns)
        new_patterns[pnum] = new_pattern

        return self.changed(patterns=new_patterns)

    #
    # Junctions; see also the class Junction
    #

    def gen_auto_junctions(self):
        for pnum in range(self.pattern_count):
            yield Junction.get_auto_junc(dim=self.dim, pnum=pnum)

    def gen_junctions_from_base(self, base_juncs):
        """Yield base junctions and their derivatives."""
        seen = set()
        to_derive = []
        for junc in base_juncs:
            yield junc
            seen.add(junc)
            to_derive.append(junc)

        while to_derive:
            junc = to_derive.pop()
            dj = self.get_derived_junction(junc)
            if dj not in seen:
                yield dj
                seen.add(dj)
                to_derive.append(dj)

    def get_base_junction(self, pnum, cnum):
        """Get base junctions if both specs are defined"""
        pattern = self.patterns[pnum]
        spec1, spec2 = pattern.specs[cnum], pattern.specs[cnum + 1]
        if spec1 is None or spec2 is None:
            raise Exception("Can't get base junction: spec is not defined!")
        delta_x = [c2j - c1j for c1j, c2j in zip(pattern.proto[cnum], pattern.proto[cnum + 1])]
        return Junction.get_junc(spec1, spec2, delta_x)

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

        cube1 = spec1.base_map.apply_cube(self.div, p1.proto[cnum1])  # spec1.base_map == id, so no need in this
        cube2 = spec2.base_map.apply_cube(self.div, p2.proto[cnum2])
        der_delta = tuple(c2j + dxj * self.div - c1j for c1j, c2j, dxj in zip(cube1, cube2, junc.delta_x))

        return Junction.get_junc(
            spec1.base_map * p1.specs[cnum1],
            spec2.base_map * p2.specs[cnum2],
            der_delta,
        )

    # словарь {junc: curve_list} кривых, приводящих к данному стыку
    def get_junctions_info(self):
        """
        Info about possible junctions.

        Returns dict {junc: curves} with specified curves that have given junction.
        """

        # config is a tuple (pnum, cnum, curve), where in curve there are specified:
        # - all begin and end specs (to get derivatives)
        # - specs at (pnum, cnum) and (pnum, cnum+1)
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

        junc_curves = {}
        for pnum, cnum, curve in configs:
            base_junc = curve.get_base_junction(pnum=pnum, cnum=cnum)
            for junc in curve.gen_junctions_from_base([base_junc]):
                junc_curves.setdefault(junc, []).append(curve)

        return junc_curves


class Curve(FuzzyCurve):
    """
    Fully-specified regular Peano curve.

    This curve defines the continuous surjective map f:[0,1]->[0,1]^d.
    """

    def get_entrance(self, pnum=None):
        """
        Entrance of a curve, i.e. point f(0).

        If pnum is set, find the entrance of pattern pnum.
        """
        if pnum is None:
            pnum = self.pnum
        start, period = self._get_cubes(pnum, 0)
        return self._get_cube_limit(start, period)

    def get_exit(self, pnum=None):
        """
        Exit of a curve, i.e. point f(1).

        If pnum is set, find the entrance of pattern pnum.
        """
        if pnum is None:
            pnum = self.pnum
        start, period = self._get_cubes(pnum, self.genus-1)
        return self._get_cube_limit(start, period)

    def _get_cubes(self, pnum, cnum):
        # we found the sequence of cubes that we obtain if we take cube #cnum in each fraction
        # returns pair (non-periodic part, periodic part)
        # возвращает пару (непериодическая часть, периодическая часть)
        cur_spec = Spec(base_map=BaseMap.id_map(self.dim), pnum=pnum)  # current curve = cur_spec * self
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
        """Get limit of sequence of nested cubes."""
        p = [0] * self.dim
        for j in range(self.dim):
            start_j = [x[j] for x in start]
            period_j = [x[j] for x in period]
            p[j] = get_periodic_sum(start_j, period_j, self.div)
        return tuple(p)

    def gen_base_junctions(self):
        """Generate base junctions."""
        seen = set()
        for pnum in range(self.pattern_count):
            for cnum in range(self.genus - 1):
                junc = self.get_base_junction(cnum=cnum, pnum=pnum)
                if junc not in seen:
                    yield junc
                    seen.add(junc)

    def gen_junctions(self):
        """Generate all regular junctions for a curve."""
        yield from self.gen_junctions_from_base(self.gen_base_junctions())

    def gen_allowed_specs(self, pnum, cnum):
        yield self.patterns[pnum].specs[cnum]

    # TODO: TEST FUZZY POLY CURVE!!!
    def forget(self, allow_time_rev=False):
        """Convert curve to a fuzzy curve, saving entrance/exit and forgetting all specs."""
        patterns = []
        patterns_symm = []

        for pnum, pattern in enumerate(self.patterns):
            entr = self.get_entrance(pnum)
            exit = self.get_exit(pnum)

            symmetries = []
            for bm in gen_constraint_cube_maps(self.dim, {entr: entr, exit: exit}):
                symmetries.append(bm)
                
            if allow_time_rev:
                for bm in gen_constraint_cube_maps(self.dim, {entr: exit, exit: entr}):
                    symmetries.append(bm.reversed_time())

            repr_specs = pattern.specs  # now old specs are only representatives
            patterns_symm.append((repr_specs, symmetries))
            patterns.append((pattern.proto, [None] * self.genus))  # forget about old specs!

        return SymmFuzzyCurve(
            dim=self.dim,
            div=self.div,
            patterns=patterns,
            patterns_symm=patterns_symm,
        )

    def get_vertex_moments(self, pnum=None):
        """Get dict {vertex: first_visit_time}."""
        if pnum is None:
            pnum = self.pnum

        # TODO: refactor, use itertools
        vertices = [[]]
        for j in range(self.dim):
            new_vertices = []
            for xj in [0,1]:
                for edge in vertices:
                    new_edge = edge + [xj]
                    new_vertices.append(new_edge)
            vertices = new_vertices

        moment = {}
        for v in vertices:
            moment[tuple(v)] = self.get_edge_touch(v, pnum)
        return moment

    def get_edge_touch(self, edge, pnum=None):
        """
        Moment of first edge touch.

        Edge is a tuple of {0,1,None} defining a set
        {(x_0,...,x_{d-1}): x_i==0 if e[i]==0, x_i==1 if e[i]==1, or arbitrary x[i] if e[i] is None.
        E.g., tuples (0,0,0) or (0,1,1) define vertices
        """
        if pnum is None:
            pnum = self.pnum

        cur_spec = Spec(base_map=BaseMap.id_map(self.dim), pnum=pnum)
        curve = cur_spec * self  # invariant
        cnums = []
        index = {}
        while True:
            cnum = curve._get_edge_cnum(edge)
            sp = curve.specs[cnum]
            cnums.append(cnum)

            index[cur_spec] = len(cnums)-1

            cur_spec = sp * cur_spec
            curve = sp * curve

            if cur_spec in index:
                period_start = index[cur_spec]
                break

        return self._get_time_limit(cnums[0:period_start], cnums[period_start:])

    def _get_edge_cnum(self, edge):
        """First cube from prototype touching edge."""
        # какой куб из прототипа первым касается грани
        N = self.div
        for cnum, cube in enumerate(self.proto):
            # проверяем, что куб касается грани
            touch = True
            for x, e in zip(cube, edge):
                if e is None:
                    continue
                elif e == 1 and x != (N-1):
                    touch = False
                    break
                elif e == 0 and x != 0:
                    touch = False
                    break
            if touch:
                return cnum

    def _get_time_limit(self, start, period):
        # задана начальная и периодическая последовательность номеров кубов, считаем время
        return get_periodic_sum(start, period, self.genus)

    def get_subdivision(self, k=1):
        """Get k-th subdivision of a curve."""
        N = self.div
        current_curve = self
        for _ in range(k):
            new_patterns = []
            for pnum, curr_pattern in enumerate(current_curve.patterns):
                new_proto = []
                new_specs = []
                for cube, spec in zip(curr_pattern.proto, curr_pattern.specs):
                    proto, specs = self.patterns[spec.pnum]  # from original curve

                    if spec.base_map.time_rev:
                        proto = reversed(proto)
                        specs = reversed(specs)

                    for c in proto:
                        nc = spec.base_map.apply_cube(N, c)
                        new_cube = [cj*N + ncj for cj, ncj in zip(cube, nc)]
                        new_proto.append(new_cube)

                    # базовые преобразования для подраздедения:
                    # пусть (cube, spec) соответствуют i-й фракции
                    # в ней мы взяли j-ю подфракцию (sp)
                    # Какое преобразование переводит кривую в j-ю фракцию внутри i-й?
                    # - сначала к исходной кривой мы применим bm, чтобы перевести её в j-ю фракцию,
                    # - потом ко всей этой картинке применяем base_map, чтобы перевести всё в i-ю фракцию (base_map)
                    # можно сделать наоборот:
                    # - сначала кривую переводим в i-ю фракцию (base_map)
                    # - применяем внутри i-й фракции преобразования для перехода в j-ю
                    #   но там оно будет сопряженное: base_map * bm * base_map^{-1}, см. apply_base_map
                    for sp in specs:
                        new_specs.append(spec.base_map * sp)

                new_patterns.append((new_proto, new_specs))

            current_curve = type(self)(
                dim=self.dim,
                div=N*current_curve.div,  # we change div so do not use ``changed''
                patterns=new_patterns,
            )

        return current_curve

    def check(self):
        """
        Check consistency of curve params.

        Main check is the continuity of the curve.
        It is equivalent to exit/entrance correspondence for all of the patterns
        (there is a subtlety for non-active patterns - we check them, too).
        """

        d, n, G, P = self.dim, self.div, self.genus, self.pattern_count

        # dummy checks
        assert all(x > 0 for x in [d, n, G, P])

        for pattern in self.patterns:
            assert len(pattern.proto) == self.genus, 'bad proto length'
            for cube in pattern.proto:
                for j in range(d):
                    assert 0 <= cube[j] < n, 'bad cube coordinates'
            assert len(set(pattern.proto)) == len(pattern.proto), 'non-unique cubes'

        # main check - continuity - for all patterns
        entr = {pnum: self.get_entrance(pnum) for pnum in range(P)}
        exit = {pnum: self.get_exit(pnum) for pnum in range(P)}
        for pnum, pattern in enumerate(self.patterns):
            curve_entr = entr[pnum]
            curve_exit = exit[pnum]

            gates = [(None, curve_entr)]  # pairs (entr, exit) in [0,1]^d

            for cube, spec in zip(pattern.proto, pattern.specs):
                bm = spec.base_map
                frac_gates_rel = [bm.apply_x_fraction(point) for point in [entr[spec.pnum], exit[spec.pnum]]]
                frac_gates = [tuple((FastFraction(c, 1) + e) * FastFraction(1, n) for c, e in zip(cube, point)) for point in frac_gates_rel]
                if bm.time_rev:
                    frac_gates.reverse()
                gates.append(frac_gates)

            gates.append((curve_exit, None))

            for i in range(len(gates)-1):
                if gates[i][1] != gates[i+1][0]:
                    msg = 'exit does not correspond to entrance at ' + str(i)
                    print(gates)
                    raise Exception(msg)


class SymmFuzzyCurve(FuzzyCurve):
    """
    Fuzzy curve with a group of symmetry acting on each fraction.

    .patterns_symm -  list of pairs (repr_specs, symmetries) for each pattern
        repr_specs -  list of left coset representatives of specs for each fraction of pattern
        symmetries -  symmetries (base_maps) subgroup
    """

    def __init__(self, *args, **kwargs):
        patterns_symm = kwargs.pop('patterns_symm')
        super().__init__(*args, **kwargs)
        self.patterns_symm = tuple(patterns_symm)

    def changed(self, *args, **kwargs):
        if 'patterns_symm' not in kwargs:
            kwargs['patterns_symm'] = self.patterns_symm
        return super().changed(*args, **kwargs)

    def reversed(self):
        new_patterns_symm = []
        for repr_specs, symmetries in self.patterns_symm:
            new_patterns_symm.append((tuple(reversed(repr_specs)), symmetries))
        return super().reversed().changed(patterns_symm=new_patterns_symm)

    def apply_cube_map(self, cube_map):
        curve = super().apply_cube_map(cube_map)

        new_patterns_symm = []
        for repr_specs, symmetries in curve.patterns_symm:
            new_repr_specs = [sp.conjugate_by(cube_map) for sp in repr_specs]
            new_symmetries = [bm.conjugate_by(cube_map) for bm in symmetries]
            new_patterns_symm.append((tuple(new_repr_specs), new_symmetries))

        return curve.changed(patterns_symm=new_patterns_symm)

    def gen_allowed_specs(self, pnum, cnum):
        pattern = self.patterns[pnum]
        if pattern.specs[cnum] is not None:
            # базовое преобразование уже определено!
            yield pattern.specs[cnum]
            return

        repr_spec = self.patterns_symm[pnum][0][cnum]
        symmetries = self.patterns_symm[repr_spec.pnum][1]  # Note! the symmetries for the other pattern
        for symm in symmetries:
            yield repr_spec * symm

    @classmethod
    def init_from_brkline(cls, dim, div, brkline, allow_time_rev):
        proto = [brk[0] for brk in brkline]

        cube_first, entr_first, _ = brkline[0]
        entr = tuple(FastFraction(cj + ej, div) for cj, ej in zip(cube_first, entr_first))

        cube_last, _, exit_last = brkline[-1]
        exit = tuple(FastFraction(cj + ej, div) for cj, ej in zip(cube_last, exit_last))

        symmetries = []
        for bm in gen_constraint_cube_maps(dim, {entr: entr, exit: exit}):
            symmetries.append(bm)
        if allow_time_rev:
            for bm in gen_constraint_cube_maps(dim, {entr: exit, exit: entr}):
                symmetries.append(bm.reversed_time())

        specs = [None] * len(proto)
        repr_specs = [None] * len(proto)

        for cnum, brk in enumerate(brkline):
            cube, rel_entr, rel_exit = brk
            rel_entr = tuple(FastFraction(xj, 1) if isinstance(xj, int) else xj for xj in rel_entr)
            rel_exit = tuple(FastFraction(xj, 1) if isinstance(xj, int) else xj for xj in rel_exit)
            repr_specs[cnum] = Spec(next(gen_constraint_cube_maps(dim, {entr: rel_entr, exit: rel_exit})))

        return cls(
            dim=dim, div=div,
            patterns=[(proto, specs)],
            patterns_symm=[(repr_specs, symmetries)],
        )


class Junction:
    """
    Junctions of two curve fractions.

    Attributes:
    .spec1:  first fraction is spec1 * curve
    .spec2:  second fractions is spec2 * curve
    .delta_x:  shift vector to get 2-nd fraction from 1-st, element of {0,1,-1}^d
    .delta_t:  time shift (=0 or 1, see below)

    Each junctions if standartized:
    - pnum1 <= pnum2
    - first spec has cube_map = id (but possibly with time_rev - this allows to get delta_t=1)

    All junctions are:
    - auto junctions (delta_t = 0)
    - regular junctions (delta_t = 1):
        * base junctions
        * derivative of base junctions

    Derivative is implemented in the curve class.
    """
    def __init__(self, spec1, spec2, delta_x, delta_t):
        self.spec1 = spec1
        self.spec2 = spec2
        self.delta_x = delta_x
        self.delta_t = delta_t

    @classmethod
    def get_junc(cls, spec1, spec2, delta_x):
        """Regular junction."""
        if spec1.pnum > spec2.pnum \
                or (spec1.pnum == spec2.pnum and spec1.base_map.time_rev and spec2.base_map.time_rev):
            # swap and reverse time
            delta_x = tuple(-dj for dj in delta_x)
            spec1, spec2 = spec2.reversed_time(), spec1.reversed_time()

        bm1_cube_inv = ~spec1.base_map.cube_map()
        return cls(
            spec1=bm1_cube_inv * spec1,  # only possible time_rev
            spec2=bm1_cube_inv * spec2,
            delta_x=bm1_cube_inv.apply_vec(delta_x),
            delta_t=1,
        )

    @classmethod
    def get_auto_junc(cls, dim, pnum=0):
        spec = Spec(base_map=BaseMap.id_map(dim), pnum=pnum)
        return cls(spec1=spec, spec2=spec, delta_x=(0,) * dim, delta_t=0)

    def _data(self):
        return (self.delta_x, self.spec1, self.spec2)

    def __eq__(self, other):
        return self._data() == other._data()

    def __hash__(self):
        return hash(self._data())

    def __repr__(self):
        return '{} | dx={}, dt={} | {}'.format(self.spec1, self.delta_x, self.delta_t, self.spec2)
