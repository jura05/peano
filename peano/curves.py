from collections import namedtuple
from functools import lru_cache
import itertools
import logging

from .fast_fractions import FastFraction
from .base_maps import BaseMap, Spec
from .utils import get_periodic_sum, gen_words
from .paths import Proto, CurvePath, Gate, PathsGenerator


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
        if last_spec is None:
            raise KeyError
        return spec.base_map * last_spec

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

    def specify(self, pnum, cnum, spec, force=False):
        """
        Check that we can set specs to spec at pnum, cnum, and return specified curve if so.

        This is the main method while dividing pairs_tree in estimators,
        so the efficiency is important here!
        """
        if not force:
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
    # Entrance/exit - defined if beg/end proto/specs are defined
    #

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
            #cur_curve = cur_spec * self
            cur_curve_proto = cur_spec.base_map * self.patterns[cur_spec.pnum].proto
            if cur_curve_proto[cnum] is None:
                raise Exception("Curve not specified enough to get cubes sequence!")
            cube = cur_curve_proto[cnum]

            cubes.append(cube)
            index[cur_spec] = len(cubes)-1

            #cur_spec = cur_curve.specs[cnum] * cur_spec
            cur_spec = self.compose_spec(cur_spec, cnum)

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

    def is_internal(self):
        """Check that curve has entrance/exit inside cube.

        Returns True is so, False if all gates are on the cube boundary, None if can't decide.
        """

        result = False
        dim = self.dim
        N = self.div
        genus = self.genus
        id_map = BaseMap.id_map(dim)
        all_edges_list = [(j, 0) for j in range(dim)] + [(j, 1) for j in range(dim)]

        for pnum in range(self.pattern_count):
            start_spec = Spec(base_map=id_map, pnum=pnum)
            for cnum in [0, genus-1]:
                if self.patterns[pnum].specs[cnum] is None:
                    continue

                # code from _get_cubes :( but faster!!!
                cur_spec = start_spec  # current curve = cur_spec * self
                seen = set()

                # limit should be on some edge, store possible edges
                possible_edges = set(all_edges_list)

                while True:
                    active_cnum = self.genus-1-cnum if cur_spec.base_map.time_rev else cnum
                    proto = self.patterns[cur_spec.pnum].proto
                    if proto[active_cnum] is None:
                        raise Exception("Curve not specified enough to get cubes sequence!")
                    cube = cur_spec.base_map.apply_cube(self.div, proto[active_cnum])

                    for j in range(dim):
                        if cube[j] != 0:
                            possible_edges.discard((j, 0))
                        if cube[j] != N-1:
                            possible_edges.discard((j, 1))

                    if not possible_edges:
                        return True

                    try:
                        cur_spec = self.compose_spec(cur_spec, cnum)
                    except KeyError:
                        result = None
                        break
                    if cur_spec in seen:
                        break  # OK
                    seen.add(cur_spec)

        return result

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

    def get_paths(self):
        # TODO: unite with check
        P = self.pattern_count
        gates = {pnum: Gate(self.get_entrance(pnum), self.get_exit(pnum)) for pnum in range(P)}
        paths = []
        for pnum, pattern in enumerate(self.patterns):
            pattern_gates = [spec.base_map * gates[spec.pnum] for spec in pattern.specs]
            paths.append(CurvePath(pattern.proto, pattern_gates))
        return paths

    def forget(self, allow_time_rev=False):
        """Convert curve to a fuzzy curve, saving entrance/exit and forgetting all specs."""
        return PathFuzzyCurve.init_from_paths(self.get_paths(), allow_time_rev=allow_time_rev)

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


class PathFuzzyCurve(FuzzyCurve):
    """
    Fuzzy curve with fixed paths.

    Additional attributes:
    Information about "standard" gates used in paths:
    .gates_symmetries  --  dict {std_gate: symmetries} - list of bms that save std_gate (doest not change at *)
    .gates_std  --  dict {std_gate: {pnum: std_map}}, such that std_map * pattern_global_gate = std_gate

    Information about patterns:
    .pattern_gates  --  array of std_gates for each fraction of each pattern;
    .pattern_reprs  --  array of reprs for each pattern; reprs[cnum] = bm that sends std_gate to gate of fraction
    """

    def __init__(self, *args, **kwargs):
        gates_symmetries = kwargs.pop('gates_symmetries')
        gates_std = kwargs.pop('gates_std')
        pattern_gates = kwargs.pop('pattern_gates')
        pattern_reprs = kwargs.pop('pattern_reprs')
        super().__init__(*args, **kwargs)
        self.gates_symmetries = gates_symmetries
        self.gates_std = gates_std
        self.pattern_gates = pattern_gates
        self.pattern_reprs = pattern_reprs

    def changed(self, *args, **kwargs):
        for add_field in ['gates_symmetries', 'gates_std', 'pattern_gates', 'pattern_reprs']:
            if add_field not in kwargs:
                kwargs[add_field] = getattr(self, add_field)
        return super().changed(*args, **kwargs)

    def reversed(self):
        # we do not change gates to keep them standard (symmetries also do not change)
        new_gates_std = {}
        for gate, data in self.gates_std.items():
            new_gates_std[gate] = {pnum: std_map.reversed_time() for pnum, std_map in data.items()}

        new_pattern_gates = [list(reversed(gates)) for gates in self.pattern_gates]
        new_pattern_reprs = []
        for reprs in self.pattern_reprs:
            new_reprs = list(reversed([bm.reversed_time() for bm in reprs]))
            new_pattern_reprs.append(new_reprs)
        return super().reversed().changed(
            gates_std=new_gates_std,
            pattern_gates=new_pattern_gates,
            pattern_reprs=new_pattern_reprs,
        )

    def apply_cube_map(self, cube_map):
        curve = super().apply_cube_map(cube_map)
        new_gates_std = {}
        for gate, data in self.gates_std.items():
            new_gates_std[gate] = {pnum: std_map * ~cube_map for pnum, std_map in data.items()}

        new_pattern_reprs = []
        for reprs in self.pattern_reprs:
            new_reprs = [cube_map * repr for repr in reprs]
            new_pattern_reprs.append(new_reprs)

        return curve.changed(
            gates_std=new_gates_std,
            pattern_reprs=new_pattern_reprs,
        )

    def gen_allowed_specs(self, pnum, cnum):
        pattern = self.patterns[pnum]
        if pattern.specs[cnum] is not None:
            yield pattern.specs[cnum]
            return

        gate = self.pattern_gates[pnum][cnum]
        repr = self.pattern_reprs[pnum][cnum]
        symmetries = self.gates_symmetries[gate]
        std = self.gates_std[gate]
        for pnum in sorted(std.keys()):
            for symm in symmetries:
                # how go we get fraction from a pattern:
                # 1) map pattern gate to std_gate
                # 2) apply symmetries for std_gate
                # 3) apply repr to get fraction gate
                yield Spec(repr * symm * std[pnum], pnum)

    @classmethod
    def init_from_paths(cls, paths, allow_time_rev=True):
        dim = paths[0].dim
        div = paths[0].div
        if allow_time_rev:
            possible_maps = list(BaseMap.gen_base_maps(dim))
        else:
            possible_maps = list(BaseMap.gen_base_maps(dim, time_rev=False))
        gates_symmetries = {}
        gates_std = {}
        for pnum, path in enumerate(paths):
            # pnum stands for path_num and also pattern_num
            std_gate = path.gate.std()
            std_map = next(bm for bm in possible_maps if bm * path.gate == std_gate)
            gates_std.setdefault(std_gate, {})[pnum] = std_map

        for gate in gates_std:
            gates_symmetries[gate] = [bm for bm in possible_maps if bm * gate == gate]

        patterns = []
        pattern_gates = []
        pattern_reprs = []
        for pnum, path in enumerate(paths):
            proto = path.proto
            specs = [None] * len(proto)  # nothing is defined
            patterns.append((proto, specs))

            reprs = []
            gates = []
            for gate in path.gates:
                repr0, gate0 = None, None
                for std_gate in gates_std:
                    allowed = [bm for bm in possible_maps if bm * std_gate == gate]
                    if allowed:
                        # should be exactly once!
                        repr0 = allowed[0]
                        gate0 = std_gate
                        break
                if repr0 is None or gate0 is None:
                    raise Exception("Can't create PathFuzzyCurve: no allowed gates")
                reprs.append(repr0)
                gates.append(gate0)
            pattern_reprs.append(reprs)
            pattern_gates.append(gates)

        return cls(
            dim=dim, div=div, patterns=patterns,
            gates_symmetries=gates_symmetries,
            gates_std=gates_std,
            pattern_gates=pattern_gates,
            pattern_reprs=pattern_reprs,
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
        self._hash = hash(self._data())

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
        return self._hash

    def __repr__(self):
        return '{} | dx={}, dt={} | {}'.format(self.spec1, self.delta_x, self.delta_t, self.spec2)


def gen_possible_gates(dim, div, pattern_count):
    """
    Generate all gates such that there is at least one path with thems.

    We consider only non-internal curves, i.e., curves with all gates on cube boundary.
    Indeed, if curve in internal, then only internal curves can "use" it,
    and non-internal curves use each other and have lower ratio.

    All objects are standartized: obj -> std(obj) = min(bm * obj for bm in base_maps).
    """

    genus = div**dim

    # curves are non-internal, so first and last cubes are on the boundary
    all_cubes = itertools.product(range(div), repeat=dim)
    edge_cubes = (cube for cube in all_cubes if any(cj == 0 or cj == div - 1 for cj in cube))
    edge_pairs = itertools.combinations(edge_cubes, 2)

    def std_pair(pair):
        def gen_all():
            for bm in BaseMap.gen_base_maps(dim, time_rev=False):
                c1 = bm.apply_cube(div, pair[0])
                c2 = bm.apply_cube(div, pair[1])
                yield (c1, c2)
                yield (c2, c1)
        return min(gen_all())

    std_pairs = sorted(set(std_pair(pair) for pair in edge_pairs))

    # as base_maps are not defined yet, we can permute pairs (i.e., assume that pairs tuple is sorted)
    # so we use combinations with replacement instead of product
    std_pairs_list = list(itertools.combinations_with_replacement(std_pairs, r=pattern_count))
    logging.info('edge_cube pairs: %d', len(std_pairs_list))

    # we can standartize all found gates, with caching for speedup
    @lru_cache(maxsize=2**16)
    def get_std_gate(gate):
        return gate.std()

    @lru_cache(maxsize=2**16)
    def get_std_gates(gates):
        return tuple(sorted(get_std_gate(gate) for gate in gates))

    seen_gates = set()

    # to find gates, we should specify:
    # * first and last cubes in proto (for each pattern) -- will use std_pairs_list for that
    # * specs on that cubes -- see below
    #
    # with fixed protos, we must test all variants of specs
    # one may store chosen specs in a "word" - tuple of length wlen and each entry in range(wbase)
    # here wlen is then number of specs to specify, and wbase is number of variants for each spec

    spec_types = [(pnum, cnum) for pnum in range(pattern_count) for cnum in [0, genus-1]]
    spec_variants = [Spec(bm, pnum) for bm in BaseMap.gen_base_maps(dim) for pnum in range(pattern_count)]
    wbase = len(spec_variants)
    wlen = len(spec_types)
    logging.info('must specify %d specs with %d variants', wlen, wbase)

    iter_no = 0
    good_count = 0
    for p_no, pairs in enumerate(std_pairs_list):
        protos = []
        for pair in pairs:
            proto = [None] * genus
            proto[0] = pair[0]
            proto[-1] = pair[-1]
            protos.append(Proto(dim, div, proto))

        # curve that given word represents
        def get_word_curve(word):
            vars = {spec_types[k]: spec_variants[wk] for k, wk in enumerate(word)}
            pattern_specs = []
            for pnum in range(pattern_count):
                specs = [None] * genus
                pattern_specs.append(specs)

            for spec_type, spec in vars.items():
                spec_pnum, spec_cnum = spec_type
                pattern_specs[spec_pnum][spec_cnum] = spec

            patterns = [(proto, specs) for proto, specs in zip(protos, pattern_specs)]
            return FuzzyCurve(dim=dim, div=div, patterns=patterns)

        def check_word(word):
            return not get_word_curve(word).is_internal()  # boundary or can't decide

        # ok, now we must check all words; we use the following optimization:
        # if some specifications (i.e., for first pattern) already forbid that curve is boundary
        # (i.e., exit of first pattern must be inside cube), then there is no need to
        # specify other specs; in terms of "words", if some subword is bad, than all words
        # that extend it are also bad; this logic is implemented in gen_words function

        for word in gen_words(base=wbase, length=wlen, check_func=check_word):
            iter_no += 1
            if iter_no % 1000 == 0:
                logging.warning(
                    'cube_pairs %d/%d, word (base %d): %s, boundary: %d, gates: seen: %d, good: %d',
                    p_no + 1, len(std_pairs_list),
                    wbase, word,
                    iter_no, len(seen_gates), good_count,
                )

            curve = get_word_curve(word)

            # the curve is checked, so gates are necessarily on the cube boundary
            gates = tuple(Gate(curve.get_entrance(pnum), curve.get_exit(pnum)) for pnum in range(pattern_count))
            std_gates = get_std_gates(gates)
            if std_gates in seen_gates:
                continue
            seen_gates.add(std_gates)

            # check that there is at least one path with given gates
            pg = PathsGenerator(dim=dim, div=div, gates=std_gates)
            if pg.get_paths_example(start_max_count=1000, finish_max_count=100000):
                logging.warning('found gates: %s', [str(g) for g in std_gates])
                good_count += 1
                yield std_gates
