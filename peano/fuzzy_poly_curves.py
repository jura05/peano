from collections import namedtuple

from .common import Junction
from .base_maps import BaseMap, Spec


# patterns - список пар (proto, specs)
# pnum - выделенный шаблон, задаёт реальную кривую [0,1]->[0,1]^d
class FuzzyPolyCurve:
    Pattern = namedtuple('Pattern', ['proto', 'specs'])

    def __init__(self, dim, div, patterns, pnum=0):
        self.dim = dim
        self.div = div

        self.patterns = []
        for proto, specs in patterns:
            proto = (tuple(cube) for cube in proto)
            specs = (sp if isinstance(sp, Spec) else Spec(sp) if sp is not None else None for sp in specs)
            pattern = FuzzyPolyCurve.Pattern(proto=tuple(proto), specs=tuple(specs))
            self.patterns.append(pattern)

        self.pattern_count = len(self.patterns)
        self.pnum = pnum

        # legacy
        self.proto = self.patterns[pnum].proto
        self.specs = self.patterns[pnum].specs
        self.base_maps = [sp.base_map if sp is not None else None for sp in self.specs]

        self.genus = div**dim

    # создать кривую с другими данным
    def changed(self, patterns=None, pnum=None, **kwargs):
        return type(self)(
            dim=self.dim,
            div=self.div,
            patterns=patterns if patterns is not None else self.patterns,
            pnum=pnum if pnum is not None else self.pnum,
            **kwargs,
        )

    def reverse(self):
        """Reverse time in a curve."""
        new_patterns = []
        for pattern in self.patterns:
            new_pattern = (
                # прототип проходится в обратном порядке
                reversed(pattern.proto),

                # базовые преобразования проходятся в обратном порядке
                # сами по себе они не меняются:
                #   - если обращения времени не было, то его и не будет
                #   - изометрия куба не меняется, т.к. время не играет роли
                reversed(pattern.specs),
            )
            new_patterns.append(new_pattern)

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

        # применяем изометрию куба

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
