from collections import namedtuple

from .common import Pattern, Junction
from .base_maps import BaseMap, Spec


# pnum - выделенный шаблон, задаёт реальную кривую [0,1]->[0,1]^d
class FuzzyPolyCurve:
    def __init__(self, dim, div, patterns, pnum=0):
        self.dim = dim
        self.div = div
        self.patterns = tuple(patterns)
        self.pattern_count = len(patterns)
        self.pnum = pnum
        self.genus = div**dim

    def reverse(self):
        """Reverse time in a curve."""
        new_patterns = []
        for pattern in self.patterns:
            new_pattern = Pattern(
                proto=tuple(reversed(pattern.proto)),
                specs=tuple(reversed(pattern.specs)),
            )
            new_patterns.append(new_pattern)

        return type(self)(
            dim=self.dim,
            div=self.div,
            patterns=new_patterns,
        )

    def __rmul__(self, other):
        """Apply base map to a fractal curve, return new curve."""
        if isinstance(other, Spec):
            new_curve = other.base_map * self
            new_curve.pnum = self.pnum
            return new_curve
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

            # базовые преобразования сопрягаются: действительно, чтобы получить
            # из преобразованной кривой её фракцию, можно сделать так:
            # - сначала вернуться к исходной кривой (inv)
            # - применить преобразование исходной кривой для перехода к фракции (bm)
            # - перейти к преобразованной кривой (base_map)
            inv = base_map.inverse()

            new_specs = [base_map * spec * inv if spec is not None else None for spec in pattern.specs]
            new_pattern = Pattern(new_proto, new_specs)
            new_patterns.append(new_pattern)

        return type(self)(
            dim=self.dim,
            div=self.div,
            patterns=new_patterns,
        )

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
