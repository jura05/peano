from collections import namedtuple

from .base_maps import BaseMap
from .common import Junction, Spec


class FuzzyPolyCurve:
    def __init__(self, dim, div, patterns):
        self.dim = dim
        self.div = div
        self.patterns = tuple(patterns)
        self.pattern_count = len(patterns)
        self.genus = div**dim

    def __getitem__(self, pnum):
        return self.patterns[pnum]

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
        pattern = self[pnum]
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
        p1 = self[spec1.pnum]
        p2 = self[spec2.pnum]

        cnum1 = 0 if spec1.base_map.time_rev else -1
        cnum2 = -1 if spec2.base_map.time_rev else 0

        if p1.specs[cnum1] is None or p2.specs[cnum2] is None:
            raise Exception("Can't get derivative: spec not defined")

        cube1 = spec1.base_map.apply_cube(self.div, p1.proto[cnum1])  # сейчас не нужно - нормализуем cube_map1 -> id
        cube2 = spec2.base_map.apply_cube(self.div, p2.proto[cnum2])
        der_delta = tuple(c2j + dxj * self.div - c1j for c1j, c2j, dxj in zip(cube1, cube2, junc.delta_x))

        bm1 = spec1.base_map * p1.specs[cnum1].base_map
        bm2 = spec2.base_map * p2.specs[cnum2].base_map

        return Junction.get_junc(
            Spec(base_map=bm1, pnum=p1.specs[cnum1].pnum),
            Spec(base_map=bm2, pnum=p2.specs[cnum2].pnum),
            der_delta,
        )
