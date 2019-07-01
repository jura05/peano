from . import fuzzy_poly_curves


class PolyCurve(fuzzy_poly_curves.FuzzyPolyCurve):

    #
    # Стыки
    #
    # Для полностью заданной кривой можно определить все стыки

    def gen_base_junctions(self):
        for pnum in range(self.pattern_count):
            for cnum in range(self.genus - 1):
                yield self.get_base_junction(cnum=cnum, pnum=pnum)

    def gen_junctions(self):
        yield from self.gen_junctions_from_base(list(self.gen_base_junctions()))

    def get_junctions(self):
        return list(self.gen_junctions())
