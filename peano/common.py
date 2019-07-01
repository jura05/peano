from collections import namedtuple

from .base_maps import BaseMap


Pattern = namedtuple('Pattern', ['proto', 'specs'])

class Spec:
    def __init__(self, base_map, pnum=0):
        self.base_map = base_map
        self.pnum = pnum

    def _data(self):
        return (self.base_map, self.pnum)

    def __eq__(self, other):
        return self._data() == other._data()

    def __hash__(self):
        return hash(self._data())


# TODO: автостыки - там только pnum ...

# стык двух фракций кривой
# всегда приводим к стандартному виду:
# - pnum1 <= pnum2
# - первая фракция стандартной пространственной ориентации, но, возможно, с обращением времени (self.time_rev)
# куб второй фракции получается из куба первого сдвигом на delta_x \in {0,1,-1}^d
class Junction:
    def __init__(self, spec1, spec2, delta_x, delta_t):
        self.spec1 = spec1
        self.spec2 = spec2
        self.delta_x = delta_x
        self.delta_t = delta_t

        # for backwards compatibility
        self.time_rev = self.spec1.base_map.time_rev
        self.base_map = self.spec2.base_map

    # обычный стык
    @classmethod
    def get_junc(cls, spec1, spec2, delta_x):
        if spec1.pnum > spec2.pnum:
            # меняем местами
            delta_x = tuple(-dj for dj in delta_x)
            spec1, spec2 = spec2, spec1

        bm1 = spec1.base_map
        bm2 = spec2.base_map

        if spec1.pnum == spec2.pnum and bm1.time_rev and bm2.time_rev:
            # обращаем время и меняем местами
            delta_x = tuple(-dj for dj in delta_x)
            bm1, bm2 = bm2.reverse_time(), bm1.reverse_time()

        bm1_cube_inv = bm1.cube_map().inverse()
        return cls(
            spec1=Spec(base_map=bm1_cube_inv * bm1, pnum=spec1.pnum),  # only possible time_rev
            spec2=Spec(base_map=bm1_cube_inv * bm2, pnum=spec2.pnum),
            delta_x=bm1_cube_inv.apply_vec(delta_x),
            delta_t=1,
        )

    # автостык
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
        return '1: ' + ('t->1-t' if self.time_rev else '') + ', 2: ' + str(self.base_map) + ', --> ' + str(self.delta_x)
