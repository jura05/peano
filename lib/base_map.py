# coding: utf-8

from fractions import Fraction
import functools
import itertools

from sympy.combinatorics.permutations import Permutation


class BaseMap:
    """Base map: isometry of cube and (possibly) time reversal.
    Immutable and hashable.
    Acts on function f:[0,1]->[0,1]^d as  Bf: [0,1]--B_t-->[0,1]--f-->[0,1]^d--B_x-->[0,1]^d
    """

    def __init__(self, perm, flip, time_rev=False):
        """Create a BaseMap instance.
        Params:
        perm, flip: params, which define isometry of cube
                      perm = [k_0,...,k_{d-1}]
                      flip = [b_0,...,b_{d-1}]
                      define the map
                      (x_0,...,x_{d-1}) --> (f(b_0;x_{k_0}), ..., f(b_{d-1};x_{k_{d-1}})),
                      where f(b;x)=x if b is False, and 1-x otherwise
        time_rev:   time reversal (boolean), default: False
        """

        self.dim = len(perm)

        # store data in tuples to make object immutable
        self.perm = tuple(perm)
        self.flip = tuple(bool(b) for b in flip)
        self.time_rev = bool(time_rev)
        self._data = (self.perm, self.flip, self.time_rev)
        self._hash = hash(self._data)

    @classmethod
    def id_map(cls, dim):
        perm = range(dim)
        flip = [False] * dim
        return cls(perm, flip)

    def cube_map(self):
        return type(self)(self.perm, self.flip)

    def __eq__(self, other):
        return self._data == other._data

    def __hash__(self):
        return self._hash

    def __repr__(self):
        if self.dim > 3:
            letters = ['x_{}'.format(i) for i in range(self.dim)]
        else:
            letters = "xyz"[:self.dim]
        s = "(" + ",".join(letters) + ")"
        s += "->("
        s += ",".join([("1-{}" if b else "{}").format(letters[k]) for k, b in zip(self.perm, self.flip)])
        s += ")"
        if self.time_rev:
            s += ",t->1-t"
        return s

    def __mul__(self, other):
        """Composition of base maps."""
        assert self.dim == other.dim
        perm = []
        flip = []
        for i in range(self.dim):
            p = self.perm[i]
            k = other.perm[p]
            b1 = self.flip[i]
            b2 = other.flip[p]
            perm.append(k)
            flip.append(b1 ^ b2)
        time_rev = self.time_rev ^ other.time_rev
        return type(self)(perm, flip, time_rev)

    def inverse(self):
        """Inverse of base map."""
        perm = [None]*self.dim
        flip = [None]*self.dim
        for i, k, b in zip(range(self.dim), self.perm, self.flip):
            perm[k] = i
            flip[k] = b
        return type(self)(perm, flip, self.time_rev)

    def reverse_time(self):
        return type(self)(self.perm, self.flip, not self.time_rev)

    def is_oriented(self):
        oriented = True
        for f in self.flip:
            if f: oriented = not oriented
        if Permutation(self.perm).signature() == -1:
            oriented = not oriented
        return oriented

    def apply_x(self, x):
        """Apply isometry to a point x."""
        return tuple(1-x[k] if b else x[k] for k, b in zip(self.perm, self.flip))

    def apply_x2(self, x, l):
        """Apply isometry to a point x of [0,l]^d."""
        return tuple(l-x[k] if b else x[k] for k, b in zip(self.perm, self.flip))

    def apply_vec(self, v):
        """Apply linear part of isometry to a vector."""
        return tuple(-v[k] if b else v[k] for k, b in zip(self.perm, self.flip))

    def apply_cube(self, div, cube):
        """Apply isometry to a sub-cube."""
        return tuple(div-cube[k]-1 if b else cube[k] for k, b in zip(self.perm, self.flip))

    def apply_cube_start(self, cube_start, cube_length):
        """Apply isometry to a cube of given length, return it's min (start) vertex."""
        return tuple(1-cube_start[k]-cube_length if b else cube_start[k] for k, b in zip(self.perm, self.flip))

    def apply_edge(self, edge):
        """Apply isometry to an edge. Not implemented!"""
        pass

    def apply_cnum(self, genus, cnum):
        return genus - 1 - cnum if self.time_rev else cnum


class PieceMap:
    """Special isometry of cube and time.
    Composition of some BaseMap, shift and scaling.
    """

    def __init__(self, base_map, shift, scale, time_shift, time_scale):
        """Create a PieceMap instance.
        Params:
            base_map
            shift
            scale       params that define space map:
                        x |--> base_map(shift + scale * x)

            time_shift  
            time_scale  params that define time map:
                        t |--> time_shift + time_scale * t
        """
        assert base_map.dim == len(shift)
        self.dim = base_map.dim
        self.base_map = base_map
        self.shift = shift
        self.scale = scale
        self.time_shift = time_shift
        self.time_scale = time_scale

    def _data(self):
        return (self.base_map, self.shift, self.scale, self.time_shift, self.time_scale)

    def __eq__(self, other):
        return self._data() == other._data()

    def __mul__(self, other):
        """Composition."""
        d = self.dim
        # B(c + lB'(c'+l'x)) = (BB')(a + ll'x), где a = lc' - B'^{-1}(0)/l + B'^{-1}(c)
        inv = other.base_map.inverse()
        inv_0 = inv.apply_x(tuple(0 for j in range(d)))
        inv_c = inv.apply_x(self.shift)
        new_shift = tuple(self.scale * other.shift[j] - self.scale * inv_0[j] + inv_c[j] for j in range(d))
        return type(self)(
            base_map=self.base_map * other.base_map,
            shift=new_shift,
            scale=self.scale * other.scale,
            time_shift=self.time_shift + self.time_scale * other.time_shift,
            time_scale=self.time_scale * other.time_scale,
        )

    def inverse(self):
        d = self.base_map.dim
        # y = B(c + lx)  <==>  x = B^{-1}(a + y/l), где a = B(-c/l)-B(0)/l
        y1 = self.base_map.apply_x(tuple(-Fraction(self.shift[j], self.scale) for j in range(d)))
        y2 = self.base_map.apply_x(tuple(0 for j in range(d)))
        new_shift = tuple(y1[j] - Fraction(y2[j], self.scale) for j in range(d))
        return type(self)(
            base_map=self.base_map.inverse(),
            shift=new_shift,
            scale=Fraction(1, self.scale),
            time_shift=-Fraction(self.time_shift, self.time_scale),
            time_scale=Fraction(1, self.time_scale),
        )

    def apply(self, x, t):
        """Apply isometry to a space point x and time t."""
        d = self.base_map.dim
        new_x = tuple(self.shift[j] + self.scale * x[j] for j in range(d))
        new_x = self.base_map.apply_x(new_x)
        new_t = self.time_shift + self.time_scale * t
        return (new_x, new_t)

    @classmethod
    def id_map(cls, d):
        return cls(
            base_map=BaseMap.id_map(d),
            shift=tuple(0 for j in range(d)),
            scale=1,
            time_shift=0,
            time_scale=1,
        )


def gen_base_maps(dim, time_rev=None):
    time_rev_variants = [True, False] if time_rev is None else [time_rev]
    for perm in itertools.permutations(range(dim)):
        for flip in itertools.product([True, False], repeat=dim):
            for time_rev in time_rev_variants:
                yield BaseMap(perm, flip, time_rev)

def gen_constraint_cube_maps(dim, points_map):
    for bm in gen_base_maps(dim, time_rev=False):
        if all(bm.apply_x(src) == dst for src, dst in points_map.items()):
            yield bm
