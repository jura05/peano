import itertools

#from sympy.combinatorics.permutations import Permutation

from .fast_fractions import FastFraction


class BaseMap:
    """
    Base map: isometry of cube and (possibly) time reversal.

    Immutable and hashable.
    Acts on function f:[0,1]->[0,1]^d as  Bf: [0,1]--B_time-->[0,1]--f-->[0,1]^d--B_cube-->[0,1]^d

    We severily use caching, because there are not so many possible base maps.
    We restrict to dim <= 6, but the limit may by increased.
    """

    basis_letters = 'ijklmn'
    var_letters = 'xyz'

    _obj_cache = {}

    def __new__(cls, coords, time_rev=False):
        """
        Create a BaseMap instance.

        Cached method. 
        Params:
            coords:     list of pairs (k, b) defining the isometry
                        B_cube(x_0,...,x_{d-1}) = (y_0,...,y_{d-1}),
                        where y_j = x_{k_j} if b_j else 1 - x_{k_j}

            time_rev:   time reversal (boolean), default: False
                        B_time(t) = 1-t if time_rev else t
        """
        
        coords = tuple((k, bool(b)) for k, b in coords)
        time_rev = bool(time_rev)

        cache = cls._obj_cache
        key = (coords, time_rev)
        if key in cache:
            return cache[key]

        obj = super().__new__(cls)
        obj.dim = len(coords)
        obj.coords = coords
        obj.time_rev = time_rev
        obj._key = key
        obj._mul_cache = {}
        obj._inv_cache = None

        cache[key] = obj
        return obj

    @classmethod
    def id_map(cls, dim):
        """Identity map."""
        coords = [(k, False) for k in range(dim)]
        return cls(coords)

    @classmethod
    def from_basis(cls, basis):
        """
        Convenient way to represent a base map.

        basis -- string, e.g., 'Ij', representing images of standard basis e_1, e_2, ...
            i=e_1, j=e_2, k=e_3, l=e_4, m=e_5, n=e_6...;
            upper-case letters correspond to negation of vectors: I = -i, J = -j, ...
            To get time reverse, place '1' at the end of the string, e.g., 'iJ1'
        """

        if basis[-1] in '01':
            time_rev = True if basis[-1] == '1' else False
            basis = basis[:-1]
        else:
            time_rev = False

        assert len(basis) <= len(cls.basis_letters)

        l2i = {l: i for i, l in enumerate(cls.basis_letters)}
        coords = []
        for l in basis:
            lk = l.lower()
            coords.append((l2i[lk], (l != lk)))

        return BaseMap(coords, time_rev)

    def repr_basis(self):
        """Concise representation of base_map (see from_basis)."""
        assert self.dim <= len(self.basis_letters)
        basis = ''
        for k, b in self.coords:
            letter = self.basis_letters[k]
            basis += (letter.upper() if b else letter)
        if self.time_rev:
            basis += '1'
        return basis

    def cube_map(self):
        """Isometry without time reversal."""
        return BaseMap(self.coords)

    def __eq__(self, other):
        return self._key == other._key

    def __hash__(self):
        return hash((self._key))

    def __repr__(self):
        if self.dim > len(self.var_letters):
            letters = ['x_{}'.format(i) for i in range(self.dim)]
        else:
            letters = self.var_letters[:self.dim]
        s = "(" + ",".join(letters) + ")"
        s += "->("
        s += ",".join([("1-{}" if b else "{}").format(letters[k]) for k, b in self.coords])
        s += ")"
        if self.time_rev:
            s += ",t->1-t"
        return s

    def __mul__(self, other):
        """
        Composition of base maps: A * B.
        """
        if not isinstance(other, BaseMap):
            return NotImplemented

        key = other._key
        val = self._mul_cache.get(key, None)
        if val is None:
            # actual multiplication
            assert self.dim == other.dim
            coords = []
            for i in range(self.dim):
                p, b1 = self.coords[i]
                k, b2 = other.coords[p]
                coords.append((k, b1 ^ b2))
            time_rev = self.time_rev ^ other.time_rev
            val = BaseMap(coords, time_rev)
            self._mul_cache[key] = val
        return val

    def __invert__(self):
        """Inverse of base map: ~B."""
        val = self._inv_cache
        if val is None:
            # actual inversion
            coords = [None] * self.dim
            for i, c in enumerate(self.coords):
                k, b = c
                coords[k] = (i, b)
            self._inv_cache = val = BaseMap(coords, self.time_rev)
        return val

    def reversed_time(self):
        return BaseMap(self.coords, not self.time_rev)

    def conjugate_by(self, other):
        """Conjugation."""
        return other * self * ~other

    def is_oriented(self):
        raise NotImplementedError
        oriented = True
        perm = []
        for k, b in self.coords:
            if b: oriented = not oriented
            perm.append(k)
        if Permutation(perm).signature() == -1:
            oriented = not oriented
        return oriented

    def apply_x(self, x, mx=1):
        """Apply isometry to a point x of [0,1]^d."""
        return tuple(mx-x[k] if b else x[k] for k, b in self.coords)

    def apply_x_fraction(self, x):
        """FastFraction version."""
        x = [FastFraction(xj, 1) if isinstance(xj, int) else xj for xj in x]
        return tuple(FastFraction(1, 1) - x[k] if b else x[k] for k, b in self.coords)

    def apply_t(self, t, mt=1):
        return mt - t if self.time_rev else t

    def apply_vec(self, v):
        """Apply linear part of isometry to a vector."""
        return tuple(-v[k] if b else v[k] for k, b in self.coords)

    def apply_cube(self, div, cube):
        """Apply isometry to a sub-cube."""
        return tuple(div-cube[k]-1 if b else cube[k] for k, b in self.coords)

    def apply_cube_start(self, cube_start, cube_length):
        """Apply isometry to a cube of given length, return it's min (start) vertex."""
        return tuple(1-cube_start[k]-cube_length if b else cube_start[k] for k, b in self.coords)

    def apply_edge(self, edge):
        """Apply isometry to an edge. Not implemented!"""
        pass

    def apply_cnum(self, genus, cnum):
        return genus - 1 - cnum if self.time_rev else cnum


def gen_base_maps(dim, time_rev=None):
    time_rev_variants = [True, False] if time_rev is None else [time_rev]
    for perm in itertools.permutations(range(dim)):
        for flip in itertools.product([True, False], repeat=dim):
            for time_rev in time_rev_variants:
                yield BaseMap(zip(perm, flip), time_rev)


def gen_constraint_cube_maps(dim, points_map):
    for bm in gen_base_maps(dim, time_rev=False):
        if all(bm.apply_x_fraction(src) == dst for src, dst in points_map.items()):
            yield bm


class Spec:
    """
    Specification of poly-fractal curve to a fraction: BaseMap and pattern choice.

    Specs act on polyfractal curves and form a semi-group.
    """

    # we do not make obj_id, to provide consistency for pickling
    _obj_cache = {}

    def __new__(cls, base_map, pnum=0):
        key = (base_map._key, pnum)
        cache = cls._obj_cache
        if key in cache:
            return cache[key]

        obj = super().__new__(cls)
        obj.base_map = base_map
        obj.pnum = pnum
        obj._key = key
        cache[key] = obj
        return obj

    def __eq__(self, other):
        return self._key == other._key

    def __hash__(self):
        return hash(self._key)

    def __repr__(self):
        return '{}[{}]'.format(self.base_map, self.pnum)

    def __mul__(self, other):
        """
        Composition of specs, i.e., (S1*S2) f = S1(S2 f).

        We also allow to combine specs and base maps,
        as they act on curves too.
        """
        if isinstance(other, BaseMap):
            other_bm = other
        elif isinstance(other, Spec):
            other_bm = other.base_map
        else:
            return NotImplemented
        return Spec(self.base_map * other_bm, self.pnum)

    def __rmul__(self, other):
        if isinstance(other, BaseMap):
            other_bm = other
            pnum = self.pnum
        elif isinstance(other, Spec):
            other_bm = other.base_map
            pnum = other.pnum
        else:
            return NotImplemented
        return Spec(other_bm * self.base_map, pnum)

    def reversed_time(self):
        return Spec(self.base_map.reversed_time(), self.pnum)

    def conjugate_by(self, other):
        """Conjugate by base map (not spec!)."""
        assert isinstance(other, BaseMap)
        return Spec(self.base_map.conjugate_by(other), self.pnum)
