# coding: utf-8

class BaseMap:
    """Base map: isometry of cube and (possibly) time reversal.
    Immutable and hashable.
    Acts on function f:[0,1]->[0,1]^d as  Bf: [0,1]--B_t-->[0,1]--f-->[0,1]^d--B_x-->[0,1]^d
    """

    def __init__(self, perm=None, flip=None, dim=None, time_rev=False):
        """Create a BaseMap instance.
        Params:
        perm, flip: params, which define isometry of cube
                      perm = [k_0,...,k_{d-1}]
                      flip = [b_0,...,b_{d-1}]
                      define the map
                      (x_0,...,x_{d-1}) --> (f(b_0;x_{k_0}), ..., f(b_{d-1};x_{k_{d-1}})),
                      where f(b;x)=x if b is False, and 1-x otherwise
          or
        dim:        dimension only (defines identity map)

        time_rev:   time reversal (boolean), default: False
        """

        if perm is None or flip is None:
            if dim is None:
                raise Exception("Can't init base_map: define perm, flip")
            perm = list(range(dim))
            flip = [False]*dim

        assert len(perm) == len(flip)
        self.dim = len(perm)

        # store data in tuples to make object immutable
        self.perm = tuple(perm)
        self.flip = tuple(bool(b) for b in flip)
        self.time_rev = bool(time_rev)

    def cube_map(self):
        return type(self)(self.perm, self.flip, False)

    def _data(self):
        return (self.perm, self.flip, self.time_rev)

    def __eq__(self, other):
        return self._data() == other._data()

    def __hash__(self):
        return hash(self._data())

    def __str__(self):
        if self.dim > 3:
            letters = ['x_{}'.format(i) for i in range(self.dim)]
        else:
            letters = "xyz"[:self.dim]
        s = "(" + ",".join(letters) + ")"
        s += "->("
        s += ",".join([("1-{}" if b else "{}").format(letters[k]) for k, b in zip(self.perm, self.flip)])
        s += ")"
        if self.time_rev:
            s += "t->-t"
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

    def apply_x(self, x):
        """Apply isometry to a point x of [0,1]^d."""
        return tuple(1-x[k] if b else x[k] for k, b in zip(self.perm, self.flip))

    def apply_vec(self, v):
        """Apply linear part of isometry to a vector."""
        return tuple(-v[k] if b else v[k] for k, b in zip(self.perm, self.flip))

    def apply_cube(self, div, cube):
        """Apply isometry to a sub-cube."""
        return tuple(div-cube[k]-1 if b else cube[k] for k, b in zip(self.perm, self.flip))

    def apply_edge(self, edge):
        """Apply isometry to an edge. Not implemented!"""
        pass

