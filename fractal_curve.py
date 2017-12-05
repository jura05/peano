# coding: utf-8

class BaseMap:
    """Base map: isometry of cube and (possibly) time reversal.
    Params:
        cube_map   isometry of cube
                      S=[(k_0,b_0),...,(k_{d-1},b_{d-1})] which defines map
                      (x_0,...,x_{d-1}) -> (f_{b_0}(x_{k_0}),...,f_{k_{d-1}}(x_{k_{d-1}})),
                      where f_b(x)=x if b is False, and f_b(x)=1-x otherwise
        time_rev   time reversal (boolean), default: False
    """
    def __init__(self, cube_map, time_rev=False):
        self.cube_map = cube_map
        self.time_rev = time_rev

    def apply(self, x):
        """Apply isometry to a point x of [0,1]^d."""
        return tuple(1-x[ki] if bi else x[ki] for ki, bi in self.cube_map)



class FractalCurve:
    """Class representing fractal peano curve in [0,1]^d.
    Params:
        div         positive integer, number of divisions of each side of the cube (characteristic of a curve)
        dim         dimension (default: 2)
        proto       curve prototype - sequence of "cubes" enumerated as (x_0,..,x_{d-1}), 0<=x_j<div
        entrance    entrance point
        exit        exit point
        base_maps   sequence of base maps (instances of BaseMap class)

    For examples see get_hilbert_curve()
    """
    def __init__(self, div, proto, base_maps, exit, entrance=(0,0), dim=2):
        self.div = div
        self.dim = dim
        self.proto = proto
        self.entrance = entrance
        self.exit = exit
        self.base_maps = base_maps

    def genus(self):
        """Fractal genus of the curve."""
        return self.div ** self.dim

    def check(self):
        """Check consistency of curve params."""
        n = self.div
        d = self.dim

        # dummy checks
        assert n > 0
        assert len(self.proto) == self.genus(), 'bad proto length'

        for cube in self.proto:
            for j in range(d):
                assert 0 <= cube[j] < n, 'bad cube coordinates'
        sqset = set(tuple(cube) for cube in self.proto)
        assert len(sqset) == len(self.proto), 'non-unique cubes'

        for j in range(d):
            assert 0 <= self.entrance[j] < n , 'bad entrance'
            assert 0 <= self.exit[j] < n, 'bad exit'

        # проверяем соответствие входов-выходов
        gates = []  # пары (вход,выход), реальные координаты в кубе, умноженные на div
        entrance_n = tuple(n*self.entrance[j] for j in range(d))  # начальный вход
        gates.append((None, entrance_n))
        for cube, bm in zip(self.proto, self.base_maps):
            entrance_pos = bm.apply(self.entrance)
            entrance_n = tuple(cube[j] + entrance_pos[j] for j in range(d))
            exit_pos = bm.apply(self.exit)
            exit_n = tuple(cube[j] + exit_pos[j] for j in range(d))
            if bm.time_rev:
                gates.append((exit_n,entrance_n))
            else:
                gates.append((entrance_n,exit_n))
        exit_n = tuple(n*self.exit[j] for j in range(d))
        gates.append((exit_n, None))

        for i in range(len(gates)-1):
            assert gates[i][1] == gates[i+1][0], 'exit does not correspond to entrance'


def get_hilbert_curve():
    """Example of fractal curve due to D.Hilbert."""
    return FractalCurve(
        dim=2, div=2,
        entrance=(0,0), exit=(1,0),
        proto=[(0,0), (0,1), (1,1), (1,0)],
        base_maps=[
            BaseMap([(1,False),(0,False)]),  # (x,y)->(y,x)
            BaseMap([(0,False),(1,False)]),  # (x,y)->(x,y)
            BaseMap([(0,False),(1,False)]),  # (x,y)->(x,y)
            BaseMap([(1,True),(0,True)]),  # (x,y)->(1-y,1-x)
        ],
    )

def get_peano_curve():
    """Example of fractal curve due to G.Peano."""
    id_map = BaseMap([(0,False),(1,False)])
    x_map  = BaseMap([(0,True),(1,False)])  # (x,y)->(1-x,y)
    y_map  = BaseMap([(0,False),(1,True)])  # (x,y)->(x,1-y)
    xy_map = BaseMap([(0,True),(1,True)])   # (x,y)->(1-x,1-y)
    return FractalCurve(
        dim=2, div=3,
        entrance=(0,0), exit=(1,1),
        proto=[(0,0),(0,1),(0,2),(1,2),(1,1),(1,0),(2,0),(2,1),(2,2)],
        base_maps=[
            id_map, x_map, id_map,
            y_map, xy_map, y_map,
            id_map, x_map, id_map,
        ],
    )

def main():
    # run some tests
    for curve in [get_hilbert_curve(), get_peano_curve()]:
        curve.check()
    print("ok!")

if __name__ == "__main__":
    main()
