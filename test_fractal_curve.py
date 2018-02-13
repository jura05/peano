#!/usr/bin/env python3
# coding: utf-8

from fractal_curve import FractalCurve, BaseMap

def get_hilbert_curve():
    """Example of fractal curve due to D.Hilbert."""
    return FractalCurve(
        dim=2, div=2,
        proto=[(0,0), (0,1), (1,1), (1,0)],
        base_maps=[
            BaseMap([1,0],[False,False]),  # (x,y)->(y,x)
            BaseMap.id_map(dim=2),  # (x,y)->(x,y)
            BaseMap.id_map(dim=2),  # (x,y)->(x,y)
            BaseMap([1,0],[True,True]),  # (x,y)->(1-y,1-x)
        ],
    )

def get_peano_curve():
    """Example of fractal curve due to G.Peano."""
    id_map = BaseMap.id_map(dim=2)
    x_map  = BaseMap([0,1],[True,False])  # (x,y)->(1-x,y)
    y_map  = BaseMap([0,1],[False,True])  # (x,y)->(x,1-y)
    xy_map = BaseMap([0,1],[True,True])   # (x,y)->(1-x,1-y)
    return FractalCurve(
        dim=2, div=3,
        proto=[(0,0),(0,1),(0,2),(1,2),(1,1),(1,0),(2,0),(2,1),(2,2)],
        base_maps=[
            id_map, x_map, id_map,
            y_map, xy_map, y_map,
            id_map, x_map, id_map,
        ],
    )

# разрывная кривая - начало в центре квадрата
def get_discontinuous_curve():
    return FractalCurve(
        dim=2, div=3,
        proto=[(1,1),(0,1),(0,0),(1,0),(2,0),(2,1),(2,2),(1,2),(0,2)],
        base_maps=[BaseMap.id_map(dim=2)]*9,
    )

def main():
    # run some tests
    for curve, name in [
        (get_hilbert_curve(), 'Hilbert'),
        (get_peano_curve(), 'Peano'),
        (get_discontinuous_curve(), 'discontinuous'),
    ]:
        print('curve:', name)
        print('entrance:', curve.entrance)
        print('exit:', curve.exit)
        try:
            id_map = BaseMap.id_map(curve.dim)
            for bm in curve.base_maps:
                # заодно проверим BaseMap
                inv = bm.inverse()
                assert bm * inv == id_map, 'base_map multiplication'
                assert inv * bm == id_map, 'base_map multiplication'
            curve.check()
            print('check ok!')
        except Exception as exc:
            print('check failed:', exc)
            continue

        for i in range(curve.genus()):
            subcurve = curve.get_subcurve(i)
            subcurve.check()
        print('subcurves are correct')

        subdiv_curve = curve.get_subdivision()
        subdiv_curve.check()
        print('subdivision is correct')

        for delta, base_map in curve.get_junctions():
            print('junction:', delta, base_map)
        print()


if __name__ == "__main__":
    main()
