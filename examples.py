# coding: utf-8

from fractal_curve import FractalCurve
from base_map import BaseMap


# some examples of curves
def get_hilbert_curve():
    """Example of fractal curve due to D.Hilbert."""
    return FractalCurve(
        dim=2, div=2,
        proto=[(0, 0), (0, 1), (1, 1), (1, 0)],
        base_maps=[
            BaseMap([1, 0], [False, False]),  # (x,y)->(y,x)
            BaseMap(dim=2),  # (x,y)->(x,y)
            BaseMap(dim=2),  # (x,y)->(x,y)
            BaseMap([1, 0], [True, True]),  # (x,y)->(1-y,1-x)
        ],
    )


def get_peano_curve():
    """Example of fractal curve due to G.Peano."""
    id_map = BaseMap(dim=2)
    x_map = BaseMap([0, 1], [True, False])  # (x,y)->(1-x,y)
    y_map = BaseMap([0, 1], [False, True])  # (x,y)->(x,1-y)
    xy_map = BaseMap([0, 1], [True, True])  # (x,y)->(1-x,1-y)
    return FractalCurve(
        dim=2, div=3,
        proto=[(0, 0), (0, 1), (0, 2), (1, 2), (1, 1), (1, 0), (2, 0), (2, 1), (2, 2)],
        base_maps=[
            id_map, x_map, id_map,
            y_map, xy_map, y_map,
            id_map, x_map, id_map,
        ],
    )


def get_tokarev_curve():
    """3-D curve."""
    dim = 3
    chain_code = 'kjKikJK'
    bases = ['jki', 'kij', 'kij', 'iJK', 'iJK', 'KIj', 'KIj', 'JkI']
    return FractalCurve(
        dim=dim, div=2,
        chain_code=chain_code,
        base_maps=[BaseMap(basis=b) for b in bases],
    )


def get_rev_curve():
    """Curve with time reversal."""
    return FractalCurve(
        dim=2,
        div=2,
        proto=[(0, 0), (0, 1), (1, 1), (1, 0)],
        base_maps=[
            BaseMap([1, 0], [False, False]),  # (x,y)->(y,x)
            BaseMap([0, 1], [True, False], True),  # (x,y)->(1-x,y), t->-t
            BaseMap(dim=2),  # (x,y)->(x,y)
            BaseMap([1, 0], [True, True]),  # (x,y)->(1-y,1-x)
        ],
    )

# def get_N_curve():
#    proto = ['j','j','i','J','J','i','j','j']
#    base_maps = ['ij','Ji','ji','iJ','JI','jI','ij','Ji','ji']


# разрывная кривая - начало в центре квадрата
def get_discontinuous_curve():
    return FractalCurve(
        dim=2, div=3,
        proto=[(1, 1), (0, 1), (0, 0), (1, 0), (2, 0), (2, 1), (2, 2), (1, 2), (0, 2)],
        base_maps=[BaseMap(dim=2)] * 9,
    )