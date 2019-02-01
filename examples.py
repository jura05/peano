# coding: utf-8

from fractal_curve import FractalCurve
from base_map import BaseMap
from utils import chain2proto, basis2base_map


# Minimal 2D monofractal curve in L_1 (9)
def get_hilbert_curve():
    """Example of fractal curve due to D.Hilbert."""
    return FractalCurve(
        proto=[(0, 0), (0, 1), (1, 1), (1, 0)],
        base_maps=[
            BaseMap([1, 0], [False, False]),  # (x,y)->(y,x)
            BaseMap(dim=2),                   # (x,y)->(x,y)
            BaseMap(dim=2),                   # (x,y)->(x,y)
            BaseMap([1, 0], [True, True]),    # (x,y)->(1-y,1-x)
        ],
    )


def get_peano_curve():
    """Example of fractal curve due to G.Peano."""
    id_map = BaseMap(dim=2)
    x_map = BaseMap([0, 1], [True, False])  # (x,y)->(1-x,y)
    y_map = BaseMap([0, 1], [False, True])  # (x,y)->(x,1-y)
    xy_map = BaseMap([0, 1], [True, True])  # (x,y)->(1-x,1-y)
    return FractalCurve(
        proto=[(0, 0), (0, 1), (0, 2), (1, 2), (1, 1), (1, 0), (2, 0), (2, 1), (2, 2)],
        base_maps=[
            id_map, x_map, id_map,
            y_map, xy_map, y_map,
            id_map, x_map, id_map,
        ],
    )


# Minimal 2D monofractal curve in L_inf (5.333) and L_2 (5.667)
def get_meurthe_curve():

    chain_code = 'jjiJJijj'
    bases = ['ij','Ji','ij','jI','JI','iJ','ji','Ji','ij']
    return FractalCurve(
        proto=chain2proto(chain_code),
        base_maps=[basis2base_map(b) for b in bases],
    )
    
    
def get_coil_curve():

    chain_code = 'jjiJJijj'
    bases = ['ji','Ji','ji','jI','JI','jI','ji','Ji','ji']
    return FractalCurve(
        proto=chain2proto(chain_code),
        base_maps=[basis2base_map(b) for b in bases],
    )


def get_serpentine_curve():

    chain_code = 'jjiJJijj'
    bases = ['ij','Ji','ji','iJ','JI','iJ','ji','Ji','ij']
    return FractalCurve(
        proto=chain2proto(chain_code),
        base_maps=[basis2base_map(b) for b in bases],
    )


def get_R_curve():

    chain_code = 'jjiiJIJi'
    bases = ['ji','ji','ij','ij','ij','IJ','JI','JI','ij']
    return FractalCurve(
        proto=chain2proto(chain_code),
        base_maps=[basis2base_map(b) for b in bases],
    )


# Minimal 3D monofractal curve with time reversal in L_inf (12.4)
def get_haverkort_curve_1():
    """3-D curve with time reversal."""
    chain_code = 'kjKikJK'
    bases = ['kji0','jik0','kIj1','iKJ0','IKJ1','KIj0','Kij1','Jki1']
    return FractalCurve(
        proto=chain2proto(chain_code),
        base_maps=[basis2base_map(b) for b in bases],
    )


# Minimal 3D monofractal curve with time reversal in L_1 (89.8) and L_2 (18.6)
def get_haverkort_curve_2():
    """3-D curve with time reversal."""
    chain_code = 'kjKikJK'
    bases = ['KIJ1','KJI1','KjI0','Jki1','jki0','kjI1','kJI0','iKJ0']
    return FractalCurve(
        proto=chain2proto(chain_code),
        base_maps=[basis2base_map(b) for b in bases],
    )


def get_tokarev_curve():
    """3-D curve."""
    chain_code = 'kjKikJK'
    bases = ['jki', 'kij', 'kij', 'iJK', 'iJK', 'KIj', 'KIj', 'JkI']
    return FractalCurve(
        proto=chain2proto(chain_code),
        base_maps=[basis2base_map(b) for b in bases],
    )


def get_rev_curve():
    """Curve with time reversal."""
    return FractalCurve(
        proto=[(0, 0), (0, 1), (1, 1), (1, 0)],
        base_maps=[
            BaseMap([1, 0], [False, False]),                # (x,y)->(y,x)
            BaseMap([0, 1], [True, False], time_rev=True),  # (x,y)->(1-x,y), t->-t
            BaseMap(dim=2),                                 # (x,y)->(x,y)
            BaseMap([1, 0], [True, True]),                  # (x,y)->(1-y,1-x)
        ],
    )


# разрывная кривая - начало в центре квадрата
def get_discontinuous_curve():
    return FractalCurve(
        proto=[(1, 1), (0, 1), (0, 0), (1, 0), (2, 0), (2, 1), (2, 2), (1, 2), (0, 2)],
        base_maps=[BaseMap(dim=2)] * 9,
    )
