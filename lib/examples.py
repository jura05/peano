# coding: utf-8

from fractal_curve import FractalCurve
from base_map import BaseMap
from utils import chain2proto, basis2base_map


# Proceedings of the Steklov Institute of Mathematics, 2008, Vol. 263, pp. 236–256.
# "Minimal Peano Curve" by E. V. Shchepin and K. E. Bauman
def get_scepin_bauman_curve():
    proto = (
        (0, 0), (0, 1), (0, 2),
        (1, 2), (1, 1), (1, 0),
        (2, 0), (2, 1), (2, 2),
    )
    base_maps = [
        BaseMap.id_map(dim=2),
        BaseMap([1,0],[True,False]),  # rot(90)
        BaseMap([1,0],[False,False]),  # (x,y)->(y,x)

        BaseMap([0,1],[False,True]),  # (x,y)->(x,1-y)
        BaseMap([1,0],[True,True]),  # (x,y)->(1-y,1-x)
        BaseMap([1,0],[False,True]),  # rot(-90)

        BaseMap.id_map(dim=2),
        BaseMap(perm=[1,0],flip=[True,False]),  # rot(90)
        BaseMap(perm=[1,0],flip=[False,False]),  # (x,y)->(y,x)
    ]
    return FractalCurve(dim=2, div=3, proto=proto, base_maps=base_maps)


# Minimal 2D monofractal curve in L_1 (9)
def get_hilbert_curve():
    """Example of fractal curve due to D.Hilbert."""
    return FractalCurve(
        proto=[(0, 0), (0, 1), (1, 1), (1, 0)],
        base_maps=[
            BaseMap([1, 0], [False, False]),  # (x,y)->(y,x)
            BaseMap.id_map(2),                # (x,y)->(x,y)
            BaseMap.id_map(2),                # (x,y)->(x,y)
            BaseMap([1, 0], [True, True]),    # (x,y)->(1-y,1-x)
        ],
    )


def get_peano_curve():
    """Example of fractal curve due to G.Peano."""
    id_map = BaseMap.id_map(2)
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

def get_peano5_curve():
    id_map = BaseMap.id_map(2)
    x_map = BaseMap([0, 1], [True, False])  # (x,y)->(1-x,y)
    y_map = BaseMap([0, 1], [False, True])  # (x,y)->(x,1-y)
    xy_map = BaseMap([0, 1], [True, True])  # (x,y)->(1-x,1-y)
    return FractalCurve(
        dim=2, div=5,
        proto=[
            (0, 0), (0, 1), (0, 2), (0, 3), (0, 4),
            (1, 4), (1, 3), (1, 2), (1, 1), (1, 0),
            (2, 0), (2, 1), (2, 2), (2, 3), (2, 4),
            (3, 4), (3, 3), (3, 2), (3, 1), (3, 0),
            (4, 0), (4, 1), (4, 2), (4, 3), (4, 4),
        ],
        base_maps=[
            id_map, x_map, id_map, x_map, id_map,
            y_map, xy_map, y_map, xy_map, y_map,
            id_map, x_map, id_map, x_map, id_map,
            y_map, xy_map, y_map, xy_map, y_map,
            id_map, x_map, id_map, x_map, id_map,
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


# Testing sample
# TODO: добавить более наглядный rev
def get_rev_curve():
    """Curve with time reversal at some middle cube."""
    return FractalCurve(
        proto=[(0, 0), (0, 1), (1, 1), (1, 0)],
        base_maps=[
            BaseMap([1, 0], [False, False]),                # (x,y)->(y,x)
            BaseMap([0, 1], [True, False], time_rev=True),  # (x,y)->(1-x,y), t->1-t
            BaseMap.id_map(2),                              # (x,y)->(x,y)
            BaseMap([1, 0], [True, True]),                  # (x,y)->(1-y,1-x)
        ],
    )

def get_rev2_curve():
    """Curve with time reversal at first cube."""
    return FractalCurve(
        proto=[(0, 0), (0, 1), (1, 1), (1, 0)],
        base_maps=[
            BaseMap([1, 0], [False, True], time_rev=True),  # (x,y)->(y,1-x), t->1-t
            BaseMap.id_map(2),                              # (x,y)->(x,y)
            BaseMap.id_map(2),                              # (x,y)->(x,y)
            BaseMap([1, 0], [True, True]),                  # (x,y)->(1-y,1-x)
        ],
    )

def get_rev3_curve():
    """Curve with time reversal at last cube."""
    return FractalCurve(
        proto=[(0, 0), (0, 1), (1, 1), (1, 0)],
        base_maps=[
            BaseMap([1, 0], [False, False]),                # (x,y)->(y,x)
            BaseMap.id_map(2),                              # (x,y)->(x,y)
            BaseMap.id_map(2),                              # (x,y)->(x,y)
            BaseMap([1, 0], [True, False], time_rev=True),  # (x,y)->(1-y,x), t->1-t
        ],
    )

# TODO: rev4


# Discontinuous curve - beginning in square center
def get_discontinuous_curve():
    return FractalCurve(
        proto=[(1, 1), (0, 1), (0, 0), (1, 0), (2, 0), (2, 1), (2, 2), (1, 2), (0, 2)],
        base_maps=[BaseMap.id_map(2)] * 9,
    )

# Discontinuous curve
def get_morton_curve():
    chain_code = 'i','Ij','i'
    bases = ['ij'] * 4
    return FractalCurve(
        proto=chain2proto(chain_code),
        base_maps=[basis2base_map(b) for b in bases],
    )

'''
# Polyfractal curve (not realized)


# Tetrafractal curve
def get_ARW_Curve(k):

    chain_code = [['i','Ij','i'],'jiJ','jiJ','jiJ'],
    bases = [['3ij0','1jI1','2Ji0','1iJ0'],
             ['3ji0','2Ij1','1ij0','1JI0'],
             ['0ji0','1Ji0','0jI0','1JI0'],
             ['0ij0','2Ji0','0jI0','3Ji1']]
    return FractalCurve(
        proto=chain2proto(chain_code),
        base_maps=[basis2base_map(b) for b in bases],
    )


# Minimal 2D monofractal curve in L_1 (9), L_2 (5), L_inf (5)
def get_beta_Omega_Curve():

    chain_code = ['jiJ','jiJ'],
    bases = [['1iJ0','1Ji0','1ji1','1IJ1'],
             ['1iJ0','1Ji0','1ji1','0jI0']]
    return FractalCurve(
        proto=chain2proto(chain_code),
        base_maps=[basis2base_map(b) for b in bases],
    )


# Minimal 3D bifractal curve in L_inf (9.45)
# Is minimal L_2 (18.3)?
def get_neptunus_curve():
    
    chain_code = ['kjKikJK','kiKjIki'],
    bases = [['0kji','1kji','1KiJ','1jKI','1ikj','1KJi','0kJI','1jKI'],
             ['0jki','1jki','1iKJ','0KiJ','1JiK','1IKj','0ikj','1ijk']]
    return FractalCurve(
        proto=chain2proto(chain_code),
        base_maps=[basis2base_map(b) for b in bases],
    )
 

# Is minimal L_1 (89.8) and L_2 (18.3)?
def get_luna_curve():
    
    chain_code = ['kjKikJK','kiKjIki'],
    bases = [['1ijk','0KJi','1KiJ','1jKI','1jik','1IKj','0kJI','1kJI'],
             ['1jik','0JKi','1iKJ','0KiJ','1KjI','1JIk','0ikj','1ikj']]
    return FractalCurve(
        proto=chain2proto(chain_code),
        base_maps=[basis2base_map(b) for b in bases],
    )
'''
