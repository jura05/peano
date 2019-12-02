# coding: utf-8

from .curves import Curve
from .base_maps import BaseMap, Spec
from .utils import chain2proto


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
        BaseMap([(1, True), (0, False)]),  # rot(90)
        BaseMap([(1, False), (0, False)]),  # (x,y)->(y,x)

        BaseMap([(0, False), (1, True)]),  # (x,y)->(x,1-y)
        BaseMap([(1, True), (0, True)]),  # (x,y)->(1-y,1-x)
        BaseMap([(1, False), (0, True)]),  # rot(-90)

        BaseMap.id_map(dim=2),
        BaseMap([(1, True), (0, False)]),  # rot(90)
        BaseMap([(1, False), (0, False)]),  # (x,y)->(y,x)
    ]
    return Curve(dim=2, div=3, patterns=[(proto, base_maps)])


# Minimal 2D monofractal curve in L_1 (9)
def get_hilbert_curve():
    """Example of fractal curve due to D.Hilbert."""
    proto = [(0, 0), (0, 1), (1, 1), (1, 0)]
    base_maps = [
        BaseMap([(1, False), (0, False)]),  # (x,y)->(y,x)
        BaseMap.id_map(2),                # (x,y)->(x,y)
        BaseMap.id_map(2),                # (x,y)->(x,y)
        BaseMap([(1, True), (0, True)]),    # (x,y)->(1-y,1-x)
    ]
    return Curve(dim=2, div=2, patterns=[(proto, base_maps)])


def get_peano_curve():
    """Example of fractal curve due to G.Peano."""
    id_map = BaseMap.id_map(2)
    x_map = BaseMap([(0, True), (1, False)])  # (x,y)->(1-x,y)
    y_map = BaseMap([(0, False), (1, True)])  # (x,y)->(x,1-y)
    xy_map = BaseMap([(0, True), (1, True)])  # (x,y)->(1-x,1-y)
    proto=[(0, 0), (0, 1), (0, 2), (1, 2), (1, 1), (1, 0), (2, 0), (2, 1), (2, 2)]
    base_maps=[
        id_map, x_map, id_map,
        y_map, xy_map, y_map,
        id_map, x_map, id_map,
    ]
    return Curve(dim=2, div=3, patterns=[(proto, base_maps)])

def get_peano5_curve():
    id_map = BaseMap.id_map(2)
    x_map = BaseMap([(0, True), (1, False)])  # (x,y)->(1-x,y)
    y_map = BaseMap([(0, False), (1, True)])  # (x,y)->(x,1-y)
    xy_map = BaseMap([(0, True), (1, True)])  # (x,y)->(1-x,1-y)
    proto = [
        (0, 0), (0, 1), (0, 2), (0, 3), (0, 4),
        (1, 4), (1, 3), (1, 2), (1, 1), (1, 0),
        (2, 0), (2, 1), (2, 2), (2, 3), (2, 4),
        (3, 4), (3, 3), (3, 2), (3, 1), (3, 0),
        (4, 0), (4, 1), (4, 2), (4, 3), (4, 4),
    ]
    base_maps = [
        id_map, x_map, id_map, x_map, id_map,
        y_map, xy_map, y_map, xy_map, y_map,
        id_map, x_map, id_map, x_map, id_map,
        y_map, xy_map, y_map, xy_map, y_map,
        id_map, x_map, id_map, x_map, id_map,
    ]
    return Curve(dim=2, div=5, patterns=[(proto, base_maps)])


# Minimal 2D monofractal curve in L_inf (5.333) and L_2 (5.667)
# equivalent to Scepin-Bauman curve
def get_meurthe_curve():
    chain_code = 'jjiJJijj'
    bases = ['ji','Ji','ij','jI','JI','iJ','ji','Ji','ij']
    return Curve(
        dim=2, div=3,
        patterns=[
            (chain2proto(chain_code), [BaseMap.from_basis(b) for b in bases]),
        ],
    )
    
    
def get_coil_curve():
    chain_code = 'jjiJJijj'
    bases = ['ji','Ji','ji','jI','JI','jI','ji','Ji','ji']
    return Curve(
        dim=2, div=3,
        patterns=[
            (chain2proto(chain_code), [BaseMap.from_basis(b) for b in bases]),
        ],
    )


def get_serpentine_curve():
    chain_code = 'jjiJJijj'
    bases = ['ij','Ji','ji','iJ','JI','iJ','ji','Ji','ij']
    return Curve(
        dim=2, div=3,
        patterns=[
            (chain2proto(chain_code), [BaseMap.from_basis(b) for b in bases]),
        ],
    )


def get_R_curve():
    chain_code = 'jjiiJIJi'
    bases = ['ji','ji','ij','ij','ij','IJ','JI','JI','ij']
    return Curve(
        dim=2, div=3,
        patterns=[
            (chain2proto(chain_code), [BaseMap.from_basis(b) for b in bases]),
        ],
    )


# Minimal 3D monofractal curve with time reversal in L_inf (12.4)
# SUSPICIOUS!!!
def get_haverkort_curve_1():
    """3-D curve with time reversal."""
    chain_code = 'kjKikJK'
    bases = ['kji0','jik0','kIj1','iKJ0','IKJ1','KIj0','Kij1','Jki1']
    return Curve(
        dim=3, div=2,
        patterns=[
            (chain2proto(chain_code), [BaseMap.from_basis(b) for b in bases]),
        ],
    )

# Minimal 3D monofractal curve with time reversal in L_inf (12.4)
# "An inventory of three-dimensional Hilbert space-filling curves", Herman Haverkort; Curve F (see p.13, p.18)
# https://arxiv.org/abs/1109.2323
# Curve A26, see p.10, p.15, p.18
def get_haverkort_curve_A26():
    """3-D curve with time reversal."""
    proto = [(0,0,0), (0,0,1), (0,1,1), (0,1,0), (1,1,0), (1,1,1), (1,0,1), (1,0,0)]
    base_maps = [
        BaseMap([(2, False), (1, False), (0, False)]),
        BaseMap([(2, False), (0, False), (1, False)]),
        BaseMap([(2, False), (0, True), (1, False)], time_rev=True),
        BaseMap([(0, False), (2, True), (1, True)]),
        BaseMap([(0, True), (2, True), (1, True)], time_rev=True),
        BaseMap([(2, True), (0, True), (1, False)]),
        BaseMap([(2, True), (0, False), (1, False)], time_rev=True),
        BaseMap([(1, True), (2, False), (0, False)], time_rev=True),
    ]
    return Curve(dim=3, div=2, patterns=[(proto, base_maps)])


# Minimal 3D monofractal curve with time reversal in L_1 (89.8) and L_2 (18.6)
# "An inventory of three-dimensional Hilbert space-filling curves", Herman Haverkort;
# https://arxiv.org/abs/1109.2323
# Curve F, see p.13, p.15, p.18
def get_haverkort_curve_2():
    """3-D curve with time reversal."""
    chain_code = 'kjKikJK'
    bases = ['KIJ1','KJI1','KjI0','Jki1','jki0','kjI1','kJI0','iKJ0']
    return Curve(
        dim=3, div=2,
        patterns=[
            (chain2proto(chain_code), [BaseMap.from_basis(b) for b in bases]),
        ],
    )


def get_tokarev_curve():
    """3-D curve."""
    chain_code = 'kjKikJK'
    bases = ['jki', 'kij', 'kij', 'iJK', 'iJK', 'KIj', 'KIj', 'JkI']
    return Curve(
        dim=3, div=2,
        patterns=[
            (chain2proto(chain_code), [BaseMap.from_basis(b) for b in bases]),
        ],
    )


# Testing sample
# TODO: добавить более наглядный rev
def get_rev_curve():
    """Curve with time reversal at some middle cube."""
    return Curve(
        dim=2, div=2,
        patterns=[
            (
                [(0, 0), (0, 1), (1, 1), (1, 0)],
                [
                    BaseMap([(1, False), (0, False)]),                # (x,y)->(y,x)
                    BaseMap([(0, True), (1, False)], time_rev=True),  # (x,y)->(1-x,y), t->1-t
                    BaseMap.id_map(2),                              # (x,y)->(x,y)
                    BaseMap([(1, True), (0, True)]),                  # (x,y)->(1-y,1-x)
                ],
            ),
        ],
    )

def get_rev2_curve():
    """Curve with time reversal at first cube."""
    proto = [(0, 0), (0, 1), (1, 1), (1, 0)]
    base_maps=[
        BaseMap([(1, False), (0, True)], time_rev=True),  # (x,y)->(y,1-x), t->1-t
        BaseMap.id_map(2),                              # (x,y)->(x,y)
        BaseMap.id_map(2),                              # (x,y)->(x,y)
        BaseMap([(1, True), (0, True)]),                  # (x,y)->(1-y,1-x)
    ]
    return Curve(
        dim=2, div=2,
        patterns=[(proto, base_maps)],
    )

def get_rev3_curve():
    """Curve with time reversal at last cube."""
    proto=[(0, 0), (0, 1), (1, 1), (1, 0)]
    base_maps=[
        BaseMap([(1, False), (0, False)]),                # (x,y)->(y,x)
        BaseMap.id_map(2),                              # (x,y)->(x,y)
        BaseMap.id_map(2),                              # (x,y)->(x,y)
        BaseMap([(1, True), (0, False)], time_rev=True),  # (x,y)->(1-y,x), t->1-t
    ]
    return Curve(dim=2, div=2, patterns=[(proto, base_maps)])


# Discontinuous curve - beginning in square center
def get_discontinuous_curve():
    proto = [(1, 1), (0, 1), (0, 0), (1, 0), (2, 0), (2, 1), (2, 2), (1, 2), (0, 2)]
    return Curve(
        dim=2, div=2,
        patterns=[
            (proto, [BaseMap.id_map(2)] * 9),
        ],
    )

# Discontinuous curve
def get_morton_curve():
    chain_code = 'i','Ij','i'
    bases = ['ij'] * 4
    return Curve(
        dim=2, div=2,
        patterns=[(chain2proto(chain_code), [BaseMap.from_basis(b) for b in bases])],
    )

#
# Polyfractal curves
#

def get_hilbert_bicurve():
    """Fake bifractal curve, actually it is Hilber curve."""
    curve = get_hilbert_curve()
    p1 = (curve.proto, [Spec(sp.base_map, pnum=cnt % 2) for cnt, sp in enumerate(curve.specs)])
    p2 = (curve.proto, [Spec(sp.base_map, pnum=(cnt + 1) % 2) for cnt, sp in enumerate(curve.specs)])
    return Curve(dim=2, div=2, patterns=[p1,p2])


def get_patterns(chain_code_list, bases_list):
    patterns = []
    for chain_code, rich_bases in zip(chain_code_list, bases_list):
        specs = []
        for rich_basis in rich_bases:
            pnum_str, basis = rich_basis[0], rich_basis[1:]
            spec = Spec(base_map=BaseMap.from_basis(basis), pnum=int(pnum_str))
            specs.append(spec)
        proto = chain2proto(chain_code)
        patterns.append((proto, specs))
    return patterns

# Tetrafractal curve
def get_ARW_Curve():
    chain_code = [['i','Ij','i'],'jiJ','jiJ','jiJ'],
    bases = [['3ij0','1jI1','2Ji0','1iJ0'],
             ['3ji0','2Ij1','1ij0','1JI0'],
             ['0ji0','1Ji0','0jI0','1JI0'],
             ['0ij0','2Ji0','0jI0','3Ji1']]
    return Curve(dim=2, div=2, patterns=get_patterns(chain_code, bases))


# Minimal 2D monofractal curve in L_1 (9), L_2 (5), L_inf (5)
def get_beta_Omega_Curve():
    chain_code_list = ['jiJ','jiJ']
    bases_list = [['1iJ0','1Ji0','1ji1','1IJ1'],
             ['1iJ0','1Ji0','1ji1','0jI0']]
    return Curve(dim=2, div=2, patterns=get_patterns(chain_code_list, bases_list))


# Minimal 3D bifractal curve in L_inf (9.45)
# Is minimal L_2 (18.3)?
def get_neptunus_curve():
    chain_code = ['kjKikJK','kiKjIki']
    bases = [['0kji','1kji','1KiJ','1jKI','1ikj','1KJi','0kJI','1jKI'],
             ['0jki','1jki','1iKJ','0KiJ','1JiK','1IKj','0ikj','1ijk']]
    return Curve(
        dim=3, div=2,
        patterns=get_patterns(chain_code, bases),
    )
 

# Is minimal L_1 (89.8) and L_2 (18.3)?
def get_luna_curve():
    chain_code = ['kjKikJK','kiKjIki']
    bases = [['1ijk','0KJi','1KiJ','1jKI','1jik','1IKj','0kJI','1kJI'],
             ['1jik','0JKi','1iKJ','0KiJ','1KjI','1JIk','0ikj','1ikj']]
    return Curve(
        dim=3, div=2,
        patterns=get_patterns(chain_code, bases),
    )
