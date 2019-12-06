# coding: utf-8

from functools import lru_cache
import math
import re

from .base_maps import BaseMap
from .fast_fractions import FastFraction


def ratio_linf(d, dv, dt):
    return FastFraction(max(abs(x) for x in dv)**d, dt)


def ratio_l1(d, dv, dt):
    return FastFraction(sum(abs(x) for x in dv)**d, dt)


def ratio_l2_squared(d, dv, dt):
    return FastFraction(sum(x**2 for x in dv)**d, dt**2)


def ratio_l2(d, dv, dt):
    assert d % 2 == 0
    d2 = d // 2
    return FastFraction(sum(x**2 for x in dv)**d2, dt)


@lru_cache(maxsize=2**20)
def get_int_cube_with_cache(dim, N, cubes):
    """Integer coordinates for sequence of embedded cubes."""
    x = [0] * dim
    Npower = 1
    for cube in reversed(cubes):
        for j in range(dim):
            x[j] += cube[j] * Npower
        Npower *= N
    return x


@lru_cache(maxsize=2**20)
def get_int_time_with_cache(dim, N, cnums):
    """Integer time for sequence of cnums of embedded cubes."""
    G = N**dim
    # multiply by G**l, l = len(cnums), i.e. depth
    # t = c0/G + c1/G**2 + ... = (c_{l-1} + c_{l-2}*G + ..) / G^l
    t = 0
    Gpower = 1
    for cnum in reversed(cnums):
        t += cnum * Gpower
        Gpower *= G
    return t


def get_periodic_sum(start, period, d):
    """
    Sum the non-periodic and periodic parts.

    The sum is:
    s_0/d + s_1/d^2 + ... + s_{k-1}/d^k (non-periodic part = start) +
       + s_k/d^{k+1} + ... + s_{k+m-1}/d^{k+m} + ... (periodic part = period)
    """
    n0 = 0
    d_power = 1
    for x in reversed(start):
        n0 += x * d_power
        d_power *= d
    t0 = FastFraction(n0, d_power)

    np = 0
    dp_power = 1
    for x in reversed(period):
        np += x * dp_power
        dp_power *= d
    tp = FastFraction(np, d_power*(dp_power - 1))

    return t0 + tp



def chain2proto(chain_code):
    """Convert chain code like 'ijK' to curve prototype."""
    dim = len(set(''.join(chain_code).lower()))

    assert dim <= 6
    letters = 'ijklmn'

    vect_dict = {}
    for k in range(dim):
        coord = [0]*dim
        coord[k] = 1
        vect_dict[letters[k]] = coord
        vect_dict[letters[k].upper()] = [-m for m in coord]
        
    def diag_coord(vector):
        arg = [vect_dict[k] for k in vector]
        coord = list(map(sum,zip(*arg)))
        return coord

    proto = [list(map(vect_dict.get,chain_code)) if len(chain_code) == 1 else diag_coord(m) for m in chain_code]
    
    proto = [[0] * dim] + proto
    for l in range(len(proto)-1):
        proto[l+1] = [c + d for c, d in zip(proto[l], proto[l+1])]

    return proto


def bmstr2base_map(bmstr):
    if '1-t' in bmstr:
        time_rev = True
    else:
        time_rev = False
    match = re.match('\(x,y\)->\((.*),(.*)\)', bmstr)
    if not match:
        print('@@@',bmstr)
        
    g1, g2 = match.groups()
    if 'x' in g1:
        perm = [0, 1]
    else:
        perm = [1, 0]
    flip = [False]*2
    flip[0] = ('1-' in g1)
    flip[1] = ('1-' in g2)
    
    return BaseMap(zip(perm, flip), time_rev)


def get_lcm(iterable):
    """Least common multiple of integer sequence."""
    lcm = 1
    for x in iterable:
        lcm = (lcm * x) // math.gcd(lcm, x)
    return lcm
