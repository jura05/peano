#!/usr/bin/env python3
# coding: UTF-8

import logging

import time
from collections import Counter
from heapq import heappop, heappush
from fractions import Fraction

from examples import *
from partial_fractal_curve import CurvePiecePosition, CurvePieceBalancedPair, forget
from fractal_curve import FractalCurve
from base_map import BaseMap, constraint_base_maps

from pysat.solvers import *

# model: [-1, -2, 3, -4, -5, -6, 7, -8, -9]
# MY: [3, 6, 6] [(1, 2), (2, 2), (2, 2)] [5, 2, 2] [(1, 0), (0, 0), (0, 0)] {2: (x,y)->(x,y), 3: (x,y)->(x,1-y), 5: (x,y)->(x,1-y), 6: (x,y)->(x,y)} 30.026835878114877 44.81251280650229
# MY:
# bm[2] id +  VAR3=True
# bm[3] -y -  VAR4=False
# bm[5] -y -  VAR6=False
# bm[6] id +  VAR7=True


def get_curve_from_model(curve, model):
    if curve.dim != 2:
        raise Exception("dim>2 not implemented!")

    base_maps = []
    for cnum, m in enumerate(model):
        allowed_maps = curve.get_allowed_maps(cnum)
        if m > 0:
            allowed_maps = [bm for bm in allowed_maps if bm.is_oriented()]
        else:
            allowed_maps = [bm for bm in allowed_maps if not bm.is_oriented()]

        assert len(allowed_maps) == 1
        base_maps.append(allowed_maps[0])

    return FractalCurve(dim=curve.dim, div=curve.div, proto=curve.proto, base_maps=base_maps)


# cnum++ т.к. в sat начинается с единицы
def get_SAT_clause(curve):
    clause = []
    for cnum, bm in enumerate(curve.base_maps):
        if bm is None:
            continue
        var_num = cnum + 1
        # just for dim=2, diagonal case!

        # Условие = отрицание запрета!
        if bm.is_oriented():
            token = -var_num
        else:
            token = +var_num
        clause.append(token)
    return clause


def my_main():
    def ratio_l2(d, dv, dt):
        return (sum(x**2 for x in dv)**d, dt**2)

    # CONFIG:
    curve = get_peano5_curve()
    THRESHOLD = 40
    rel_tol = 1.1

    #curve = get_peano_curve()
    #THRESHOLD = 20
    #rel_tol = 1.1

    THRESHOLD_LO = THRESHOLD
    THRESHOLD_UP = THRESHOLD * rel_tol
    print('thresholds:', THRESHOLD_LO, THRESHOLD_UP)

    res = curve.estimate_ratio(ratio_l2, rel_tol=0.01, junctions=[None], verbose=1)
    print(res)
    print('ORIGINAL upper bound:', str(res['upper_bound']))

    pcurve = forget(curve)

    stats = Counter()

    active_pairs = []
    iter_count = 0

    bads = []
    all_SAT = []

    def add_pair(pair):
        bms = pair.curve.base_maps
        #logging.debug(
        #    'add_pair %s %s %s',
        #    {cnum: bm for cnum, bm in enumerate(bms) if bm is not None},
        #    pair.pos1.cnums,
        #    pair.pos2.cnums,
        #)

        nonlocal iter_count
        iter_count += 1
        up = pair.upper_bound(ratio_l2)

        if float(up) < THRESHOLD_UP:
            # Not interesed
            logging.debug(' => good')
            stats['good'] += 1
            return

        lo = pair.lower_bound(ratio_l2)
        if float(lo) > THRESHOLD_LO:
            bads.append(pair)
            # found zapret
            logging.debug(' => bad')
            stats['good'] += 1
            stats['bad'] += 1
            return

        logging.debug(' => active')

        priority = -float(up)
        heappush(active_pairs, (priority, iter_count, {'pair': pair, 'up': up, 'lo': lo}))

    for cnum1 in range(pcurve.genus()):
        for cnum2 in range(cnum1 + 2, pcurve.genus()):
            pos1 = pcurve.get_piece_position(cnum1)
            pos2 = pcurve.get_piece_position(cnum2)
            pair = CurvePieceBalancedPair(pcurve, pos1, pos2)
            add_pair(pair)

    N_ITER = 1000000
    for it in range(N_ITER):
        worst_item = heappop(active_pairs)[-1]
        worst_pair = worst_item['pair']

        if it % 1000 == 0 or it == N_ITER-1:
            print({
                'iter': it, 'pairs': len(active_pairs),
                'up': float(worst_item['up']),
                'lo': float(worst_item['lo']),
                'depth': (worst_pair.pos1.depth(), worst_pair.pos2.depth()),
            })
            print(stats)

        for subpair in worst_pair.divide():
            add_pair(subpair)

        if 0 and (it % 10000 == 0 or it == N_ITER-1 or not active_pairs):
            SAT = [get_SAT_clause(pair.curve) for pair in bads]
            SAT = [clause for clause in SAT if clause]
            all_SAT += SAT
            bads = []
            print('all_SAT:', len(all_SAT))
            if not all_SAT:
                continue

            solver = Glucose3()
            solver.append_formula(all_SAT)
            res = solver.solve()
            if not res:
                print('NO SAT MODEL!')
                break
            print('res:', res)
            model = solver.get_model()
            print('model:', model)

            found_curve = get_curve_from_model(pcurve, model)
            found_curve.check()
            res2 = found_curve.estimate_ratio(ratio_l2, rel_tol=0.01, junctions=[None], verbose=1, find_argmax=True)
            print(res2)
            print('FOUND upper bound:', str(res2['upper_bound']))

        if not active_pairs:
            break

def time_to_cnums(genus, t, max_length=10):
    cnums = []
    for _ in range(max_length):
        t = t * genus
        cnum = int(t)
        cnums.append(cnum)
        t = t - cnum
    return cnums

def point_to_cubes(div, x, max_length=10):
    cubes = []
    for _ in range(max_length):
        x = [xj * div for xj in x]
        cube = tuple(int(xj) for xj in x)
        cubes.append(cube)
        x = [xj - cubej for xj, cubej in zip(x, cube)]
    return cubes


if __name__ == "__main__":
    #print(time_to_cnums(9, Fraction(5,12)))
    #print(time_to_cnums(9, Fraction(7,12)))
    #print(point_to_cubes(3, (Fraction(2,3),Fraction(1,1))))
    #print(point_to_cubes(3, (Fraction(1,3),Fraction(0,1))))
    #test0()
    my_main()
