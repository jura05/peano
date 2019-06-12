#!/usr/bin/env python3
# coding: UTF-8

import sys
import logging
import time
from collections import Counter
from heapq import heappop, heappush
from fractions import Fraction

from examples import *
from partial_fractal_curve import CurvePiecePosition, CurvePieceBalancedPair, forget, PartialFractalCurve,  get_int_cube_with_cache
from fractal_curve import FractalCurve
from base_map import BaseMap, constraint_base_maps, list_base_maps
from curve_sat_adapter import CurveSATAdapter
from gen_curve import CurveGenerator

#logging.basicConfig(level=logging.DEBUG)

global_timer = {}

def ratio_l1(d, dv, dt):
    return (sum(abs(x) for x in dv)**d, dt)

def ratio_linf(d, dv, dt):
    return (max(abs(x) for x in dv)**d, dt)

def ratio_l2(d, dv, dt):
    return (sum(x**2 for x in dv)**d, dt**2)

ratio = ratio_l2


def my_main(pcurve, threshold, rel_tol, sat_pack=10):

    # CONFIG:

    adapter = CurveSATAdapter(dim=pcurve.dim)

    threshold_up = threshold * (1 + rel_tol)
    threshold_lo = threshold
    print('thresholds:', threshold_lo, threshold_up)

    for junc, curves in pcurve.get_junctions_info().items():
        adapter.make_junc_var(junc, curves)

    stats = Counter()

    active_pairs = []
    iter_count = 0

    bads = []
    all_SAT = []

    def add_pair(pair):
        nonlocal iter_count

        iter_count += 1
        up = pair.upper_bound(ratio)

        if float(up) < threshold_up:
            # Not interesed
            stats['good'] += 1
            return

        lo = pair.lower_bound(ratio)
        if float(lo) > threshold_lo:
            bads.append(pair)
            # found zapret
            stats['good'] += 1
            stats['bad'] += 1
            return

        priority = -float(up)
        heappush(active_pairs, (priority, iter_count, {'pair': pair, 'up': up, 'lo': lo}))

    G = pcurve.genus()
    for cnum1 in range(G):
        for cnum2 in range(cnum1 + 2, G):
            pos1 = pcurve.get_piece_position(cnum1)
            pos2 = pcurve.get_piece_position(cnum2)
            pair = CurvePieceBalancedPair(pcurve, None, pos1, pos2)
            add_pair(pair)

    for junc in pcurve.get_junctions_info().keys():
        for cnum1 in range(G):
            for cnum2 in range(G):
                if (cnum1, cnum2) == (G - 1, 0):
                    continue
                pos1 = pcurve.get_piece_position(cnum1)
                pos2 = pcurve.get_piece_position(cnum2)
                pair = CurvePieceBalancedPair(pcurve, junc, pos1, pos2)
                add_pair(pair)

    no_model = None
    N_ITER = 1000000
    it = 0
    while active_pairs and it < N_ITER:
        it += 1
        worst_item = heappop(active_pairs)[-1]
        worst_pair = worst_item['pair']

        if it % 10 == 0 or it == N_ITER-1:
            print({
                'iter': it, 'pairs': len(active_pairs),
                'up': float(worst_item['up']),
                'lo': float(worst_item['lo']),
                'depth': (worst_pair.pos1.depth, worst_pair.pos2.depth),
            })
            print(stats)

        for subpair in worst_pair.divide():
            add_pair(subpair)

        if it % sat_pack == 0 and bads:
            start_time = time.time()
            print('it:', it, 'bads:', len(bads))
            for bad_pair in bads:
                adapter.add_forbid_clause(bad_pair.junc, bad_pair.curve)
            bads = []
            res = adapter.solve()

            elapsed = time.time() - start_time
            global_timer['SAT'] = global_timer.get('SAT', 0) + elapsed

            if res:
                print('SAT model exists!')
            else:
                print('NO SAT MODEL!')
                no_model = True
                break

    if no_model:
        return None

    for bad_pair in bads:
        adapter.add_forbid_clause(bad_pair.junc, bad_pair.curve)

    if not adapter.solve():
        print('NO SAT MODEL!')
        return


    print('ADAPTER STATS:', adapter.stats())
    #return True

    # это если попросят модель
    model = adapter.get_model()
    found_curves = adapter.get_curves_from_model(pcurve, model)
    for found_curve in found_curves:
        print('found curve:', found_curve.proto, found_curve.base_maps)
        res2 = found_curve.estimate_ratio(ratio, rel_tol=0.002, verbose=1)
        print(res2)
        print('FOUND upper bound:', str(res2['upper_bound']))

    return True


def get_pcurve_for_brkline(dim, div, brkline):
    proto = []
    gates = []
    for cube, entrance, exit in brkline:
        proto.append(cube)
        gates.append((entrance, exit))
    base_maps = [None] * len(proto)
    return PartialFractalCurve(dim=dim, div=div, proto=proto, base_maps=base_maps, gates=gates)

# time=16s {'threshold': 80,  'rel_tol': 0.01, 'sat_pack': 10, 'max_iter': 20},; Counter({'not_found': 11, 'found': 10})
# l40: found: 2
def perebor():
    dim = 2
    div = 7
    configs = [
        #{'threshold': 60, 'rel_tol': 0.1, 'sat_pack': 50},
        #{'threshold': 40, 'rel_tol': 0.01, 'sat_pack': 10, 'max_iter': None},
        {'threshold': 40, 'rel_tol': 0.01, 'sat_pack': 10, 'max_iter': None},
    ]

    curve_gen = CurveGenerator(dim=dim, div=div, hdist=2, verbose=1, max_cdist=1)
    brklines = list(curve_gen.generate_brklines())

    for conf in configs:
        print('CONF:', conf)
        print('BRKLINES:', len(brklines))
        good_brklines = []
        stats = Counter()
        for it, brkline in enumerate(brklines):
            print('ITER:', it, stats)
            pcurve = get_pcurve_for_brkline(dim, div, brkline)
            has_good = my_main(pcurve, conf['threshold'], conf['rel_tol'], conf['sat_pack'])
            if has_good:
                stats['found'] += 1
                good_brklines.append(brkline)
            else:
                stats['not_found'] += 1
            if conf.get('max_iter') and it >= conf['max_iter']: break
        print(stats)
        brklines = good_brklines

    print('BRKLINES:', len(brklines))
    #print('get_int_cube_with_cache:', get_int_cube_with_cache.cache_info())
    #print('get_int_time_with_cache:', get_int_time_with_cache.cache_info())

def bauman():
    proto = ((0, 0), (0, 1), (0, 2), (1, 2), (1, 1), (1, 0), (2, 0), (2, 1), (2, 2))
    base_maps = [
        BaseMap(perm=[1,0],flip=[False,False]),
        BaseMap(perm=[1,0],flip=[True,False]),
        BaseMap(dim=2),
        BaseMap(perm=[1,0],flip=[False,True]),
        BaseMap(perm=[1,0],flip=[True,True]),
        BaseMap(perm=[0,1],flip=[False,True]),
        BaseMap(perm=[1,0],flip=[False,False]),
        BaseMap(perm=[1,0],flip=[True,False]),
        BaseMap(dim=2),
    ]
    good = FractalCurve(dim=2, div=3, proto=proto, base_maps=base_maps)
    #good.estimate_ratio(ratio, rel_tol=0.002, verbose=1)
    good = good.get_subdivision(2)
    #good.estimate_ratio(ratio, rel_tol=0.002, verbose=1)
    my_main(forget(good), 32.1, 0.01, 200)

if __name__ == "__main__":
    perebor()
    #bauman()
