#!/usr/bin/env python3
# coding: UTF-8

import sys
import logging
import itertools
import time
from collections import Counter
from fractions import Fraction

from utils import bmstr2base_map
from examples import *
from partial_fractal_curve import CurvePiecePosition, CurvePieceBalancedPair, PartialFractalCurve,  get_int_cube_with_cache
from fractal_curve import FractalCurve
from base_map import BaseMap
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



def get_pcurve_for_brkline(dim, div, brkline, allow_time_rev):
    proto = []
    gates = []
    for cube, entrance, exit in brkline:
        proto.append(cube)
        gates.append((entrance, exit))
    base_maps = [None] * len(proto)
    return PartialFractalCurve(dim=dim, div=div, proto=proto, base_maps=base_maps, gates=gates, allow_time_rev=allow_time_rev)

# l40: found: 2
def perebor(conf):
    curve_gen = CurveGenerator(dim=conf['dim'], div=conf['div'], hdist=conf['hdist'], max_cdist=conf['max_cdist'], verbose=1)
    brklines = list(curve_gen.generate_brklines())

    print('CONF:', conf)
    print('BRKLINES:', len(brklines))
    stats = Counter()
    for it, brkline in enumerate(brklines):
        print('ITER:', it, stats)
        pcurve = get_pcurve_for_brkline(conf['dim'], conf['div'], brkline, allow_time_rev=conf['allow_time_rev'])
        has_good = pcurve.estimate_ratio(
            ratio_l2,
            lower_bound=conf['lower_bound'],
            upper_bound=conf['upper_bound'],
            sat_pack=conf['sat_pack'],
            find_model=conf.get('find_model'),
        )
        if has_good:
            stats['found'] += 1
        else:
            stats['not_found'] += 1
        if conf.get('max_iter') and it >= conf['max_iter']: break
    print(stats)

def pristalno():
    proto = [
        (0, 0), (1, 0), (1, 1), (0, 1), (0, 2),
        (1, 2), (1, 3), (0, 3), (0, 4), (1, 4),
        (2, 4), (2, 3), (3, 3), (3, 4), (4, 4),
        (4, 3), (4, 2), (3, 2), (2, 2), (2, 1),
        (2, 0), (3, 0), (3, 1), (4, 1), (4, 0),
    ]
    bmstrs = [
        '(x,y)->(x,y)', '(x,y)->(y,1-x),t->1-t', '(x,y)->(y,1-x),t->1-t', '(x,y)->(1-x,1-y)', '(x,y)->(x,y)',
        '(x,y)->(y,1-x),t->1-t', '(x,y)->(y,1-x),t->1-t', '(x,y)->(1-x,1-y)', '(x,y)->(x,y)', '(x,y)->(x,y)',
        '(x,y)->(x,y)', '(x,y)->(1-y,x),t->1-t', '(x,y)->(y,1-x),t->1-t', '(x,y)->(x,y)', '(x,y)->(x,y)',
        '(x,y)->(1-y,x),t->1-t', '(x,y)->(1-x,1-y)', '(x,y)->(1-x,1-y)', '(x,y)->(1-y,x),t->1-t', '(x,y)->(1-y,x),t->1-t',
        '(x,y)->(1-y,x),t->1-t', '(x,y)->(y,1-x),t->1-t', '(x,y)->(x,y)', '(x,y)->(x,y)', '(x,y)->(1-y,x),t->1-t',
    ]
    base_maps = [bmstr2base_map(x) for x in bmstrs]

    curve = FractalCurve(dim=2, div=5, proto=proto, base_maps=base_maps)
    curve.check()
    print(curve.get_junctions())
    print('entr:', curve.get_entrance())
    print('exit:', curve.get_exit())

    curve.estimate_ratio_new(ratio_l2, rel_tol=0.001, max_iter=1000)

    pcurve = curve.forget()
    pcurve.allow_time_rev = True
    result = pcurve.estimate_ratio(ratio_l2, lower_bound=31.2, upper_bound=31.3, find_model=True)
    if result:
        print('found:', result.proto, result.base_maps)
    else:
        print('NOT FOUND')


def bauman():
    proto = ((0, 0), (0, 1), (0, 2), (1, 2), (1, 1), (1, 0), (2, 0), (2, 1), (2, 2))
    base_maps = [
        BaseMap(perm=[1,0],flip=[False,False]),
        BaseMap(perm=[1,0],flip=[True,False]),
        BaseMap.id_map(dim=2),
        BaseMap(perm=[1,0],flip=[False,True]),
        BaseMap(perm=[1,0],flip=[True,True]),
        BaseMap(perm=[0,1],flip=[False,True]),
        BaseMap(perm=[1,0],flip=[False,False]),
        BaseMap(perm=[1,0],flip=[True,False]),
        BaseMap.id_map(dim=2),
    ]
    bauman = FractalCurve(dim=2, div=3, proto=proto, base_maps=base_maps)
    bauman.estimate_ratio_vertex_brkline(ratio_l2, 3)

    #good.estimate_ratio(ratio, rel_tol=0.002, verbose=1)
    bauman = bauman.get_subdivision(1)
    pcurve = bauman.forget()
    pcurve = pcurve.changed(allow_time_rev=True)
    pcurve.estimate_ratio(ratio_l2, lower_bound=32.2, upper_bound=32.1, sat_pack=1000, find_model=True)

# Counter({'not_found': 11, 'found': 10})
# time=16s
def timings1():
    conf = {
        'dim': 2,
        'div': 5,
        'hdist': 2,
        'max_cdist': 2,
        'lower_bound': 80,
        'upper_bound': 80.8,
        'sat_pack': 10,
        'max_iter': 20,
        'allow_time_rev': False,
        'find_model': False,
    }
    perebor(conf)

# Counter({'not_found': 76, 'found': 25})
# 53s
def timings2():
    conf = {
        'dim': 2,
        'div': 5,
        'hdist': 2,
        'max_cdist': 2,
        'lower_bound': 80,
        'upper_bound': 80.8,
        'sat_pack': 10,
        'max_iter': 100,
        'allow_time_rev': False,
        'find_model': False,
    }
    perebor(conf)
    print(get_int_cube_with_cache.cache_info())

if __name__ == "__main__":
    conf = {
        'dim': 2,
        'div': 5,
        'hdist': 1,
        'max_cdist': 1,
        'lower_bound': 32.05,
        'upper_bound': 32.1,
        'sat_pack': 100,
        #'max_iter': 20,
        'allow_time_rev': True,
    }
    timings2()
    #perebor(conf)
    #bauman()
    #pristalno()
