#!/usr/bin/env python3
# coding: UTF-8

import pickle
import sys
import logging
import itertools
import time
from collections import Counter
from fractions import Fraction

from utils import bmstr2base_map
from examples import *
from partial_fractal_curves import PartialFractalCurve,  get_int_cube_with_cache
from fractal_curves import FractalCurve
from base_maps import BaseMap, gen_constraint_cube_maps
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
    proto = [brk[0] for brk in brkline]

    cube_first, entr_first, _ = brkline[0]
    entr = tuple(Fraction(cj + ej, div) for cj, ej in zip(cube_first, entr_first))

    cube_last, _, exit_last = brkline[-1]
    exit = tuple(Fraction(cj + ej, div) for cj, ej in zip(cube_last, exit_last))

    symmetries = []
    for bm in gen_constraint_cube_maps(dim, {entr: entr, exit: exit}):
        symmetries.append(bm)
    if allow_time_rev:
        for bm in gen_constraint_cube_maps(dim, {entr: exit, exit: entr}):
            symmetries.append(bm.reverse_time())

    base_maps = [None] * len(proto)
    repr_maps = [None] * len(proto)

    for cnum, brk in enumerate(brkline):
        cube, rel_entr, rel_exit = brk
        for bm in gen_constraint_cube_maps(dim, {entr: rel_entr, exit: rel_exit}):
            repr_maps[cnum] = bm
            break

    return PartialFractalCurve(
        dim=dim, div=div,
        proto=proto,
        base_maps=base_maps,
        repr_maps=repr_maps,
        symmetries=symmetries,
    )

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
            verbose=True,
        )
        if has_good:
            stats['found'] += 1
        else:
            stats['not_found'] += 1
        if conf.get('max_iter') and it >= conf['max_iter']: break
    print(stats)

def clever_perebor(conf):
    curve_gen = CurveGenerator(dim=conf['dim'], div=conf['div'], hdist=conf['hdist'], max_cdist=conf['max_cdist'], verbose=1)
    brklines = list(curve_gen.generate_brklines())
    best_up = 40  # АПРИОРИ
    ratio_func = ratio_l2

    print('CONF:', conf)
    stats = Counter()
    for it, brkline in enumerate(brklines):
        print('BRKLINE: {} of {} ({:.2f} %)'.format(it + 1, len(brklines), 100 * (it + 1) / len(brklines)), flush=True)
        pcurve = get_pcurve_for_brkline(conf['dim'], conf['div'], brkline, allow_time_rev=conf['allow_time_rev'])

        if best_up is not None:
            # проверка, что тут есть шансы
            has_good = pcurve.estimate_ratio(
                ratio_func,
                lower_bound=best_up * (1 - conf['rel_tol']),
                upper_bound=best_up,
                sat_pack=conf['sat_pack'],
                find_model=False,
                verbose=True,
            )
            if not has_good:
                print('no chances')
                continue

        # начинается поиск!
        curr_lo = 0
        curr_up = best_up
        # инвариант: лучшая кривая в данном классе лежит в [curr_lo, curr_up]

        curr_threshold = curr_up * 0.6 + curr_lo * 0.4
        # делаем bisect для поиска хорошей кривой
        while curr_up * (1 - 2 * conf['rel_tol']) > curr_lo:
            print('best in ', curr_lo, curr_up, 'seek with threshold:', curr_threshold)
            has_good = pcurve.estimate_ratio(
                ratio_func,
                lower_bound=curr_threshold * (1 - conf['rel_tol']),
                upper_bound=curr_threshold,
                sat_pack=conf['sat_pack'],
                find_model=False,
                verbose=True,
            )
            if has_good:
                curr_up = curr_threshold
                print('NEW BEST UP:', curr_up)
                best_up = curr_up
            else:
                curr_lo = curr_threshold * (1 - conf['rel_tol'])

            curr_threshold = curr_up * 0.8 + curr_lo * 0.2

    print('BEST UP:', best_up)

def triple_perebor(conf):
    curve_gen = CurveGenerator(dim=conf['dim'], div=conf['div'], hdist=conf['hdist'], max_cdist=conf['max_cdist'], verbose=1)
    brklines = list(curve_gen.generate_brklines())
    best_up = 36
    apriori_up = 80
    ratio_func = ratio_l2

    print('CONF:', conf)
    stats = Counter()
    for it, brkline in enumerate(brklines):
        print('BRKLINE: {} of {} ({:.2f} %)'.format(it + 1, len(brklines), 100 * (it + 1) / len(brklines)), flush=True)
        pcurve = get_pcurve_for_brkline(conf['dim'], conf['div'], brkline, allow_time_rev=conf['allow_time_rev'])

        if it + 1 != 49:
            continue

        # начинается поиск!
        curr_lo = 0
        curr_up = apriori_up
        # инвариант: лучшая кривая в данном классе лежит в [curr_lo, curr_up]

        # делаем bisect для поиска хорошей кривой
        while curr_up * (1 - conf['rel_tol']) > curr_lo:
            new_lo = 2/3 * curr_lo + 1/3 * curr_up
            new_up = 1/3 * curr_lo + 2/3 * curr_up
            print('best in ', curr_lo, curr_up, 'seek with thresholds:', new_lo, new_up)
            has_good = pcurve.estimate_ratio(
                ratio_func,
                lower_bound=new_lo,
                upper_bound=new_up,
                sat_pack=conf['sat_pack'],
                find_model=False,
                verbose=True,
            )
            if has_good:
                curr_up = new_up
            else:
                curr_lo = new_lo

            if curr_up < best_up:
                print('NEW BEST UP:', curr_up)
                best_up = curr_up

            if curr_lo > best_up:
                print('NO CHANCES!')
                break

        print('ratio in:', curr_lo, curr_up)
        data = pcurve.estimate_ratio(
            ratio_func,
            lower_bound=best_up * 1.0001,
            upper_bound=best_up * 1.0002,
            sat_pack=conf['sat_pack'],
            find_model=True,
            verbose=True,
        )
        curve = data['curve']
        pickle.dump(curve, open('best_curve.pickle', 'wb'))
        res = curve.estimate_ratio_new(ratio_l2, rel_tol=0.000001)
        print(res)
        print(curve.proto)
        print(curve.base_maps)

    print('BEST UP:', best_up)


def pristalno2():
    curve_gen = CurveGenerator(dim=2, div=5, hdist=1, max_cdist=1, verbose=1)
    brklines = list(curve_gen.generate_brklines())
    brkline = brklines[48]
    for cnum, data in enumerate(brkline):
        print('cnum:', cnum, 'cube:', data[0], 'entr:', data[1], 'exit:', data[2])
    pcurve = get_pcurve_for_brkline(2, 5, brkline, allow_time_rev=True)

    res1 = pcurve.estimate_ratio(
        ratio_l2,
        lower_bound=5.55556**2,
        upper_bound=5.55557**2,
        sat_pack=100,
        find_model=True,
        verbose=True,
    )
    res2 = res1['curve'].estimate_ratio_new(ratio_l2, rel_tol=0.00001)
    allres = {
        "model": res1["model"],
        "curve": res1["curve"],
        'log1': res1['pairs_tree'].LOG,
        'log2': res2['pairs_tree'].LOG,
    }
    #pickle.dump(allres, open("50_9.pickle", "wb"))



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


def test_perebor():
    conf = {
        'dim': 2,
        'div': 5,
        'hdist': 2,
        'max_cdist': 1,
        'lower_bound': 40,
        'upper_bound': 60,
        'sat_pack': 10,
        'allow_time_rev': True,
        'find_model': True,
        'max_iter': 50,
    }
    perebor(conf)

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
# 40s
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
        'verbose': True,
    }
    perebor(conf)
    print(get_int_cube_with_cache.cache_info())

# Counter({'not_found': 70, 'found': 14})
# 25s -> 14s
def timings3():
    conf = {
        'dim': 3,
        'div': 2,
        'hdist': 1,
        'max_cdist': 2,
        'lower_bound': 700,
        'upper_bound': 750,
        'sat_pack': 50,
        'max_iter': 100,
        'allow_time_rev': False,
        'find_model': False,
    }
    perebor(conf)

# Counter({'not_found': 84, 'found': 2})
# real    0m42.613s
# -> 32s
def timings4():
    conf = {
        'dim': 2,
        'div': 5,
        'hdist': 1,
        'max_cdist': 1,
        'lower_bound': 31.5,
        'upper_bound': 32,
        'sat_pack': 10,
        'max_iter': 10000,
        'allow_time_rev': True,
        'find_model': True,
    }
    perebor(conf)

if __name__ == "__main__":
    conf = {
        'dim': 2,
        'div': 5,
        'hdist': 1,
        'max_cdist': 1,
        'rel_tol': 0.00001,
        'sat_pack': 100,
        'allow_time_rev': True,
    }
    curve = get_bauman_curve()
    print(curve.estimate_ratio_new(ratio_l2, rel_tol=0.0001))

    #test_perebor()
    #triple_perebor(conf),
    #timings4()
    #perebor(conf)
    #bauman()
    #pristalno2()
