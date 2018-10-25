#!/usr/bin/env python3

import sys
from math import sqrt
import logging

import gen_curve

logging.basicConfig(level=0, stream=sys.stdout, format='[%(process)d] %(asctime)s %(message)s)')

def ratio_linf(d, dv, dt):
    return (max(abs(x) for x in dv)**d, dt)

def ratio_l1(d, dv, dt):
    return (sum(abs(x) for x in dv)**d, dt)

def ratio_l2(d, dv, dt):
    return (sum(x**2 for x in dv)**d, dt**2)

def main():
    transit = True
    oriented = False
    div = 5
    curve_type = 'edge'
    setups = [
        {'div': div, 'type': curve_type, 'transit': transit, 'oriented': oriented, 'ratio_func': ratio_l1},
        {'div': div, 'type': curve_type, 'transit': transit, 'oriented': oriented, 'ratio_func': ratio_l2},
        {'div': div, 'type': curve_type, 'transit': transit, 'oriented': oriented, 'ratio_func': ratio_linf},
    ]
    for setup in setups:
        go(setup)

def go(setup):
    total_iter = 0
    seen_curves = 0
    exit = (1,0) if setup['type'] == 'edge' else (1,1)
    generator = gen_curve.CurveGenerator(
        div=setup['div'], exit=exit, allow_vertex_transit=setup.get('transit'), oriented=setup.get('oriented'),
    )
    best_upper_bound = None
    curves = list(generator.generate_curves())
    for curve in curves:
        seen_curves += 1
        logging.info('process curve %d of %d, current upper_bound: %s', seen_curves, len(curves), best_upper_bound)
        res = curve.estimate_ratio(setup['ratio_func'], rel_tol=0.01, upper_bound=best_upper_bound, verbose=0)
        total_iter += res['stats'].get('iter', 0)
        if best_upper_bound is None or res['upper_bound'] < best_upper_bound:
            best_upper_bound = res['upper_bound']
            logging.info('new upper bound found: %f', best_upper_bound)

    real_upper_bound = float(best_upper_bound)
    if setup['ratio_func'] == ratio_l2:
        real_upper_bound = sqrt(real_upper_bound)

    print('BEST for setup', setup)
    print('ratio:', real_upper_bound, 'total_iter:', total_iter)

if __name__ == "__main__":
    main()
