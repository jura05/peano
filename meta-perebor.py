#!/usr/bin/env python3

import logging
import argparse

import peano.utils as utils
from peano.paths import PathsGenerator
from peano.curves import gen_possible_gates, PathFuzzyCurve
from peano.ratio import Estimator


def meta_perebor(dim, div, pattern_count, ratio_func, rel_tol_inv, rel_tol_inv_mult):
    estimator = Estimator(ratio_func, cache_max_size=2**16)
    for gates in gen_possible_gates(dim=dim, div=div, pattern_count=pattern_count):
        paths_gen = PathsGenerator(dim=dim, div=div, gates=gates)
        paths_list = list(paths_gen.generate_paths(uniq=True))
        logging.warning('gates: %s', [str(g) for g in gates])
        logging.warning('paths: %d', len(paths_list))
        pcurves = (PathFuzzyCurve.init_from_paths(paths) for paths in paths_list)
        res = estimator.estimate_ratio_sequence(
            pcurves,
            rel_tol_inv=rel_tol_inv,
            rel_tol_inv_mult=rel_tol_inv_mult,
            sat_strategy={'type': 'geometric', 'multiplier': 1.3},
        )
        print('Result for gates:', [str(g) for g in gates])
        print('paths:', len(paths_list))
        print('lower bound:', float(res['lo']))
        print('upper bound:', float(res['up']))


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument('--dim', type=int, required=True)
    argparser.add_argument('--div', type=int, required=True)
    argparser.add_argument('--pattern-count', type=int, required=True)
    argparser.add_argument('--metric', type=str, required=True, choices=['l1','l2','l2_squared','linf'])
    argparser.add_argument('--rel_tol_inv', type=int, default=1000)
    argparser.add_argument('--rel_tol_inv_mult', type=int, default=5)
    args = argparser.parse_args()
    funcs = {
        'l1': utils.ratio_l1,
        'l2': utils.ratio_l2,
        'l2_squared': utils.ratio_l2_squared,
        'linf': utils.ratio_linf,
    }
    meta_perebor(
        dim=args.dim, div=args.div, pattern_count=args.pattern_count, ratio_func=funcs[args.metric],
        rel_tol_inv=args.rel_tol_inv, rel_tol_inv_mult=args.rel_tol_inv_mult,
    )
