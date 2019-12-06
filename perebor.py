#!/usr/bin/env python3
# coding: UTF-8

import argparse
import json
import pickle
import logging

from peano import utils
from peano.paths import PathsGenerator, gen_uniq
from peano.ratio import Estimator
from peano.curves import SymmFuzzyCurve


def perebor(conf):
    logging.warning('CONF: %s', conf)
    funcs = {
        'l1': utils.ratio_l1,
        'l2': utils.ratio_l2,
        'l2_squared': utils.ratio_l2_squared,
        'linf': utils.ratio_linf,
    }
    ratio_func = funcs[conf['ratio_func_name']]
    dim, div = conf['dim'], conf['div']
    paths_gen = PathsGenerator(dim=dim, div=div, hdist=conf['hdist'], max_cdist=conf['max_cdist'])
    paths = list(paths_gen.generate_paths())
    logging.warning('paths: %d', len(paths))
    uniq_paths = list(gen_uniq(dim, paths))
    logging.warning('uniq paths: %d', len(uniq_paths))

    pcurves = (SymmFuzzyCurve.init_from_path(dim, div, path, allow_time_rev=True) for path in uniq_paths)
    estimator = Estimator(ratio_func)
    res = estimator.estimate_ratio_sequence(
        pcurves,
        rel_tol_inv=conf['rel_tol_inv'],
        rel_tol_inv_mult=conf.get('rel_tol_inv_mult', 2),
        sat_strategy={'type': 'geometric', 'multiplier': 1.3},
    )
    logging.warning('FINAL BOUNDS: %.6f <= r <= %.6f', res['lo'], res['up'])
    res['uniq_paths'] = len(uniq_paths)
    print(res)


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument('config')
    args = argparser.parse_args()

    with open(args.config) as fh:
        config = json.load(fh)

    perebor(config)
