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


#logging.basicConfig(level=logging.DEBUG)


def perebor(conf):
    #logging.warning('CONF: %s', conf)
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

    estimator = Estimator(ratio_func, cache_max_size=2**16)

    logging.warning('sort by paths ratio ...')
    path_ratio = {path: estimator.estimate_path(path) for path in uniq_paths}
    uniq_paths.sort(key=lambda path: path_ratio[path])

    pcurves = (SymmFuzzyCurve.init_from_path(dim, div, path, allow_time_rev=True) for path in uniq_paths)
    res = estimator.estimate_ratio_sequence(
        pcurves,
        rel_tol_inv=conf['rel_tol_inv'],
        rel_tol_inv_mult=conf.get('rel_tol_inv_mult', 2),
        sat_strategy={'type': 'geometric', 'multiplier': 1.3},
    )
    logging.warning('CONF: %s', conf)
    logging.warning('FINAL BOUNDS: %.6f <= r <= %.6f', res['lo'], res['up'])
    res['uniq_paths'] = len(uniq_paths)
    print(res)


def perebor_basic():
    configs = []
    for div in range(4, 7):
        hdists = [1] if div % 2 == 0 else [1, 2]
        if div <= 3:
            max_cdists = [2]
        elif div == 4:
            max_cdists = [1, 2]
        else:
            max_cdists = [1]
        for hdist in hdists:
            for max_cdist in max_cdists:
                for name in ['l1', 'l2', 'linf']:
                    conf = {
                        "dim": 2,
                        "div": div,
                        "hdist": hdist,
                        "max_cdist": max_cdist,
                        "ratio_func_name": name,
                        "rel_tol_inv": 1000,
                        "rel_tol_inv_mult": 10
                    }
                    configs.append(conf)

    for conf in configs:
        perebor(conf)
    logging.warning('======')


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument('config')
    args = argparser.parse_args()
    with open(args.config) as fh:
        config = json.load(fh)

    perebor(config)
