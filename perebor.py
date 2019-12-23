#!/usr/bin/env python3
# coding: UTF-8

import argparse
import json
import logging

from peano import utils
from peano.paths import PathsGenerator, Gate
from peano.ratio import Estimator
from peano.curves import PathFuzzyCurve
from peano.fast_fractions import FastFraction


logging.basicConfig(level=logging.DEBUG)


def get_pt(pt_str):
    xjs = []
    for xj in pt_str.strip('()').split(','):
        if '/' in xj:
            n, d = xj.split('/')
            xjf = FastFraction(int(n), int(d))
        else:
            xjf = FastFraction(int(xj), 1)
        xjs.append(xjf)
    return tuple(xjs)


def get_gate(gate_str):
    entr_str, exit_str = gate_str.split('->')
    return Gate(entrance=get_pt(entr_str), exit=get_pt(exit_str))


def perebor(conf, idx=None):
    logging.warning('CONF: %s', conf)
    funcs = {
        'l1': utils.ratio_l1,
        'l2': utils.ratio_l2,
        'l2_squared': utils.ratio_l2_squared,
        'linf': utils.ratio_linf,
    }
    ratio_func = funcs[conf['ratio_func_name']]
    dim, div = conf['dim'], conf['div']
    gates = None
    if 'gates' in conf:
        gates = [get_gate(g) for g in conf['gates']]
    paths_gen = PathsGenerator(
        dim=dim, div=div, hdist=conf.get('hdist'), gates=gates, max_cdist=conf.get('max_cdist'),
    )
    paths_list = list(paths_gen.generate_paths(uniq=True))
    logging.warning('paths: %d', len(paths_list))

    if idx is not None:
        logging.warning('uniq paths idx: %d:%d', idx[0], idx[1])
        paths = paths[idx[0]:idx[1]]

    estimator = Estimator(ratio_func, cache_max_size=2**16)

    #logging.warning('sort by paths ratio ...')
    #path_ratio = {path: estimator.estimate_path(path) for path in paths}
    #paths.sort(key=lambda path: path_ratio[path])

    pcurves = (PathFuzzyCurve.init_from_paths(paths) for paths in paths_list)
    res = estimator.estimate_ratio_sequence(
        pcurves,
        rel_tol_inv=conf['rel_tol_inv'],
        rel_tol_inv_mult=conf.get('rel_tol_inv_mult', 2),
        sat_strategy={'type': 'geometric', 'multiplier': 1.3},
        upper_bound=conf.get('upper_bound'),
    )

    logging.warning('CONF: %s', conf)
    logging.warning('FINAL BOUNDS: %.6f <= r <= %.6f', res['lo'], res['up'])
    res['paths_count'] = len(paths_list)
    print(res)


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument('config')
    #argparser.add_argument('--begin', type=int)
    #argparser.add_argument('--end', type=int)
    args = argparser.parse_args()
    with open(args.config) as fh:
        config = json.load(fh)

    #config['upper_bound'] = FastFraction(559, 100)
    perebor(config)
