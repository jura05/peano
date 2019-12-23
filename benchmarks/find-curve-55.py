#!/usr/bin/env python3
"""
Benchmark: find best 5*5 curve.
Run on 2.6Ghz CPU: 1m35s, bisect: 152, get_bounds: 3.9m hit + 560k miss
"""

import os
import psutil

import sys
sys.path.append('.')

from peano import paths
from peano.ratio import Estimator
from peano.utils import ratio_l2_squared
from peano.curves import PathFuzzyCurve

import logging

logging.basicConfig(level=0, stream=sys.stdout, format='[%(process)d] %(asctime)s %(message)s')


def main():
    curve_gen = paths.PathsGenerator(dim=2, div=5, hdist=1, max_cdist=1, verbose=1)
    pcurves = list(PathFuzzyCurve.init_from_paths([path]) for path in curve_gen.generate_paths(uniq=True))
    estimator = Estimator(ratio_l2_squared)
    res = estimator.estimate_ratio_sequence(
        pcurves, rel_tol_inv=1000000, rel_tol_inv_mult=2, sat_strategy={'type': 'geometric', 'multiplier': 1.5},
    )
    print(res)
    print(estimator.stats)
    process = psutil.Process(os.getpid())
    print('RSS:', process.memory_info().rss)  # in bytes


if __name__ == "__main__":
    main()
