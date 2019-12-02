#!/usr/bin/env python3
"""
Benchmark: find best 5*5 curve.
Run on 2.6Ghz CPU: 3m14s, bisect: 291, get_bounds: 8.8m
"""

import sys
sys.path.append('.')

from peano import gen_curve
from peano.ratio import Estimator
from peano.utils import ratio_l2_squared
from peano.curves import SymmFuzzyCurve

import logging

logging.basicConfig(level=0, stream=sys.stdout, format='[%(process)d] %(asctime)s %(message)s')


def main():
    curve_gen = gen_curve.CurveGenerator(dim=2, div=5, hdist=1, max_cdist=1, verbose=1)
    pcurves = []
    for brkline in curve_gen.generate_brklines():
        pcurves.append(SymmFuzzyCurve.init_from_brkline(2, 5, brkline, allow_time_rev=True))
    estimator = Estimator(ratio_l2_squared)
    res = estimator.estimate_ratio_sequence(pcurves, rel_tol_inv=1000000)
    print(res)

if __name__ == "__main__":
    main()
