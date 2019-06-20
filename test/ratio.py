#!/usr/bin/env python3
# coding: utf-8

import unittest

# run script from peano directory
import sys
import os
sys.path.append(os.path.dirname(sys.argv[0]) + '/../lib')

from examples import *
import utils


class TestCurve(unittest.TestCase):
    def setUp(self):
        self.data = [
            {'curve': get_hilbert_curve(), 'ratio': {'l2': 6, 'l1': 9, 'linf': 6}},
            {'curve': get_peano_curve(), 'ratio': {'l2': 8, 'l1': 32/3, 'linf': 8}},
            {'curve': get_haverkort_curve_1(), 'ratio': {'l1': (99 + 5/9)}},
        ]

    def test_ratio(self):
        for data in self.data:
            curve = data['curve']
            for metric, ratio in data['ratio'].items():
                if metric == 'l2':
                    func = utils.ratio_l2_squared
                    ratio = ratio**2
                elif metric == 'l1':
                    func = utils.ratio_l1
                elif metric == 'linf':
                    func = utils.ratio_linf

                res = curve.estimate_ratio_new(func, rel_tol=0.0001)
                assert res['up'] <= ratio * 1.001
                assert res['lo'] >= ratio * 0.999

if __name__ == "__main__":
    unittest.main()
