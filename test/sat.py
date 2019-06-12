#!/usr/bin/env python3
# coding: utf-8

import unittest

# run script from peano directory
import sys
import os
sys.path.append(os.path.dirname(sys.argv[0]) + '/..')

from examples import *
from fractal_curve import FractalCurve
from curve_sat_adapter import CurveSATAdapter

class TestCurve(unittest.TestCase):
    def setUp(self):
        self.curves = [
            get_peano_curve().forget(),
            get_meurthe_curve().forget(),
        ]

    def test_curves(self):
        for pcurve in self.curves:
            adapter = CurveSATAdapter(dim=pcurve.dim)
            junc_info = pcurve.get_junctions_info()

            for junc, curves in junc_info.items():
                adapter.make_junc_var(junc, curves)

            for curve in FractalCurve.get_possible_curves(pcurve):
                juncs = curve.get_junctions()
                model = adapter.get_model_from_curve(curve)
                for junc in juncs:
                    junc_var = ('junc', junc)
                    if not model[junc_var]:
                        raise Exception("Bad junc_var: False for existent junc!")
                for junc in junc_info:
                    if junc not in juncs:
                        junc_var = ('junc', junc)
                        if model[junc_var]:
                            raise Exception("Bad junc_var: True for non-existent junc!")
                print('.', end='', flush=True)

            print('*', flush=True)




if __name__ == "__main__":
    unittest.main()
