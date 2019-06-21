#!/usr/bin/env python3
# coding: utf-8

import unittest

# run script from peano directory
import sys
import os
sys.path.append(os.path.dirname(sys.argv[0]) + '/../lib')

from examples import *
from fractal_curves import FractalCurve
from sat_adapters import CurveSATAdapter

class TestCurve(unittest.TestCase):
    def setUp(self):
        self.curves = [
            get_peano_curve().forget(),
            get_meurthe_curve().forget(),
        ]

    def test_curves(self):
        for pcurve in self.curves:
            for curve in FractalCurve.gen_possible_curves(pcurve):
                adapter = CurveSATAdapter(dim=pcurve.dim)
                adapter.init_curve(pcurve)
                model = adapter.get_model_from_curve(curve)
                juncs = curve.get_junctions()
                for junc in juncs:
                    junc_var = adapter.get_junc_var(junc)
                    if not model[junc_var]:
                        raise Exception("Bad junc_var: False for existent junc!")
#                for junc in junc_info:
#                    if junc not in juncs:
#                        junc_var = adapter.get_junc_var(junc)
#                        if model[junc_var]:
#                            raise Exception("Bad junc_var: True for non-existent junc!")
                print('.', end='', flush=True)

            print('*', flush=True)




if __name__ == "__main__":
    unittest.main()

