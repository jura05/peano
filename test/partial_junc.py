#!/usr/bin/env python3
# coding: utf-8

import unittest

# run script from peano directory
import sys
import os
sys.path.append(os.path.dirname(sys.argv[0]) + '/..')

from examples import *
from fractal_curve import FractalCurve

class TestCurve(unittest.TestCase):
    def setUp(self):
        self.curves = [
            get_peano_curve().forget(),
            get_meurthe_curve().forget(),
        ]

    def test_curves(self):
        for pcurve in self.curves:
            junc_info = pcurve.get_junctions_info()
            for curve in FractalCurve.get_possible_curves(pcurve):
                juncs = curve.get_junctions()
                # проверяем, что для каждого найденного стыка есть порождающая кривая
                for junc in juncs:
                    if junc not in junc_info:
                        raise Exception("Unknown junc!")
                    if not any(curve.is_specialization(tmpl) for tmpl in junc_info[junc]):
                        raise Exception("Can't found consistent curve!")

                # проверяем, что для каждого не найденного стыка нет порождающих кривых
                for junc, curves in junc_info.items():
                    if junc in juncs:
                        continue
                    if any(curve.is_specialization(tmpl) for tmpl in junc_info[junc]):
                        raise Exception("Found consistent curve for wrong junc!")

if __name__ == "__main__":
    unittest.main()
