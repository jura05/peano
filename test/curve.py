#!/usr/bin/env python3
# coding: utf-8

import unittest

# run script from peano directory
import sys
import os
sys.path.append(os.path.dirname(sys.argv[0]) + '/..')

from examples import *

class TestCurve(unittest.TestCase):
    def setUp(self):
        self.curves = [
            get_hilbert_curve(),
            get_peano_curve(),
            get_peano5_curve(),
            get_tokarev_curve(),
            get_meurthe_curve(),
            get_coil_curve(),
            get_serpentine_curve(),
            get_R_curve(),
            get_haverkort_curve_1(),
            get_haverkort_curve_2(),
            get_rev_curve(),
            get_rev2_curve(),
            get_rev3_curve(),
        ]

    def test_curves(self):
        for curve in self.curves:
            rrcurve = curve.reverse().reverse()
            self.assertEqual(curve.proto, rrcurve.proto)
            self.assertEqual(curve.base_maps, rrcurve.base_maps)

            curve.check()
            rev_curve = curve.reverse()
            rev_curve.check()

    def test_fractions(self):
        for curve in self.curves:
            for i in range(curve.genus):
                fraction = curve.get_fraction(i)
                fraction.check()

    def test_subdivision(self):
        for curve in self.curves:
            curve.get_subdivision().check()

    def test_subsubdivision(self):
        for curve in self.curves:
            curve.get_subdivision(2).check()

    def test_juncs(self):
        for curve in self.curves:
            juncs = curve.get_junctions()
            print('junctions count:', len(juncs))

    def test_apply_base_map(self):
        # проверяем, что если применить к кривой последовательно ряд преобразований,
        # произведение которых равно единичному, то мы получим исходную кривую
        bm_rev = BaseMap.id_map(dim=2).reverse_time()
        self.assertEqual(
            get_hilbert_curve().apply_base_map(bm_rev).proto,
            ((1, 0), (1, 1), (0, 1), (0, 0)),
        )

    def test_apply_base_map_unit(self):
        # проверяем, что если применить к кривой последовательно ряд преобразований,
        # произведение которых равно единичному, то мы получим исходную кривую
        rot_90 = BaseMap([1,0],[True,False])  # поворот на 90 градусов
        rot_90_t = BaseMap([1,0],[True,False], time_rev=True)
        bms_list = [
            [rot_90] * 3,
            [rot_90_t] * 3,
            [
                BaseMap([1,0],[False,False], time_rev=True),
                BaseMap([0,1],[True,False]),
            ],
            [
                BaseMap([0,2,1],[True,False,True], time_rev=True),
                BaseMap([1,0,2],[False,True,True], time_rev=True),
            ],
        ]
        for bms in bms_list:
            # will make bms[0] * bms[1] * ... * bms[-1] * last_map = id <=> last_map = bms[-1]^{-1} * ... * bms[0]^{-1}
            last_map = BaseMap.id_map(dim=bms[0].dim)
            for bm in bms:
                last_map = bm.inverse() * last_map

            for curve in self.curves:
                if curve.dim != bms[0].dim:
                    continue
                orig = curve
                current = curve
                for bm in reversed(bms + [last_map]):
                    current = current.apply_base_map(bm)
                self.assertEqual(orig.proto, current.proto)
                self.assertEqual(orig.base_maps, current.base_maps)


if __name__ == "__main__":
    unittest.main()
