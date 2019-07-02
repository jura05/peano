#!/usr/bin/env python3
# coding: utf-8

import unittest

# run script from peano directory
import sys
import os
sys.path.append(os.path.dirname(sys.argv[0]) + '/..')

from peano.examples import *

class TestCurve(unittest.TestCase):
    def setUp(self):
        self.curves = [
            get_hilbert_bicurve(),
            #get_ARW_Curve(),
            get_neptunus_curve(),
            get_luna_curve(),
        ]

    def test_junc(self):
        for curve in self.curves:
            for jcnt, junc in enumerate(curve.gen_junctions()):
                if jcnt > 100:
                    raise Exception("Too many juncs!")

    def test_apply_base_map_unit(self):
        # проверяем, что если применить к кривой последовательно ряд преобразований,
        # произведение которых равно единичному, то мы получим исходную кривую
        rev_t  = BaseMap([0,1],[False,False], time_rev=True)
        rot_90 = BaseMap([1,0],[True,False])  # поворот на 90 градусов
        rot_90_t = BaseMap([1,0],[True,False], time_rev=True)
        bms_list = [
            [rev_t],
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
                    current = bm * current
                self.assertEqual(orig.patterns[0].proto, current.patterns[0].proto)
                self.assertEqual(orig.patterns[0].specs, current.patterns[0].specs)


if __name__ == "__main__":
    unittest.main()
