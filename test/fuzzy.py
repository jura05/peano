#!/usr/bin/env python3
# coding: utf-8

import unittest

# run script from peano directory
import sys
import os
sys.path.append(os.path.dirname(sys.argv[0]) + '/..')

from peano.examples import *
from peano.base_maps import Spec

class TestCurve(unittest.TestCase):
    def setUp(self):
        self.curves = [
            get_hilbert_curve(),
            get_peano_curve(),
            get_tokarev_curve(),
            get_rev_curve(),
        ]

    def test_curves(self):
        for curve in self.curves:
            pcurve = curve.forget()
            rrcurve = pcurve.reverse().reverse()
            self.assertEqual(pcurve.proto, rrcurve.proto)
            self.assertEqual(pcurve.sp_info(), rrcurve.sp_info())

        for curve in self.curves:
            pcurve = curve.forget(allow_time_rev=True)
            cnum = 0
            for bm in pcurve.gen_allowed_maps(cnum):
                scurve = pcurve.specify(pnum=0, cnum=cnum, spec=Spec(bm))
                piece = scurve.get_fraction(cnum)

        curve0 = self.curves[0].forget()
        for c in curve0.gen_possible_curves():
            c.check()

        curve0 = self.curves[0].forget(allow_time_rev=True)
        for c in curve0.gen_possible_curves():
            c.check()

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
                    current = bm * current
                self.assertEqual(orig.proto, current.proto)
                self.assertEqual(orig.base_maps, current.base_maps)

    def test_junc(self):
        for num, curve in enumerate(self.curves):
            if num == 0:
                pcurve = curve.forget(allow_time_rev=True)
            else:
                pcurve = curve.forget()

            junc_info = pcurve.get_junctions_info()
            for junc in junc_info:
                if any(dx not in set([0,1,-1]) for dx in junc.delta_x):
                    raise Exception("Bad delta!")

            for curve in pcurve.gen_possible_curves():
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
