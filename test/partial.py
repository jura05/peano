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
            get_hilbert_curve().forget(),
            get_peano_curve().forget(),
            get_tokarev_curve().forget(),
        ]

    def test_curves(self):
        for curve in self.curves:
            cnum = 0
            bms = curve.get_allowed_maps(cnum)
            for bm in bms:
                scurve = curve.specify(cnum, bm)
                piece = scurve.get_fraction(cnum)


if __name__ == "__main__":
    unittest.main()
