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
            get_tokarev_curve(),
        ]

    def test_curves(self):
        for curve in self.curves:
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


if __name__ == "__main__":
    unittest.main()
