# coding: utf-8

import unittest

# run script from peano directory
import sys
import os
sys.path.append(os.path.dirname(sys.argv[0]) + '/..')

from gen_curve import *

class TestGen(unittest.TestCase):
    def test_brkline1(self):
        gen = CurveGenerator(div=4, exit=(1,0), allow_vertex_transit=True)
        res = gen.generate_brkline(start_width=5, finish_width=8)
        assert len(res) == 298

    def test_brkline2(self):
        gen = CurveGenerator(div=5, exit=(1,0), allow_vertex_transit=True)
        res = gen.generate_brkline(start_width=6, finish_width=11)
        assert len(res) == 49700

    def test_brkline3(self):
        gen = CurveGenerator(div=5, exit=(1,1), allow_vertex_transit=True)
        res = gen.generate_brkline(start_width=5, finish_width=10)
        assert len(res) == 2592

if __name__ == "__main__":
    unittest.main()
