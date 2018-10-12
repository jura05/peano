# coding: utf-8

import unittest

# run script from peano directory
import sys
import os
sys.path.append(os.path.dirname(sys.argv[0]) + '/..')

import logging
logging.basicConfig(level=0, stream=sys.stdout)

from gen_curve import *

class TestGen(unittest.TestCase):
    def test_brkline1(self):
        gen = CurveGenerator(div=4, exit=(1,0), allow_vertex_transit=True)
        res = gen.generate_brklines(start_max_count=1, finish_max_count=10**6)
        assert len(res) == 298

    def test_brkline2(self):
        gen = CurveGenerator(div=5, exit=(1,0), allow_vertex_transit=True)
        res = gen.generate_brklines(start_max_count=1, finish_max_count=10**6)
        assert len(res) == 49700

    def test_brkline3(self):
        gen = CurveGenerator(div=5, exit=(1,1), allow_vertex_transit=True)
        res = gen.generate_brklines(start_max_count=1, finish_max_count=10**6)
        assert len(res) == 2592

if __name__ == "__main__":
    unittest.main()
