#!/usr/bin/env python3
# coding: utf-8

import unittest

# run script from peano directory
import sys
import os
sys.path.append(os.path.dirname(sys.argv[0]) + '/../lib')

import logging
logging.basicConfig(level=0, stream=sys.stdout)

from gen_curve import *

class TestGen(unittest.TestCase):
    def test_brkline1(self):
        gen = CurveGenerator(dim=2, div=4, hdist=1)
        res = list(gen.generate_brklines(start_max_count=1, finish_max_count=10**6))
        assert len(res) == 298

    def test_brkline2(self):
        gen = CurveGenerator(dim=2, div=5, hdist=1)
        res = list(gen.generate_brklines(start_max_count=1, finish_max_count=10**6))
        assert len(res) == 49700

    def test_brkline3(self):
        gen = CurveGenerator(dim=2, div=5, hdist=2)
        res = list(gen.generate_brklines(start_max_count=1, finish_max_count=10**6))
        assert len(res) == 2592

if __name__ == "__main__":
    unittest.main()
