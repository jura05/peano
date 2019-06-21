#!/usr/bin/env python3

import unittest

import sys
import os

sys.path.append(os.path.dirname(sys.argv[0]) + '/../lib')
from fast_fractions import FastFraction

class TestFF(unittest.TestCase):

    def test_mul(self):
        x = FastFraction(2, 3)
        y = FastFraction(-1, 2)
        z = FastFraction(-1, 3)
        assert x * y == z

    def test_add(self):
        x = FastFraction(1, 3)
        y = FastFraction(1, 6)
        z = FastFraction(1, 2)
        assert x + y == z


if __name__ == "__main__":
    unittest.main()
