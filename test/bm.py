#!/usr/bin/env python3

import unittest

import sys
import os
from fractions import Fraction

sys.path.append(os.path.dirname(sys.argv[0]) + '/..')
from peano.base_maps import BaseMap

class TestBaseMap(unittest.TestCase):

    def test_mul(self):
        bm1 = BaseMap([0,1],[True,False])  # (1-x,y)
        bm2 = BaseMap([1,0],[False,False])  # (y,x)
        self.assertEqual(bm1 * bm2, BaseMap([1,0],[True,False]))
        self.assertEqual(bm2 * bm1, BaseMap([1,0],[False,True]))

    def test_inv(self):
        bm = BaseMap(perm=[3,2,1,0], flip=[True,True,True,True])
        self.assertEqual(bm * bm.inverse(), BaseMap.id_map(dim=4))
        self.assertEqual(bm, bm.reverse_time().reverse_time())


if __name__ == "__main__":
    unittest.main()
