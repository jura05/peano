import unittest

import sys
import os
from fractions import Fraction

sys.path.append(os.path.dirname(sys.argv[0]) + '/..')
from base_map import BaseMap, PieceMap

class TestBaseMap(unittest.TestCase):

    def test_mul(self):
        bm1 = BaseMap([0,1],[True,False])  # (1-x,y)
        bm2 = BaseMap([1,0],[False,False])  # (y,x)
        self.assertEqual(bm1 * bm2, BaseMap([1,0],[True,False]))
        self.assertEqual(bm2 * bm1, BaseMap([1,0],[False,True]))

    def test_inv(self):
        bm = BaseMap(perm=[3,2,1,0], flip=[True,True,True,True])
        self.assertEqual(bm * bm.inverse(), BaseMap(dim=4))


class TestPieceMap(unittest.TestCase):

    def setUp(self):
        self.pm_examples = [
            PieceMap(
                base_map=BaseMap([3,2,1,0],[True,False,False,False]),
                shift=(-1,-2,-3,-4),
                scale=3,
                time_shift=-10,
                time_scale=1,
            ),
            PieceMap(
                base_map=BaseMap([4,0,1,3,2],[True,False,False,True,True]),
                shift=(-Fraction(1, 2), Fraction(1, 2), Fraction(1, 3), Fraction(1, 5), Fraction(1, 7)),
                scale=Fraction(1, 11),
                time_shift=0,
                time_scale=1,
            ),
        ]

    def test_id_map(self):
        id_map = PieceMap.id_map(3)
        self.assertEqual(
            id_map.apply((0,1,2), 3),
            ((0,1,2), 3),
        )

    def test_inv(self):
        for pm in self.pm_examples:
            id_map = PieceMap.id_map(pm.dim)
            self.assertEqual(pm * pm.inverse(), id_map)
            self.assertEqual(pm.inverse() * pm, id_map)

    def test_apply(self):
        for pm in self.pm_examples:
            x = tuple(j for j in range(pm.dim))
            t = 1
            x1, t1 = pm.apply(x, t)
            x2, t2 = pm.inverse().apply(x1, t1)
            self.assertEqual([x,t], [x2,t2])


if __name__ == "__main__":
    unittest.main()
