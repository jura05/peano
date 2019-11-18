import unittest

from peano.base_maps import BaseMap


class TestBaseMap(unittest.TestCase):

    def test_mul(self):
        bm1 = BaseMap([(0, True), (1, False)])  # (1-x,y)
        bm2 = BaseMap([(1, False), (0, False)])  # (y,x)
        self.assertEqual(bm1 * bm2, BaseMap([(1, True), (0, False)]))
        self.assertEqual(bm2 * bm1, BaseMap([(1, False), (0, True)]))

    def test_inv(self):
        bm = BaseMap([(3, True), (2, True), (1, True), (0, True)])
        self.assertEqual(bm * ~bm, BaseMap.id_map(dim=4))
        self.assertEqual(~bm * bm, BaseMap.id_map(dim=4))
        self.assertEqual(bm, bm.reverse_time().reverse_time())
