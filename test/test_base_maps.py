import unittest

from peano.base_maps import BaseMap


class TestBaseMap(unittest.TestCase):
    def setUp(self):
        self.base_maps = []
        for dim in range(2, 6):
            self.base_maps += list(BaseMap.gen_base_maps(dim))

    def test_mul(self):
        bm1 = BaseMap.from_basis('Ij')
        bm2 = BaseMap.from_basis('ji')
        self.assertEqual(bm1 * bm2, BaseMap.from_basis('Ji'))
        self.assertEqual(bm2 * bm1, BaseMap.from_basis('jI'))

    def test_inv(self):
        for bm in self.base_maps:
            self.assertEqual(bm * ~bm, BaseMap.id_map(dim=bm.dim))
            self.assertEqual(~bm * bm, BaseMap.id_map(dim=bm.dim))

    def test_conj(self):
        dim3 = [bm for bm in self.base_maps if bm.dim == 3]
        for bm1 in dim3:
            for bm2 in dim3:
                self.assertEqual(bm1.conjugate_by(bm2), bm2 * bm1 * ~bm2)
