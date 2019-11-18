import unittest

from peano.fast_fractions import FastFraction


class TestFF(unittest.TestCase):

    def test_mul(self):
        assert FastFraction(2, 3) * FastFraction(-1, 2) == FastFraction(-1, 3)
        assert FastFraction(1, 5) * FastFraction(5, -1) == FastFraction(-1, 1)

    def test_add(self):
        assert FastFraction(1, 3) + FastFraction(1, 6) == FastFraction(1, 2)
        assert FastFraction(1, 5) + FastFraction(1, -5) == FastFraction(0, 1)

    def test_sub(self):
        assert FastFraction(1, 2) - FastFraction(1, 6) == FastFraction(1, 3)

    def test_neg(self):
        assert -FastFraction(3, 5) == FastFraction(-3, 5)
        assert -FastFraction(0, 1) == FastFraction(0, -1)