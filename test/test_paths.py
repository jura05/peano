import unittest

import sys

from peano.paths import *


logging.basicConfig(level=0, stream=sys.stdout)


class TestGen(unittest.TestCase):
    def test_p1(self):
        gen = PathsGenerator(dim=2, div=4, hdist=1)
        res = list(gen.generate_paths(start_max_count=1, finish_max_count=10 ** 6))
        assert len(res) == 298
        assert len(list(gen_uniq(2, res))) == 162

    def test_p2(self):
        gen = PathsGenerator(dim=2, div=5, hdist=1)
        res = list(gen.generate_paths(start_max_count=1, finish_max_count=10 ** 6))
        assert len(res) == 49700
        assert len(list(gen_uniq(2, res))) == 24850

    def test_p3(self):
        gen = PathsGenerator(dim=2, div=5, hdist=2)
        res = list(gen.generate_paths(start_max_count=1, finish_max_count=10 ** 6))
        assert len(res) == 2592
        assert len(list(gen_uniq(2, res))) == 659
