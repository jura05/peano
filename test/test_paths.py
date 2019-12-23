import unittest

import logging
import sys

from peano.paths import PathsGenerator


class TestGen(unittest.TestCase):
    def setUp(self):
        logging.basicConfig(level=0, stream=sys.stdout)
        self.kws = {'start_max_count': 1, 'finish_max_count': 10 ** 6}

    def test_p1(self):
        gen = PathsGenerator(dim=2, div=4, hdist=1)
        assert len(list(gen.generate_paths(**self.kws))) == 298
        assert len(list(gen.generate_paths(uniq=True, **self.kws))) == 162

    def test_p2(self):
        gen = PathsGenerator(dim=2, div=5, hdist=1)
        assert len(list(gen.generate_paths(**self.kws))) == 49700
        assert len(list(gen.generate_paths(uniq=True, **self.kws))) == 24850

    def test_p3(self):
        gen = PathsGenerator(dim=2, div=5, hdist=2)
        assert len(list(gen.generate_paths(**self.kws))) == 2592
        assert len(list(gen.generate_paths(uniq=True, **self.kws))) == 659
