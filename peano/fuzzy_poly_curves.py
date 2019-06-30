from collections import namedtuple

from .base_maps import BaseMap


Pattern = namedtuple('Pattern', ['proto', 'specs'])

class FuzzyPolyCurve:
    def __init__(self, dim, div, patterns):
        self.dim = dim
        self.div = div
        self.patterns = patterns

    def __getitem__(self, pnum):
        return self.patterns[pnum]
