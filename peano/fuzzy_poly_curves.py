from collections import namedtuple

from .base_maps import BaseMap


Spec = namedtuple('Spec', ['pnum', 'base_map'])
Pattern = namedtuple('Pattern', ['proto', 'specs'])

class FuzzyPolyCurve:
    pass
