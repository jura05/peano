import unittest

from peano.examples import *
from peano.sat_adapters import CurveSATAdapter


class TestSAT(unittest.TestCase):
    def setUp(self):
        self.curves = [
            get_peano_curve().forget(),
            get_meurthe_curve().forget(),
        ]

    def test_sat(self):
        for pcurve in self.curves:
            for curve in pcurve.gen_possible_curves():
                adapter = CurveSATAdapter(dim=pcurve.dim)
                adapter.init_curve(pcurve)
                model = adapter.get_model_from_curve(curve)
                juncs = list(curve.gen_junctions())
                for junc in juncs:
                    junc_var = adapter.get_junc_var(junc)
                    if not model[junc_var]:
                        raise Exception("Bad junc_var: False for existent junc!")
#                for junc in junc_info:
#                    if junc not in juncs:
#                        junc_var = adapter.get_junc_var(junc)
#                        if model[junc_var]:
#                            raise Exception("Bad junc_var: True for non-existent junc!")
                print('.', end='', flush=True)

            print('*', flush=True)
