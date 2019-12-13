import unittest

from peano.fast_fractions import FastFraction
from peano.examples import *
from peano import utils
from peano.paths import PathsGenerator
from peano.curves import PathFuzzyCurve, Proto
from peano.ratio import Estimator


class TestCurve(unittest.TestCase):

    def test_curve_ratio(self):
        known_bounds = [
            {
                'curve': get_hilbert_curve(),
                'ratio': { 'l2': 6, 'l1': 9, 'linf': 6},
            },
            {
                'curve': get_peano_curve(),
                'ratio': {'l2': 8, 'l1': 32/3, 'linf': 8},
            },
            {
                'curve': get_haverkort_curve_1(),
                'ratio': {'l1': (99 + 5/9)},
            },
            {
                'curve': get_tokarev_curve(),
                'ratio': {'l1': [98.2, 98.4], 'l2': [26.1, 26.3], 'linf': [24.1, 24.3]},
            },
            {
                'curve': get_scepin_bauman_curve(),
                'ratio': {'l1': (10 + 2/3), 'l2': (5 + 2/3), 'linf': (5 + 1/3)},
            },
            {
                'curve': get_meurthe_curve(),
                'ratio': {'l1': (10 + 2/3), 'l2': (5 + 2/3), 'linf': (5 + 1/3)},
            },
            {
                'curve': get_serpentine_curve(),
                'ratio': {'l1': 10, 'l2': 6.25, 'linf': 5.625},
            },
            {
                'curve': get_coil_curve(),
                'ratio': {'l1': (10 + 2/3), 'l2': (6 + 2/3), 'linf': (6 + 2/3)},
            },
            {
                'curve': get_R_curve(),
                'ratio': {'l1': (10 + 2/3), 'l2': (6 + 2/3), 'linf': (6 + 2/3)},
            },
            {   
                'curve': get_haverkort_curve_A26(),
                'ratio': {'l1': (99 + 5/9), 'l2': [22.7,22.9], 'linf': (12 + 4/9)},
            },
            {   
                'curve': get_haverkort_curve_2(),
                'ratio': {'l1': [89.7, 89.8], 'l2': [18,19], 'linf': 14},
            },
        ]
        for data in known_bounds:
            curve = data['curve']
            ratio_dict = data['ratio']
            for metric in sorted(ratio_dict.keys()):
                ratio = ratio_dict[metric]
                if isinstance(ratio, list):
                    ratio_lo, ratio_up = ratio
                else:
                    ratio_lo = ratio * 0.999
                    ratio_up = ratio * 1.001

                if metric == 'l2':
                    func = utils.ratio_l2_squared
                    ratio_lo, ratio_up = ratio_lo**2, ratio_up**2
                elif metric == 'l1':
                    func = utils.ratio_l1
                elif metric == 'linf':
                    func = utils.ratio_linf

                res = Estimator(func).estimate_ratio(curve, rel_tol_inv=10**5, verbose=False)
                print(res)
                assert float(res['up']) <= ratio_up, 'metric {} up failed: {} > {}'.format(metric, res['up'], ratio_up)
                assert float(res['lo']) >= ratio_lo, 'metric {} lo failed: {} < {}'.format(metric, res['lo'], ratio_lo)

    def test_polycurve_ratio(self):
        known_bounds = [
            {
                'curve': get_beta_Omega_Curve(),
                'ratio': { 'l2': 5, 'l1': 9, 'linf': 5},
            },
            {
                'curve': get_neptunus_curve(),
                'ratio': { 'l2': [18.2, 18.4], 'linf': [9.44, 9.46]},
            },
        ]
        for data in known_bounds:
            curve = data['curve']
            for metric, ratio in data['ratio'].items():
                if isinstance(ratio, list):
                    ratio_lo, ratio_up = ratio
                else:
                    ratio_lo = ratio * 0.999
                    ratio_up = ratio * 1.001

                if metric == 'l2':
                    func = utils.ratio_l2_squared
                    ratio_lo, ratio_up = ratio_lo**2, ratio_up**2
                elif metric == 'l1':
                    func = utils.ratio_l1
                elif metric == 'linf':
                    func = utils.ratio_linf

                res = Estimator(func).estimate_ratio(curve, rel_tol_inv=10**4)
                assert float(res['up']) <= ratio_up, 'metric {} up failed: {} > {}'.format(metric, res['up'], ratio_up)
                assert float(res['lo']) >= ratio_lo, 'metric {} lo failed: {} < {}'.format(metric, res['lo'], ratio_lo)

    #def test_pcurve_ratio(self):
    #    pcurve = get_peano5_curve().forget(allow_time_rev=True)
    #    assert pcurve.test_ratio(utils.ratio_l2_squared, lower_bound=FastFraction(90, 1), upper_bound=FastFraction(100, 1))

    def test_55_ratio(self):
        good_proto = Proto(dim=2, div=5, cubes=[
            (0, 0), (0, 1), (1, 1), (1, 0), (2, 0),
            (2, 1), (2, 2), (1, 2), (0, 2), (0, 3),
            (0, 4), (1, 4), (1, 3), (2, 3), (2, 4),
            (3, 4), (4, 4), (4, 3), (3, 3), (3, 2),
            (4, 2), (4, 1), (3, 1), (3, 0), (4, 0),
        ])
        # in new version we have (0,0)->(0,1) gate
        good_proto = BaseMap.from_basis('ji') * good_proto

        paths_gen = PathsGenerator(dim=2, div=5, hdist=1, max_cdist=1, verbose=1)
        for path in paths_gen.generate_paths():
            if path.proto == good_proto:
                path0 = path
                break

        pcurve = PathFuzzyCurve.init_from_paths([path0])
        estimator = Estimator(utils.ratio_l2_squared)
        curve = estimator.estimate_ratio(pcurve, rel_tol_inv=10000, verbose=False)['curve']
        ratio = estimator.estimate_ratio(curve, rel_tol_inv=10000, use_vertex_brkline=True, verbose=False)

        assert ratio['lo'] == (FastFraction(408, 73)**2)
