# coding: utf-8

# run script from peano directory
import sys
import os
sys.path.append(os.path.dirname(sys.argv[0]) + '/..')

from examples import *

def main():
    # run some tests
    for curve, name in [
        (get_hilbert_curve(), 'Hilbert'),
        (get_peano_curve(), 'Peano'),
        (get_tokarev_curve(), 'Tokarev'),
        #(get_rev_curve(), 'with_rev'),
        #(get_discontinuous_curve(), 'discontinuous'),
    ]:
        print('curve:', name)
        print('entrance:', curve.get_entrance())
        print('exit:', curve.get_exit())
        try:
            id_map = BaseMap(dim=curve.dim)
            for bm in curve.base_maps:
                # заодно проверим BaseMap
                inv = bm.inverse()
                assert bm * inv == id_map, 'base_map multiplication'
                assert inv * bm == id_map, 'base_map multiplication'
            curve.check()
            print('check ok!')
        except Exception as exc:
            print('check failed:', exc)
            continue

        rev_curve = curve.reverse()
        rev_curve.check()
        print('check reversed ok!')

        for i in range(curve.genus()):
            fraction = curve.get_fraction(i)
            fraction.check()
        print('fractions are correct')

        subdiv_curve = curve.get_subdivision()
        subdiv_curve.check()
        print('sub-division is correct')

        subsub_curve = curve.get_subdivision(2)
        subsub_curve.check()
        print('sub-division(2) is correct')
        print('sub-division(2) genus:', subsub_curve.genus())

        for delta, base_map in curve.get_junctions():
            print('junction:', delta, base_map)
        print()


if __name__ == "__main__":
    main()