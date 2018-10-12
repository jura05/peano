#!/usr/bin/env python3

from fractions import Fraction

import sys
import gen_curve

from multiprocessing import Process, Queue

import logging

logging.basicConfig(level=0, stream=sys.stdout, format='[%(process)d] %(asctime)s %(message)s')

def ratio_linf(d, dv, dt):
    return Fraction(max(abs(x) for x in dv)**d, dt)

def ratio_l1(d, dv, dt):
    return Fraction(sum(abs(x) for x in dv)**d, dt)

def ratio_l2(d, dv, dt):
    return Fraction(sum(x**2 for x in dv)**d, dt**2)


def submain(tasks_queue, results_queue):
    while True:
        item = tasks_queue.get()
        if item['type'] == 'task':
            curve = item['curve']
            res = curve.estimate_ratio(ratio_l2, rel_tol=0.01, upper_bound=item['upper_bound'], verbose=0)
            logging.info('processed curve, result: %s', res)
            results_queue.put(res)
        elif item['type'] == 'STOP':
            break

def main():
    total_iter = 0
    n_proc = 2
    procs = []
    tasks_queue = Queue(maxsize=n_proc)
    results_queue = Queue(maxsize=n_proc * 2)
    for j in range(n_proc):
        proc = Process(target=submain, args=(tasks_queue, results_queue))
        proc.start()

    seen_curves = processed_curves = 0
    generator = gen_curve.CurveGenerator(div=4, exit=(1,0), allow_vertex_transit=False)
    best_upper_bound = None
    for curve in generator.generate_curves():
        seen_curves += 1
        logging.info('dispatch curve %d', seen_curves)
        if best_upper_bound is not None:
            logging.info('BEST_UPPER_BOUND: %f', best_upper_bound)
        tasks_queue.put({'type': 'task', 'curve': curve, 'upper_bound': best_upper_bound})
        while not results_queue.empty():
            res = results_queue.get()
            processed_curves += 1
            total_iter += res['stats']['iter']
            if best_upper_bound is None or res['upper_bound'] < best_upper_bound:
                best_upper_bound = res['upper_bound']
                logging.info('new upper bound found: %f', best_upper_bound)

    while processed_curves < seen_curves:
        res = results_queue.get()
        processed_curves += 1
        if best_upper_bound is None or res['upper_bound'] < best_upper_bound:
            best_upper_bound = res['upper_bound']
            logging.info('new upper bound found: %f', best_upper_bound)

    for j in range(n_proc):
        tasks_queue.put({'type': 'STOP'})

    print('BEST:', best_upper_bound, 'total_iter:', total_iter)

if __name__ == "__main__":
    main()
