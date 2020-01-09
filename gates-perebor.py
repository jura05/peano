import json
import logging

from peano.paths import Gate, PathsGenerator
from peano.fast_fractions import FastFraction

from perebor import perebor

logging.basicConfig(level=logging.DEBUG)

#input = 'gates-322.txt.sorted'
input = 'gates-223.txt.sorted'



def get_pt(pt_str):
    xjs = []
    for xj in pt_str.strip('()').split(','):
        if '/' in xj:
            n, d = xj.split('/')
            xjf = FastFraction(int(n), int(d))
        else:
            xjf = FastFraction(int(xj), 1)
        xjs.append(xjf)
    return tuple(xjs)


def get_gate(gate_str):
    entr_str, exit_str = gate_str.split('->')
    return Gate(entrance=get_pt(entr_str), exit=get_pt(exit_str))


gates_list = []
with open(input) as ifh:
    for line in ifh:
        gate_strs = json.loads(line.replace("'", '"'))
        gates = [get_gate(gs) for gs in gate_strs]
        gates_list.append(gates)

pattern_count = len(gates_list[0])

def is_vertex(pt):
    if any(FastFraction(0, 1) < xj < FastFraction(1, 1) for xj in pt):
        return False
    return True

def cnt_bnd_coords(pt):
    return len(list(xj for xj in pt if xj == FastFraction(0, 1) or xj == FastFraction(1, 1)))

for idx, gates in enumerate(gates_list):
    pts = []
    for g in gates:
        pts += [g.entrance, g.exit]
    vertex_count = len(list(pt for pt in pts if is_vertex(pt)))
    if vertex_count != 0:
        continue

    total_bnd_coords = sum(cnt_bnd_coords(gate.entrance) + cnt_bnd_coords(gate.exit) for gate in gates)
    add_bnd_coords = total_bnd_coords - 2 * pattern_count
    if add_bnd_coords != 0:
        continue

    conf = {
        'dim': 2,
        'div': 2,
        'ratio_func_name': 'linf',
        #'ratio_func_name': 'l2',
        'gates': gates,
        'rel_tol_inv': 10000,
        'upper_bound': FastFraction(5, 1),
        #'upper_bound': FastFraction(51, 10),
    }
    print('GATE:', idx, [str(g) for g in gates])
    perebor(conf)
    print('===')
    print('')
