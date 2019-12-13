import time
import resource
from collections import namedtuple, Counter
import logging
import itertools

from .base_maps import BaseMap
from .fast_fractions import FastFraction


class Proto:
    """
    Curve prototype -- sequence of cubes.

    We allow None-s for some of cubes, to support usage of get_entrance/get_exit methods.
    """

    def __init__(self, dim, div, cubes):
        self.dim = dim
        self.div = div
        self._cubes = tuple(tuple(cube) if cube is not None else None for cube in cubes)

    def __iter__(self):
        return iter(self._cubes)

    def __getitem__(self, idx):
        return self._cubes[idx]

    def __len__(self):
        return len(self._cubes)

    def __eq__(self, other):
        return self._cubes == other._cubes

    def __hash__(self):
        return hash(self._cubes)

    def __rmul__(self, base_map):
        if not isinstance(base_map, BaseMap):
            return NotImplemented
        cubes = [base_map.apply_cube(self.div, cube) if cube is not None else None for cube in self._cubes]
        if base_map.time_rev:
            cubes = reversed(cubes)
        return type(self)(self.dim, self.div, cubes)


class Gate(namedtuple('Gate', ['entrance', 'exit'])):

    def __rmul__(self, base_map):
        ne, nx = base_map.apply_x_fraction(self.entrance), base_map.apply_x_fraction(self.exit)
        if base_map.time_rev:
            ne, nx = nx, ne
        return type(self)(ne, nx)

    def reversed(self):
        return type(self)(self.exit, self.entrance)

    def std(self):
        dim = len(self.entrance)
        return min(bm * self for bm in BaseMap.gen_base_maps(dim))

    def __str__(self):
        entr_str = '(' + ','.join([str(pj) for pj in self.entrance]) + ')'
        exit_str = '(' + ','.join([str(pj) for pj in self.exit]) + ')'
        return entr_str + '-->' + exit_str


class CurvePath:
    """Prototype + gates."""

    def __init__(self, proto, gates):
        self.proto = proto
        self.dim = proto.dim
        self.div = proto.div
        self.gates = tuple(gates)
        div_inv = FastFraction(1, proto.div)
        self.entrance = tuple((FastFraction(cj, 1) + ej) * div_inv for cj, ej in zip(proto[0], gates[0].entrance))
        self.exit = tuple((FastFraction(cj, 1) + ej) * div_inv for cj, ej in zip(proto[-1], gates[-1].exit))
        self.gate = Gate(self.entrance, self.exit)

    @classmethod
    def from_data(cls, dim, div, data):
        """Create path from triples (cube, entr, exit)."""
        cubes = (d[0] for d in data)
        proto = Proto(dim, div, cubes)
        gates = tuple(d[1] for d in data)
        return cls(proto, gates)

    def __eq__(self, other):
        return self.proto == other.proto and self.gates == other.gates

    def __hash__(self):
        return hash((self.proto, self.gates))

    def __rmul__(self, base_map):
        new_proto = base_map * self.proto
        new_gates = [base_map * gate for gate in self.gates]
        if base_map.time_rev:
            new_gates.reverse()
        return CurvePath(new_proto, new_gates)

    def reversed(self):
        return BaseMap.id_map(self.proto.dim).reversed_time() * self


def gen_uniq(dim, paths):
    seen = set()
    for path in paths:
        entr, exit = path.gate
        bms = list(BaseMap.gen_constraint_cube_maps(dim, {entr: entr, exit: exit}))
        bms += [bm.reversed_time() for bm in BaseMap.gen_constraint_cube_maps(dim, {entr: exit, exit: entr})]
        if any(bm * path in seen for bm in bms):
            continue
        seen.add(path)
        yield path


class PathNode(namedtuple('PathNode', ['head', 'gate', 'prev', 'len'])):
    """
    Node of the tree representing partial paths.

    Path is a curve prototype enriched with entr-exit gates.
    Here we allow partial paths, to be used during generation.

    head  --  end of path cube
    gate  --  gate for head cube
    prev  --  link to tail path (without head)
    len  --  path length
    """

    def get_data(self):
        data = []
        node = self
        while node is not None:
            data.append((node.head, node.gate))
            node = node.prev
        data.reverse()
        return data

    def flat(self):
        length = self.len
        path = self
        cubes = [None] * length
        for i in range(length):
            cubes[i] = path.head
            path = path.prev
        return cubes

    def support(self):
        """Set of cubes."""
        return set(self.flat())

    # состояние, которое определяет продолжимость пути
    # entrance нужен лишь для удобства пересечения с финишем
    def state(self):
        return (frozenset(self.flat()[1:]), self.head, self.gate)


class PathsGenerator:
    """Generate vertex curves of given div, with entrance at 0 and exit (1,1,..,1,0,0,..,0)."""

    def __init__(self, dim, div, gates=None, hdist=None, max_cdist=None, verbose=0):
        """
        Init paths generator.

        dim, div  --  subj
        hdist  --  k means curve with entr (0,0,..,0), exit (1,1,..,1,0,..,0), with k "1"-s
        gate  --  TODO
        max_cdist  --  maximum l1-distance between cubes
        verbose  --  subj
        """

        self.dim = dim
        self.div = div
        self.verbose = verbose
        self.stats = Counter()
        N = self.div

        if gates is None:
            gates = [Gate(
                entrance=(FastFraction(0, 1),) * dim,
                exit=(FastFraction(0, 1),) * (dim - hdist) + (FastFraction(1, 1),) * hdist,
            )]

        self.gates = gates
        self.next_dict = self.get_next_dict(gates, max_cdist)


    # point is in R^d, get integer cubes for it
    @staticmethod
    def gen_cubes(point):
        # ranges for cubej
        minjs = []
        maxjs = []
        for xj in point:
            ij = int(xj)
            if xj == FastFraction(ij, 1):
                minj, maxj = ij-1, ij
            else:
                minj, maxj = ij, ij
            minjs.append(minj)
            maxjs.append(maxj)
        ranges = [range(minj, maxj + 1) for minj, maxj in zip(minjs, maxjs)]
        yield from itertools.product(*ranges)

    def get_next_dict(self, gates, max_cdist):
        """
        Get a dict: exit => [(cube_delta, new_gate), ...]

        Params:
            dim:    subj
            gate:   TODO
            max_cdist: do not allow cube changes greater than it
        """

        dim = self.dim

        gate_set = set()
        for bm in BaseMap.gen_base_maps(dim):
            for gate in gates:
                gate_set.add(bm * gate)
        gate_list = sorted(gate_set)

        # start positions - all exits
        result = {}
        for pt in [g.exit for g in gate_list]:
            if pt in result:
                continue
            new_pos = []
            for cube in self.gen_cubes(pt):
                if cube == (0,) * dim:
                    continue
                if max_cdist is not None:
                    if sum(abs(cj) for cj in cube) > max_cdist:
                        continue
                rel_pos = tuple(pj - FastFraction(cj, 1) for pj, cj in zip(pt, cube))
                for g in gate_list:
                    if g.entrance == rel_pos:
                        new_pos.append((cube, g))
            result[pt] = new_pos

        return result

    def continue_path(self, path, cubeset=None, finish_pids=None):
        """
        Add one edge to a path, yields paths.

        Params:
            path -      subj
            cubeset -   support of path (optimization!)
        """

        self.stats['continue'] += 1

        head = path.head
        N = self.div
        d = self.dim
        if cubeset is None:
            cubeset = path.support()

        for delta, gate in self.next_dict[path.gate.exit]:
            new_head = tuple(head[j] + delta[j] for j in range(d))
            if any(nj < 0 or nj == N for nj in new_head) or new_head in cubeset:
                continue
            if finish_pids is not None:
                pid = hash((frozenset(cubeset), new_head, gate))
                if pid not in finish_pids:
                    continue
            yield PathNode(head=new_head, gate=gate, prev=path, len=path.len + 1)

    def depth_search(self, paths, max_len, finish_pids):
        """Depth-first search for paths.
        Params:
            paths -     starting paths
            max_len -   subj
            finish_pids -    finish path ids
        """
        todo = [(path, path.support()) for path in paths]
        while todo:
            path, cubeset = todo.pop()
            if path.len == (max_len - 1):
                # путь почти закончен, нужно проверять на соответствие финишу
                yield from self.continue_path(path, cubeset, finish_pids=finish_pids)
            else:
                for np in self.continue_path(path, cubeset):
                    new_cubeset = cubeset.copy()
                    new_cubeset.add(np.head)
                    todo.append((np, new_cubeset))

    # в худшем случае мы удваиваем работу:(
    def width_search(self, paths, max_steps, max_count):
        """Width-search for given starting paths."""
        curr_paths = paths
        for i in range(max_steps):
            new_paths = []
            for path in curr_paths:
                new_paths += self.continue_path(path)
                if len(new_paths) > max_count:
                    break
            if len(new_paths) > max_count:
                break
            curr_paths = new_paths
            logging.info('width_search: step %d, expanded to %d' % (i+1, len(curr_paths)))
        return curr_paths

    def generate_paths(self, **kwargs):
        #TEMPORORAY
        yield from self.generate_paths_generic(self.gates[0], **kwargs)
        return

        path_lists = [list(self.generate_paths_generic(gate, **kwargs)) for gate in self.gates]
        yield from itertools.product(*path_lists)

    def generate_paths_generic(self, gate, start_max_count=100, finish_max_count=10 ** 6):
        """Generate entrance-exit broken line."""
        N = self.div
        d = self.dim
        max_steps = (N**d // 2) - 1

        # COPY-PASTE :( TODO: use self.next_dict
        gate_set = set()
        for bm in BaseMap.gen_base_maps(d):
            for g in self.gates:
                gate_set.add(bm * g)
        gate_list = sorted(gate_set)

        start_init = []
        finish_init = []
        for name, point in [('entrance', gate.entrance), ('exit', gate.exit)]:
            scaled = tuple(xj * FastFraction(N, 1) for xj in point)
            for cube in self.gen_cubes(scaled):
                if any(cj < 0 or cj >= N for cj in cube):
                    continue
                rel = tuple(sj - FastFraction(cj, 1) for sj, cj in zip(scaled, cube))
                for g in gate_list:
                    # here we use that we allow time-rev bms for gate_list
                    if g.entrance == rel:
                        path = PathNode(head=cube, prev=None, len=1, gate=g)
                        if name == 'entrance':
                            start_init.append(path)
                        else:
                            finish_init.append(path)


        start_st = time.time()
        start_paths = self.width_search(start_init, max_steps=max_steps, max_count=start_max_count)
        if not start_paths:
            return

        start = {}
        for path in start_paths:
            pid = hash(path.state())
            start.setdefault(pid, []).append(path)
        logging.info('start: width: %d, configurations: %d', start_paths[0].len - 1, len(start))

        finish_st = time.time()
        finish_paths = self.width_search(finish_init, max_steps=max_steps, max_count=finish_max_count)
        if not finish_paths:
            return
        
        finish = {}
        all_cubes = set(itertools.product(range(N), repeat=d))
        for path in finish_paths:
            complement_cubeset = all_cubes - path.support()
            complement_state = (frozenset(complement_cubeset), path.head, path.gate.reversed())
            pid = hash(complement_state)
            finish.setdefault(pid, []).append(path)
        finish_width = finish_paths[0].len - 1
        logging.info('finish: width %d, configurations: %d', finish_width, len(finish))

        finish_pids = frozenset(finish.keys())

        depth_st = time.time()
        max_len = N**d - finish_width
        found_paths = 0

        for i, pid in enumerate(sorted(start.keys())):
            # вообще говоря, могут быть коллизии из-за hash()
            # но делать ключом state а не hash(state) накладно по памяти
            states = []  # stabilize sort
            state_paths = {}
            for path in start[pid]:
                state = path.state()
                if state not in state_paths:
                    state_paths[state] = []
                    states.append(state)
                state_paths[state].append(path)

            # важна лишь группировка путей по одинаковому state
            for state in states:
                paths = state_paths[state]

                # теперь пытаемся всё склеить в один путь
                for mid_path in self.depth_search([paths[0]], max_len, finish_pids):
                    mid_pid = hash(mid_path.state())
                    for start_path in paths:
                        for fin_path in finish[mid_pid]:
                            path_data = self.glue_paths(start_path, mid_path, fin_path)
                            if path_data is not None:
                                found_paths += 1
                                yield CurvePath.from_data(self.dim, self.div, path_data)

            if self.verbose:
                logging.info(
                    '(%f) %d of %d (%.1f %%): total found: %d',
                    time.time(), i+1, len(start.keys()), 100 * (i+1)/len(start.keys()), found_paths,
                )

        if self.verbose:
            print('max_len:', max_len)
            print('mem:', resource.getrusage(resource.RUSAGE_SELF)[2])
            print('global:', self.stats['continue'])
            print('start time:', finish_st - start_st)
            print('finish time:', depth_st - finish_st)
            print('depth time:', time.time() - depth_st)

    @staticmethod
    def glue_paths(start, mid, fin):
        # сначала проверим, что fin корректно соотносится с mid
        if mid.head != fin.head or mid.gate != fin.gate.reversed():
            return
        fin_past = fin.support() - set([fin.head])
        mid_past = mid.support() - set([mid.head])
        if len(fin_past & mid_past):
            return

        start_data = start.get_data()
        mid_data = mid.get_data()
        fin_data = fin.get_data()
        rev_fin_data = []
        for cube, gate in reversed(fin_data):
            rev_fin_data.append((cube, gate.reversed()))

        return start_data + mid_data[start.len:] + rev_fin_data[1:]
