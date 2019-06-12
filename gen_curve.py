# coding: utf-8

import time
import resource
from collections import namedtuple
import logging
import itertools

import base_map
from fractal_curve import FractalCurve


# head - квадрат конца пути
# entrance - вход в head (относительный, т.е. точка {0,1}^d)
# exit - выход из head (относительный)
# prev - ссылка на путь без последнего квадрата
# len - длина пути (для оптимизации)
class Path(namedtuple('Path', ['head','entrance','exit','prev','len'])):
    def get_data(self):
        data = []
        node = self
        while node is not None:
            data.append((node.head, node.entrance, node.exit))
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

    # состояние, которое определяет продолжимость пути
    # entrance нужен лишь для удобства пересечения с финишем
    def state(self):
        return (frozenset(self.flat()[1:]), self.head, self.entrance, self.exit)


def get_next_dict(dim, hdist, max_cdist=None):
    """Get a dict: exit => [(delta, new_entrance, new_exit), ...]
    Params:
        dim:    subj
        hdist:  arrow changes hdist coordinates
        max_cdist: do not allow cube changes greater than it
    """
    result = {}

    # стартовые позиции - все вершины куба
    start_positions = itertools.product((0, 1), repeat=dim)

    # сдвиги - вектора из {0,1,-1} с ||v||_1 = hdist
    deltas = []
    for delta in itertools.product((0,1,-1), repeat=dim):
        l1norm = sum(abs(dj) for dj in delta)
        if l1norm == hdist:
            deltas.append(delta)

    for start_pos in start_positions:
        for delta in deltas:
            new_pos = tuple(start_pos[j] + delta[j] for j in range(dim))
            # ищем кубы, содержащие и start_pos, и new_pos
            # куб: {(x0,..): cube[j] <= x[j] <= cube[j]+1}
            x_max = [min(start_pos[j], new_pos[j]) for j in range(dim)]
            x_min = [max(start_pos[j], new_pos[j]) - 1 for j in range(dim)]
            ranges = [range(x_min[j], x_max[j] + 1) for j in range(dim)]
            new_cubes = set()
            for cj in itertools.product(*ranges):
                if cj != (0,) * dim:
                    new_cubes.add(cj)

            for new_cube in new_cubes:
                start_rel_pos = tuple(start_pos[j] - new_cube[j] for j in range(dim))
                new_rel_pos = tuple(new_pos[j] - new_cube[j] for j in range(dim))
                if start_pos not in result:
                    result[start_pos] = []
                result[start_pos].append((new_cube, start_rel_pos, new_rel_pos))

    if max_cdist is not None:
        for start_pos, data in result.items():
            result[start_pos] = [(d, entr, exit) for d, entr, exit in data if sum(abs(dj) for dj in d) <= max_cdist]

    return result


class CurveGenerator:
    """Generate vertex curves of given div, with entrance at 0 and exit (1,1,..,1,0,0,..,0)."""
    def __init__(self, dim, div, hdist, verbose=0, max_cdist=None, oriented=False):
        self.dim = dim
        self.div = div
        self.verbose = verbose
        N = self.div

        # вершинную кривую всегда можно привести к виду, когда exit=(1,1,..,1,0,..,0)
        self.exit = (1,) * hdist + (0,) * (dim - hdist)

        # здесь бы использовать continue_path !!!
        start_init = []
        finish_init = []
        for arr_coords in itertools.combinations(range(dim), hdist):
            # положение меняется на координатах arr_coords
            st_arr_0 = (0,) * dim
            st_arr_1 = tuple(1 if j in arr_coords else 0 for j in range(dim))
            st_path = Path(head=(0,) * dim, prev=None, len=1, entrance=st_arr_0, exit=st_arr_1)
            start_init.append(st_path)

            fn_head = (N - 1,) * hdist + (0,) * (dim - hdist)
            fn_arr_0 = (1,) * hdist + (0,) * (dim - hdist)
            fn_arr_1 = tuple(1 - fn_arr_0[j] if j in arr_coords else fn_arr_0[j] for j in range(dim))
            fn_path = Path(head=fn_head, prev=None, len=1, entrance=fn_arr_0, exit=fn_arr_1)
            finish_init.append(fn_path)

        self.start_init = start_init
        self.finish_init = finish_init

        self.next_dict = get_next_dict(dim, hdist, max_cdist)
        self.oriented = oriented

    def get_non_cube(self, reverse=False):
        N = self.div
        d = self.dim
        if reverse:
            return (0,) * d
        else:
            return tuple((N-1) * ej for ej in self.exit)

    def continue_path(self, path, cubeset=None, non_cube=None, finish_pids=None):
        """Add one edge to a path.
        Params:
            path -      subj
            cubeset -   support of path (optimization!)
            non_cube -  avoid this cube
        """
        self.global_counter += 1

        head = path.head
        N = self.div
        d = self.dim
        result = []
        if cubeset is None:
            cubeset = set(path.flat())

        cubeset_id = None
        for delta, entrance, exit in self.next_dict[path.exit]:
            new_head = tuple(head[j] + delta[j] for j in range(d))
            if any(nj < 0 or nj == N for nj in new_head) or new_head in cubeset or new_head == non_cube:
                continue
            do_add = True
            if finish_pids is not None:
                pid = hash((frozenset(cubeset), new_head, entrance, exit))
                if pid not in finish_pids:
                    do_add = False
            if do_add:
                new_path = Path(head=new_head, entrance=entrance, exit=exit, prev=path, len=path.len + 1)
                result.append(new_path)

        return result


    def depth_search(self, paths, max_len, finish_pids):
        """Depth-first search for paths.
        Params:
            paths -     starting paths
            max_len -   subj
            finish_pids -    finish path ids
        """
        todo = paths.copy()
        cubesets = [set(p.flat()) for p in todo]  # кэшируем множество кубов в "активном" пути
        result = []
        non_cube = self.get_non_cube()
        while todo:
            path = todo.pop()
            cubeset = cubesets.pop()
            if path.len == (max_len - 1):
                # путь почти закончен, нужно проверять на соответствие финишу
                rr = self.continue_path(path, cubeset, non_cube, finish_pids=finish_pids)
                result += rr
            else:
                new_paths = self.continue_path(path, cubeset, non_cube)
                for np in new_paths:
                    todo.append(np)  # важно дописывать в конец
                    new_cubeset = cubeset.copy()
                    new_cubeset.add(np.head)
                    cubesets.append(new_cubeset)
        return result

    # в худшем случае мы удваиваем работу:(
    def width_search(self, paths, max_steps, max_count, reverse=False):
        """Width-search for given starting paths."""
        curr_paths = paths
        non_cube = self.get_non_cube(reverse)
        for i in range(max_steps):
            new_paths = []
            for path in curr_paths:
                new_paths += self.continue_path(path, non_cube=non_cube)
                if len(new_paths) > max_count:
                    break
            if len(new_paths) > max_count:
                break
            curr_paths = new_paths
            logging.info('width_search: step %d, expanded to %d' % (i+1, len(curr_paths)))
        return curr_paths

    def generate_brklines(self, start_max_count=100, finish_max_count=10**6):
        """Generate entrance-exit broken line."""
        self.global_counter = 0
        N = self.div
        d = self.dim
        max_steps = (N**d // 2) - 1;
        start_st = time.time()
        start_paths = self.width_search(self.start_init, max_steps=max_steps, max_count=start_max_count)
        if not start_paths:
            return

        start = {}
        for path in start_paths:
            pid = hash(path.state())
            start.setdefault(pid, []).append(path)
        logging.info('start: width: %d, configurations: %d', start_paths[0].len - 1, len(start))

        finish_st = time.time()
        finish_paths = self.width_search(self.finish_init, max_steps=max_steps, max_count=finish_max_count, reverse=True)
        if not finish_paths:
            return
        
        finish = {}
        all_cubes = set(itertools.product(range(N), repeat=d))
        for path in finish_paths:
            complement_cubeset = all_cubes - set(path.flat())
            complement_state = (frozenset(complement_cubeset), path.head, path.exit, path.entrance)
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
                mid_paths = self.depth_search([paths[0]], max_len, finish_pids)

                # теперь пытаемся всё склеить в один путь
                for start_path in paths:
                    for mid_path in mid_paths:
                        mid_pid = hash(mid_path.state())
                        for fin_path in finish[mid_pid]:
                            found_paths += 1
                            yield glue_paths(start_path, mid_path, fin_path)

            if self.verbose:
                logging.info(
                    '(%f) %d of %d (%.1f %%): found %d, total: %d',
                    time.time(), i+1, len(start.keys()), 100 * (i+1)/len(start.keys()), len(mid_paths), found_paths,
                )

        if self.verbose:
            print('max_len:', max_len)
            #print('result:', len(result_paths))
            print('mem:', resource.getrusage(resource.RUSAGE_SELF)[2])
            print('global:', self.global_counter)
            print('start time:', finish_st - start_st)
            print('finish time:', depth_st - finish_st)
            print('depth time:', time.time() - depth_st)

    def generate_curves(self, **kwargs):
        for brkline in self.generate_brklines(**kwargs):
            for curve in self.generate_curves_for_brkline(brkline):
                yield curve

    def generate_curves_for_brkline(self, brkline):
        global_entrance = (0,) * self.dim
        global_exit = self.exit
        bms_variants = []
        proto = []
        for cube, entrance, exit in brkline:
            constr = {
                global_entrance: entrance,
                global_exit: exit,
            }
            bms_for_cube = base_map.constraint_base_maps(self.dim, constr)
            if self.oriented:
                bms_for_cube = [bm for bm in bms_for_cube if bm.is_oriented()]
            bms_variants.append(bms_for_cube)
            proto.append(cube)

        for bms in itertools.product(*bms_variants):
            curve = FractalCurve(dim=self.dim, div=self.div, proto=proto, base_maps=bms)
            yield curve
            

def glue_paths(start, mid, fin):
    # сначала проверим, что fin корректно соотносится с mid
    if mid.head != fin.head or mid.entrance != fin.exit or mid.exit != fin.entrance:
        return []
    fin_past = set(fin.flat()[1:])
    mid_past = set(mid.flat()[1:])
    if len(fin_past & mid_past):
        return []

    start_data = start.get_data()
    mid_data = mid.get_data()
    fin_data = fin.get_data()
    rev_fin_data = []
    for cube, entrance, exit in reversed(fin_data):
        rev_fin_data.append((cube, exit, entrance))

    return start_data + mid_data[start.len:] + rev_fin_data[1:]
