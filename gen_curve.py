# coding: utf-8

import time
import resource
from collections import namedtuple
import logging

# head - квадрат конца пути
# arr - положение стрелки от входа до выхода в квадрате (A=00, B=10, C=11, D=01, напр.: arr='AB')
# prev - ссылка на путь без последнего квадрата
# len - длина пути (для оптимизации)
class Path(namedtuple('Path', ['head','prev','len','arr'])):

    def get_data(self):
        data = []
        node = self
        while node is not None:
            data.append((node.head, node.arr))
            node = node.prev
        data.reverse()
        return data

    def flat(self):
        length = self.len
        path = self
        squares = [None] * length
        for i in range(length):
            squares[i] = path.head
            path = path.prev
        return squares

    # состояние, которое определяет продолжимость пути
    # arr[0] нужен лишь для удобства пересечения с финишем
    def state(self):
        return (frozenset(self.flat()[1:]), self.head, self.arr)

    def pprint(self):
        path = self
        length = self.len
        nodes = []
        for i in range(length):
            nodes.append(path)
            path = path.prev
        nodes.reverse()
        print(length, ' '.join([str(n.head)+'['+n.arr+']' for n in nodes]))


# словарь arr[1] => [(delta_square, new_arr), ...]
def get_next_dict(arr_type, allow_vertex_transit):
    letter_to_position = {'A': (0, 0), 'B': (1, 0), 'C': (1, 1), 'D': (0, 1)}
    position_to_letter = {v: k for k, v in letter_to_position.items()}
    result = {}

    if arr_type == 'edge':
        deltas = [(1,0), (0,1), (-1,0), (0,-1)]
    else:
        deltas = [(1,1), (-1,1), (-1,-1), (1,-1)]
    for letter in 'ABCD':
        start_pos = letter_to_position[letter]
        for d in deltas:
            new_pos = (start_pos[0] + d[0], start_pos[1] + d[1])
            # ищем квадраты, содержащие и start_pos, и new_pos
            # квадрат sq: {(x,y): sq[0] <= x <= sq[0]+1,  sq[1] <= y <= sq[1]+1}
            sx_up = min(start_pos[0], new_pos[0])
            sx_low = max(start_pos[0], new_pos[0]) - 1

            sy_up = min(start_pos[1], new_pos[1])
            sy_low = max(start_pos[1], new_pos[1]) - 1

            new_squares = set()
            for sx in range(sx_low, sx_up + 1):
                for sy in range(sy_low, sy_up + 1):
                    if (sx, sy) != (0, 0):
                        new_squares.add((sx, sy))

            for new_square in new_squares:
                start_rel_pos = (start_pos[0] - new_square[0], start_pos[1] - new_square[1])
                new_rel_pos = (new_pos[0] - new_square[0], new_pos[1] - new_square[1])
                arr = position_to_letter[start_rel_pos] + position_to_letter[new_rel_pos]
                if letter not in result:
                    result[letter] = []
                result[letter].append((new_square, arr))

    if not allow_vertex_transit:
        for letter in 'ABCD':
            result[letter] = [(d, a) for d, a in result[letter] if abs(d[0]) + abs(d[1]) <= 1]

    return result


class CurveGenerator:
    """Generate 2D curves of given div, with entrance at (0,0) and given exit (1,1) or (1,0)."""
    def __init__(self, div, exit, verbose=0, allow_vertex_transit=False):
        self.div = div
        N = self.div
        self.verbose = verbose

        assert exit == (1, 1) or exit == (1, 0)

        if exit == (1, 0):
            self.is_edge = True
            curve_type = 'edge'
            self.start_init = [
                Path(head=(0,0), prev=None, len=1, arr='AB'),
                Path(head=(0,0), prev=None, len=1, arr='AD'),
            ]
            self.finish_init = [
                Path(head=(N-1,0), prev=None, len=1, arr='BC'),
                Path(head=(N-1,0), prev=None, len=1, arr='BA'),
            ]
        else:
            self.is_edge = False
            curve_type = 'diag'
            self.start_init = [
                Path(head=(0,0), prev=None, len=1, arr='AC'),
            ]
            self.finish_init = [
                Path(head=(N-1,N-1), prev=None, len=1, arr='CA'),
            ]

        self.next_dict = get_next_dict(curve_type, allow_vertex_transit)


    def get_non_sq(self, reverse=False):
        N = self.div
        if reverse:
            return (0,0)
        elif self.is_edge:
            return (N-1,0)
        else:
            return (N-1,N-1)

    def continue_path(self, path, sqset=None, non_sq=None, finish_pids=None):
        """Add one edge to a path.
        Params:
            path -      subj
            sqset -     support of path (optimization!)
            non_sq -    avoid this square
        """
        self.global_counter += 1

        sq = path.head
        N = self.div
        result = []
        if sqset is None:
            sqset = set(path.flat())

        sqset_id = None
        for d, a in self.next_dict[path.arr[1]]:
            n = (sq[0] + d[0], sq[1] + d[1])  # new head
            if n[0] < 0 or n[1] < 0 or n[0] == N or n[1] == N or n in sqset or n == non_sq:
                continue
            do_add = True
            if finish_pids is not None:
                pid = hash((frozenset(sqset), n, a))
                if pid not in finish_pids:
                    do_add = False
            if do_add:
                new_path = Path(head=n, arr=a, prev=path, len=path.len + 1)
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
        sqsets = [set(p.flat()) for p in todo]  # кэшируем множество квадратов в "активном" пути
        result = []
        non_sq = self.get_non_sq()
        while todo:
            path = todo.pop()
            sqset = sqsets.pop()
            if path.len == (max_len - 1):
                # путь почти закончен, нужно проверять на соответствие финишу
                rr = self.continue_path(path, sqset, non_sq, finish_pids=finish_pids)
                result += rr
            else:
                new_paths = self.continue_path(path, sqset, non_sq)
                for np in new_paths:
                    todo.append(np)  # важно дописывать в конец
                    new_sqset = sqset.copy()
                    new_sqset.add(np.head)
                    sqsets.append(new_sqset)
        return result

    def width_search(self, paths, steps, reverse=False):
        """Width-search for given starting paths."""
        curr_paths = paths
        non_sq = self.get_non_sq(reverse)
        for i in range(steps):
            new_paths = []
            for path in curr_paths:
                new_paths += self.continue_path(path, non_sq=non_sq)
            curr_paths = new_paths
            logging.warning('width_search: step %d, expanded to %d' % (i+1, len(curr_paths)))
        return curr_paths

    def generate_brkline(self, start_width, finish_width):
        """Generate entrance-exit broken line."""
        self.global_counter = 0
        start_st = time.time()
        start_paths = self.width_search(self.start_init, start_width)
        start = {}
        for path in start_paths:
            pid = hash(path.state())
            start.setdefault(pid, []).append(path)
        logging.warning('start: width: %d, configurations: %d', start_width, len(start))

        finish_st = time.time()
        finish_paths = self.width_search(self.finish_init, finish_width, reverse=True)
        
        finish = {}
        N = self.div
        all_sq = set([(i,j) for i in range(N) for j in range(N)])
        for path in finish_paths:
            rev_arr = path.arr[1] + path.arr[0]
            complement_sqset = all_sq - set(path.flat())
            complement_state = (frozenset(complement_sqset), path.head, rev_arr)
            pid = hash(complement_state)
            finish.setdefault(pid, []).append(path)
        logging.info('finish: width %d, configurations: %d', finish_width, len(finish))

        finish_pids = frozenset(finish.keys())

        result_paths = []
        depth_st = time.time()
        max_len = N**2 - finish_width
        print('max_len:', max_len)

        for i, pid in enumerate(start.keys()):
            # вообще говоря, могут быть коллизии из-за hash()
            # но делать ключом state а не hash(state) накладно по памяти
            state_paths = {}
            for path in start[pid]:
                state_paths.setdefault(path.state(), []).append(path)

            # важна лишь группировка путей по одинаковому state
            for state, paths in state_paths.items():
                mid_paths = self.depth_search([paths[0]], max_len, finish_pids)

                # теперь пытаемся всё склеить в один путь
                for start_path in paths:
                    for mid_path in mid_paths:
                        mid_pid = hash(mid_path.state())
                        for fin_path in finish[mid_pid]:
                            result_paths.append(glue_paths(start_path, mid_path, fin_path))

            if self.verbose:
                logging.info(
                    '(%f) %d of %d (%.1f %%): found %d, total: %d',
                    time.time(), i+1, len(start.keys()), 100 * (i+1)/len(start.keys()), len(mid_paths), len(result_paths),
                )

        if self.verbose:
            print('result:', len(result_paths))
            print('mem:', resource.getrusage(resource.RUSAGE_SELF)[2])
            print('global:', self.global_counter)
            print('start time:', finish_st - start_st)
            print('finish time:', depth_st - finish_st)
            print('depth time:', time.time() - depth_st)

        return result_paths


def glue_paths(start, mid, fin):
    # сначала проверим, что fin корректно соотносится с mid
    if mid.head != fin.head:
        return []
    if mid.arr != fin.arr[1] + fin.arr[0]:
        return []
    fin_past = set(fin.flat()[1:])
    mid_past = set(mid.flat()[1:])
    if len(fin_past & mid_past):
        return []

    start_data = start.get_data()
    mid_data = mid.get_data()
    fin_data = fin.get_data()

    return start_data + mid_data[start.len:] + fin_data[1:]


# assert: 4 edge(?) => 298
# assert: 5 edge => 49700
# 5-6-11  global: 1.8m
#   dict: mem: 500, time: 0+5+12
#   namedtuple: mem: 363, time: 0+6+15 
#   namedtuple + hash: mem 170
#   -globalcounter: time: 0+5+12
#   func->dict: 0+4.5+10.5
#   all_todo_sqsets: 0+4.5+8.5

# diag: 5 => 2592
