import operator
from collections import Counter, namedtuple
from collections.abc import Sized
from heapq import heappop, heappush
import logging

from .utils import get_lcm, get_int_cube_with_cache, get_int_time_with_cache

from .base_maps import BaseMap, Spec
from .fast_fractions import FastFraction
from . import sat_adapters
from .curves import Curve


class CurvePiecePosition:
    """
    Information about time/space position of a fraction in a curve.

    A support of the fraction is the last cube in the nested sequence of cubes.
    """

    def __init__(self, dim, div, cnums, cubes):
        """
        Init CurvePiecePosition instance.
        Params:
        dim, div  --  the same as in curve
        cnums  --  sequence of cnums of nested sequence
        cubes  --  sequence of integer coords of nested sequence

        Each cnum corresponds to the curve in the corresponding fraction.
        """
        self.dim = dim
        self.div = div
        self.cnums = tuple(cnums)
        self.cubes = tuple(cubes)

        self.depth = len(self.cnums)
        self.sub_div = div**self.depth
        self.sub_genus = self.sub_div**dim

    def specify(self, cnum, cube):
        """Create sub-fraction position."""
        return type(self)(
            dim=self.dim,
            div=self.div,
            cnums=self.cnums + (cnum,),
            cubes=self.cubes + (cube,),
        )

    def get_int_coords(self):
        """
        Natural integer coordinates - time and lower-left cube vertex
        returns l, x, t:
        l - depth
        x - int cube: cj <= xj <= cj+1, where abs cube: cj/N^l <= xj <= (cj+1)/N^l
        t - int time: [t, t+1], abs time: [t/G^l, (t+1)/G^l]

        These may be viewed as coordinates in curve [0,G^l] -> [0,N^l]^d
        """
        return (
            self.depth,
            get_int_cube_with_cache(self.dim, self.div, self.cubes),
            get_int_time_with_cache(self.dim, self.div, self.cnums),
        )


class CurvePiece:
    """
    Fraction of a curve.

    Defined as triple (curve, pattern number, position).
    The specs in curve must be specified in all cubes of the position, except last one.
    """

    def __init__(self, curve, pnum, pos):
        self.curve = curve
        self.pnum = pnum
        self.pos = pos

    # only for full fractal curves - ?? later
    def get_last_map(self):
        curve = self.curve
        dim, G = curve.dim, curve.genus
        last_map = BaseMap.id_map(dim)
        for cnum in self.pos.cnums:
            cnum = last_map.apply_cnum(G, cnum)
            last_map = last_map * curve.specs[cnum].base_map  # именно в таком порядке!
        return last_map

    def divide(self):
        """
        Divide the fraction in all possible ways.
        """
        curve = self.curve
        dim, div, G = curve.dim, curve.div, curve.genus

        # define orientation of last_but_one fraction of a curve
        # in the last fraction we do not know the orientation yet!
        prev_spec = Spec(BaseMap.id_map(dim), self.pnum)
        for cnum in self.pos.cnums[:-1]:
            prev_spec = self.curve.compose_spec(prev_spec, cnum)

        prev_curve = prev_spec * self.curve

        active_pnum = prev_curve.pnum
        active_cnum = self.pos.cnums[-1]  # the cube that will be divided
        spec_cnum = prev_spec.base_map.apply_cnum(G, active_cnum)

        for sp in prev_curve.gen_allowed_specs(active_pnum, active_cnum):
            # see apply_cube_map in curves.py:
            # prev_curve.specs[active_cnum] = bm  =>  orig_curve.specs[spec_cnum] = ...
            new_spec = sp.conjugate_by(~prev_spec.base_map)
            specified_curve = curve.specify(active_pnum, spec_cnum, new_spec)

            # last_curve = sp * prev_curve
            # last_curve_proto = (sp * prev_curve).proto, we do not want to use multiplication
            last_curve_proto = sp.base_map * prev_curve.patterns[sp.pnum].proto
            for cnum, cube in enumerate(last_curve_proto):
                new_pos = self.pos.specify(cnum, cube)
                new_piece = CurvePiece(specified_curve, self.pnum, new_pos)
                yield new_piece


class CurvePieceBalancedPair:
    """
    Pair or curve fractions with almost equal depth.

    We ensure that depth1 == depth2 or depth1 == depth2 + 1.
    """

    def __init__(self, curve, junc, pos1, pos2):
        """
        Attributes:
        curve, junc  --  subj
        pos1, pos2  --  position in patterns of curve (not rotated!)
        """
        self.curve = curve
        self.junc = junc
        self.pos1 = pos1
        self.pos2 = pos2
        self.piece1 = CurvePiece(self.curve, junc.spec1.pnum, self.pos1)
        self.piece2 = CurvePiece(self.curve, junc.spec2.pnum, self.pos2)

    def divide(self):
        # use curve from divided piece because it has specified curve
        if self.pos1.depth > self.pos2.depth:
            for subpiece in self.piece2.divide():
                yield type(self)(subpiece.curve, self.junc, self.pos1, subpiece.pos)
        else:
            for subpiece in self.piece1.divide():
                yield type(self)(subpiece.curve, self.junc, subpiece.pos, self.pos2)


class Threshold:
    """
    Generic thresholding.

    Examples:
    * Threshold('<', 2):  thr.check(x) <=> x < 2
    """

    def __init__(self, sign, threshold):
        self.threshold = threshold
        if sign == '<':
            op = operator.lt
        elif sign == '<=':
            op = operator.le
        elif sign == '>':
            op = operator.gt
        elif sign == '>=':
            op = operator.ge
        else:
            raise Exception("Unknown comparison method")
        self.operator = op

    def check(self, x):
        """Check if x is compliant with threshold."""
        t = self.threshold
        return self.operator(x, t)


class RichPairsTree:
    """
    RichPairsTree - structure which holds pairs of curve fractions.

    Main element occurring in algorithms for bounding curve ratio is a pair or curve fractions.
    The pairs are represented by the class CurvePieceBalancedPair.
    Here we hold pairs with their up/lo bounds and keep some invariants.

    Attributes:
    .nodes  -  active nodes (rich pairs), stored in a heap (priority is given externally)
    .new_pairs  -  pairs waiting to be processed, with unknown up/lo
    .bad_pairs  -  temporary storage for "bad" pairs with bad ratio
    .good_threshold  -  if up < good, pair is considered "good" (hence insignificant); class Threshold
    .bad_threshold  -  if lo > bad, pair is considered "bad"; class Threshold
    .max_lo_node  -  node with maximal seen lo
    """

    RichPair = namedtuple('RichPair', [
        'sort_field',  # first field, will used for comparison by heapq
        'pair',  # CurvePieceBalancedPair
        'up',  # upper_bound
        'lo',  # lower_bound
        'data',  # some arbitrary data
    ])

    def __init__(self, pairs):
        self.nodes = []  # use heapq algorithm for plain list
        self.new_pairs = pairs
        self.bad_pairs = []
        self.good_threshold = None
        self.bad_threshold = None
        self.max_lo_node = None
        self.stats = Counter()
        self._inc = 0

    def set_good_threshold(self, threshold):
        self.good_threshold = threshold

    def set_bad_threshold(self, threshold):
        self.bad_threshold = threshold

    def add_pair(self, pair, lo, up, priority, data):
        """Add node keeping the invariants."""
        self.stats['adds'] += 1

        # the order of checks may be important; as in current algorithms
        # we add bad pairs to SAT clauses, it may be more effective to
        # check for goodness first ...

        if self.good_threshold is not None and self.good_threshold.check(up):
            self.stats['good'] += 1
            return

        if self.bad_threshold is not None and self.bad_threshold.check(lo):
            self.bad_pairs.append(pair)
            self.stats['bad'] += 1
            return

        self._inc += 1  # to restrict comparison in heap to 'sort_field' of RichPair
        node = self.RichPair(sort_field=(-priority, self._inc), pair=pair, lo=lo, up=up, data=data)

        if self.max_lo_node is None or lo > self.max_lo_node.lo:
            self.max_lo_node = node

        heappush(self.nodes, node)

    def divide(self):
        """
        Divide node of the tree with highest priority.
        """
        if self.nodes:
            worst_node = heappop(self.nodes)
            self.new_pairs += worst_node.pair.divide()

    def copy(self):
        assert not self.bad_pairs
        assert self.good_threshold is None
        assert self.bad_threshold is None
        new_tree = RichPairsTree([])
        new_tree.new_pairs = self.new_pairs.copy()
        new_tree.nodes = self.nodes.copy()
        new_tree.max_lo_node = self.max_lo_node
        new_tree._inc = self._inc
        return new_tree


class RunOutOfIterationsException(Exception):
    pass

class Estimator:
    """
    Estimator - estimates curve maximal ratio.

    This is main class for curve ratio estimation.
    It orchestrates work of helper classes: RichPairsTree, CurveBalancedPairs and others.
    """

    def __init__(self, ratio_func, cache_max_size=2**20):
        """
        Init Estimator instance and set some basic properties.

        Params:
        ratio_func  --  function (dim, dx, dt) -> FastFraction, it is assumed to be d-uniform
        cache_max_size  --  subj for pairs bounds cache
        """

        self.ratio_func = ratio_func
        self.stats = Counter()
        self._get_bounds_cache = {}
        self._cache_max_size = cache_max_size

    def get_bounds(self, pair, brkline=None):
        """
        Get lower and upper bounds for max ratio of given fractions pair.

        brkline  --  instance of IntegerBrokenLine class
        Return triple (lo, up, argmax), argmax only for brkline
        """
        dim = pair.curve.dim
        N = pair.curve.div
        G = pair.curve.genus

        pos1, pos2 = pair.pos1, pair.pos2

        use_cache = (brkline is None)
        if use_cache:
            # do not use caching for brklines, maybe later
            cache = self._get_bounds_cache
            cache_key = (dim, N, pair.junc, pos1.cnums, pos1.cubes, pos2.cnums, pos2.cubes)
            if cache_key in cache:
                self.stats['get_bounds_cache_hit'] += 1
                return cache[cache_key]
            else:
                self.stats['get_bounds_cache_miss'] += 1

        # these are integer positions in original curve patterns
        l1, x1, t1 = pos1.get_int_coords()
        l2, x2, t2 = pos2.get_int_coords()

        use_brkline = (brkline is not None)
        if use_brkline:
            # ломаные нужно поворачивать:(
            brk1_bm = pair.piece1.get_last_map()
            brk2_bm = pair.piece2.get_last_map()

        junc = pair.junc

        # junc: apply base_maps to coordinates
        x1 = junc.spec1.base_map.apply_cube(pos1.sub_div, x1)
        t1 = junc.spec1.base_map.apply_cnum(pos1.sub_genus, t1)
        if use_brkline: brk1_bm = junc.spec1.base_map * brk1_bm

        x2 = junc.spec2.base_map.apply_cube(pos2.sub_div, x2)
        t2 = junc.spec2.base_map.apply_cnum(pos2.sub_genus, t2)
        if use_brkline: brk2_bm = junc.spec2.base_map * brk2_bm

        # common scale
        if l1 == l2:
            mx2 = mt2 = 1
        elif l1 == l2 + 1:
            x2 = [x2j * N for x2j in x2]
            t2 *= G
            mx2 = N
            mt2 = G
        else:
            raise Exception("Bad coordinates!")

        mx = pos1.sub_div
        mt = pos1.sub_genus

        # now we have the following integer coordinates:
        # cube1: x1j <= xj <= x1j + 1    -- cube inside [0, mx]^d
        # cube2: x2j <= xj <= x2j + mx2  -- cube inside [0, mx2 * mx]^d, + shift junc_dx * mx
        #
        # time1: t1 <= t <= t1 + 1
        # time2: t2 <= t <= t2 + mt2,  after shift: t2 + junc_dt * mt <= t <= t2 + mt2 + junc_dt * mt

        # junc: shifts
        t2 += junc.delta_t * mt
        x2 = [x2j + junc_dxj * mx for x2j, junc_dxj in zip(x2, junc.delta_x)]

        max_dx = [max(abs(x1j - x2j + 1), abs(x1j - x2j - mx2)) for x1j, x2j in zip(x1, x2)]

        max_dt = t2 + mt2 - t1  # max(t_2 - t_1)
        min_dt = t2 - (t1 + 1)  # min(t_2 - t_1)

        lo = self.ratio_func(dim, max_dx, max_dt)
        up = self.ratio_func(dim, max_dx, min_dt)

        argmax = None
        if use_brkline:
            brk_mx = brkline.mx
            brk_mt = brkline.mt
            brkline1 = [(brk1_bm.apply_x(x, mx=brk_mx), brk1_bm.apply_t(t, mt=brk_mt)) for x, t in brkline.points]
            brkline2 = [(brk2_bm.apply_x(x, mx=brk_mx), brk2_bm.apply_t(t, mt=brk_mt)) for x, t in brkline.points]

            t1 *= brk_mt
            t2 *= brk_mt
            x1 = [xj * brk_mx for xj in x1]
            x2 = [xj * brk_mx for xj in x2]
            for x1rel, t1rel in brkline1:
                t1_point = t1 + t1rel
                x1_point = [x1j + x1relj for x1j, x1relj in zip(x1, x1rel)]
                for x2rel, t2rel in brkline2:
                    t2_point = t2 + t2rel * mt2
                    x2_point = [x2j + x2relj * mx2 for x2j, x2relj in zip(x2, x2rel)]

                    dx = [x1j - x2j for x1j, x2j in zip(x1_point, x2_point)]
                    dt = t2_point - t1_point

                    lo_point = self.ratio_func(dim, dx, dt)
                    if lo_point > lo:
                        lo = lo_point
                        x1_real = [FastFraction(x1j, mx * brk_mx) for x1j in x1_point]
                        x2_real = [FastFraction(x2j, mx * brk_mx) for x2j in x2_point]
                        t1_real = FastFraction(t1_point, mt * brk_mt)
                        t2_real = FastFraction(t2_point, mt * brk_mt)
                        argmax = {'x1': x1_real, 't1': t1_real, 'x2': x2_real, 't2': t2_real, 'junc': junc}

        if use_cache:
            if len(cache) == self._cache_max_size:
                # poor man's LRU cache :(
                cache.clear()
                self.stats['get_bounds_cache_cleanup'] += 1
            cache[cache_key] = (lo, up, None)

        return lo, up, argmax

    def bound_new_pairs(self, tree, brkline=None):
        """
        Process new pairs of the tree by estimating their ratio.

        Here we get lo/up bounds for the ratio and define priority for a pair.
        """

        for pair in tree.new_pairs:
            lo, up, argmax = self.get_bounds(pair, brkline=brkline)
            priority = up
            tree.add_pair(pair, priority=priority, lo=lo, up=up, data={'argmax': argmax})
        tree.new_pairs = []

    def estimate_ratio(self, curve, *args, **kwargs):
        if isinstance(curve, Curve):
            return self.estimate_ratio_regular(curve, *args, **kwargs)
        else:
            return self.estimate_ratio_fuzzy(curve, *args, **kwargs)

    @staticmethod
    def get_piece_position(curve, pnum, cnum):
        return CurvePiecePosition(
            dim=curve.dim,
            div=curve.div,
            cnums=[cnum],
            cubes=[curve.patterns[pnum].proto[cnum]],
        )

    @classmethod
    def init_pairs_tree(cls, curve):
        """Create initial pairs tree from a curve."""
        pairs = []
        G = curve.genus
        for junc in curve.gen_auto_junctions():
            for cnum1 in range(G):
                for cnum2 in range(cnum1 + 2, G):
                    pos1 = cls.get_piece_position(curve, junc.spec1.pnum, cnum1)
                    pos2 = cls.get_piece_position(curve, junc.spec2.pnum, cnum2)
                    pairs.append(CurvePieceBalancedPair(curve, junc, pos1, pos2))

        # all possible juncs -- coincides with gen_junctions for regular curves
        for junc in curve.get_junctions_info().keys():
            last_cnum1 = 0 if junc.spec1.base_map.time_rev else G - 1
            first_cnum2 = G - 1 if junc.spec2.base_map.time_rev else 0
            for cnum1 in range(G):
                for cnum2 in range(G):
                    if (cnum1, cnum2) == (last_cnum1, first_cnum2):
                        continue
                    pos1 = cls.get_piece_position(curve, pnum=junc.spec1.pnum, cnum=cnum1)
                    pos2 = cls.get_piece_position(curve, pnum=junc.spec2.pnum, cnum=cnum2)
                    pairs.append(CurvePieceBalancedPair(curve, junc, pos1, pos2))

        return RichPairsTree(pairs)

    def estimate_path(self, path):
        points = []
        entr = tuple(cj + gj for cj, gj in zip(path.proto[0], path.gates[0][0]))
        points.append((0, entr))
        for cnum, cube, gate in zip(range(len(path.proto)), path.proto, path.gates):
            exit = tuple(cj + gj for cj, gj in zip(cube, gate[1]))
            points.append((cnum + 1, exit))

        maxr = None
        for t1, x1 in points:
            for t2, x2 in points:
                if t2 <= t1:
                    continue
                dt = t2 - t1
                dx = tuple(x2j - x1j for x1j, x2j in zip(x1, x2))
                r = self.ratio_func(path.proto.dim, dx, dt)
                if maxr is None or r > maxr:
                    maxr = r
        return maxr

    def estimate_ratio_regular(self, curve, rel_tol_inv=100, max_iter=None, use_vertex_brkline=False, verbose=False):
        """
        Estimate maximal ratio for a regular peano curve (class Curve).

        We maintain a tree of pairs of all non-adjacent curve fractions, for all junctions.
        For each pair we get lower and upper bounds for maximal ratio.
        At each iteration we divide the worst pair (with max upper bound).

        Return value is the dict with keys:
        'lo':  lower bound for the curve: ratio(curve) >= lo
        'up':  upper bound for the curve: ratio(curve) <= up
        'argmax':  pair of points where lo is achieved (if use_vertex_brkline is set)
        """
        if use_vertex_brkline:
            if curve.pattern_count > 1:
                raise NotImplementedError("Brklines for multiple patterns not implemented!")
            vertex_brkline = [(x, t) for x, t in curve.get_vertex_moments().items()]
            brkline = IntegerBrokenLine(curve.dim, vertex_brkline)
        else:
            brkline = None

        pairs_tree = self.init_pairs_tree(curve)
        self.bound_new_pairs(pairs_tree, brkline=brkline)

        # invariant: ratio of the curve is in [curr_lo, curr_up]

        # if the was a pair with lower bound lo, then lo is the bound for the whole curve
        # and we do not need to consider pairs with less-or-equal upper bound
        curr_lo = pairs_tree.max_lo_node.lo
        argmax = pairs_tree.max_lo_node.data['argmax']
        pairs_tree.set_good_threshold(Threshold('<=', curr_lo))

        # since the bound of the curve is attended at one of the active pairs,
        # upper bound for worst pair gives us upper bound for the curve
        # we use here that:
        # * pairs_tree.new_pairs is empty after bound_new_pairs, so all pairs are bounded
        # * priority is "up" and nodes are stored in a heap, so the worst node is [0]
        curr_up = pairs_tree.nodes[0].up
        if verbose:
            print('start bounds: ', curr_lo, curr_up)

        tolerance = FastFraction(rel_tol_inv + 1, rel_tol_inv)
        iter_no = 0
        while curr_up > curr_lo * tolerance:
            iter_no += 1
            if not pairs_tree.nodes:
                break
            if max_iter is not None and iter_no > max_iter:
                # ok, leave current estimates
                break

            pairs_tree.divide()
            self.bound_new_pairs(pairs_tree, brkline=brkline)

            node = pairs_tree.max_lo_node
            if node.lo > curr_lo:
                curr_lo = node.lo
                argmax = node.data['argmax']
                if verbose:
                    print('new lower bound: ', curr_lo, curr_up)
            pairs_tree.set_good_threshold(Threshold('<=', curr_lo))

            new_up = pairs_tree.nodes[0].up
            if new_up < curr_up:
                if verbose:
                    print('new upper bound: ', curr_lo, new_up)
                curr_up = new_up

        if verbose:
            print('Pairs tree stats:', pairs_tree.stats)
        res = {'up': curr_up, 'lo': curr_lo}
        if argmax is not None:
            res['argmax'] = argmax

        return res

    def estimate_ratio_fuzzy(self, curve, rel_tol_inv=1000, upper_bound=None,
                             start_lower_bound=None, start_upper_bound=None, start_curve=None, start_pairs_tree=None,
                             max_iter=None, sat_strategy=None, verbose=False):
        """
        Estimate minimal ratio of a fuzzy curve.

        Arguments:
        curve  --  fuzzy curve
        rel_tol_inv  --  inverted relative tolerance, integer
        upper_bound  --  do not proceed if there is ratio is higher than this
        start_lower_bound  --  known lo bound for ratio, start bisection with it
        start_upper_bound  --  known up bound for ratio, start bisection with it
        start_curve  --  known curve with start_upper_bound
        start_pairs_tree  --  just init_pairs_tree
        max_iter  --  subj
        sat_strategy  --  passed to test_ratio_fuzzy, see there (TODO: use kwargs for that)
        find_model  --  bool, will find Curve if set True
        verbose  --  subj

        Return dict with keys:
        'lo' -- lower_bound
        'up' -- upper_bound
        'curve' -- curve_example with ratio in [lo, up]
        """

        # this method is simply "bisection" algorithm based on test_ratio_fuzzy

        stats = Counter()
        # start lower bound: it would be profitable to use good theoretical
        # bounds like 5**2 for ratio_l2_squared, dim=2, pattern_count=1 (?)
        if start_lower_bound is None:
            curr_lo = FastFraction(0, 1)
        else:
            curr_lo = start_lower_bound

        if start_upper_bound is None:
            # got from first regular curve
            curve0 = next(curve.gen_possible_curves())
            curr_up = self.estimate_ratio_regular(curve0, rel_tol_inv=rel_tol_inv)['up']
            curr_curve = curve0
        else:
            curr_up = start_upper_bound
            curr_curve = start_curve

        if start_pairs_tree is None:
            pairs_tree = self.init_pairs_tree(curve)
        else:
            pairs_tree = start_pairs_tree.copy()

        # invariant: best curve in the class is in [curr_lo, curr_up]
        # curr_curve ratio also in [curr_lo, curr_up]
        tolerance = FastFraction(rel_tol_inv + 1, rel_tol_inv)
        while curr_up > curr_lo * tolerance:
            if curr_lo == FastFraction(0, 1):
                # optimization ?
                new_lo = FastFraction(1, 2) * curr_up
                new_up = FastFraction(2, 3) * curr_up
            else:
                new_lo = FastFraction(2, 3) * curr_lo + FastFraction(1, 3) * curr_up
                new_up = FastFraction(1, 3) * curr_lo + FastFraction(2, 3) * curr_up
            stats['bisect_iter'] += 1
            logging.warning(
                '#%d. best in: [%.5f, %.5f]; seek with thresholds: [%.5f, %.5f]', stats['bisect_iter'],
                curr_lo, curr_up, new_lo, new_up,
            )
            logging.debug('precise test thresholds: %s, %s', new_lo, new_up)
            try:
                # we always ask for a model, to get actual curve with guarantees
                test_result = self.test_ratio_fuzzy(
                    curve,
                    bad_threshold=Threshold('>=', new_lo),
                    good_threshold=Threshold('<=', new_up),
                    max_iter=max_iter,
                    sat_strategy=sat_strategy,
                    start_pairs_tree=pairs_tree,
                    verbose=verbose,
                )
                stats.update(test_result['stats'])
            except RunOutOfIterationsException:
                # run out of iterations
                break

            if test_result.get('curve'):
                curr_curve = test_result['curve']
                # try to estimate ratio for curr_curve ???
                curr_up = new_up
            else:
                curr_lo = new_lo

            if upper_bound is not None and curr_lo > upper_bound:
                break

        return {
            'curve': curr_curve,
            'pairs_tree': pairs_tree,
            'lo': curr_lo,
            'up': curr_up,
            'stats': stats,
        }

    def test_ratio_fuzzy(self, curve, bad_threshold, good_threshold,
                         max_iter=None, sat_strategy=None, verbose=False,
                         start_pairs_tree=None):
        """
        Test if there is a "good" curve.

        If the method returns no curve, then ratio is "bad" for all regular curves for given fuzzy curve.
        If the method returns a curve, then it is good one.

        Arguments:
        curve  --  FuzzyCurve instance
        bad_threshold  --  Threshold, bad_threshold.check(ratio)==True means that curve is "bad",
                            so it is Threshold('>', lo) or Threshold('>=', lo)
        good_threshold  --  Threshold, good_threshold.check(ratio)==True means that curve is "good",
                            so it is Threshold('<', up) or Threshold('<=', up)
                            (It is required that bad_threshold is less than good_threshold, so
                            it is possible that there are good curves, but all curves are bad.
                            It this case the return value is not specified.)
        max_iter  --  max number of divisions; raise Exception if maximum iterations reached
        sat_strategy  --    when do we call sat solver:
                            strategy['type'] == 'equal': call every strategy['count'] divisions
                            strategy['type'] == 'geometric': call on strategy['multiplier']**k iterations
        find_model  --  subj
        start_pairs_tree  --  just cache init_pairs_tree
        """

        adapter = sat_adapters.CurveSATAdapter(dim=curve.dim)
        adapter.init_curve(curve)

        # how often should we call sat solver? default is equidistant strategy
        if sat_strategy is None:
            sat_strategy = {'type': 'equal', 'count': 100}
        if sat_strategy['type'] == 'geometric':
            sat_current_iter = 1

        # okay, we allow FastFractions as thresholds
        if not isinstance(good_threshold, Threshold):
            good_threshold = Threshold('<=', good_threshold)
        if not isinstance(bad_threshold, Threshold):
            bad_threshold = Threshold('>=', bad_threshold)

        # all possible regular curves from given curve are encoded using
        # boolean variables in sat adapter
        #
        # we grow the pairs tree with fixed good and bad thresholds
        # all pairs below good threshold are forgotten
        # all pairs above bad threshold are added to the list of forbidden configurations,
        # i.e. a list of boolean clauses in sat adapter
        #
        # if we can't find a model, then all curves are necessarily bad
        # if there is a model, it avoids bad pairs and all other pairs are good,
        # because we grow the tree till we can
        # (see also estimate_ratio_regular for more explanations)

        if start_pairs_tree is None:
            pairs_tree = self.init_pairs_tree(curve)
        else:
            pairs_tree = start_pairs_tree.copy()
        pairs_tree.set_good_threshold(good_threshold)
        pairs_tree.set_bad_threshold(bad_threshold)
        self.bound_new_pairs(pairs_tree)
        while pairs_tree.bad_pairs:
            bad_pair = pairs_tree.bad_pairs.pop()
            adapter.add_forbid_clause(bad_pair.junc, bad_pair.curve)

        no_model = None
        stats = Counter()
        result = {'stats': stats}

        while pairs_tree.nodes:
            stats['divide_iter'] += 1
            if max_iter is not None and stats['divide_iter'] > max_iter:
                raise RunOutOfIterationsException()

            pairs_tree.divide()
            self.bound_new_pairs(pairs_tree)
            while pairs_tree.bad_pairs:
                bad_pair = pairs_tree.bad_pairs.pop()
                adapter.add_forbid_clause(bad_pair.junc, bad_pair.curve)

            try_sat = False
            if sat_strategy['type'] == 'equal':
                if stats['divide_iter'] % sat_strategy['count'] == 0:
                    try_sat = True
            elif sat_strategy['type'] == 'geometric':
                if stats['divide_iter'] >= sat_current_iter:
                    try_sat = True
                    sat_current_iter = int(sat_current_iter * sat_strategy['multiplier']) + 1

            if try_sat:
                logging.info('current stats: %s', result['stats'])
                if not adapter.solve():
                    no_model = True
                    break

        result['stats'].update({'ptree_' + k: v for k, v in pairs_tree.stats.items()})
        result['stats'].update({'sat_' + k: v for k, v in adapter.get_stats().items()})

        if no_model or not adapter.solve():
            return result

        model = adapter.get_model()
        result['model'] = model
        result['curve'] = adapter.get_curve_from_model(curve, model)
        return result

    def estimate_ratio_sequence(self, curves, rel_tol_inv, rel_tol_inv_mult=2, **kwargs):
        """
        Estimate minimal curve ratio for sequence of fuzzy curves.

        This method relies totally on estimate_ratio_fuzzy.
        Params:
        rel_tol_inv  --  subj
        rel_tol_inv_mult  --  current rel_tol_inv is multiplied by this every epoch

        Additional kwargs are passed as is to estimate_ratio_fuzzy.
        Returns lo, up, and list of curve candidates.
        """

        CurveItem = namedtuple('CurveItem', ['priority', 'lo', 'up', 'curve', 'example', 'pairs_tree'])
        _inc = 0

        def get_item(curve, lo=None, up=None, example=None, pairs_tree=None):
            nonlocal _inc
            _inc += 1
            priority = ((-lo if lo is not None else None), _inc)
            return CurveItem(priority, lo, up, curve, example, pairs_tree)

        curr_lo = FastFraction(0, 1)
        curr_up = None

        active = (get_item(curve) for curve in curves)
        if isinstance(curves, Sized):
            active = list(active)
        curr_rel_tol_inv = 1
        epoch = 0
        stats = Counter()
        tolerance = FastFraction(rel_tol_inv + 1, rel_tol_inv)
        while curr_up is None or curr_up > curr_lo * tolerance:
            curr_rel_tol_inv *= rel_tol_inv_mult
            epoch += 1
            total = len(active) if isinstance(active, list) else -1
            new_active = []  # heap of CurveItem
            for cnt, item in enumerate(active):
                logging.warning('E%d, curve %d / %d', epoch, cnt + 1, total)
                res = self.estimate_ratio_fuzzy(
                    item.curve, rel_tol_inv=curr_rel_tol_inv, upper_bound=curr_up,
                    start_lower_bound=item.lo, start_upper_bound=item.up,
                    start_curve=item.example, start_pairs_tree=item.pairs_tree,
                    **kwargs,
                )
                if curr_up is None or res['up'] < curr_up:
                    curr_up = res['up']
                    logging.warning('new upper bound: %.3f', curr_up)

                if res['lo'] <= curr_up:
                    # have a chance to be the best
                    if total < 0 or total > 200:
                        # avoid memory leak
                        res['pairs_tree'] = None
                    new_item = get_item(curve=item.curve, lo=res['lo'], up=res['up'], example=res['curve'], pairs_tree=res['pairs_tree'])
                    heappush(new_active, new_item)
                    logging.warning('added new active item!')

                while new_active[0].lo > curr_up:  # priority = -lo
                    heappop(new_active)

                logging.warning('current active: %d, stats: %s', len(new_active), res['stats'])
                stats.update(res['stats'])

            active = sorted(new_active, key=lambda item: item.up)  # better to start with good curves
            curr_lo = min(item.lo for item in active)
            logging.warning('current bounds: [%.5f, %.5f]', curr_lo, curr_up)

        return {'lo': curr_lo, 'up': curr_up, 'curves': [d.curve for d in active], 'stats': stats}


class IntegerBrokenLine:
    def __init__(self, dim, brkline):
        denoms = set()
        for x, t in brkline:
            if isinstance(t, FastFraction):
                denoms.add(t.d)
            for xj in x:
                if isinstance(xj, FastFraction):
                    denoms.add(xj.d)
        lcm = get_lcm(denoms)
        mx = lcm
        mt = lcm**dim
        points = []
        for x, t in brkline:
            xp = tuple(int(FastFraction.convert(xj) * FastFraction(mx, 1)) for xj in x)
            tp = int(FastFraction.convert(t) * FastFraction(mt, 1))
            points.append((xp, tp))
        self.points = points
        self.mx = mx
        self.mt = mt
