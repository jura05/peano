import logging
import itertools

from pysat.solvers import *

from . import curves
from .base_maps import gen_base_maps
from .common import Junction

# clause = {var: True|False}
# var - hashable, в числа переводится непосредственно перед вызовом солвера
# token - пара (var, True|False), предполагаем что каждый sp задаётся токеном
class CurveSATAdapter:
    def __init__(self, dim):
        self.dim = dim
        self.int_clauses = []
        self.curve_vars = set()  # переменные кривых, условия для которых уже записаны
        self.var_no = {}

    def get_sp_var(self, pnum, cnum, sp):
        return ('spec', pnum, cnum, sp)

    def get_junc_var(self, junc):
        return ('junc', junc)

    # добавляем условие, что истинна ровно одна из переменных
    # (var1 or var2 or ... or varN) and (!var_i or !var_j)
    # здесь есть возможность использовать add_atmost
    def make_only(self, only_vars):
        self.append_clause({var: True for var in only_vars})
        for var_pair in itertools.combinations(only_vars, 2):
            self.append_clause({var_pair[0]: False, var_pair[1]: False})

    # инициализация условий для заданной Partial-кривой
    def init_curve(self, curve):
        # возможные base_map-ы:
        for pnum in range(curve.pattern_count):
            for cnum in range(curve.genus):
                self.make_only([self.get_sp_var(pnum, cnum, sp) for sp in curve.gen_allowed_specs(pnum, cnum)])

        # автостыки - есть у каждой кривой
        for junc in curve.gen_auto_junctions():
            self.append_clause({self.get_junc_var(junc): True})

        # стыки
        for junc, curves in curve.get_junctions_info().items():
            self.make_junc_var(junc, curves)


    # создать переменную Z, такую что Z=True <=> кривая согласуется с данной
    # добавляем условия эквивалентности
    def make_curve_var(self, curve):
        curve_info = tuple(curve.sp_info())
        Z = ('curve', curve_info)
        if Z in self.curve_vars:
            return Z
        # Z <-> curve  <=>  (Z->curve) and (curve->Z)
        # Z->curve  <=>  !Z or curve  <=>  (!Z or bm1) and (!Z or bm2) and ... (!Z or bmk)
        for pnum, cnum, sp in curve_info:
            self.append_clause({Z: False, self.get_sp_var(pnum, cnum, sp): True})

        # curve->Z  <=>  !curve or Z  <=>  !bm1 or !bm2 or ... or !bmk or Z
        clause_rev = {self.get_sp_var(pnum, cnum, sp): False for pnum, cnum, sp in curve_info}
        clause_rev[Z] = True
        self.append_clause(clause_rev)

        self.curve_vars.add(Z)
        return Z
        
    # создать переменную J, такую что J=True <=> в кривой есть стык J
    # методу нужно передать стык и список кривых, в которых он возникает
    def make_junc_var(self, junc, curves):
        J = self.get_junc_var(junc)

        curve_vars = [self.make_curve_var(curve) for curve in curves]
        # J <-> (c1 or c2 or .. or ck)

        # J->(c1 or c2 or .. or ck)  <=>  !J or c1 or c2 .. or ck
        clause = {cv: True for cv in curve_vars}
        clause[J] = False
        self.append_clause(clause)

        # (c1 or c2 .. or ck)->J  <=>  !(c1 or .. or ck) or J  <=>  (!c1 and ... !ck) or J  <=> (!c1 or J) and ... (!ck or J)
        for cv in curve_vars:
            self.append_clause({cv: False, J: True})

        return J

    # !(J and bm1 and bm2 .. and bmk) = !J or !bm1 or !bm1 ..
    def add_forbid_clause(self, junc, curve):
        junc_var = self.get_junc_var(junc)
        clause = {self.get_sp_var(pnum, cnum, sp): False for pnum, cnum, sp in curve.sp_info()}
        clause[junc_var] = False
        self.append_clause(clause)

    # переводим var в натуральные числа здесь
    def append_clause(self, clause):
        int_clause = []
        for var, val in clause.items():
            if var not in self.var_no:
                max_var_no = 1 + len(self.var_no)
                self.var_no[var] = max_var_no
            var_no = self.var_no[var]
            token = var_no if val else -var_no
            int_clause.append(token)
        self.int_clauses.append(int_clause)

    def append_formula(self, clauses):
        for clause in clauses:
            self.append_clause(clause)

    def stats(self):
        return {
            'clauses_count': len(self.int_clauses),
            'variables_count': len(self.var_no),
        }

    def solve(self):
        self.solver = Glucose3()
        self.solver.append_formula(self.int_clauses)
        return self.solver.solve()

    def get_model(self):
        int_model = self.solver.get_model()
        return self.get_model_from_int_model(int_model)

    def get_model_from_int_model(self, int_model):
        model = {}
        no_var = {var_no: var for var, var_no in self.var_no.items()}
        for int_tok in int_model:
            var = no_var[abs(int_tok)]
            model[var] = (int_tok > 0)
        return model

    def get_curve_from_model(self, curve, model):
        base_maps = []
        allowed_by_model_variants = []
        for pnum in range(curve.pattern_count):
            for cnum in range(curve.genus):
                allowed_by_model = []
                for sp in curve.gen_allowed_specs(pnum, cnum):
                    sp_var = self.get_sp_var(pnum, cnum, sp)
                    if (sp_var not in model) or model[sp_var]:
                        allowed_by_model.append(sp)
                if not allowed_by_model:
                    return
                allowed_by_model_variants.append(allowed_by_model)

        for all_specs in itertools.product(*allowed_by_model_variants):
            # опять надо нарезать по pnum :( 
            patterns = []
            for pnum, pattern in enumerate(curve.patterns):
                specs = all_specs[pnum * curve.genus : (pnum + 1) * curve.genus]
                patterns.append((pattern.proto, specs))
            full_curve = curves.Curve(dim=curve.dim, div=curve.div, patterns=patterns)
            has_bad_juncs = False
            for junc in full_curve.gen_junctions():
                junc_var = self.get_junc_var(junc)
                if junc_var in model and not model[junc_var]:
                    # bad curve
                    has_bad_juncs = True
                    break
            # структура задачи такова, что достаточно проверить отсутствие плохих стыков!
            # в будущем можно также учесть стыки, не входящие в кривую
            if not has_bad_juncs:
                return full_curve

    # для отладки
    def get_model_from_curve(self, curve):
        for cnum, pnum, sp in curve.sp_info():
            sp_var = self.get_sp_var(pnum, cnum, sp)
            self.append_clause({sp_var: True})

        if not self.solve():
            raise Exception("Can't get model, no such curve!")
        return self.get_model()
