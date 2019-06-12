import logging
import itertools

from pysat.solvers import *

from fractal_curve import FractalCurve

# clause = {var: True|False}
# var - hashable, в числа переводится непосредственно перед вызовом солвера
# token - пара (var, True|False), предполагаем что каждый bm задаётся токеном
class CurveSATAdapter:
    def __init__(self, dim):
        self.dim = dim
        self.int_clauses = []
        self.curve_vars = set()  # переменные кривых, условия для которых уже записаны
        self.var_no = {}

        self.append_formula([
            {('junc', None): True},  # у каждой кривой есть автостык
        ])

    # токен (переменная или её отрицание), истинный тттк в кривой base_maps[cnum] = bm
    # предполагаем, таким образом, что каждый bm задаётся одним токеном
    def get_bm_token(self, cnum, bm):
        if self.dim != 2:
            raise Exception("dim>2 not implemented!")
        bm_var = ('base_map', cnum)
        bm_val = bm.is_oriented()
        return (bm_var, bm_val)

    # создать переменную Z, такую что Z=True <=> кривая согласуется с данной
    # добавляем условия эквивалентности
    def make_curve_var(self, curve):
        curve_info = tuple((cnum, bm) for cnum, bm in enumerate(curve.base_maps) if bm is not None)
        Z = ('curve', curve_info)
        if Z in self.curve_vars:
            return Z
        # Z <-> curve  <=>  (Z->curve) and (curve->Z)
        # Z->curve  <=>  !Z or curve  <=>  (!Z or bm1) and (!Z or bm2) and ... (!Z or bmk)
        clauses = []
        for cnum, bm in curve_info:
            bm_var, bm_val = self.get_bm_token(cnum, bm)
            clause = {Z: False, bm_var: bm_val}
            clauses.append(clause)
        self.append_formula(clauses)

        # curve->Z  <=>  !curve or Z  <=>  !bm1 or !bm2 or ... or !bmk or Z
        clause_rev = {Z: True}
        for cnum, bm in curve_info:
            bm_var, bm_val = self.get_bm_token(cnum, bm)
            clause_rev[bm_var] = not bm_val

        self.append_formula([clause_rev])
        self.curve_vars.add(Z)
        return Z
        
    # создать переменную J, такую что J=True <=> в кривой есть стык J
    # методу нужно передать стык и список кривых, в которых он возникает
    def make_junc_var(self, junc, curves):
        J = ('junc', junc)

        curve_vars = [self.make_curve_var(curve) for curve in curves]
        # J <-> (c1 or c2 or .. or ck)

        # J->(c1 or c2 or .. or ck)  <=>  !J or c1 or c2 .. or ck
        clause1 = {cv: True for cv in curve_vars}
        clause1[J] = False
        self.append_formula([clause1])

        # (c1 or c2 .. or ck)->J  <=>  !(c1 or .. or ck) or J  <=>  (!c1 and ... !ck) or J  <=> (!c1 or J) and ... (!ck or J)
        clauses = []
        for cv in curve_vars:
            clause = {cv: False, J: True}
            clauses.append(clause)
        self.append_formula(clauses)

        return J

    # !(J and bm1 and bm2 .. and bmk) = !J or !bm1 or !bm1 ..
    def add_forbid_clause(self, junc, curve):
        clause = {}
        junc_var = ('junc', junc)
        clause[junc_var] = False
        for cnum, bm in enumerate(curve.base_maps):
            if bm is None:
                continue
            bm_var, bm_val = self.get_bm_token(cnum, bm)
            clause[bm_var] = not bm_val

        self.append_formula([clause])

    # переводим var в натуральные числа здесь
    def append_formula(self, clauses):
        for clause in clauses:
            for var in clause.keys():
                if var not in self.var_no:
                    max_var_no = 1 + len(self.var_no)
                    self.var_no[var] = max_var_no
                    
        int_clauses = []
        for clause in clauses:
            int_clause = []
            for var, val in clause.items():
                var_no = self.var_no[var]
                token = var_no if val else -var_no
                int_clause.append(token)
            int_clauses.append(int_clause)

        self.int_clauses += int_clauses

    def stats(self):
        return {
            'clauses_count', len(self.int_clauses),
            'variables_count', len(self.var_no),
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

    def get_curves_from_model(self, curve, model):
        base_maps = []
        allowed_by_model_variants = []
        for cnum in range(curve.genus):
            allowed_maps = curve.get_allowed_maps(cnum)
            allowed_by_model = []
            for bm in allowed_maps:
                bm_var, bm_val = self.get_bm_token(cnum, bm)
                if (bm_var not in model) or (model[bm_var] == bm_val):
                    allowed_by_model.append(bm)
            if not allowed_by_model:
                return None
            allowed_by_model_variants.append(allowed_by_model)

        full_curves = []
        for base_maps in itertools.product(*allowed_by_model_variants):
            full_curve = FractalCurve(dim=curve.dim, div=curve.div, proto=curve.proto, base_maps=base_maps)
            juncs = full_curve.get_junctions()
            has_bad_juncs = False
            for junc in juncs:
                junc_var = ('junc', junc)
                if junc_var in model and not model[junc_var]:
                    # bad curve
                    has_bad_juncs = True
                    break
            if has_bad_juncs:
                continue
            # структура задачи такова, что достаточно проверить отсутствие плохих стыков!
            # в будущем можно также учесть стыки, не входящие в кривую
            full_curves.append(full_curve)

        return full_curves

    # для отладки
#    def get_model_from_curve(self, curve):
#        curve_clauses = []
#        for cnum, bm in enumerate(curve.base_maps):
#            if bm is not None:
#                bm_var, bm_val = self.get_bm_token(cnum, bm)
#                curve_clauses.append({bm_var: bm_val})
#
#        clauses = self.common_clauses + self.def_clauses + curve_clauses
#        int_clauses = self.get_int_clauses(clauses)
#        self.solver = Glucose3()
#        self.solver.append_formula(int_clauses)
#        self.solver.solve()
#        int_model = self.solver.get_model()
#        return self.get_model_from_int_model(int_model)
