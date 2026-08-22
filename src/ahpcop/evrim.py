from __future__ import annotations

import time
from dataclasses import dataclass, replace
from itertools import combinations

import gurobipy as gp
import numpy as np
from gurobipy import GRB

from .metrics import violation_score
from .pcm import LOG_SCALE, geometric_consistency_index, matrix_from_upper_logs, softmax, upper_pairs
from .priority import (
    GurobiSettings,
    _add_total_preorder_cuts,
    _configure,
    _raw_relation_lookup,
    _relation_groups,
    _status_name,
)


@dataclass
class EVRIMResult:
    variant: str
    status: str
    solved: bool
    revised_matrix: np.ndarray | None
    weights: np.ndarray | None
    nrp: int | None
    aoc: float | None
    gci: float | None
    nv: float | None
    runtime: float
    stage1_runtime: float
    stage2_runtime: float
    gap: float | None
    value_protected: list[tuple[int, int]]
    direction_protected: list[tuple[int, int]]
    early_stop: bool = False
    screening_status: str | None = None
    screening_solved: bool | None = None
    screening_nrp: int | None = None
    screening_aoc: float | None = None
    screening_gci: float | None = None
    screening_gap: float | None = None
    screening_runtime: float = 0.0
    bnc_runtime: float = 0.0
    callback_mipsol_checks: int = 0
    callback_mipnode_checks: int = 0
    lazy_cuts: int = 0
    user_cuts: int = 0
    maximum_gci_excess: float = 0.0
    gci_threshold: float | None = None
    gci_constraint_enforced: bool = True
    stage1_status: str | None = None
    stage2_status: str | None = None
    stage1_gap: float | None = None
    stage2_gap: float | None = None
    stage1_objective: float | None = None
    stage1_bound: float | None = None
    stage2_objective: float | None = None
    stage2_bound: float | None = None
    node_count: float = 0.0
    solution_count: int = 0
    work: float | None = None
    max_violation: float | None = None
    num_variables: int = 0
    num_binary_variables: int = 0
    num_constraints: int = 0
    num_quadratic_constraints: int = 0
    num_general_constraints: int = 0


def gci_threshold(n: int) -> float:
    return 0.31 if n == 3 else (0.35 if n == 4 else 0.37)


def _original_scale_indices(a: np.ndarray, pairs: list[tuple[int, int]]) -> list[int]:
    return [int(np.argmin(np.abs(LOG_SCALE - np.log(a[i, j])))) for i, j in pairs]


def _normalize_pair(pair: tuple[int, int]) -> tuple[int, int]:
    i, j = pair
    return (i, j) if i < j else (j, i)


def default_protected_sets(a: np.ndarray) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    pairs = upper_pairs(a.shape[0])
    ordered = sorted(pairs, key=lambda p: abs(np.log(a[p])), reverse=True)
    value = [ordered[0]] if ordered else []
    direction = [ordered[1]] if len(ordered) > 1 else []
    return value, direction


def solve_evrim(
    a: np.ndarray,
    settings: GurobiSettings,
    threshold: float | None = None,
    value_protected: list[tuple[int, int]] | None = None,
    direction_protected: list[tuple[int, int]] | None = None,
    variant: str = "EVRIM",
    backend: str = "direct",
    enforce_gci: bool = True,
) -> EVRIMResult:
    if backend not in {"direct", "oa_callback"}:
        raise ValueError("backend must be 'direct' or 'oa_callback'")
    n = a.shape[0]
    deadline = time.perf_counter() + settings.time_limit
    threshold = gci_threshold(n) if threshold is None else threshold
    value_protected = sorted({_normalize_pair(p) for p in (value_protected or [])})
    direction_protected = sorted({_normalize_pair(p) for p in (direction_protected or [])})
    value_set, direction_set = set(value_protected), set(direction_protected)
    pairs = upper_pairs(n)
    pair_index = {pair: p for p, pair in enumerate(pairs)}
    original = _original_scale_indices(a, pairs)

    model = gp.Model(f"{variant}_{backend}_lex")
    _configure(model, settings)
    model.Params.TimeLimit = settings.time_limit
    y = model.addVars(n, lb=-settings.y_bound, ub=settings.y_bound, name="y")
    model.addConstr(gp.quicksum(y[i] for i in range(n)) == 0.0, name="gauge")
    choose = model.addVars(len(pairs), len(LOG_SCALE), vtype=GRB.BINARY, name="choose")
    x = model.addVars(len(pairs), lb=float(LOG_SCALE[0]), ub=float(LOG_SCALE[-1]), name="x")

    for p, pair in enumerate(pairs):
        model.addConstr(gp.quicksum(choose[p, m] for m in range(len(LOG_SCALE))) == 1)
        model.addConstr(x[p] == gp.quicksum(float(LOG_SCALE[m]) * choose[p, m] for m in range(len(LOG_SCALE))))
        if pair in value_set:
            model.addConstr(choose[p, original[p]] == 1, name=f"value_protect[{p}]")
        if pair in direction_set:
            sign = np.sign(np.log(a[pair]))
            if sign > 0:
                model.addConstr(gp.quicksum(choose[p, m] for m, v in enumerate(LOG_SCALE) if v >= -1e-12) == 1)
            elif sign < 0:
                model.addConstr(gp.quicksum(choose[p, m] for m, v in enumerate(LOG_SCALE) if v <= 1e-12) == 1)

    items = pairs + [(0, 0)]
    neutral_index = len(pairs)
    relation_groups = _relation_groups(a)
    relation_vars: list[tuple[gp.Var, gp.Var, gp.Var]] = []
    for r, (coefficients, _n_positive, _n_negative, _n_equal) in enumerate(relation_groups):
        dy = gp.quicksum(coefficients[i] * y[i] for i in range(n))
        greater = model.addVar(vtype=GRB.BINARY, name=f"g[{r}]")
        equal = model.addVar(vtype=GRB.BINARY, name=f"e[{r}]")
        less = model.addVar(vtype=GRB.BINARY, name=f"l[{r}]")
        model.addConstr(greater + equal + less == 1)
        model.addGenConstrIndicator(greater, True, dy, GRB.GREATER_EQUAL, settings.epsilon)
        model.addGenConstrIndicator(equal, True, dy, GRB.EQUAL, 0.0)
        model.addGenConstrIndicator(less, True, dy, GRB.LESS_EQUAL, -settings.epsilon)
        relation_vars.append((greater, equal, less))

    _add_total_preorder_cuts(model, relation_vars, relation_groups, n)
    _item_count, relation_lookup = _raw_relation_lookup(n, [group[0] for group in relation_groups])
    for r, (p, q) in enumerate(combinations(range(len(items)), 2)):
        xp = 0.0 if p == neutral_index else x[p]
        xq = 0.0 if q == neutral_index else x[q]
        dx = xp - xq
        group, orientation = relation_lookup[p, q]
        greater, equal, less = relation_vars[group]
        raw_greater, raw_less = (greater, less) if orientation > 0 else (less, greater)
        model.addGenConstrIndicator(raw_greater, True, dx, GRB.GREATER_EQUAL, settings.epsilon)
        model.addGenConstrIndicator(equal, True, dx, GRB.EQUAL, 0.0)
        model.addGenConstrIndicator(raw_less, True, dx, GRB.LESS_EQUAL, -settings.epsilon)

    callback_stats = {
        "mipsol_checks": 0,
        "mipnode_checks": 0,
        "lazy_cuts": 0,
        "user_cuts": 0,
        "maximum_gci_excess": 0.0,
    }
    callback = None
    if enforce_gci:
        # The POIP/COP certificate y and the LLSM vector u play different roles.
        # For a complete reciprocal PCM, u_i is the row mean of the revised log PCM.
        u_bound = (n - 1) * float(LOG_SCALE[-1]) / n
        u = model.addVars(n, lb=-u_bound, ub=u_bound, name="u")
        for i in range(n):
            row_sum = gp.LinExpr()
            for j in range(n):
                if i < j:
                    row_sum += x[pair_index[i, j]]
                elif j < i:
                    row_sum -= x[pair_index[j, i]]
            model.addConstr(u[i] == row_sum / n, name=f"llsm_u[{i}]")

        c = 2.0 / ((n - 1) * (n - 2))
        residuals = [x[p] - u[i] + u[j] for p, (i, j) in enumerate(pairs)]
        if backend == "direct":
            model.addQConstr(c * gp.quicksum(r * r for r in residuals) <= threshold, name="gci")
        else:
            model.Params.LazyConstraints = 1
            model.Params.PreCrush = 1

            def oa_callback(callback_model: gp.Model, where: int) -> None:
                if where == GRB.Callback.MIPSOL:
                    callback_stats["mipsol_checks"] += 1
                    x_hat = np.array([callback_model.cbGetSolution(x[p]) for p in range(len(pairs))])
                    u_hat = np.array([callback_model.cbGetSolution(u[i]) for i in range(n)])
                    add_cut = callback_model.cbLazy
                    cut_counter = "lazy_cuts"
                elif where == GRB.Callback.MIPNODE:
                    if callback_model.cbGet(GRB.Callback.MIPNODE_STATUS) != GRB.OPTIMAL:
                        return
                    callback_stats["mipnode_checks"] += 1
                    x_hat = np.array([callback_model.cbGetNodeRel(x[p]) for p in range(len(pairs))])
                    u_hat = np.array([callback_model.cbGetNodeRel(u[i]) for i in range(n)])
                    add_cut = callback_model.cbCut
                    cut_counter = "user_cuts"
                else:
                    return
                r_hat = np.array([
                    x_hat[p] - u_hat[i] + u_hat[j]
                    for p, (i, j) in enumerate(pairs)
                ])
                g_hat = c * float(np.dot(r_hat, r_hat))
                callback_stats["maximum_gci_excess"] = max(
                    callback_stats["maximum_gci_excess"],
                    g_hat - threshold,
                )
                if g_hat <= threshold + settings.gci_tolerance:
                    return
                gradient_x = 2.0 * c * r_hat
                gradient_u = np.zeros(n)
                for p, (i, j) in enumerate(pairs):
                    gradient_u[i] -= 2.0 * c * r_hat[p]
                    gradient_u[j] += 2.0 * c * r_hat[p]
                tangent = gp.LinExpr(g_hat)
                tangent += gp.quicksum(
                    float(gradient_x[p]) * (x[p] - float(x_hat[p])) for p in range(len(pairs))
                )
                tangent += gp.quicksum(
                    float(gradient_u[i]) * (u[i] - float(u_hat[i])) for i in range(n)
                )
                add_cut(tangent <= threshold)
                callback_stats[cut_counter] += 1

            callback = oa_callback

    def optimize() -> None:
        remaining = max(0.001, deadline - time.perf_counter())
        model.Params.TimeLimit = remaining
        if callback is None:
            model.optimize()
        else:
            model.optimize(callback)

    nrp_expr = gp.quicksum(1 - choose[p, original[p]] for p in range(len(pairs)))
    aoc_expr = gp.quicksum(
        abs(float(LOG_SCALE[m] - LOG_SCALE[original[p]])) * choose[p, m]
        for p in range(len(pairs)) for m in range(len(LOG_SCALE))
    )

    model.setObjective(nrp_expr, GRB.MINIMIZE)
    optimize()
    status1 = _status_name(model.Status)
    runtime1 = float(model.Runtime)
    gap1 = float(model.MIPGap) if model.SolCount else None
    stage1_objective = float(model.ObjVal) if model.SolCount else None
    stage1_bound = float(model.ObjBound) if model.IsMIP else stage1_objective
    certified1 = model.Status == GRB.OPTIMAL or (
        gap1 is not None and gap1 <= settings.mip_gap
    )
    if not certified1:
        node_count = float(model.NodeCount)
        solution_count = int(model.SolCount)
        work = float(model.Work) if hasattr(model, "Work") else None
        max_violation = float(model.MaxVio) if model.SolCount else None
        num_variables = int(model.NumVars)
        num_binary_variables = int(model.NumBinVars)
        num_constraints = int(model.NumConstrs)
        num_quadratic_constraints = int(model.NumQConstrs)
        num_general_constraints = int(model.NumGenConstrs)
        model.dispose()
        return EVRIMResult(
            variant, status1, False, None, None, None, None, None, None,
            runtime1, runtime1, 0.0, gap1, value_protected, direction_protected,
            callback_mipsol_checks=int(callback_stats["mipsol_checks"]),
            callback_mipnode_checks=int(callback_stats["mipnode_checks"]),
            lazy_cuts=int(callback_stats["lazy_cuts"]),
            user_cuts=int(callback_stats["user_cuts"]),
            maximum_gci_excess=float(callback_stats["maximum_gci_excess"]),
            gci_threshold=float(threshold),
            gci_constraint_enforced=enforce_gci,
            stage1_status=status1, stage1_gap=gap1,
            stage1_objective=stage1_objective, stage1_bound=stage1_bound,
            node_count=node_count, solution_count=solution_count,
            work=work, max_violation=max_violation,
            num_variables=num_variables, num_binary_variables=num_binary_variables,
            num_constraints=num_constraints, num_quadratic_constraints=num_quadratic_constraints,
            num_general_constraints=num_general_constraints,
        )

    nrp_star = int(round(model.ObjVal))
    model.addConstr(nrp_expr == nrp_star, name="fix_nrp")
    model.setObjective(aoc_expr, GRB.MINIMIZE)
    optimize()
    status2 = _status_name(model.Status)
    runtime2 = float(model.Runtime)
    has_solution = model.SolCount > 0
    gap2 = float(model.MIPGap) if has_solution else None
    stage2_objective = float(model.ObjVal) if has_solution else None
    stage2_bound = float(model.ObjBound) if model.IsMIP else stage2_objective
    solved = model.Status == GRB.OPTIMAL or (
        gap2 is not None and gap2 <= settings.mip_gap
    )
    if has_solution:
        selected_logs: dict[tuple[int, int], float] = {}
        for p, pair in enumerate(pairs):
            selected_logs[pair] = float(x[p].X)
        revised = matrix_from_upper_logs(n, selected_logs)
        weights = softmax(np.array([y[i].X for i in range(n)]))
        aoc = float(sum(abs(selected_logs[pair] - np.log(a[pair])) for pair in pairs))
        gci = geometric_consistency_index(revised)
        nv = violation_score(revised, weights)
        gap = gap2
        if enforce_gci and gci > threshold + settings.gci_tolerance:
            solved = False
    else:
        revised = weights = None
        aoc = gci = nv = gap = None
    node_count = float(model.NodeCount)
    solution_count = int(model.SolCount)
    work = float(model.Work) if hasattr(model, "Work") else None
    max_violation = float(model.MaxVio) if has_solution else None
    num_variables = int(model.NumVars)
    num_binary_variables = int(model.NumBinVars)
    num_constraints = int(model.NumConstrs)
    num_quadratic_constraints = int(model.NumQConstrs)
    num_general_constraints = int(model.NumGenConstrs)
    model.dispose()
    return EVRIMResult(
        variant=variant,
        status=status2,
        solved=solved,
        revised_matrix=revised,
        weights=weights,
        nrp=nrp_star if has_solution else None,
        aoc=aoc,
        gci=gci,
        nv=nv,
        runtime=runtime1 + runtime2,
        stage1_runtime=runtime1,
        stage2_runtime=runtime2,
        gap=gap,
        value_protected=value_protected,
        direction_protected=direction_protected,
        callback_mipsol_checks=int(callback_stats["mipsol_checks"]),
        callback_mipnode_checks=int(callback_stats["mipnode_checks"]),
        lazy_cuts=int(callback_stats["lazy_cuts"]),
        user_cuts=int(callback_stats["user_cuts"]),
        maximum_gci_excess=float(callback_stats["maximum_gci_excess"]),
        gci_threshold=float(threshold),
        gci_constraint_enforced=enforce_gci,
        stage1_status=status1, stage2_status=status2,
        stage1_gap=gap1, stage2_gap=gap2,
        stage1_objective=stage1_objective, stage1_bound=stage1_bound,
        stage2_objective=stage2_objective, stage2_bound=stage2_bound,
        node_count=node_count, solution_count=solution_count,
        work=work, max_violation=max_violation,
        num_variables=num_variables, num_binary_variables=num_binary_variables,
        num_constraints=num_constraints, num_quadratic_constraints=num_quadratic_constraints,
        num_general_constraints=num_general_constraints,
    )


def solve_evrim_check_first(
    a: np.ndarray,
    settings: GurobiSettings,
    threshold: float | None = None,
    value_protected: list[tuple[int, int]] | None = None,
    direction_protected: list[tuple[int, int]] | None = None,
    variant: str = "EVRIM-check-first",
) -> EVRIMResult:
    """Solve EVRIM without GCI first, then run fresh B&C only if needed.

    No MIP start is passed to the second model.  The two phases share one
    wall-clock budget.  Early termination is exact only when the GCI-free
    lexicographic optimum is certified and independently satisfies GCI.
    """
    a = np.asarray(a, dtype=float)
    actual_threshold = gci_threshold(a.shape[0]) if threshold is None else float(threshold)
    started = time.perf_counter()
    screening = solve_evrim(
        a=a,
        settings=settings,
        threshold=actual_threshold,
        value_protected=value_protected,
        direction_protected=direction_protected,
        variant=f"{variant}-screen",
        backend="direct",
        enforce_gci=False,
    )
    screening_runtime = time.perf_counter() - started
    screening_passes = bool(
        screening.solved
        and screening.gci is not None
        and screening.gci <= actual_threshold + settings.gci_tolerance
    )
    common = {
        "variant": variant,
        "screening_status": screening.status,
        "screening_solved": screening.solved,
        "screening_nrp": screening.nrp,
        "screening_aoc": screening.aoc,
        "screening_gci": screening.gci,
        "screening_gap": screening.gap,
        "screening_runtime": screening_runtime,
        "gci_threshold": actual_threshold,
    }
    if screening_passes:
        return replace(
            screening,
            **common,
            early_stop=True,
            runtime=screening_runtime,
            bnc_runtime=0.0,
            gci_constraint_enforced=False,
        )

    remaining = settings.time_limit - screening_runtime
    if remaining <= 0.01:
        return replace(
            screening,
            **common,
            status="TIME_LIMIT_BEFORE_BNC",
            solved=False,
            early_stop=False,
            runtime=screening_runtime,
            bnc_runtime=0.0,
        )
    bnc_settings = replace(settings, time_limit=remaining)
    bnc_started = time.perf_counter()
    bnc = solve_evrim(
        a=a,
        settings=bnc_settings,
        threshold=actual_threshold,
        value_protected=value_protected,
        direction_protected=direction_protected,
        variant=f"{variant}-bnc",
        backend="oa_callback",
        enforce_gci=True,
    )
    bnc_runtime = time.perf_counter() - bnc_started
    return replace(
        bnc,
        **common,
        early_stop=False,
        runtime=screening_runtime + bnc_runtime,
        bnc_runtime=bnc_runtime,
    )


HOUSE_PCM = np.array([
    [1, 5, 3, 7, 6, 6, 1/3, 1/4],
    [1/5, 1, 1/3, 5, 3, 3, 1/5, 1/7],
    [1/3, 3, 1, 6, 3, 4, 6, 1/5],
    [1/7, 1/5, 1/6, 1, 1/3, 1/4, 1/7, 1/8],
    [1/6, 1/3, 1/3, 3, 1, 1/2, 1/5, 1/6],
    [1/6, 1/3, 1/4, 4, 2, 1, 1/5, 1/6],
    [3, 5, 1/6, 7, 5, 5, 1, 1/2],
    [4, 7, 5, 8, 6, 6, 2, 1],
], dtype=float)


def house_cases(settings: GurobiSettings) -> list[EVRIMResult]:
    # Manuscript indices are one-based: a_37 corresponds to Python pair (2, 6).
    return [
        solve_evrim(HOUSE_PCM, settings, variant="House-A-unprotected"),
        solve_evrim(HOUSE_PCM, settings, value_protected=[(2, 6)], variant="House-B-protect-a37"),
        solve_evrim(HOUSE_PCM, settings, value_protected=[(0, 2), (2, 6)], variant="House-C-protect-a13-a37"),
    ]
