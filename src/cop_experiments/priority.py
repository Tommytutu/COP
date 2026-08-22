from __future__ import annotations

import time
from dataclasses import dataclass
from itertools import combinations

import gurobipy as gp
import numpy as np
from gurobipy import GRB

from .pcm import em_weights, llsm_weights, softmax, upper_pairs


@dataclass
class GurobiSettings:
    epsilon: float = 1e-4
    y_bound: float = 10.0
    time_limit: float = 60.0
    mip_gap: float = 1e-5
    mip_gap_abs: float = 1e-5
    feasibility_tol: float = 1e-5
    optimality_tol: float = 1e-5
    integer_feasibility_tol: float = 1e-5
    barrier_convergence_tol: float = 1e-5
    barrier_qcp_convergence_tol: float = 1e-5
    gci_tolerance: float = 1e-5
    numeric_focus: int = 3
    integrality_focus: int = 1
    scale_flag: int = 2
    presolve: int = -1
    threads: int = 1
    seed: int = 20260815
    output_flag: int = 0


@dataclass
class Stage1Result:
    status: str
    solved: bool
    nv: float | None
    nv2: int | None
    runtime: float
    gap: float | None
    node_count: float | None
    signs: list[int] | None
    y: np.ndarray | None
    model_variant: str


@dataclass
class PriorityResult:
    method: str
    status: str
    solved: bool
    weights: np.ndarray | None
    runtime: float
    stage1_runtime: float
    stage2_runtime: float
    nv_star: float | None
    objective: float | None
    gap: float | None


@dataclass
class AlphaMNVDMResult:
    """Complete result of one blended alpha-MNVDM solve."""

    alpha: float
    status: str
    status_code: int
    solved: bool
    certified: bool
    has_solution: bool
    weights: np.ndarray | None
    y: np.ndarray | None
    relation_states: list[int] | None
    relation_details: list[dict] | None
    nv: float | None
    nv2: int | None
    n_order_relations: int
    nvr: float | None
    gci_deviation: float | None
    gci_normalizer: float
    weighted_objective: float | None
    solver_objective: float | None
    objective_bound: float | None
    gap: float | None
    runtime: float
    work: float | None
    node_count: float
    iteration_count: float
    barrier_iteration_count: int
    solution_count: int
    max_violation: float | None
    max_integrality_violation: float | None
    min_strict_slack: float | None
    max_equality_residual: float | None
    warm_start_source: str
    num_variables: int
    num_binary_variables: int
    num_constraints: int
    num_quadratic_constraints: int
    num_general_constraints: int
    endpoint_state_rejections: int
    endpoint_refinement_runtime: float
    endpoint_fixed_state_status: str | None
    endpoint_fallback_used: bool


def _status_name(code: int) -> str:
    return {
        GRB.OPTIMAL: "OPTIMAL",
        GRB.TIME_LIMIT: "TIME_LIMIT",
        GRB.INFEASIBLE: "INFEASIBLE",
        GRB.INF_OR_UNBD: "INF_OR_UNBD",
        GRB.UNBOUNDED: "UNBOUNDED",
        GRB.INTERRUPTED: "INTERRUPTED",
        GRB.NUMERIC: "NUMERIC",
    }.get(code, f"STATUS_{code}")


def _configure(model: gp.Model, settings: GurobiSettings) -> None:
    model.Params.OutputFlag = settings.output_flag
    model.Params.Threads = settings.threads
    model.Params.Seed = settings.seed
    model.Params.TimeLimit = settings.time_limit
    model.Params.MIPGap = settings.mip_gap
    model.Params.MIPGapAbs = settings.mip_gap_abs
    model.Params.NumericFocus = settings.numeric_focus
    model.Params.IntegralityFocus = settings.integrality_focus
    model.Params.ScaleFlag = settings.scale_flag
    model.Params.Presolve = settings.presolve
    model.Params.FeasibilityTol = settings.feasibility_tol
    model.Params.IntFeasTol = settings.integer_feasibility_tol
    model.Params.OptimalityTol = settings.optimality_tol
    model.Params.BarConvTol = settings.barrier_convergence_tol
    model.Params.BarQCPConvTol = settings.barrier_qcp_convergence_tol


def _relation_groups(a: np.ndarray) -> list[tuple[tuple[int, ...], int, int, int]]:
    """Aggregate order relations that have the same difference expression.

    Each returned tuple contains ``(coefficients, n_positive, n_negative,
    n_equal)``.  Reversed expressions are put in the same canonical group.
    This is an exact algebraic presolve: for example,
    ``(y_i-y_j)-(y_k-y_l)`` and ``(y_i-y_k)-(y_j-y_l)`` are identical.
    The neutral item represents the intensity 1 and therefore also includes
    preference-order-preservation relations in the same audit.
    """
    n = a.shape[0]
    items = upper_pairs(n) + [(0, 0)]
    x = np.array([np.log(a[i, j]) for i, j in items])
    item_coefficients: list[np.ndarray] = []
    for i, j in items:
        coefficient = np.zeros(n, dtype=int)
        coefficient[i] += 1
        coefficient[j] -= 1
        item_coefficients.append(coefficient)

    grouped: dict[tuple[int, ...], list[int]] = {}
    for p, q in combinations(range(len(items)), 2):
        coefficient = item_coefficients[p] - item_coefficients[q]
        delta = float(x[p] - x[q])
        # Inputs reconstructed from log-scale values can differ by a few ulps
        # even when they denote the same Saaty intensity.  Equality is a
        # theoretical case, so classify it with a scale-independent numerical
        # tolerance rather than the bitwise sign of a floating-point residual.
        desired_sign = 0 if abs(delta) <= 1e-10 else (1 if delta > 0 else -1)
        first_nonzero = int(np.flatnonzero(coefficient)[0])
        if coefficient[first_nonzero] < 0:
            coefficient = -coefficient
            desired_sign = -desired_sign
        key = tuple(int(value) for value in coefficient)
        counts = grouped.setdefault(key, [0, 0, 0])
        if desired_sign > 0:
            counts[0] += 1
        elif desired_sign < 0:
            counts[1] += 1
        else:
            counts[2] += 1
    return [(key, counts[0], counts[1], counts[2]) for key, counts in grouped.items()]


def _raw_relation_lookup(
    n: int,
    group_keys: list[tuple[int, ...]],
) -> tuple[int, dict[tuple[int, int], tuple[int, int]]]:
    """Map every raw intensity comparison to its canonical relation group."""
    items = upper_pairs(n) + [(0, 0)]
    item_coefficients: list[np.ndarray] = []
    for i, j in items:
        coefficient = np.zeros(n, dtype=int)
        coefficient[i] += 1
        coefficient[j] -= 1
        item_coefficients.append(coefficient)
    group_index = {key: index for index, key in enumerate(group_keys)}
    lookup: dict[tuple[int, int], tuple[int, int]] = {}
    for p, q in combinations(range(len(items)), 2):
        coefficient = item_coefficients[p] - item_coefficients[q]
        first_nonzero = int(np.flatnonzero(coefficient)[0])
        orientation = 1
        if coefficient[first_nonzero] < 0:
            coefficient = -coefficient
            orientation = -1
        lookup[p, q] = (group_index[tuple(int(value) for value in coefficient)], orientation)
    return len(items), lookup


def _add_total_preorder_cuts(
    model: gp.Model,
    relation_vars: list[tuple[gp.Var, ...]],
    relation_groups: list[tuple[tuple[int, ...], int, int, int]],
    n: int,
) -> None:
    """Add valid transitivity inequalities for the induced scalar preorder."""
    item_count, lookup = _raw_relation_lookup(n, [group[0] for group in relation_groups])

    def weak_states(p: int, q: int) -> tuple[gp.LinExpr, gp.LinExpr]:
        group, orientation = lookup[p, q]
        greater, equal, less = relation_vars[group]
        if orientation > 0:
            return greater + equal, less + equal
        return less + equal, greater + equal

    for p, q, r in combinations(range(item_count), 3):
        ge_pq, le_pq = weak_states(p, q)
        ge_qr, le_qr = weak_states(q, r)
        ge_pr, le_pr = weak_states(p, r)
        model.addConstr(ge_pq + ge_qr - 1 <= ge_pr, name=f"trans_ge[{p},{q},{r}]")
        model.addConstr(le_pq + le_qr - 1 <= le_pr, name=f"trans_le[{p},{q},{r}]")


def _add_strong_order_system(
    model: gp.Model,
    y: gp.tupledict,
    a: np.ndarray,
    settings: GurobiSettings,
) -> gp.LinExpr:
    """Add the exact trichotomy system and return the doubled NV expression."""
    n = a.shape[0]
    relation_groups = _relation_groups(a)
    relation_vars: list[tuple[gp.Var, gp.Var, gp.Var]] = []
    nv2_terms: list[gp.LinExpr] = []
    for r, (coefficients, n_positive, n_negative, n_equal) in enumerate(relation_groups):
        d = gp.quicksum(coefficients[i] * y[i] for i in range(n))
        greater = model.addVar(vtype=GRB.BINARY, name=f"g[{r}]")
        equal = model.addVar(vtype=GRB.BINARY, name=f"e[{r}]")
        less = model.addVar(vtype=GRB.BINARY, name=f"l[{r}]")
        model.addConstr(greater + equal + less == 1, name=f"trichotomy[{r}]")
        # With |y_i| <= y_bound, this particular difference has the valid
        # relation-specific bounds [lower, upper].  Explicit linear disjunctions
        # avoid the looser feasibility handling observed for indicator equalities
        # near the strict/equality boundary, while retaining a theoretically
        # valid (non-heuristic) M for every relation.
        radius = settings.y_bound * float(sum(abs(value) for value in coefficients))
        lower, upper = -radius, radius
        model.addConstr(
            d >= settings.epsilon - (settings.epsilon - lower) * (1 - greater),
            name=f"greater_lb[{r}]",
        )
        model.addConstr(d <= upper * greater, name=f"greater_ub[{r}]")
        model.addConstr(
            d <= upper * (1 - equal),
            name=f"equal_ub[{r}]",
        )
        model.addConstr(
            d >= lower * (1 - equal),
            name=f"equal_lb[{r}]",
        )
        model.addConstr(
            d <= -settings.epsilon + (upper + settings.epsilon) * (1 - less),
            name=f"less_ub[{r}]",
        )
        model.addConstr(d >= lower * less, name=f"less_lb[{r}]")
        nv2_terms.append(
            n_positive * (2 * less + equal)
            + n_negative * (2 * greater + equal)
            + n_equal * (greater + less)
        )
        relation_vars.append((greater, equal, less))
    _add_total_preorder_cuts(model, relation_vars, relation_groups, n)
    return gp.quicksum(nv2_terms)


def alpha_grid(step: float = 0.01) -> list[float]:
    """Return 1, 1-step, ..., 0 with both endpoints included."""
    if step <= 0 or step > 1:
        raise ValueError("alpha step must be in (0, 1]")
    count = int(round(1.0 / step))
    if not np.isclose(count * step, 1.0, atol=1e-12):
        raise ValueError("alpha step must divide one exactly")
    return [round(1.0 - index * step, 12) for index in range(count + 1)]


def _solve_fixed_state_gci(
    a: np.ndarray,
    relation_groups: list[tuple[tuple[int, ...], int, int, int]],
    states: list[int],
    settings: GurobiSettings,
    time_limit: float,
) -> tuple[str, np.ndarray | None, float | None, float, float | None]:
    """Check one discrete order state without big-M and select its best GCI y."""
    n = a.shape[0]
    model = gp.Model("alpha1_fixed_state_check")
    _configure(model, settings)
    model.Params.OutputFlag = 0
    model.Params.TimeLimit = max(1e-3, float(time_limit))
    # The fixed-state model has no disjunction and can safely use a tighter
    # internal feasibility tolerance to certify the requested epsilon boundary.
    model.Params.FeasibilityTol = min(settings.feasibility_tol, 1e-9)
    y = model.addVars(n, lb=-settings.y_bound, ub=settings.y_bound, name="fixed_y")
    model.addConstr(gp.quicksum(y[i] for i in range(n)) == 0.0, name="fixed_gauge")
    for r, ((coefficients, _positive, _negative, _equal), state) in enumerate(
        zip(relation_groups, states, strict=True)
    ):
        difference = gp.quicksum(coefficients[i] * y[i] for i in range(n))
        if state > 0:
            model.addConstr(difference >= settings.epsilon, name=f"fixed_greater[{r}]")
        elif state < 0:
            model.addConstr(difference <= -settings.epsilon, name=f"fixed_less[{r}]")
        else:
            model.addConstr(difference == 0.0, name=f"fixed_equal[{r}]")
    coefficient = 2.0 / ((n - 1) * (n - 2))
    objective = coefficient * gp.quicksum(
        (float(np.log(a[i, j])) - y[i] + y[j])
        * (float(np.log(a[i, j])) - y[i] + y[j])
        for i, j in upper_pairs(n)
    )
    model.setObjective(objective, GRB.MINIMIZE)
    model.optimize()
    status = _status_name(int(model.Status))
    has_solution = model.SolCount > 0
    values = np.array([y[i].X for i in range(n)], dtype=float) if has_solution else None
    objective_value = float(model.ObjVal) if has_solution else None
    runtime = float(model.Runtime)
    max_violation = float(model.MaxVio) if has_solution else None
    model.dispose()
    return status, values, objective_value, runtime, max_violation


def solve_alpha_mnvdm(
    a: np.ndarray,
    alpha: float,
    settings: GurobiSettings,
    warm_start_y: np.ndarray | None = None,
    warm_start_states: list[int] | None = None,
    warm_start_source: str = "LLSM",
    gci_normalizer: float = 1.0,
    log_file: str | None = None,
) -> AlphaMNVDMResult:
    """Solve the single-objective alpha-MNVDM formulation.

    The blended objective is alpha*NVR + (1-alpha)*(D_GCI/gci_normalizer).
    D_GCI uses the manuscript's GCI-normalized squared log residual.  The
    routine records both solver attributes and independently recomputed
    feasibility diagnostics so time-limited incumbents remain auditable.
    """
    if not 0.0 <= alpha <= 1.0:
        raise ValueError("alpha must lie in [0, 1]")
    if gci_normalizer <= 0:
        raise ValueError("gci_normalizer must be positive")
    a = np.asarray(a, dtype=float)
    n = a.shape[0]
    relation_groups = _relation_groups(a)
    n_order_relations = int(sum(pos + neg + eq for _, pos, neg, eq in relation_groups))

    model = gp.Model(f"alpha_mnvdm_{n}_{alpha:.2f}")
    _configure(model, settings)
    # The endpoint is highly degenerate; retain the configured presolve and
    # rely on the independent residual audit plus the fixed-state refinement.
    # Disabling presolve here can materially weaken the numerical search and
    # has produced falsely inferior endpoint incumbents on otherwise easy
    # instances.
    if log_file:
        model.Params.OutputFlag = 1
        model.Params.LogToConsole = 0
        model.Params.LogFile = str(log_file)

    y = model.addVars(n, lb=-settings.y_bound, ub=settings.y_bound, name="y")
    model.addConstr(gp.quicksum(y[i] for i in range(n)) == 0.0, name="gauge")
    relation_vars: list[tuple[gp.Var, gp.Var, gp.Var]] = []
    nv2_terms: list[gp.LinExpr] = []
    for r, (coefficients, n_positive, n_negative, n_equal) in enumerate(relation_groups):
        d = gp.quicksum(coefficients[i] * y[i] for i in range(n))
        greater = model.addVar(vtype=GRB.BINARY, name=f"g[{r}]")
        equal = model.addVar(vtype=GRB.BINARY, name=f"e[{r}]")
        less = model.addVar(vtype=GRB.BINARY, name=f"l[{r}]")
        model.addConstr(greater + equal + less == 1, name=f"trichotomy[{r}]")
        radius = settings.y_bound * float(sum(abs(value) for value in coefficients))
        lower, upper = -radius, radius
        model.addConstr(
            d >= settings.epsilon - (settings.epsilon - lower) * (1 - greater),
            name=f"greater_lb[{r}]",
        )
        model.addConstr(d <= upper * (1 - equal), name=f"equal_ub[{r}]")
        model.addConstr(d >= lower * (1 - equal), name=f"equal_lb[{r}]")
        model.addConstr(
            d <= -settings.epsilon + (upper + settings.epsilon) * (1 - less),
            name=f"less_ub[{r}]",
        )
        nv2_terms.append(
            n_positive * (2 * less + equal)
            + n_negative * (2 * greater + equal)
            + n_equal * (greater + less)
        )
        relation_vars.append((greater, equal, less))
    _add_total_preorder_cuts(model, relation_vars, relation_groups, n)

    nv2_expr = gp.quicksum(nv2_terms)
    nvr_expr = nv2_expr / (2.0 * n_order_relations)
    c = 2.0 / ((n - 1) * (n - 2))
    gci_expr = c * gp.quicksum(
        (float(np.log(a[i, j])) - y[i] + y[j])
        * (float(np.log(a[i, j])) - y[i] + y[j])
        for i, j in upper_pairs(n)
    )
    model.setObjective(
        float(alpha) * nvr_expr + (1.0 - float(alpha)) * gci_expr / gci_normalizer,
        GRB.MINIMIZE,
    )

    if warm_start_y is None:
        initial = np.log(llsm_weights(a))
        initial -= initial.mean()
    else:
        initial = np.asarray(warm_start_y, dtype=float).copy()
        initial -= initial.mean()
    initial = np.clip(initial, -settings.y_bound, settings.y_bound)
    for i in range(n):
        y[i].Start = float(initial[i])
    if warm_start_states is not None and len(warm_start_states) == len(relation_vars):
        for state, (greater, equal, less) in zip(warm_start_states, relation_vars, strict=True):
            greater.Start = float(state > 0)
            equal.Start = float(state == 0)
            less.Start = float(state < 0)

    solve_started = time.perf_counter()
    model.optimize()
    endpoint_state_rejections = 0
    endpoint_refinement_runtime = 0.0
    endpoint_fixed_state_status: str | None = None
    endpoint_values_y: np.ndarray | None = None
    endpoint_gci_value: float | None = None
    endpoint_states_override: list[int] | None = None
    endpoint_fallback_y: np.ndarray | None = None
    endpoint_fallback_states: list[int] | None = None
    endpoint_fallback_gci: float | None = None
    endpoint_fallback_used = False
    endpoint_validated = not np.isclose(alpha, 1.0, atol=1e-12)
    if not endpoint_validated:
        if warm_start_states is not None and len(warm_start_states) == len(relation_groups):
            fallback_y = initial.copy()
            fallback_ok = True
            for (coefficients, _positive, _negative, _equal), state in zip(
                relation_groups, warm_start_states, strict=True
            ):
                difference = float(np.dot(np.asarray(coefficients, dtype=float), fallback_y))
                if state > 0 and difference < settings.epsilon - settings.feasibility_tol:
                    fallback_ok = False
                elif state < 0 and difference > -settings.epsilon + settings.feasibility_tol:
                    fallback_ok = False
                elif state == 0 and abs(difference) > settings.feasibility_tol:
                    fallback_ok = False
            if fallback_ok:
                fallback_residuals = np.array([
                    float(np.log(a[i, j])) - fallback_y[i] + fallback_y[j]
                    for i, j in upper_pairs(n)
                ])
                endpoint_fallback_y = fallback_y
                endpoint_fallback_states = list(warm_start_states)
                endpoint_fallback_gci = c * float(np.dot(fallback_residuals, fallback_residuals))
        while model.SolCount > 0:
            candidate_states = []
            selected_variables = []
            for greater, equal, less in relation_vars:
                if greater.X > 0.5:
                    candidate_states.append(1)
                    selected_variables.append(greater)
                elif equal.X > 0.5:
                    candidate_states.append(0)
                    selected_variables.append(equal)
                else:
                    candidate_states.append(-1)
                    selected_variables.append(less)
            remaining = settings.time_limit - (time.perf_counter() - solve_started)
            if remaining <= 1e-3:
                endpoint_fixed_state_status = "TIME_LIMIT"
                break
            fixed_status, fixed_y, fixed_gci, fixed_runtime, _fixed_max_vio = _solve_fixed_state_gci(
                a, relation_groups, candidate_states, settings, remaining
            )
            endpoint_refinement_runtime += fixed_runtime
            endpoint_fixed_state_status = fixed_status
            if fixed_y is not None and fixed_status == "OPTIMAL":
                endpoint_values_y = fixed_y
                endpoint_gci_value = fixed_gci
                endpoint_states_override = candidate_states
                endpoint_validated = True
                break
            endpoint_state_rejections += 1
            model.addConstr(
                gp.quicksum(selected_variables) <= len(selected_variables) - 1,
                name=f"reject_infeasible_endpoint_state[{endpoint_state_rejections}]",
            )
            remaining = settings.time_limit - (time.perf_counter() - solve_started)
            if remaining <= 1e-3:
                endpoint_fixed_state_status = "TIME_LIMIT"
                break
            model.Params.TimeLimit = remaining
            model.optimize()
        if not endpoint_validated and endpoint_fallback_y is not None:
            endpoint_values_y = endpoint_fallback_y
            endpoint_gci_value = endpoint_fallback_gci
            endpoint_states_override = endpoint_fallback_states
            endpoint_fallback_used = True
            endpoint_validated = True
    status_code = int(model.Status)
    status = _status_name(status_code)
    has_solution = model.SolCount > 0 and endpoint_validated
    gap = float(model.MIPGap) if has_solution and model.IsMIP else None
    certified = has_solution and not endpoint_fallback_used and (
        status_code == GRB.OPTIMAL
        or (gap is not None and gap <= settings.mip_gap)
    )

    values_y: np.ndarray | None = None
    weights: np.ndarray | None = None
    states: list[int] | None = None
    details: list[dict] | None = None
    nv2: int | None = None
    nv: float | None = None
    nvr: float | None = None
    gci_value: float | None = None
    weighted_value: float | None = None
    min_strict_slack: float | None = None
    max_equality_residual: float | None = None
    if has_solution:
        values_y = (
            endpoint_values_y.copy()
            if endpoint_values_y is not None
            else np.array([y[i].X for i in range(n)], dtype=float)
        )
        weights = softmax(values_y)
        details = []
        states = []
        nv2_value = 0
        strict_slacks: list[float] = []
        equality_residuals: list[float] = []
        for r, ((coefficients, n_positive, n_negative, n_equal), variables) in enumerate(
            zip(relation_groups, relation_vars, strict=True)
        ):
            greater, equal, less = variables
            state = (
                endpoint_states_override[r]
                if endpoint_states_override is not None
                else (1 if greater.X > 0.5 else (0 if equal.X > 0.5 else -1))
            )
            difference = float(np.dot(np.asarray(coefficients, dtype=float), values_y))
            states.append(state)
            nv2_contribution = int(
                n_positive * (2 * (state < 0) + (state == 0))
                + n_negative * (2 * (state > 0) + (state == 0))
                + n_equal * ((state > 0) + (state < 0))
            )
            nv2_value += nv2_contribution
            if state > 0:
                strict_slacks.append(difference - settings.epsilon)
            elif state < 0:
                strict_slacks.append(-difference - settings.epsilon)
            else:
                equality_residuals.append(abs(difference))
            details.append({
                "group": r,
                "coefficients": list(coefficients),
                "n_positive": n_positive,
                "n_negative": n_negative,
                "n_equal": n_equal,
                "state": state,
                "difference": difference,
                "nv2_contribution": nv2_contribution,
            })
        nv2 = int(nv2_value)
        nv = nv2 / 2.0
        nvr = nv / n_order_relations
        residuals = np.array([
            float(np.log(a[i, j])) - values_y[i] + values_y[j]
            for i, j in upper_pairs(n)
        ])
        gci_value = (
            float(endpoint_gci_value)
            if endpoint_gci_value is not None
            else c * float(np.dot(residuals, residuals))
        )
        weighted_value = float(alpha) * nvr + (1.0 - float(alpha)) * gci_value / gci_normalizer
        min_strict_slack = min(strict_slacks) if strict_slacks else None
        max_equality_residual = max(equality_residuals) if equality_residuals else 0.0

    solver_objective = float(model.ObjVal) if has_solution else None
    objective_bound = float(model.ObjBound) if model.IsMIP else solver_objective
    result = AlphaMNVDMResult(
        alpha=float(alpha), status=status, status_code=status_code,
        solved=certified, certified=certified, has_solution=has_solution,
        weights=weights, y=values_y, relation_states=states, relation_details=details,
        nv=nv, nv2=nv2, n_order_relations=n_order_relations, nvr=nvr,
        gci_deviation=gci_value, gci_normalizer=float(gci_normalizer),
        weighted_objective=weighted_value, solver_objective=solver_objective,
        objective_bound=objective_bound, gap=gap,
        runtime=float(time.perf_counter() - solve_started),
        work=float(model.Work) if hasattr(model, "Work") else None,
        node_count=float(model.NodeCount), iteration_count=float(model.IterCount),
        barrier_iteration_count=int(model.BarIterCount), solution_count=int(model.SolCount),
        max_violation=float(model.MaxVio) if has_solution else None,
        max_integrality_violation=float(model.IntVio) if has_solution and model.IsMIP else None,
        min_strict_slack=min_strict_slack, max_equality_residual=max_equality_residual,
        warm_start_source=warm_start_source, num_variables=int(model.NumVars),
        num_binary_variables=int(model.NumBinVars), num_constraints=int(model.NumConstrs),
        num_quadratic_constraints=int(model.NumQConstrs),
        num_general_constraints=int(model.NumGenConstrs),
        endpoint_state_rejections=endpoint_state_rejections,
        endpoint_refinement_runtime=endpoint_refinement_runtime,
        endpoint_fixed_state_status=endpoint_fixed_state_status,
        endpoint_fallback_used=endpoint_fallback_used,
    )
    model.dispose()
    return result


def solve_stage1(
    a: np.ndarray,
    settings: GurobiSettings,
    variant: str = "indicator",
    warm_start_y: np.ndarray | None = None,
    warm_start_signs: list[int] | None = None,
) -> Stage1Result:
    n = a.shape[0]
    relation_groups = _relation_groups(a)
    model = gp.Model(f"mnvdm_stage1_{variant}")
    _configure(model, settings)
    y = model.addVars(n, lb=-settings.y_bound, ub=settings.y_bound, name="y")
    model.addConstr(gp.quicksum(y[i] for i in range(n)) == 0.0, name="gauge")
    nv2_terms: list[gp.LinExpr] = []
    relation_vars: list[tuple[gp.Var, ...]] = []

    if variant == "indicator":
        for r, (coefficients, n_positive, n_negative, n_equal) in enumerate(relation_groups):
            d = gp.quicksum(coefficients[i] * y[i] for i in range(n))
            greater = model.addVar(vtype=GRB.BINARY, name=f"g[{r}]")
            equal = model.addVar(vtype=GRB.BINARY, name=f"e[{r}]")
            less = model.addVar(vtype=GRB.BINARY, name=f"l[{r}]")
            model.addConstr(greater + equal + less == 1, name=f"trichotomy[{r}]")
            model.addGenConstrIndicator(greater, True, d, GRB.GREATER_EQUAL, settings.epsilon)
            model.addGenConstrIndicator(equal, True, d, GRB.EQUAL, 0.0)
            model.addGenConstrIndicator(less, True, d, GRB.LESS_EQUAL, -settings.epsilon)
            nv2_terms.append(
                n_positive * (2 * less + equal)
                + n_negative * (2 * greater + equal)
                + n_equal * (greater + less)
            )
            relation_vars.append((greater, equal, less))
        _add_total_preorder_cuts(model, relation_vars, relation_groups, n)
    elif variant in {"basic", "strong"}:
        for r, (coefficients, n_positive, n_negative, n_equal) in enumerate(relation_groups):
            d = gp.quicksum(coefficients[i] * y[i] for i in range(n))
            # Since |y_i| <= y_bound, this is the tight bound implied by the
            # declared variable box for this particular difference expression.
            d_bound = settings.y_bound * sum(abs(value) for value in coefficients)
            big_m = d_bound + settings.epsilon
            greater = model.addVar(vtype=GRB.BINARY, name=f"g[{r}]")
            equal = model.addVar(vtype=GRB.BINARY, name=f"e[{r}]")
            # Basic forward implication system for d > 0 and d = 0.
            model.addConstr(d >= settings.epsilon - big_m * (1 - greater), name=f"basic_g_lb[{r}]")
            model.addConstr(d <= d_bound * greater, name=f"basic_g_ub[{r}]")
            model.addConstr(-d <= d_bound * (1 - equal), name=f"basic_e_lb[{r}]")
            model.addConstr(d <= d_bound * (1 - equal), name=f"basic_e_ub[{r}]")
            model.addConstr(d + settings.epsilon <= big_m * (greater + equal), name=f"basic_l[{r}]")
            model.addConstr(greater + equal <= 1, name=f"basic_partition[{r}]")
            if variant == "basic":
                nv2_terms.append(
                    n_positive * (2 - 2 * greater - equal)
                    + n_negative * (2 * greater + equal)
                    + n_equal * (1 - equal)
                )
                relation_vars.append((greater, equal))
            else:
                # The strong model retains every basic inequality, introduces
                # the explicit reverse state, and adds the valid trichotomy
                # equality.  Thus its continuous relaxation is a subset of the
                # basic relaxation after projection onto (y,g,e).
                less = model.addVar(vtype=GRB.BINARY, name=f"l[{r}]")
                model.addConstr(-d >= settings.epsilon - big_m * (1 - less), name=f"reverse_l_lb[{r}]")
                model.addConstr(-d <= d_bound * less, name=f"reverse_l_ub[{r}]")
                model.addConstr(-d + settings.epsilon <= big_m * (less + equal), name=f"reverse_g[{r}]")
                model.addConstr(less + equal <= 1, name=f"reverse_partition[{r}]")
                model.addConstr(greater + equal + less == 1, name=f"trichotomy[{r}]")
                nv2_terms.append(
                    n_positive * (2 * less + equal)
                    + n_negative * (2 * greater + equal)
                    + n_equal * (greater + less)
                )
                relation_vars.append((greater, equal, less))
        if variant == "strong":
            _add_total_preorder_cuts(model, relation_vars, relation_groups, n)
    else:
        raise ValueError(f"Unknown Stage-1 variant: {variant}")

    nv2_expr = gp.quicksum(nv2_terms)
    model.setObjective(nv2_expr, GRB.MINIMIZE)
    if warm_start_y is not None:
        initial = np.asarray(warm_start_y, dtype=float).copy()
        initial -= initial.mean()
        initial = np.clip(initial, -settings.y_bound, settings.y_bound)
        for i in range(n):
            y[i].Start = float(initial[i])
    if warm_start_signs is not None and len(warm_start_signs) == len(relation_vars):
        for sign, variables in zip(warm_start_signs, relation_vars, strict=True):
            if variant in {"indicator", "strong"}:
                greater, equal, less = variables
                greater.Start = float(sign > 0)
                equal.Start = float(sign == 0)
                less.Start = float(sign < 0)
            else:
                greater, equal = variables
                greater.Start = float(sign > 0)
                equal.Start = float(sign == 0)
    model.optimize()
    status = _status_name(model.Status)
    has_solution = model.SolCount > 0
    gap = float(model.MIPGap) if has_solution and model.IsMIP else None
    solved = model.Status == GRB.OPTIMAL or (
        has_solution and gap is not None and gap <= settings.mip_gap
    )
    if has_solution:
        nv2 = int(round(model.ObjVal))
        values_y = np.array([y[i].X for i in range(n)])
        signs: list[int] = []
        for variables in relation_vars:
            if variant in {"indicator", "strong"}:
                greater, equal, _less = variables
                signs.append(1 if greater.X > 0.5 else (0 if equal.X > 0.5 else -1))
            else:
                greater, equal = variables
                signs.append(1 if greater.X > 0.5 else (0 if equal.X > 0.5 else -1))
    else:
        nv2, values_y, signs = None, None, None
    nodes = float(model.NodeCount) if model.IsMIP else 0.0
    runtime = float(model.Runtime)
    model.dispose()
    return Stage1Result(
        status=status,
        solved=solved,
        nv=None if nv2 is None else nv2 / 2.0,
        nv2=nv2,
        runtime=runtime,
        gap=gap,
        node_count=nodes,
        signs=signs,
        y=values_y,
        model_variant=variant,
    )


def _add_fixed_order_constraints(
    model: gp.Model,
    y: gp.tupledict,
    a: np.ndarray,
    signs: list[int],
    epsilon: float,
) -> None:
    relation_groups = _relation_groups(a)
    for r, ((coefficients, _n_positive, _n_negative, _n_equal), sign) in enumerate(
        zip(relation_groups, signs, strict=True)
    ):
        d = gp.quicksum(coefficients[i] * y[i] for i in range(a.shape[0]))
        if sign > 0:
            model.addConstr(d >= epsilon, name=f"fixed_g[{r}]")
        elif sign < 0:
            model.addConstr(d <= -epsilon, name=f"fixed_l[{r}]")
        else:
            model.addConstr(d == 0.0, name=f"fixed_e[{r}]")


def solve_mnvllsm(a: np.ndarray, stage1: Stage1Result, settings: GurobiSettings) -> PriorityResult:
    if not stage1.solved or stage1.nv is None or stage1.nv2 is None:
        return PriorityResult(
            "MNVLLSM", stage1.status, False, None, stage1.runtime,
            stage1.runtime, 0.0, stage1.nv, None, stage1.gap,
        )
    n = a.shape[0]
    model = gp.Model("mnvllsm_stage2")
    _configure(model, settings)
    y = model.addVars(n, lb=-settings.y_bound, ub=settings.y_bound, name="y")
    model.addConstr(gp.quicksum(y[i] for i in range(n)) == 0.0)
    nv2_expr = _add_strong_order_system(model, y, a, settings)
    model.addConstr(nv2_expr == stage1.nv2, name="fix_optimal_nv")
    c = 2.0 / ((n - 1) * (n - 2))
    objective = c * gp.quicksum(
        (float(np.log(a[i, j])) - y[i] + y[j]) *
        (float(np.log(a[i, j])) - y[i] + y[j])
        for i, j in upper_pairs(n)
    )
    model.setObjective(objective, GRB.MINIMIZE)
    model.optimize()
    status = _status_name(model.Status)
    has_solution = model.SolCount > 0
    gap2 = float(model.MIPGap) if has_solution and model.IsMIP else None
    solved = model.Status == GRB.OPTIMAL or (
        has_solution and gap2 is not None and gap2 <= settings.mip_gap
    )
    weights = softmax(np.array([y[i].X for i in range(n)])) if has_solution else None
    obj = float(model.ObjVal) if has_solution else None
    runtime2 = float(model.Runtime)
    model.dispose()
    return PriorityResult(
        "MNVLLSM", status, solved, weights, stage1.runtime + runtime2,
        stage1.runtime, runtime2, stage1.nv, obj, gap2,
    )


def solve_mnvem(a: np.ndarray, stage1: Stage1Result, settings: GurobiSettings) -> PriorityResult:
    if not stage1.solved or stage1.nv is None or stage1.nv2 is None:
        return PriorityResult(
            "MNVEM", stage1.status, False, None, stage1.runtime,
            stage1.runtime, 0.0, stage1.nv, None, stage1.gap,
        )
    n = a.shape[0]
    model = gp.Model("mnvem_stage2")
    _configure(model, settings)
    model.Params.FuncNonlinear = 1
    y = model.addVars(n, lb=-settings.y_bound, ub=settings.y_bound, name="y")
    model.addConstr(gp.quicksum(y[i] for i in range(n)) == 0.0)
    nv2_expr = _add_strong_order_system(model, y, a, settings)
    model.addConstr(nv2_expr == stage1.nv2, name="fix_optimal_nv")
    lam = model.addVar(lb=0.0, name="lambda")
    exp_terms: dict[tuple[int, int], gp.Var] = {}
    for i in range(n):
        for j in range(n):
            u = model.addVar(lb=-2 * settings.y_bound - np.log(9),
                             ub=2 * settings.y_bound + np.log(9), name=f"u[{i},{j}]")
            t = model.addVar(lb=0.0, name=f"exp[{i},{j}]")
            model.addConstr(u == float(np.log(a[i, j])) + y[j] - y[i])
            model.addGenConstrExp(u, t, name=f"expcon[{i},{j}]")
            exp_terms[i, j] = t
    for i in range(n):
        model.addConstr(gp.quicksum(exp_terms[i, j] for j in range(n)) <= lam)
    model.setObjective(lam, GRB.MINIMIZE)
    model.optimize()
    status = _status_name(model.Status)
    has_solution = model.SolCount > 0
    gap2 = float(model.MIPGap) if has_solution and model.IsMIP else None
    solved = model.Status == GRB.OPTIMAL or (
        has_solution and gap2 is not None and gap2 <= settings.mip_gap
    )
    weights = softmax(np.array([y[i].X for i in range(n)])) if has_solution else None
    obj = float(model.ObjVal) if has_solution else None
    runtime2 = float(model.Runtime)
    model.dispose()
    return PriorityResult(
        "MNVEM", status, solved, weights, stage1.runtime + runtime2,
        stage1.runtime, runtime2, stage1.nv, obj, gap2,
    )


def classical_priorities(a: np.ndarray) -> list[PriorityResult]:
    output: list[PriorityResult] = []
    for name, function in (("EM", em_weights), ("LLSM", llsm_weights)):
        start = time.perf_counter()
        weights = function(a)
        elapsed = time.perf_counter() - start
        output.append(PriorityResult(name, "OPTIMAL", True, weights, elapsed, 0.0, elapsed, None, None, 0.0))
    return output


def cop_llsm_from_mnvllsm(result: PriorityResult) -> PriorityResult:
    if result.nv_star == 0:
        if result.weights is not None:
            return PriorityResult(
                "COP-LLSM", result.status, result.solved, result.weights.copy(),
                result.stage2_runtime, 0.0, result.stage2_runtime, 0.0, result.objective, result.gap,
            )
        return PriorityResult(
            "COP-LLSM", result.status, False, None, result.stage2_runtime,
            0.0, result.stage2_runtime, 0.0, None, result.gap,
        )
    return PriorityResult("COP-LLSM", "INFEASIBLE", False, None, 0.0, 0.0, 0.0, result.nv_star, None, None)
