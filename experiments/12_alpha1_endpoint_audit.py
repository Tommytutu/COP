from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import gurobipy as gp
from gurobipy import GRB


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ahpcop.metrics import recovery_metrics  # noqa: E402
from ahpcop.pcm import load_record_matrix, load_record_weights, softmax, upper_pairs  # noqa: E402
from ahpcop.priority import (  # noqa: E402
    GurobiSettings,
    _add_total_preorder_cuts,
    _configure,
    _relation_groups,
    solve_stage1,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Certify and GCI-tie-break the degenerate alpha=1 endpoint."
    )
    parser.add_argument(
        "--config", type=Path, default=ROOT / "config_weighted_pareto_20260819.json"
    )
    return parser.parse_args()


def settings_from_config(config: dict) -> GurobiSettings:
    return GurobiSettings(
        epsilon=float(config["epsilon"]),
        y_bound=float(config["y_bound"]),
        time_limit=float(config["time_limit_seconds"]),
        mip_gap=float(config["mip_gap"]),
        mip_gap_abs=float(config["mip_gap_abs"]),
        feasibility_tol=float(config["feasibility_tol"]),
        optimality_tol=float(config["optimality_tol"]),
        integer_feasibility_tol=float(config["integer_feasibility_tol"]),
        barrier_convergence_tol=float(config["barrier_convergence_tol"]),
        barrier_qcp_convergence_tol=float(config["barrier_qcp_convergence_tol"]),
        gci_tolerance=float(config["gci_tolerance"]),
        numeric_focus=int(config["numeric_focus"]),
        integrality_focus=int(config["integrality_focus"]),
        scale_flag=int(config["scale_flag"]),
        threads=int(config["threads"]),
        seed=int(config["gurobi_seed"]),
        output_flag=int(config["output_flag"]),
    )


def append_csv(path: Path, row: dict) -> None:
    pd.DataFrame([row]).to_csv(path, mode="a", header=not path.exists(), index=False)


def gci_log_deviation(a: np.ndarray, weights: np.ndarray) -> float:
    n = a.shape[0]
    y = np.log(np.asarray(weights, dtype=float))
    squared_residual = sum(
        (float(np.log(a[i, j])) - y[i] + y[j]) ** 2
        for i in range(n) for j in range(i + 1, n)
    )
    return 2.0 * squared_residual / ((n - 1) * (n - 2))


def solve_indicator_gci_tiebreak(
    a: np.ndarray,
    nv2_star: int,
    settings: GurobiSettings,
    warm_y: np.ndarray | None,
    warm_signs: list[int] | None,
) -> tuple[str, bool, np.ndarray | None, float | None, float, float | None]:
    n = a.shape[0]
    relation_groups = _relation_groups(a)
    model = gp.Model("alpha1_indicator_gci_tiebreak")
    _configure(model, settings)
    y = model.addVars(n, lb=-settings.y_bound, ub=settings.y_bound, name="y")
    model.addConstr(gp.quicksum(y[i] for i in range(n)) == 0.0, name="gauge")
    relation_vars = []
    nv2_terms = []
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
    model.addConstr(gp.quicksum(nv2_terms) == int(nv2_star), name="fix_minimum_nv2")
    c = 2.0 / ((n - 1) * (n - 2))
    gci = c * gp.quicksum(
        (float(np.log(a[i, j])) - y[i] + y[j]) ** 2 for i, j in upper_pairs(n)
    )
    model.setObjective(gci, GRB.MINIMIZE)
    if warm_y is not None:
        initial = np.asarray(warm_y, dtype=float).copy()
        initial -= initial.mean()
        for i in range(n):
            y[i].Start = float(initial[i])
    if warm_signs is not None and len(warm_signs) == len(relation_vars):
        for sign, (greater, equal, less) in zip(warm_signs, relation_vars, strict=True):
            greater.Start = float(sign > 0)
            equal.Start = float(sign == 0)
            less.Start = float(sign < 0)
    model.optimize()
    status = "OPTIMAL" if model.Status == GRB.OPTIMAL else (
        "TIME_LIMIT" if model.Status == GRB.TIME_LIMIT else f"STATUS_{model.Status}"
    )
    has_solution = model.SolCount > 0
    gap = float(model.MIPGap) if has_solution else None
    certified = bool(model.Status == GRB.OPTIMAL or (gap is not None and gap <= settings.mip_gap))
    values = np.array([y[i].X for i in range(n)], dtype=float) if has_solution else None
    objective = float(model.ObjVal) if has_solution else None
    runtime = float(model.Runtime)
    model.dispose()
    return status, certified, values, objective, runtime, gap


def main() -> None:
    args = parse_args()
    config = json.loads(args.config.resolve().read_text(encoding="utf-8"))
    settings = settings_from_config(config)
    data = pd.read_csv(ROOT / str(config["source_dataset"]))
    data = data[
        (data["n"] <= 7) | ((data["n"].isin([8, 9])) & (data["replicate"] < 3))
    ].sort_values(["n", "regime", "replicate"])
    results = ROOT / str(config["results_dir"])
    raw_dir = results / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    output = raw_dir / "alpha1_indicator_tight_lex_endpoint_audit.csv"
    path_results = pd.read_csv(raw_dir / "alpha_mnvdm_results.csv", low_memory=False)
    completed = set()
    if output.exists():
        completed = set(pd.read_csv(output)["instance_id"].astype(str))
    total = len(data)
    for index, (_, row) in enumerate(data.iterrows(), start=1):
        instance_id = str(row["instance_id"])
        if instance_id in completed:
            continue
        a = load_record_matrix(row["matrix"])
        w0 = load_record_weights(row["latent_weights"])
        candidates = path_results[
            (path_results["instance_id"].astype(str) == instance_id)
            & path_results["python_feasible"].astype(bool)
        ].sort_values(["nvr", "gci_deviation"])
        warm_y = None
        warm_signs = None
        if not candidates.empty:
            best = candidates.iloc[0]
            warm_y = np.asarray(json.loads(str(best["y"])), dtype=float)
            warm_signs = [int(value) for value in json.loads(str(best["relation_states"]))]
        endpoint_settings = replace(
            settings,
            feasibility_tol=min(settings.feasibility_tol, 1e-7),
        )
        stage1 = solve_stage1(
            a,
            endpoint_settings,
            variant="indicator",
            warm_start_y=warm_y,
            warm_start_signs=warm_signs,
        )
        tiebreak_status = None
        tiebreak_certified = False
        tiebreak_runtime = 0.0
        tiebreak_gci = None
        tiebreak_gap = None
        if stage1.nv2 is not None:
            (
                tiebreak_status,
                tiebreak_certified,
                tiebreak_y,
                tiebreak_gci,
                tiebreak_runtime,
                tiebreak_gap,
            ) = solve_indicator_gci_tiebreak(
                a,
                stage1.nv2,
                endpoint_settings,
                stage1.y,
                stage1.signs,
            )
        else:
            tiebreak_y = None
        if tiebreak_y is not None:
            weights = softmax(tiebreak_y)
            endpoint_source = "indicator_min_nv_gci_tiebreak"
        elif stage1.y is not None:
            weights = softmax(stage1.y)
            endpoint_source = "indicator_stage1_incumbent"
        else:
            weights = None
            endpoint_source = "no_incumbent"
        metrics = recovery_metrics(a, weights, w0) if weights is not None else {}
        n_order_relations = metrics.get("n_order_relations")
        audited_nv = metrics.get("nv")
        audited_nvr = metrics.get("nvr")
        stage1_audit_pass = bool(
            stage1.nv is not None
            and audited_nv is not None
            and abs(float(stage1.nv) - float(audited_nv)) <= 1e-8
        )
        output_row = {
            "instance_id": instance_id,
            "n": int(row["n"]),
            "regime": str(row["regime"]),
            "replicate": int(row["replicate"]),
            "alpha": 1.0,
            "endpoint_source": endpoint_source,
            "stage1_status": stage1.status,
            "stage1_certified": bool(stage1.solved),
            "stage1_audit_pass": stage1_audit_pass,
            "stage1_runtime": stage1.runtime,
            "stage1_gap": stage1.gap,
            "stage1_node_count": stage1.node_count,
            "tiebreak_status": tiebreak_status,
            "tiebreak_certified": tiebreak_certified,
            "tiebreak_runtime": tiebreak_runtime,
            "tiebreak_gap": tiebreak_gap,
            "total_runtime": stage1.runtime + tiebreak_runtime,
            "model_nv": stage1.nv,
            "audited_nv": audited_nv,
            "nv2": stage1.nv2,
            "n_order_relations": n_order_relations,
            "nvr": audited_nvr,
            "gci_deviation": (
                float(tiebreak_gci)
                if tiebreak_gci is not None
                else (None if weights is None else gci_log_deviation(a, weights))
            ),
            "kendall_tau_b": metrics.get("kendall_tau_b"),
            "best_choice_accuracy": metrics.get("best_choice_accuracy"),
            "lrmse": metrics.get("lrmse"),
            "weights": "" if weights is None else json.dumps(weights.tolist(), separators=(",", ":")),
        }
        append_csv(output, output_row)
        completed.add(instance_id)
        if index % 25 == 0 or not stage1.solved or not stage1_audit_pass:
            print(
                f"[{index}/{total}] {instance_id} indicator={stage1.status} "
                f"audit={stage1_audit_pass} NVR={audited_nvr}",
                flush=True,
            )
    print(f"results -> {output}")


if __name__ == "__main__":
    main()
