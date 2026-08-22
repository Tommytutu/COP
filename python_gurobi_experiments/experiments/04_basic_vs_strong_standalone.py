"""Standalone Basic-versus-Strong Stage-1 Gurobi experiment.

This script deliberately does not call the project's pipeline or Stage-1 MIP
wrappers.  It reuses the existing ``results/raw/datasets.csv`` instances,
selects the first two replicates in every (n, regime) cell, and writes the
Basic and Strong Stage-1 Gurobi formulations explicitly in this file.

Each completed optimization is appended to CSV immediately and followed by
``flush()`` + ``os.fsync()`` so an interrupted run keeps every finished solve.

Default experiment size
-----------------------
7 matrix sizes x 4 regimes x 2 PCMs per cell = 56 PCMs
56 PCMs x 2 formulations = 112 Stage-1 solves
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from itertools import combinations
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "config.json"
DEFAULT_DATASET = ROOT / "results" / "raw" / "datasets.csv"
DEFAULT_OUTPUT = ROOT / "results" / "raw" / "basic_vs_strong_standalone.csv"

CSV_FIELDS = [
    "instance_id",
    "n",
    "regime",
    "replicate",
    "seed",
    "formulation",
    "status",
    "has_solution",
    "certified",
    "nv_star",
    "nv2_star",
    "runtime_seconds",
    "mip_gap",
    "obj_bound",
    "node_count",
    "raw_order_relations",
    "relation_groups",
    "binary_variables",
    "continuous_variables",
    "linear_constraints",
    "general_constraints",
    "constraints",
    "epsilon",
    "y_bound",
    "time_limit_seconds",
    "threads",
    "gurobi_seed",
]


def parse_matrix(value: str) -> np.ndarray:
    """Parse one JSON-encoded PCM from ``datasets.csv``."""
    matrix = np.asarray(json.loads(value), dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"matrix must be square, got shape {matrix.shape}")
    return matrix


def upper_pairs(n: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def select_instances(
    data: pd.DataFrame,
    regimes: Iterable[str],
    samples_per_cell: int,
) -> pd.DataFrame:
    """Select the first replicas per (n, regime) in deterministic config order."""
    regime_list = list(regimes)
    regime_rank = {name: rank for rank, name in enumerate(regime_list)}
    selected = data[data["regime"].isin(regime_list)].copy()
    selected["_regime_rank"] = selected["regime"].map(regime_rank)
    selected = selected.sort_values(
        ["n", "_regime_rank", "replicate", "instance_id"], kind="stable"
    )
    selected = (
        selected.groupby(["n", "regime"], sort=False, group_keys=False)
        .head(int(samples_per_cell))
        .sort_values(["n", "_regime_rank", "replicate", "instance_id"], kind="stable")
        .drop(columns=["_regime_rank"])
        .reset_index(drop=True)
    )
    return selected


def relation_groups(a: np.ndarray) -> list[tuple[tuple[int, ...], int, int, int]]:
    """Exact algebraic aggregation used by the original Stage-1 implementation.

    Returns ``(coefficients, n_positive, n_negative, n_equal)`` for every
    distinct canonical log-priority difference expression.  The audit items are
    all upper-triangular elicited judgments plus one neutral intensity 1.
    """
    n = int(a.shape[0])
    items = upper_pairs(n) + [(0, 0)]
    log_intensity = np.array([math.log(float(a[i, j])) for i, j in items])

    item_coefficients: list[np.ndarray] = []
    for i, j in items:
        coefficient = np.zeros(n, dtype=int)
        coefficient[i] += 1
        coefficient[j] -= 1
        item_coefficients.append(coefficient)

    grouped: dict[tuple[int, ...], list[int]] = {}
    for p, q in combinations(range(len(items)), 2):
        coefficient = item_coefficients[p] - item_coefficients[q]
        desired_sign = int(np.sign(log_intensity[p] - log_intensity[q]))
        nonzero = np.flatnonzero(coefficient)
        if len(nonzero) == 0:
            raise RuntimeError("unexpected zero relation expression")
        first_nonzero = int(nonzero[0])
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


def raw_relation_lookup(
    n: int,
    group_keys: list[tuple[int, ...]],
) -> tuple[int, dict[tuple[int, int], tuple[int, int]]]:
    """Map each raw item comparison to ``(group_index, orientation)``."""
    items = upper_pairs(n) + [(0, 0)]
    coefficients: list[np.ndarray] = []
    for i, j in items:
        c = np.zeros(n, dtype=int)
        c[i] += 1
        c[j] -= 1
        coefficients.append(c)

    group_index = {key: index for index, key in enumerate(group_keys)}
    lookup: dict[tuple[int, int], tuple[int, int]] = {}
    for p, q in combinations(range(len(items)), 2):
        c = coefficients[p] - coefficients[q]
        nonzero = np.flatnonzero(c)
        if len(nonzero) == 0:
            raise RuntimeError("unexpected zero raw relation expression")
        orientation = 1
        if c[int(nonzero[0])] < 0:
            c = -c
            orientation = -1
        key = tuple(int(value) for value in c)
        lookup[p, q] = (group_index[key], orientation)
    return len(items), lookup


def status_name(code: int) -> str:
    """Return readable Gurobi optimization status without importing Gurobi here."""
    return {
        2: "OPTIMAL",
        3: "INFEASIBLE",
        4: "INF_OR_UNBD",
        5: "UNBOUNDED",
        6: "CUTOFF",
        7: "ITERATION_LIMIT",
        8: "NODE_LIMIT",
        9: "TIME_LIMIT",
        10: "SOLUTION_LIMIT",
        11: "INTERRUPTED",
        12: "NUMERIC",
        13: "SUBOPTIMAL",
        14: "INPROGRESS",
        15: "USER_OBJ_LIMIT",
        16: "WORK_LIMIT",
        17: "MEM_LIMIT",
    }.get(int(code), f"STATUS_{code}")


def append_result_durably(path: Path, row: dict, fieldnames: list[str]) -> None:
    """Append exactly one row, then flush and fsync it before returning."""
    path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not path.exists() or path.stat().st_size == 0
    with path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        if write_header:
            writer.writeheader()
        writer.writerow(row)
        handle.flush()
        os.fsync(handle.fileno())


def load_completed_pairs(path: Path) -> set[tuple[str, str]]:
    """Return completed ``(instance_id, formulation)`` keys for resume."""
    if not path.exists() or path.stat().st_size == 0:
        return set()
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return set()
        required = {"instance_id", "formulation"}
        if not required.issubset(reader.fieldnames):
            raise ValueError(
                f"existing output {path} is missing required columns {sorted(required)}"
            )
        return {
            (str(row["instance_id"]), str(row["formulation"]))
            for row in reader
            if row.get("instance_id") and row.get("formulation")
        }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Standalone Basic-vs-Strong Stage-1 Gurobi experiment"
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--dataset", type=Path, default=DEFAULT_DATASET)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Delete the output CSV first and rerun completed solves.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional number of selected PCM instances (smoke testing only).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate selection and print the target solve count without Gurobi.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config_path = args.config.resolve()
    dataset_path = args.dataset.resolve()
    output_path = args.output.resolve()

    config = json.loads(config_path.read_text(encoding="utf-8"))
    data = pd.read_csv(dataset_path)
    selected = select_instances(
        data,
        config["regimes"],
        int(config["formulation_samples_per_cell"]),
    )

    expected_cells = len(config["n_values"]) * len(config["regimes"])
    expected_instances = expected_cells * int(config["formulation_samples_per_cell"])
    if len(selected) != expected_instances:
        raise RuntimeError(
            f"expected {expected_instances} selected PCMs, found {len(selected)}"
        )
    if selected["instance_id"].nunique() != expected_instances:
        raise RuntimeError("selected dataset contains duplicate instance_id values")

    if args.limit is not None:
        if args.limit < 1:
            raise ValueError("--limit must be a positive integer")
        selected = selected.head(args.limit).copy()

    target_solves = len(selected) * 2
    if args.dry_run:
        print(
            f"selected_instances={len(selected)} target_solves={target_solves} "
            f"dataset={dataset_path}",
            flush=True,
        )
        return

    if args.overwrite and output_path.exists():
        output_path.unlink()
    completed_pairs = load_completed_pairs(output_path)

    try:
        import gurobipy as gp
        from gurobipy import GRB
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "gurobipy is required to run the optimization. Install Gurobi Python "
            "bindings in the environment, then rerun this script."
        ) from exc

    epsilon = float(config["epsilon"])
    y_bound = float(config["y_bound"])
    time_limit = float(config["formulation_time_limit_seconds"])
    mip_gap_target = float(config["mip_gap"])
    threads = int(config["threads"])
    gurobi_seed = int(config["gurobi_seed"])
    output_flag = int(config["output_flag"])

    total_full_experiment = expected_instances * 2
    already_in_target = sum(
        (str(row.instance_id), formulation) in completed_pairs
        for row in selected.itertuples(index=False)
        for formulation in ("basic", "strong")
    )
    completed_counter = already_in_target

    for row in selected.itertuples(index=False):
        instance_id = str(row.instance_id)
        a = parse_matrix(str(row.matrix))
        n = int(row.n)
        groups = relation_groups(a)
        group_keys = [group[0] for group in groups]
        item_count = len(upper_pairs(n)) + 1
        raw_order_relations = math.comb(item_count, 2)

        # ================================================================
        # BASIC STAGE-1 FORMULATION
        # ================================================================
        if (instance_id, "basic") not in completed_pairs:
            model = gp.Model(f"mnvdm_stage1_basic_{instance_id}")
            model.Params.OutputFlag = output_flag
            model.Params.Threads = threads
            model.Params.Seed = gurobi_seed
            model.Params.TimeLimit = time_limit
            model.Params.MIPGap = mip_gap_target
            model.Params.NumericFocus = 3
            model.Params.FeasibilityTol = 1e-9
            model.Params.IntFeasTol = 1e-9
            model.Params.OptimalityTol = 1e-9

            y = model.addVars(n, lb=-y_bound, ub=y_bound, name="y")
            model.addConstr(gp.quicksum(y[i] for i in range(n)) == 0.0, name="gauge")

            basic_nv2_terms = []
            basic_big_m = 4.0 * y_bound + epsilon
            for r, (coefficients, n_positive, n_negative, n_equal) in enumerate(groups):
                d = gp.quicksum(coefficients[i] * y[i] for i in range(n))
                greater = model.addVar(vtype=GRB.BINARY, name=f"g[{r}]")
                equal = model.addVar(vtype=GRB.BINARY, name=f"e[{r}]")

                model.addConstr(d >= epsilon - basic_big_m * (1 - greater))
                model.addConstr(d <= basic_big_m * greater)
                model.addConstr(-d <= basic_big_m * (1 - equal))
                model.addConstr(d <= basic_big_m * (1 - equal))
                model.addConstr(d + epsilon <= basic_big_m * (greater + equal))
                model.addConstr(greater + equal <= 1)

                basic_nv2_terms.append(
                    n_positive * (2 - 2 * greater - equal)
                    + n_negative * (2 * greater + equal)
                    + n_equal * (1 - equal)
                )

            basic_nv2 = gp.quicksum(basic_nv2_terms)
            model.setObjective(basic_nv2, GRB.MINIMIZE)
            model.optimize()
            model.update()

            basic_has_solution = model.SolCount > 0
            basic_nv2_star = int(round(model.ObjVal)) if basic_has_solution else None
            basic_row = {
                "instance_id": instance_id,
                "n": n,
                "regime": str(row.regime),
                "replicate": int(row.replicate),
                "seed": int(row.seed),
                "formulation": "basic",
                "status": status_name(model.Status),
                "has_solution": bool(basic_has_solution),
                "certified": bool(model.Status == GRB.OPTIMAL),
                "nv_star": None if basic_nv2_star is None else basic_nv2_star / 2.0,
                "nv2_star": basic_nv2_star,
                "runtime_seconds": float(model.Runtime),
                "mip_gap": float(model.MIPGap) if basic_has_solution and model.IsMIP else None,
                "obj_bound": float(model.ObjBound) if model.IsMIP else None,
                "node_count": float(model.NodeCount) if model.IsMIP else 0.0,
                "raw_order_relations": raw_order_relations,
                "relation_groups": len(groups),
                "binary_variables": 2 * len(groups),
                "continuous_variables": n,
                "linear_constraints": int(model.NumConstrs),
                "general_constraints": int(model.NumGenConstrs),
                "constraints": int(model.NumConstrs + model.NumGenConstrs),
                "epsilon": epsilon,
                "y_bound": y_bound,
                "time_limit_seconds": time_limit,
                "threads": threads,
                "gurobi_seed": gurobi_seed,
            }
            model.dispose()

            # Persist this solve before constructing the Strong model.
            append_result_durably(output_path, basic_row, CSV_FIELDS)
            completed_pairs.add((instance_id, "basic"))
            completed_counter += 1
            print(
                f"{completed_counter}/{total_full_experiment} {instance_id} basic "
                f"{basic_row['status']} NV*={basic_row['nv_star']} "
                f"time={basic_row['runtime_seconds']:.6f}s -> {output_path}",
                flush=True,
            )

        # ================================================================
        # STRONG STAGE-1 FORMULATION
        # ================================================================
        if (instance_id, "strong") not in completed_pairs:
            model = gp.Model(f"mnvdm_stage1_strong_{instance_id}")
            model.Params.OutputFlag = output_flag
            model.Params.Threads = threads
            model.Params.Seed = gurobi_seed
            model.Params.TimeLimit = time_limit
            model.Params.MIPGap = mip_gap_target
            model.Params.NumericFocus = 3
            model.Params.FeasibilityTol = 1e-9
            model.Params.IntFeasTol = 1e-9
            model.Params.OptimalityTol = 1e-9

            y = model.addVars(n, lb=-y_bound, ub=y_bound, name="y")
            model.addConstr(gp.quicksum(y[i] for i in range(n)) == 0.0, name="gauge")

            strong_nv2_terms = []
            relation_vars = []
            for r, (coefficients, n_positive, n_negative, n_equal) in enumerate(groups):
                d = gp.quicksum(coefficients[i] * y[i] for i in range(n))
                greater = model.addVar(vtype=GRB.BINARY, name=f"g[{r}]")
                equal = model.addVar(vtype=GRB.BINARY, name=f"e[{r}]")
                less = model.addVar(vtype=GRB.BINARY, name=f"l[{r}]")

                model.addConstr(greater + equal + less == 1, name=f"trichotomy[{r}]")
                model.addGenConstrIndicator(
                    greater, True, d, GRB.GREATER_EQUAL, epsilon,
                    name=f"indicator_g[{r}]",
                )
                model.addGenConstrIndicator(
                    equal, True, d, GRB.EQUAL, 0.0,
                    name=f"indicator_e[{r}]",
                )
                model.addGenConstrIndicator(
                    less, True, d, GRB.LESS_EQUAL, -epsilon,
                    name=f"indicator_l[{r}]",
                )

                strong_nv2_terms.append(
                    n_positive * (2 * less + equal)
                    + n_negative * (2 * greater + equal)
                    + n_equal * (greater + less)
                )
                relation_vars.append((greater, equal, less))

            # Valid total-preorder / transitivity cuts.  The raw comparison
            # pairs are mapped back to the exact aggregated relation states.
            strong_item_count, lookup = raw_relation_lookup(n, group_keys)
            for p, q, r_item in combinations(range(strong_item_count), 3):
                group_pq, orientation_pq = lookup[p, q]
                group_qr, orientation_qr = lookup[q, r_item]
                group_pr, orientation_pr = lookup[p, r_item]

                g_pq, e_pq, l_pq = relation_vars[group_pq]
                g_qr, e_qr, l_qr = relation_vars[group_qr]
                g_pr, e_pr, l_pr = relation_vars[group_pr]

                if orientation_pq > 0:
                    ge_pq, le_pq = g_pq + e_pq, l_pq + e_pq
                else:
                    ge_pq, le_pq = l_pq + e_pq, g_pq + e_pq

                if orientation_qr > 0:
                    ge_qr, le_qr = g_qr + e_qr, l_qr + e_qr
                else:
                    ge_qr, le_qr = l_qr + e_qr, g_qr + e_qr

                if orientation_pr > 0:
                    ge_pr, le_pr = g_pr + e_pr, l_pr + e_pr
                else:
                    ge_pr, le_pr = l_pr + e_pr, g_pr + e_pr

                model.addConstr(
                    ge_pq + ge_qr - 1 <= ge_pr,
                    name=f"trans_ge[{p},{q},{r_item}]",
                )
                model.addConstr(
                    le_pq + le_qr - 1 <= le_pr,
                    name=f"trans_le[{p},{q},{r_item}]",
                )

            strong_nv2 = gp.quicksum(strong_nv2_terms)
            model.setObjective(strong_nv2, GRB.MINIMIZE)
            model.optimize()
            model.update()

            strong_has_solution = model.SolCount > 0
            strong_nv2_star = int(round(model.ObjVal)) if strong_has_solution else None
            strong_row = {
                "instance_id": instance_id,
                "n": n,
                "regime": str(row.regime),
                "replicate": int(row.replicate),
                "seed": int(row.seed),
                "formulation": "strong",
                "status": status_name(model.Status),
                "has_solution": bool(strong_has_solution),
                "certified": bool(model.Status == GRB.OPTIMAL),
                "nv_star": None if strong_nv2_star is None else strong_nv2_star / 2.0,
                "nv2_star": strong_nv2_star,
                "runtime_seconds": float(model.Runtime),
                "mip_gap": float(model.MIPGap) if strong_has_solution and model.IsMIP else None,
                "obj_bound": float(model.ObjBound) if model.IsMIP else None,
                "node_count": float(model.NodeCount) if model.IsMIP else 0.0,
                "raw_order_relations": raw_order_relations,
                "relation_groups": len(groups),
                "binary_variables": 3 * len(groups),
                "continuous_variables": n,
                "linear_constraints": int(model.NumConstrs),
                "general_constraints": int(model.NumGenConstrs),
                "constraints": int(model.NumConstrs + model.NumGenConstrs),
                "epsilon": epsilon,
                "y_bound": y_bound,
                "time_limit_seconds": time_limit,
                "threads": threads,
                "gurobi_seed": gurobi_seed,
            }
            model.dispose()

            # Persist Strong immediately as well.
            append_result_durably(output_path, strong_row, CSV_FIELDS)
            completed_pairs.add((instance_id, "strong"))
            completed_counter += 1
            print(
                f"{completed_counter}/{total_full_experiment} {instance_id} strong "
                f"{strong_row['status']} NV*={strong_row['nv_star']} "
                f"time={strong_row['runtime_seconds']:.6f}s -> {output_path}",
                flush=True,
            )


if __name__ == "__main__":
    main()
