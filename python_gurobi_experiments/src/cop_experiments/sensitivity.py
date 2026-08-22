"""Exact small-instance epsilon check reported in the manuscript."""

from __future__ import annotations

from itertools import combinations, product
from pathlib import Path

import gurobipy as gp
import numpy as np
import pandas as pd
from gurobipy import GRB

from .metrics import violation_score
from .pcm import upper_pairs
from .priority import GurobiSettings, solve_mnvllsm, solve_stage1


EXAMPLE_PCM = np.array([
    [1, 2, 4, 9],
    [1 / 2, 1, 3, 7],
    [1 / 4, 1 / 3, 1, 5],
    [1 / 9, 1 / 7, 1 / 5, 1],
], dtype=float)


def _ordered_weak_orders(item_count: int):
    """Yield all ordered set partitions as contiguous integer rank vectors."""
    for level_count in range(1, item_count + 1):
        required = set(range(level_count))
        for ranks in product(range(level_count), repeat=item_count):
            if set(ranks) == required:
                yield ranks


def enumerate_feasible_weak_orders(
    matrix: np.ndarray = EXAMPLE_PCM,
    epsilon: float = 1e-4,
    y_bound: float = 10.0,
) -> dict[str, float | int]:
    """Enumerate the 4,683 weak orders and test each with a Gurobi LP."""
    n = matrix.shape[0]
    pairs = upper_pairs(n)
    if len(pairs) != 6:
        raise ValueError("The exhaustive manuscript check is defined for its 4x4 example.")
    input_logs = np.array([np.log(matrix[pair]) for pair in pairs])
    coefficients: list[np.ndarray] = []
    for i, j in pairs:
        coefficient = np.zeros(n)
        coefficient[i] = 1.0
        coefficient[j] = -1.0
        coefficients.append(coefficient)

    enumerated = feasible = 0
    minimum_nv = float("inf")
    environment = gp.Env(empty=True)
    environment.setParam("OutputFlag", 0)
    environment.start()
    try:
        for ranks in _ordered_weak_orders(len(pairs)):
            enumerated += 1
            model = gp.Model("weak_order_feasibility", env=environment)
            model.Params.OutputFlag = 0
            y = model.addVars(n, lb=-y_bound, ub=y_bound, name="y")
            model.addConstr(gp.quicksum(y[i] for i in range(n)) == 0.0)
            differences = [
                gp.quicksum(float(coefficients[p][i]) * y[i] for i in range(n))
                for p in range(len(pairs))
            ]
            blocks = [[p for p, rank in enumerate(ranks) if rank == level]
                      for level in range(max(ranks) + 1)]
            for block in blocks:
                representative = block[0]
                for p in block[1:]:
                    model.addConstr(differences[p] == differences[representative])
            for lower, upper in zip(blocks[:-1], blocks[1:], strict=True):
                model.addConstr(differences[upper[0]] - differences[lower[0]] >= epsilon)
            # Every upper-triangular judgment in this example exceeds the neutral value 1.
            model.addConstr(differences[blocks[0][0]] >= epsilon)
            model.setObjective(0.0, GRB.MINIMIZE)
            model.optimize()
            if model.Status == GRB.OPTIMAL:
                feasible += 1
                nv = 0.0
                for p, q in combinations(range(len(pairs)), 2):
                    input_sign = int(np.sign(input_logs[p] - input_logs[q]))
                    represented_sign = int(np.sign(ranks[p] - ranks[q]))
                    if input_sign != represented_sign:
                        nv += 0.5 if input_sign == 0 or represented_sign == 0 else 1.0
                minimum_nv = min(minimum_nv, nv)
            model.dispose()
    finally:
        environment.dispose()
    return {
        "enumerated_weak_orders": enumerated,
        "feasible_weak_orders": feasible,
        "enumerated_min_nv": minimum_nv,
    }


def run_example_epsilon_sensitivity(
    epsilons: list[float],
    base_settings: GurobiSettings,
    output_path: Path,
) -> pd.DataFrame:
    """Run the exhaustive check once and MNVLLSM for every numerical epsilon."""
    audit = enumerate_feasible_weak_orders(
        EXAMPLE_PCM, epsilon=min(epsilons), y_bound=base_settings.y_bound
    )
    records = []
    for epsilon in epsilons:
        settings = GurobiSettings(
            epsilon=float(epsilon),
            y_bound=base_settings.y_bound,
            time_limit=base_settings.time_limit,
            mip_gap=base_settings.mip_gap,
            threads=base_settings.threads,
            seed=base_settings.seed,
            output_flag=base_settings.output_flag,
        )
        stage1 = solve_stage1(EXAMPLE_PCM, settings, "indicator")
        result = solve_mnvllsm(EXAMPLE_PCM, stage1, settings)
        ranking = None if result.weights is None else (np.argsort(-result.weights) + 1).tolist()
        records.append({
            "epsilon": float(epsilon),
            "status": result.status,
            "solved": result.solved,
            "nv_star": result.nv_star,
            "returned_nv": None if result.weights is None else violation_score(EXAMPLE_PCM, result.weights),
            "ranking": None if ranking is None else ">".join(map(str, ranking)),
            "llsm_deviation": result.objective,
            "runtime": result.runtime,
            **audit,
        })
    frame = pd.DataFrame(records)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_path, index=False)
    return frame
