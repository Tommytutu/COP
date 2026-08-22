from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from cop_experiments.evrim import HOUSE_PCM, solve_evrim_check_first  # noqa: E402
from cop_experiments.metrics import violation_score  # noqa: E402
from cop_experiments.pcm import (  # noqa: E402
    LOG_SCALE,
    consistency_ratio,
    geometric_consistency_index,
    llsm_weights,
    upper_pairs,
)
from cop_experiments.priority import GurobiSettings  # noqa: E402


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
        threads=int(config["threads"]),
        seed=int(config["gurobi_seed"]),
        output_flag=int(config["output_flag"]),
    )


def json_value(value: object) -> str:
    if value is None:
        return ""
    return json.dumps(np.asarray(value).tolist(), separators=(",", ":"))


def protected_values_hold(
    original: np.ndarray,
    revised: np.ndarray,
    protected: list[tuple[int, int]],
) -> bool:
    return all(np.isclose(original[pair], revised[pair], atol=1e-10, rtol=1e-10) for pair in protected)


def scale_feasible(a: np.ndarray) -> bool:
    return all(
        float(np.min(np.abs(LOG_SCALE - np.log(a[i, j])))) <= 1e-9
        for i, j in upper_pairs(a.shape[0])
    )


def main() -> None:
    config_path = ROOT / "config_weighted_pareto_20260819.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    settings = settings_from_config(config)
    results = ROOT / str(config["results_dir"])
    raw_dir = results / "raw"
    summary_dir = results / "summary"
    raw_dir.mkdir(parents=True, exist_ok=True)
    summary_dir.mkdir(parents=True, exist_ok=True)
    full_path = raw_dir / "house_evrim_check_first.csv"
    if full_path.exists():
        existing = pd.read_csv(full_path)
        completed = set(existing["case"].astype(str))
    else:
        completed = set()
    cases = [
        ("No protection", [], "None"),
        ("Protect a37", [(2, 6)], r"a37=6"),
        ("Protect a13,a37", [(0, 2), (2, 6)], r"a13=3; a37=6"),
    ]
    for case, protected, protected_label in cases:
        if case in completed:
            continue
        result = solve_evrim_check_first(
            HOUSE_PCM,
            settings,
            value_protected=protected,
            variant=f"House-{case}",
        )
        revised = result.revised_matrix
        if revised is not None:
            cr = consistency_ratio(revised)
            gci = geometric_consistency_index(revised)
            llsm = llsm_weights(revised)
            llsm_nv = violation_score(revised, llsm)
            protected_ok = protected_values_hold(HOUSE_PCM, revised, protected)
            scale_ok = scale_feasible(revised)
        else:
            cr = gci = llsm_nv = None
            llsm = None
            protected_ok = scale_ok = False
        row = {
            "case": case,
            "protected_value": protected_label,
            "status": result.status,
            "solved": result.solved,
            "early_stop": result.early_stop,
            "nrp": result.nrp,
            "aoc": result.aoc,
            "cr": cr,
            "gci": gci,
            "gci_threshold": result.gci_threshold,
            "gci_excess": None if gci is None else gci - float(result.gci_threshold),
            "nv_certificate": result.nv,
            "nv_llsm": llsm_nv,
            "runtime": result.runtime,
            "screening_runtime": result.screening_runtime,
            "bnc_runtime": result.bnc_runtime,
            "screening_status": result.screening_status,
            "screening_solved": result.screening_solved,
            "screening_nrp": result.screening_nrp,
            "screening_aoc": result.screening_aoc,
            "screening_gci": result.screening_gci,
            "screening_gap": result.screening_gap,
            "stage1_status": result.stage1_status,
            "stage2_status": result.stage2_status,
            "stage1_gap": result.stage1_gap,
            "stage2_gap": result.stage2_gap,
            "stage1_objective": result.stage1_objective,
            "stage1_bound": result.stage1_bound,
            "stage2_objective": result.stage2_objective,
            "stage2_bound": result.stage2_bound,
            "gap": result.gap,
            "callback_mipsol_checks": result.callback_mipsol_checks,
            "callback_mipnode_checks": result.callback_mipnode_checks,
            "lazy_cuts": result.lazy_cuts,
            "user_cuts": result.user_cuts,
            "maximum_gci_excess": result.maximum_gci_excess,
            "node_count": result.node_count,
            "solution_count": result.solution_count,
            "work": result.work,
            "max_violation": result.max_violation,
            "num_variables": result.num_variables,
            "num_binary_variables": result.num_binary_variables,
            "num_constraints": result.num_constraints,
            "num_quadratic_constraints": result.num_quadratic_constraints,
            "num_general_constraints": result.num_general_constraints,
            "protected_values_feasible": protected_ok,
            "scale_feasible": scale_ok,
            "gci_feasible": bool(gci is not None and gci <= float(result.gci_threshold) + settings.gci_tolerance),
            "revised_matrix": json_value(revised),
            "certificate_weights": json_value(result.weights),
            "llsm_weights": json_value(llsm),
            "ranking": json_value(None if llsm is None else np.argsort(-llsm)),
        }
        pd.DataFrame([row]).to_csv(full_path, mode="a", header=not full_path.exists(), index=False)
        print(f"{case}: {result.status}, early_stop={result.early_stop}, runtime={result.runtime:.3f}s", flush=True)

    full = pd.read_csv(full_path)
    table = full[["case", "protected_value", "nrp", "aoc", "cr", "gci", "runtime", "status"]].copy()
    table.columns = ["Case", "Protected value", "NRP", "AOC", "CR", "GCI", "Time", "Status"]
    table.to_csv(summary_dir / "house_table6.csv", index=False)
    print(table.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
