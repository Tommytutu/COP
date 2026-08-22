from __future__ import annotations

import json
import platform
import sys
from pathlib import Path

import gurobipy as gp
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ahpcop.evrim import (  # noqa: E402
    default_protected_sets,
    solve_evrim_check_first,
)
from ahpcop.pcm import (  # noqa: E402
    consistency_ratio,
    load_record_matrix,
)
from ahpcop.priority import GurobiSettings  # noqa: E402


RESULTS_DIR = ROOT / "results_table7_check_first_20260821"
RAW_PATH = RESULTS_DIR / "raw" / "table7_check_first_rows.csv"
SUMMARY_PATH = RESULTS_DIR / "summary" / "table7_updated.csv"
FULL_STATS_PATH = RESULTS_DIR / "summary" / "table7_full_runtime_statistics.csv"
DIRECT_SOURCE = ROOT / "results_model_aligned_20260819" / "raw" / "evrim_results.csv"
DATASET_SOURCE = ROOT / "results_model_aligned_20260819" / "raw" / "datasets.csv"
CONFIG_PATH = ROOT / "config_weighted_pareto_20260819.json"
TARGET_SIZES = (7, 8, 9)


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
        scale_flag=int(config.get("scale_flag", 2)),
        threads=int(config["threads"]),
        seed=int(config["gurobi_seed"]),
        output_flag=int(config["output_flag"]),
    )


def json_value(value: object) -> str:
    if value is None:
        return ""
    return json.dumps(np.asarray(value).tolist(), separators=(",", ":"))


def append_row(path: Path, row: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([row]).to_csv(path, mode="a", header=not path.exists(), index=False)


def completed_keys(path: Path) -> set[tuple[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return set()
    frame = pd.read_csv(path)
    return set(zip(frame["instance_id"].astype(str), frame["variant"].astype(str)))


def certified_percentage(values: pd.Series) -> float:
    return 100.0 * float(values.astype(bool).mean())


def summarize_runtime(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for (n, variant), group in frame.groupby(["n", "variant"], sort=False):
        runtime = pd.to_numeric(group["runtime"], errors="coerce").dropna()
        certified = group["certified"].astype(bool)
        rows.append(
            {
                "n": int(n),
                "Variant": variant,
                "Instances": int(len(group)),
                "Mean": float(runtime.mean()),
                "Median": float(runtime.median()),
                "Std": float(runtime.std(ddof=1)),
                "P90": float(runtime.quantile(0.90)),
                "P95": float(runtime.quantile(0.95)),
                "Max": float(runtime.max()),
                "Certified (%)": certified_percentage(certified),
                "Timeout (%)": 100.0 - certified_percentage(certified),
            }
        )
    return pd.DataFrame(rows)


def build_direct_rows(source: pd.DataFrame, target_ids: set[str]) -> pd.DataFrame:
    direct = source[
        source["instance_id"].astype(str).isin(target_ids)
        & source["variant"].eq("EVRIM-Direct")
    ].copy()
    direct["variant"] = "Direct"
    direct["certified"] = direct["solved"].astype(str).str.lower().eq("true")
    direct["early_stop"] = False
    direct["screening_runtime"] = 0.0
    direct["bnc_runtime"] = 0.0
    direct["callback_mipsol_checks"] = 0
    direct["lazy_cuts"] = 0
    return direct


def main() -> None:
    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    settings = settings_from_config(config)
    RESULTS_DIR.joinpath("raw").mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.joinpath("summary").mkdir(parents=True, exist_ok=True)

    direct_source = pd.read_csv(DIRECT_SOURCE)
    direct_target = direct_source[
        direct_source["variant"].eq("EVRIM-Direct")
        & direct_source["n"].astype(int).isin(TARGET_SIZES)
        & direct_source["regime"].isin(["high", "cyclic"])
    ].copy()
    counts = direct_target.groupby(["n", "regime"])["instance_id"].nunique()
    if len(direct_target) != 18 or not (counts == 3).all():
        raise RuntimeError("The archived Table 7 instance set is incomplete or has changed.")

    datasets = pd.read_csv(DATASET_SOURCE).set_index("instance_id", drop=False)
    target_ids = set(direct_target["instance_id"].astype(str))
    missing = sorted(target_ids.difference(datasets.index.astype(str)))
    if missing:
        raise RuntimeError(f"Missing Table 7 matrices: {missing}")

    completed = completed_keys(RAW_PATH)
    for instance_id in sorted(target_ids, key=lambda key: (int(datasets.loc[key, "n"]), key)):
        record = datasets.loc[instance_id]
        matrix = load_record_matrix(record["matrix"])
        value_set, direction_set = default_protected_sets(matrix)
        variants = [
            ("Check first OA", [], []),
            ("Check first OA with protected judgments", value_set, direction_set),
        ]
        for variant, protected_values, protected_directions in variants:
            if (instance_id, variant) in completed:
                continue
            result = solve_evrim_check_first(
                matrix,
                settings,
                value_protected=protected_values,
                direction_protected=protected_directions,
                variant=variant,
            )
            revised = result.revised_matrix
            row = {
                "instance_id": instance_id,
                "n": int(record["n"]),
                "regime": str(record["regime"]),
                "replicate": int(record["replicate"]),
                "variant": variant,
                "status": result.status,
                "certified": bool(result.solved),
                "early_stop": bool(result.early_stop),
                "nrp": result.nrp,
                "aoc": result.aoc,
                "cr": None if revised is None else consistency_ratio(revised),
                "gci": result.gci,
                "gci_threshold": result.gci_threshold,
                "nv": result.nv,
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
                "stage1_runtime": result.stage1_runtime,
                "stage2_runtime": result.stage2_runtime,
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
                "node_count": result.node_count,
                "solution_count": result.solution_count,
                "work": result.work,
                "max_violation": result.max_violation,
                "value_protected": json.dumps(result.value_protected),
                "direction_protected": json.dumps(result.direction_protected),
                "revised_matrix": json_value(revised),
                "certificate_weights": json_value(result.weights),
            }
            append_row(RAW_PATH, row)
            completed.add((instance_id, variant))
            print(
                f"{instance_id} | {variant} | {result.status} | "
                f"early_stop={result.early_stop} | {result.runtime:.3f}s",
                flush=True,
            )

    rerun = pd.read_csv(RAW_PATH)
    rerun["certified"] = rerun["certified"].astype(str).str.lower().eq("true")
    direct = build_direct_rows(direct_source, target_ids)
    combined = pd.concat([direct, rerun], ignore_index=True, sort=False)
    order = {
        "Direct": 0,
        "Check first OA": 1,
        "Check first OA with protected judgments": 2,
    }
    combined["_variant_order"] = combined["variant"].map(order)
    combined = combined.sort_values(["n", "_variant_order", "instance_id"])
    stats = summarize_runtime(combined)
    stats["_variant_order"] = stats["Variant"].map(order)
    stats = stats.sort_values(["n", "_variant_order"]).drop(columns="_variant_order")
    stats.to_csv(FULL_STATS_PATH, index=False)
    stats[["n", "Variant", "Mean", "Median", "Max", "Certified (%)"]].to_csv(
        SUMMARY_PATH, index=False
    )

    environment = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "gurobi": ".".join(map(str, gp.gurobi.version())),
        "config": config,
        "direct_source": str(DIRECT_SOURCE),
        "dataset_source": str(DATASET_SOURCE),
        "target_sizes": list(TARGET_SIZES),
        "algorithm": "certified GCI-free lexicographic solve followed by B&C/OA only if required",
    }
    (RESULTS_DIR / "environment.json").write_text(
        json.dumps(environment, indent=2), encoding="utf-8"
    )
    print("\nUpdated Table 7 summary", flush=True)
    print(pd.read_csv(SUMMARY_PATH).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
