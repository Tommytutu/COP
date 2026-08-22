from __future__ import annotations

import argparse
import json
import platform
import shutil
import subprocess
import sys
import time
from dataclasses import replace
from pathlib import Path

import gurobipy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from cop_experiments.metrics import recovery_metrics  # noqa: E402
from cop_experiments.pcm import (  # noqa: E402
    llsm_weights,
    load_record_matrix,
    load_record_weights,
)
from cop_experiments.priority import (  # noqa: E402
    GurobiSettings,
    alpha_grid,
    solve_alpha_mnvdm,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the 101-point alpha-MNVDM Pareto experiment.")
    parser.add_argument(
        "--config",
        type=Path,
        default=ROOT / "config_weighted_pareto_20260819.json",
    )
    parser.add_argument(
        "--mode",
        choices=("pilot", "n3_n6_full", "n7_full", "selected", "full"),
        default="pilot",
        help=(
            "pilot: one PCM for n=7,8,9; n3_n6_full: all 320 PCMs for n=3,...,6; "
            "n7_full: all 80 n=7 PCMs; "
            "selected: all n=7 plus three PCMs per regime for n=8 and n=9; "
            "full: all 560 PCMs"
        ),
    )
    parser.add_argument("--limit-alpha", type=int, default=None, help="Testing only: keep the first k alphas.")
    parser.add_argument(
        "--deciles-only",
        action="store_true",
        help="Run only the eleven reported values alpha=1,0.9,...,0.",
    )
    parser.add_argument(
        "--rerun-alpha1",
        action="store_true",
        help="Archive and recompute alpha=1 rows for the instances selected by --mode.",
    )
    return parser.parse_args()


def append_csv(path: Path, row: dict) -> None:
    pd.DataFrame([row]).to_csv(path, mode="a", header=not path.exists(), index=False)


def load_completed(path: Path) -> dict[tuple[str, str], dict]:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    frame = pd.read_csv(path)
    # CSV readers normalize values such as ``1.00`` to ``1.0``. Rebuild the
    # checkpoint key from numeric alpha so resumed runs skip completed points.
    frame["alpha_key"] = frame["alpha"].map(lambda value: f"{float(value):.2f}")
    return {
        (str(row["instance_id"]), str(row["alpha_key"])): row.to_dict()
        for _, row in frame.iterrows()
    }


def array_json(value: np.ndarray | list[int] | None) -> str:
    if value is None:
        return ""
    return json.dumps(np.asarray(value).tolist(), separators=(",", ":"))


def optional_array(value: object) -> np.ndarray | None:
    if value is None or (isinstance(value, float) and np.isnan(value)) or value == "":
        return None
    return np.asarray(json.loads(str(value)), dtype=float)


def optional_states(value: object) -> list[int] | None:
    if value is None or (isinstance(value, float) and np.isnan(value)) or value == "":
        return None
    return [int(item) for item in json.loads(str(value))]


def gci_log_deviation(a: np.ndarray, weights: np.ndarray) -> float:
    """Recompute the manuscript's GCI-normalized log-fit deviation."""
    n = a.shape[0]
    y = np.log(np.asarray(weights, dtype=float))
    squared_residual = sum(
        (float(np.log(a[i, j])) - y[i] + y[j]) ** 2
        for i in range(n) for j in range(i + 1, n)
    )
    return 2.0 * squared_residual / ((n - 1) * (n - 2))


def write_method_comparison_tables(
    results: Path,
    flagged: pd.DataFrame,
    summary_dir: Path,
) -> None:
    """Build EM and decile-alpha MNVDM tables for every available size."""
    priority_path = ROOT / "results_model_aligned_20260819" / "raw" / "priority_results.csv"
    instances_path = results / "raw" / "instances.csv"
    if not priority_path.exists() or not instances_path.exists():
        return

    selected_ids = set(flagged["instance_id"].astype(str))
    instances = pd.read_csv(instances_path)
    matrices = {
        str(row["instance_id"]): load_record_matrix(row["matrix"])
        for _, row in instances[instances["instance_id"].astype(str).isin(selected_ids)].iterrows()
    }
    priority = pd.read_csv(priority_path)
    priority = priority[
        priority["instance_id"].astype(str).isin(selected_ids)
        & priority["method"].isin(["EM", "LLSM"])
        & priority["solved"].astype(bool)
    ].copy()
    if len(priority) != 2 * len(selected_ids):
        raise RuntimeError(
            f"Expected {2 * len(selected_ids)} solved EM/LLSM rows, found {len(priority)}"
        )

    baseline_rows: list[dict] = []
    for _, row in priority.iterrows():
        instance_id = str(row["instance_id"])
        weights = np.asarray(json.loads(str(row["weights"])), dtype=float)
        gci_value = gci_log_deviation(matrices[instance_id], weights)
        baseline_rows.append({
            "instance_id": instance_id,
            "n": int(row["n"]),
            "regime": str(row["regime"]),
            "method": str(row["method"]),
            "nvr": float(row["nvr"]),
            "gci_deviation": gci_value,
            "reported_gci_input": float(row["gci_input"]),
        })
    baseline_frame = pd.DataFrame(baseline_rows).sort_values(
        ["regime", "instance_id", "method"]
    )
    baseline_frame.to_csv(summary_dir / "alpha_mnvdm_baseline_rows_all_n.csv", index=False)

    llsm_check = baseline_frame[baseline_frame["method"] == "LLSM"]
    max_llsm_residual = float(
        (llsm_check["gci_deviation"] - llsm_check["reported_gci_input"]).abs().max()
    )
    if max_llsm_residual > 1e-10:
        raise RuntimeError(f"LLSM/GCI definition audit failed: {max_llsm_residual}")

    baseline_summary = baseline_frame.groupby(["n", "regime", "method"], as_index=False).agg(
        instances=("instance_id", "nunique"),
        mean_nvr=("nvr", "mean"),
        std_nvr=("nvr", "std"),
        mean_gci=("gci_deviation", "mean"),
        std_gci=("gci_deviation", "std"),
    )
    baseline_summary["setting"] = baseline_summary["method"]

    selected_values = [round(value / 10.0, 1) for value in range(10, -1, -1)]
    endpoint_audit_path = results / "raw" / "alpha1_indicator_tight_lex_endpoint_audit.csv"
    if not endpoint_audit_path.exists():
        raise RuntimeError("alpha=1 indicator endpoint audit is missing")
    endpoint_audit = pd.read_csv(endpoint_audit_path)
    endpoint_audit = endpoint_audit[
        endpoint_audit["instance_id"].astype(str).isin(selected_ids)
    ].copy()
    if endpoint_audit["instance_id"].nunique() != len(selected_ids):
        raise RuntimeError("alpha=1 indicator endpoint audit is incomplete")
    canonical_endpoint_rows: list[dict] = []
    for _, endpoint in endpoint_audit.iterrows():
        instance_id = str(endpoint["instance_id"])
        minimum_nvr = float(endpoint["nvr"])
        candidates = flagged[
            (flagged["instance_id"].astype(str) == instance_id)
            & flagged["python_feasible"].astype(bool)
            & np.isclose(flagged["nvr"].astype(float), minimum_nvr, atol=1e-9)
        ].sort_values("gci_deviation")
        audited_gci = float(endpoint["gci_deviation"])
        if candidates.empty or audited_gci <= float(candidates.iloc[0]["gci_deviation"]):
            representative_gci = audited_gci
            representative_source = "indicator_min_nv_gci_tiebreak"
        else:
            representative_gci = float(candidates.iloc[0]["gci_deviation"])
            representative_source = f"best_path_tie_alpha_{float(candidates.iloc[0]['alpha']):.2f}"
        canonical_endpoint_rows.append({
            "instance_id": instance_id,
            "n": int(endpoint["n"]),
            "regime": str(endpoint["regime"]),
            "replicate": int(endpoint["replicate"]),
            "alpha": 1.0,
            "nvr": minimum_nvr,
            "gci_deviation": representative_gci,
            "certified": bool(endpoint["stage1_certified"] and endpoint["stage1_audit_pass"]),
            "python_feasible": True,
            "status": str(endpoint["stage1_status"]),
            "runtime": float(endpoint["total_runtime"]),
            "endpoint_representative_source": representative_source,
        })
    canonical_endpoint = pd.DataFrame(canonical_endpoint_rows)
    canonical_endpoint.to_csv(
        summary_dir / "alpha_mnvdm_alpha1_canonical_rows.csv", index=False
    )

    alpha_rows = []
    for (n, regime), point_rows in canonical_endpoint.groupby(["n", "regime"]):
        alpha_rows.append({
            "n": int(n),
            "regime": str(regime),
            "method": "MNVDM",
            "setting": "MNVDM alpha=1.0",
            "instances": int(point_rows["instance_id"].nunique()),
            "mean_nvr": float(point_rows["nvr"].mean()),
            "std_nvr": float(point_rows["nvr"].std()),
            "mean_gci": float(point_rows["gci_deviation"].mean()),
            "std_gci": float(point_rows["gci_deviation"].std()),
        })
    for (n, regime), group in flagged.groupby(["n", "regime"]):
        for alpha in selected_values[1:-1]:
            point_rows = group[np.isclose(group["alpha"], alpha)]
            if point_rows.empty:
                continue
            alpha_rows.append({
                "n": int(n),
                "regime": str(regime),
                "method": "MNVDM",
                "setting": f"MNVDM alpha={alpha:.1f}",
                "instances": int(point_rows["instance_id"].nunique()),
                "mean_nvr": float(point_rows["nvr"].mean()),
                "std_nvr": float(point_rows["nvr"].std()),
                "mean_gci": float(point_rows["gci_deviation"].mean()),
                "std_gci": float(point_rows["gci_deviation"].std()),
            })

    # At alpha=0 the scalarized problem is exactly the unconstrained LLSM
    # log-least-squares problem.  Use the independently computed LLSM row so
    # that NVR is evaluated from the LLSM weights rather than from arbitrary
    # zero-cost binary relation states in the degenerate MIQP endpoint.
    llsm_endpoint = baseline_summary[baseline_summary["method"] == "LLSM"].copy()
    llsm_endpoint["method"] = "MNVDM"
    llsm_endpoint["setting"] = "MNVDM alpha=0.0 (LLSM)"
    alpha_frame = pd.concat([pd.DataFrame(alpha_rows), llsm_endpoint], ignore_index=True)
    comparison = pd.concat([
        baseline_summary[baseline_summary["method"] == "EM"],
        alpha_frame,
    ], ignore_index=True)

    llsm_reference = baseline_summary[baseline_summary["method"] == "LLSM"].set_index(["n", "regime"])
    comparison["delta_nvr_vs_llsm"] = comparison.apply(
        lambda row: float(row["mean_nvr"] - llsm_reference.at[(row["n"], row["regime"]), "mean_nvr"]), axis=1
    )
    comparison["delta_gci_vs_llsm"] = comparison.apply(
        lambda row: float(row["mean_gci"] - llsm_reference.at[(row["n"], row["regime"]), "mean_gci"]), axis=1
    )
    setting_order = {"EM": 0}
    setting_order.update({f"MNVDM alpha={alpha:.1f}": idx + 1 for idx, alpha in enumerate(selected_values[:-1])})
    setting_order["MNVDM alpha=0.0 (LLSM)"] = 11
    comparison["_order"] = comparison["setting"].map(setting_order)
    comparison = comparison.sort_values(["n", "regime", "_order"]).drop(columns="_order")
    comparison.to_csv(summary_dir / "alpha_mnvdm_decile_comparison_all_n.csv", index=False)

    reported_mask = flagged["alpha"].map(
        lambda value: any(np.isclose(float(value), target) for target in selected_values[1:-1])
    )
    llsm_points = baseline_frame[baseline_frame["method"] == "LLSM"].copy()
    llsm_points["alpha"] = 0.0
    llsm_points["certified"] = True
    llsm_points["python_feasible"] = True
    llsm_points["status"] = "OPTIMAL"
    llsm_points["runtime"] = 0.0
    reported_columns = [
        "instance_id", "n", "regime", "alpha", "nvr", "gci_deviation",
        "certified", "python_feasible", "status", "runtime",
    ]
    reported_base = flagged.loc[reported_mask, reported_columns].copy()
    reported_flagged = add_pareto_flags(pd.concat([
        canonical_endpoint[reported_columns],
        reported_base,
        llsm_points[reported_columns],
    ], ignore_index=True))
    reported_flagged.to_csv(
        summary_dir / "alpha_mnvdm_reported_decile_points.csv", index=False
    )
    supported = reported_flagged.groupby(["n", "regime", "instance_id"], as_index=False).agg(
        supported_points=("supported_pareto_flag", "sum")
    ).groupby(["n", "regime"], as_index=False).agg(
        mean_supported_points=("supported_points", "mean")
    )
    wide_rows = []
    for (n, regime), group in comparison.groupby(["n", "regime"]):
        row = {
            "n": int(n),
            "regime": str(regime),
            "instances": int(group["instances"].max()),
        }
        for _, item in group.iterrows():
            setting = str(item["setting"])
            if setting == "EM":
                key = "em"
            else:
                alpha_text = setting.split("alpha=", 1)[1].split()[0]
                key = "alpha_" + alpha_text.replace(".", "p")
            row[f"{key}_nvr"] = float(item["mean_nvr"])
            row[f"{key}_gci"] = float(item["mean_gci"])
        match = supported[(supported["n"] == n) & (supported["regime"] == regime)]
        row["mean_supported_points"] = float(match.iloc[0]["mean_supported_points"])
        wide_rows.append(row)
    pd.DataFrame(wide_rows).sort_values(["n", "regime"]).to_csv(
        summary_dir / "alpha_mnvdm_decile_tables_wide.csv", index=False
    )

    sample_counts = baseline_frame.groupby(["n", "regime", "method"]).size()
    (summary_dir / "alpha_mnvdm_baseline_audit_all_n.json").write_text(
        json.dumps({
            "instances": len(selected_ids),
            "baseline_rows": len(baseline_frame),
            "minimum_rows_per_size_regime_method": int(sample_counts.min()),
            "maximum_rows_per_size_regime_method": int(sample_counts.max()),
            "max_llsm_gci_definition_residual": max_llsm_residual,
            "reported_alpha_values": selected_values,
            "alpha_zero_source": "independently computed LLSM weights",
            "alpha_one_source": "warm-started indicator Stage-1 with internal FeasibilityTol=1e-7 and independent NV audit; indicator GCI tie-break and path candidates supply the best feasible minimum-NVR representative",
            "supported_point_grid": "alpha=1,0.9,...,0",
        }, indent=2),
        encoding="utf-8",
    )

    # Preserve the historical n=7 filename for downstream users while giving
    # it the new requested EM + alpha-decile content.
    comparison[comparison["n"] == 7].to_csv(
        summary_dir / "alpha_mnvdm_table5_n7.csv", index=False
    )


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


def write_environment(results: Path, config: dict, mode: str) -> None:
    try:
        processor = subprocess.check_output(
            ["powershell", "-NoProfile", "-Command", "(Get-CimInstance Win32_Processor).Name"],
            text=True,
            timeout=10,
        ).strip()
    except (OSError, subprocess.SubprocessError):
        processor = platform.processor()
    payload = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "mode": mode,
        "platform": platform.platform(),
        "python": platform.python_version(),
        "processor": processor,
        "cpu_count": __import__("os").cpu_count(),
        "gurobi": ".".join(map(str, gp.gurobi.version())),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "config": config,
    }
    (results / "environment.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")


def add_pareto_flags(frame: pd.DataFrame) -> pd.DataFrame:
    output = frame.copy()
    output["dominated_flag"] = False
    output["duplicate_point_flag"] = False
    output["supported_pareto_flag"] = False
    tolerance = 1e-9
    for _instance_id, indices in output.groupby("instance_id").groups.items():
        certified_indices = [
            idx for idx in indices
            if bool(output.at[idx, "certified"])
            and pd.notna(output.at[idx, "nvr"])
            and pd.notna(output.at[idx, "gci_deviation"])
        ]
        seen: list[tuple[float, float]] = []
        for idx in certified_indices:
            point = (float(output.at[idx, "nvr"]), float(output.at[idx, "gci_deviation"]))
            duplicate = any(abs(point[0] - old[0]) <= tolerance and abs(point[1] - old[1]) <= tolerance for old in seen)
            output.at[idx, "duplicate_point_flag"] = duplicate
            if not duplicate:
                seen.append(point)
            dominated = any(
                float(output.at[other, "nvr"]) <= point[0] + tolerance
                and float(output.at[other, "gci_deviation"]) <= point[1] + tolerance
                and (
                    float(output.at[other, "nvr"]) < point[0] - tolerance
                    or float(output.at[other, "gci_deviation"]) < point[1] - tolerance
                )
                for other in certified_indices if other != idx
            )
            output.at[idx, "dominated_flag"] = dominated
            output.at[idx, "supported_pareto_flag"] = not dominated and not duplicate
    return output


def summarize(results: Path, raw_path: Path) -> None:
    if not raw_path.exists():
        return
    summary_dir = results / "summary"
    figures_dir = results / "figures"
    summary_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)
    frame = pd.read_csv(raw_path, dtype={"alpha_key": str}).sort_values(
        ["n", "regime", "replicate", "alpha"], ascending=[True, True, True, False]
    ).reset_index(drop=True)
    endpoint_audit_path = results / "raw" / "alpha1_indicator_tight_lex_endpoint_audit.csv"
    if endpoint_audit_path.exists():
        endpoint_audit = pd.read_csv(endpoint_audit_path).set_index("instance_id")
        for instance_id, group in frame.groupby("instance_id"):
            if instance_id not in endpoint_audit.index:
                continue
            endpoint = endpoint_audit.loc[instance_id]
            if not bool(endpoint["stage1_certified"] and endpoint["stage1_audit_pass"]):
                continue
            best_path_nvr = float(group.loc[group["python_feasible"].astype(bool), "nvr"].min())
            if float(endpoint["nvr"]) > best_path_nvr + 1e-9:
                raise RuntimeError(
                    f"indicator alpha=1 endpoint failed the NVR audit for {instance_id}: "
                    f"{float(endpoint['nvr'])} > {best_path_nvr}"
                )
    flagged = add_pareto_flags(frame).sort_values(
        ["n", "regime", "replicate", "alpha"], ascending=[True, True, True, False]
    )
    flagged.to_csv(summary_dir / "alpha_mnvdm_all_points.csv", index=False)
    flagged[flagged["supported_pareto_flag"]].to_csv(
        summary_dir / "alpha_mnvdm_supported_pareto.csv", index=False
    )
    runtime = flagged.groupby(["instance_id", "n"], as_index=False).agg(
        runs=("alpha", "size"),
        certified=("certified", "sum"),
        solver_incumbents=("has_solution", "sum"),
        python_feasible_incumbents=("python_feasible", "sum"),
        mean_runtime=("runtime", "mean"),
        median_runtime=("runtime", "median"),
        max_runtime=("runtime", "max"),
        total_runtime=("runtime", "sum"),
        supported_points=("supported_pareto_flag", "sum"),
    )
    runtime.to_csv(summary_dir / "alpha_mnvdm_runtime_summary.csv", index=False)

    runtime_distribution = flagged.groupby(["n", "regime"], as_index=False).agg(
        runs=("alpha", "size"),
        instances=("instance_id", "nunique"),
        mean_runtime=("runtime", "mean"),
        median_runtime=("runtime", "median"),
        std_runtime=("runtime", "std"),
        p90_runtime=("runtime", lambda values: values.quantile(0.90)),
        p95_runtime=("runtime", lambda values: values.quantile(0.95)),
        max_runtime=("runtime", "max"),
        certified_rate=("certified", "mean"),
        timeout_rate=("status", lambda values: (values == "TIME_LIMIT").mean()),
        feasible_incumbent_rate=("python_feasible", "mean"),
    )
    runtime_distribution.to_csv(
        summary_dir / "alpha_mnvdm_runtime_by_n_regime.csv", index=False
    )
    pareto_counts = flagged.groupby(
        ["instance_id", "n", "regime", "replicate"], as_index=False
    ).agg(
        runs=("alpha", "size"),
        certified=("certified", "sum"),
        supported_points=("supported_pareto_flag", "sum"),
        dominated_points=("dominated_flag", "sum"),
        duplicate_points=("duplicate_point_flag", "sum"),
    )
    pareto_counts.to_csv(summary_dir / "alpha_mnvdm_pareto_counts.csv", index=False)
    n7_alpha = flagged[flagged["n"] == 7].groupby(
        ["regime", "alpha"], as_index=False
    ).agg(
        instances=("instance_id", "nunique"),
        certified_rate=("certified", "mean"),
        mean_nvr=("nvr", "mean"),
        median_nvr=("nvr", "median"),
        mean_gci=("gci_deviation", "mean"),
        median_gci=("gci_deviation", "median"),
        mean_runtime=("runtime", "mean"),
    )
    n7_alpha.to_csv(summary_dir / "alpha_mnvdm_alpha_summary_n7.csv", index=False)
    write_method_comparison_tables(results, flagged, summary_dir)

    # Treat movements smaller than the declared optimization tolerance as
    # numerical ties.  A stricter audit would flag harmless differences among
    # alternative optima (for example, GCI changes of order 1e-7).
    monotonicity_tolerance = 1e-5
    monotonic_violations: list[dict] = []
    for instance_id, group in flagged.groupby("instance_id"):
        certified = group[
            group["certified"].astype(bool) & group["python_feasible"].astype(bool)
        ].sort_values("alpha", ascending=False)
        previous = None
        for _, point in certified.iterrows():
            if previous is not None and float(previous["alpha"]) < 0.999999:
                if float(point["nvr"]) < float(previous["nvr"]) - monotonicity_tolerance:
                    monotonic_violations.append({
                        "instance_id": instance_id,
                        "criterion": "NVR",
                        "alpha_high": float(previous["alpha"]),
                        "alpha_low": float(point["alpha"]),
                    })
                if float(point["gci_deviation"]) > float(previous["gci_deviation"]) + monotonicity_tolerance:
                    monotonic_violations.append({
                        "instance_id": instance_id,
                        "criterion": "GCI",
                        "alpha_high": float(previous["alpha"]),
                        "alpha_low": float(point["alpha"]),
                    })
            previous = point
    pd.DataFrame(monotonic_violations).to_csv(
        summary_dir / "alpha_mnvdm_monotonicity_audit.csv", index=False
    )
    audit = {
        "rows": int(len(flagged)),
        "unique_instance_alpha": int(flagged[["instance_id", "alpha"]].drop_duplicates().shape[0]),
        "instances": int(flagged["instance_id"].nunique()),
        "certified": int(flagged["certified"].sum()),
        "python_feasible": int(flagged["python_feasible"].sum()),
        "timeouts": int((flagged["status"] == "TIME_LIMIT").sum()),
        "endpoint_fallbacks": int(flagged["endpoint_fallback_used"].fillna(False).sum()),
        "monotonicity_tolerance": monotonicity_tolerance,
        "monotonicity_violations": int(len(monotonic_violations)),
        "max_certified_objective_residual": float(
            flagged.loc[flagged["certified"].astype(bool), "objective_residual"].dropna().max()
        ),
    }
    (summary_dir / "alpha_mnvdm_audit_summary.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )

    figure, axis = plt.subplots(figsize=(7.0, 4.8))
    representative_ids = set(
        flagged.sort_values(["n", "regime", "replicate"])
        .groupby(["n", "regime"], as_index=False).first()["instance_id"]
    )
    for instance_id, group in flagged[flagged["instance_id"].isin(representative_ids)].groupby("instance_id"):
        valid = group[group["python_feasible"]].sort_values("nvr")
        first = valid.iloc[0]
        legend_label = f"{str(first['regime']).capitalize()}, $n={int(first['n'])}$"
        line = axis.plot(
            valid["nvr"], valid["gci_deviation"], marker="o",
            markersize=2.5, linewidth=1.0, label=legend_label,
        )[0]
        uncertified = valid[~valid["certified"].astype(bool)]
        if not uncertified.empty:
            axis.scatter(
                uncertified["nvr"], uncertified["gci_deviation"],
                marker="x", s=30, linewidths=1.2, color=line.get_color(), zorder=4,
            )
    axis.set_xlabel("Normalized violation rate (NVR)")
    axis.set_ylabel("GCI-based log deviation")
    axis.grid(alpha=0.25)
    axis.legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(figures_dir / "alpha_mnvdm_pareto_selected.pdf")
    figure.savefig(figures_dir / "alpha_mnvdm_pareto_selected.png", dpi=220)
    plt.close(figure)


def main() -> None:
    args = parse_args()
    config_path = args.config.resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    settings = settings_from_config(config)
    results = ROOT / str(config["results_dir"])
    raw_dir = results / "raw"
    detail_dir = results / "details"
    log_dir = results / "logs"
    staging_log_dir = Path(str(config.get(
        "solver_log_staging_dir",
        ROOT / "solver_log_staging",
    )))
    for directory in (results, raw_dir, detail_dir, log_dir):
        directory.mkdir(parents=True, exist_ok=True)
    staging_log_dir.mkdir(parents=True, exist_ok=True)
    write_environment(results, config, args.mode)

    dataset_path = ROOT / str(config["source_dataset"])
    data = pd.read_csv(dataset_path)
    if args.mode == "pilot":
        requested = list(config["pilot_instances"])
        data = data[data["instance_id"].isin(requested)].copy()
        missing = sorted(set(requested) - set(data["instance_id"]))
        if missing:
            raise RuntimeError(f"pilot instances missing from source dataset: {missing}")
        data["_order"] = data["instance_id"].map({name: idx for idx, name in enumerate(requested)})
        data = data.sort_values("_order").drop(columns="_order")
    elif args.mode == "n3_n6_full":
        data = data[data["n"].between(3, 6)].sort_values(["n", "regime", "replicate"])
    elif args.mode == "n7_full":
        data = data[data["n"] == 7].sort_values(["regime", "replicate"])
    elif args.mode == "selected":
        requested_large = data[(data["n"].isin([8, 9])) & (data["replicate"] < 3)]["instance_id"]
        data = data[(data["n"] == 7) | data["instance_id"].isin(requested_large)].copy()
        data = data.sort_values(["n", "regime", "replicate"])
    else:
        data = data.sort_values(["n", "regime", "replicate"])

    alpha_spec = config["alpha_values"]
    grid_step = 0.1 if args.deciles_only else float(alpha_spec["step"])
    alphas = [1.0] if args.rerun_alpha1 else alpha_grid(grid_step)
    if args.limit_alpha is not None:
        alphas = alphas[: args.limit_alpha]
    expected_alpha_count = 1 if args.rerun_alpha1 else (11 if args.deciles_only else 101)
    if len(alphas) != expected_alpha_count and args.limit_alpha is None:
        raise RuntimeError(f"expected {expected_alpha_count} alpha values, found {len(alphas)}")
    # The alpha=1 endpoint has no metric-fit coefficient and is highly
    # degenerate. Solve alpha=0.99 first, use its certified state as the
    # endpoint start, and sort the saved tables back to descending alpha.
    if not args.rerun_alpha1 and args.limit_alpha is None and len(alphas) >= 2:
        alphas = [alphas[1], alphas[0], *alphas[2:]]

    raw_path = raw_dir / "alpha_mnvdm_results.csv"
    instances_path = raw_dir / "instances.csv"
    if instances_path.exists():
        existing_instances = pd.read_csv(instances_path)
        pd.concat([existing_instances, data], ignore_index=True).drop_duplicates(
            subset=["instance_id"], keep="first"
        ).sort_values(["n", "regime", "replicate"]).to_csv(instances_path, index=False)
    else:
        data.to_csv(instances_path, index=False)
    if args.rerun_alpha1 and raw_path.exists():
        current_raw = pd.read_csv(raw_path)
        selected_ids = set(data["instance_id"].astype(str))
        replace_mask = (
            current_raw["instance_id"].astype(str).isin(selected_ids)
            & np.isclose(current_raw["alpha"].astype(float), 1.0)
        )
        archived = current_raw[replace_mask].copy()
        archive_path = results / "logs" / f"alpha1_before_strong_fix_{args.mode}.csv"
        archive_index = 1
        while archive_path.exists():
            archive_path = results / "logs" / (
                f"alpha1_before_strong_fix_{args.mode}_{archive_index}.csv"
            )
            archive_index += 1
        archived.to_csv(archive_path, index=False)
        current_raw[~replace_mask].to_csv(raw_path, index=False)
    completed = load_completed(raw_path)
    detail_path = detail_dir / "relation_details.jsonl"
    total = len(data) * len(alphas)
    completed_this_run = 0
    for _, row in data.iterrows():
        instance_id = str(row["instance_id"])
        a = load_record_matrix(row["matrix"])
        w0 = load_record_weights(row["latent_weights"])
        llsm = llsm_weights(a)
        previous_y: np.ndarray | None = None
        previous_states: list[int] | None = None
        previous_alpha: float | None = None
        instance_log_dir = log_dir / instance_id
        instance_log_dir.mkdir(exist_ok=True)
        for alpha in alphas:
            alpha_key = f"{alpha:.2f}"
            key = (instance_id, alpha_key)
            if key in completed:
                cached = completed[key]
                if bool(cached.get("has_solution", False)):
                    previous_y = optional_array(cached.get("y"))
                    previous_states = optional_states(cached.get("relation_states"))
                    previous_alpha = alpha
                continue
            warm_source = "LLSM" if previous_y is None else f"alpha_{previous_alpha:.2f}"
            log_file = None
            final_log_file = instance_log_dir / f"alpha_{alpha_key}.log"
            if bool(config.get("save_solver_logs", False)):
                staged_instance_dir = staging_log_dir / instance_id
                staged_instance_dir.mkdir(parents=True, exist_ok=True)
                staged_log_file = staged_instance_dir / f"alpha_{alpha_key}.log"
                log_file = str(staged_log_file)
            result = solve_alpha_mnvdm(
                a=a,
                alpha=alpha,
                settings=settings,
                warm_start_y=previous_y,
                warm_start_states=previous_states,
                warm_start_source=warm_source,
                gci_normalizer=float(config.get("gci_normalizer", 1.0)),
                log_file=log_file,
            )
            preliminary_feasible = bool(
                result.has_solution
                and (result.min_strict_slack is None or result.min_strict_slack >= -settings.feasibility_tol)
                and (result.max_equality_residual is None or result.max_equality_residual <= settings.feasibility_tol)
            )
            if result.certified and not preliminary_feasible:
                # Gurobi can occasionally return OPTIMAL with a scaled-row
                # residual above the declared independent 1e-5 audit.  Do not
                # relax the audit: re-solve that point with a tighter internal
                # feasibility tolerance and record the retry in the existing
                # warm-start provenance field.
                retry_settings = replace(settings, feasibility_tol=min(settings.feasibility_tol, 1e-7))
                result = solve_alpha_mnvdm(
                    a=a,
                    alpha=alpha,
                    settings=retry_settings,
                    warm_start_y=result.y,
                    warm_start_states=result.relation_states,
                    warm_start_source=f"audit_retry({warm_source})",
                    gci_normalizer=float(config.get("gci_normalizer", 1.0)),
                    log_file=None,
                )
            if log_file and Path(log_file).exists():
                shutil.copy2(log_file, final_log_file)
            metrics = recovery_metrics(a, result.weights, w0) if result.weights is not None else {}
            objective_residual = None
            if result.weighted_objective is not None and result.solver_objective is not None:
                objective_residual = abs(result.weighted_objective - result.solver_objective)
            python_feasible = bool(
                result.has_solution
                and (result.min_strict_slack is None or result.min_strict_slack >= -settings.feasibility_tol)
                and (result.max_equality_residual is None or result.max_equality_residual <= settings.feasibility_tol)
            )
            if result.num_general_constraints != 0:
                raise RuntimeError(
                    "Production alpha-MNVDM unexpectedly contains general/indicator constraints"
                )
            if result.certified and not python_feasible:
                raise RuntimeError(
                    f"Certified incumbent failed independent feasibility audit: "
                    f"{instance_id}, alpha={alpha_key}"
                )
            output_row = {
                "instance_id": instance_id,
                "n": int(row["n"]),
                "regime": row["regime"],
                "replicate": int(row["replicate"]),
                "seed": int(row["seed"]),
                "cr_input": float(row["cr"]),
                "gci_input": float(row["gci"]),
                "alpha": result.alpha,
                "alpha_key": alpha_key,
                "order_formulation": "relation_specific_bigm",
                "status": result.status,
                "status_code": result.status_code,
                "solved": result.solved,
                "certified": result.certified,
                "has_solution": result.has_solution,
                "nv": result.nv,
                "nv2": result.nv2,
                "n_order_relations": result.n_order_relations,
                "nvr": result.nvr,
                "gci_deviation": result.gci_deviation,
                "gci_normalizer": result.gci_normalizer,
                "weighted_objective": result.weighted_objective,
                "solver_objective": result.solver_objective,
                "objective_residual": objective_residual,
                "objective_bound": result.objective_bound,
                "gap": result.gap,
                "runtime": result.runtime,
                "work": result.work,
                "node_count": result.node_count,
                "iteration_count": result.iteration_count,
                "barrier_iteration_count": result.barrier_iteration_count,
                "solution_count": result.solution_count,
                "max_violation": result.max_violation,
                "max_integrality_violation": result.max_integrality_violation,
                "min_strict_slack": result.min_strict_slack,
                "max_equality_residual": result.max_equality_residual,
                "python_feasible": python_feasible,
                "usable_incumbent": python_feasible,
                "warm_start_source": result.warm_start_source,
                "num_variables": result.num_variables,
                "num_binary_variables": result.num_binary_variables,
                "num_constraints": result.num_constraints,
                "num_quadratic_constraints": result.num_quadratic_constraints,
                "num_general_constraints": result.num_general_constraints,
                "endpoint_state_rejections": result.endpoint_state_rejections,
                "endpoint_refinement_runtime": result.endpoint_refinement_runtime,
                "endpoint_fixed_state_status": result.endpoint_fixed_state_status,
                "endpoint_fallback_used": result.endpoint_fallback_used,
                "kendall_tau_b": metrics.get("kendall_tau_b"),
                "best_choice_accuracy": metrics.get("best_choice_accuracy"),
                "lrmse": metrics.get("lrmse"),
                "weights": array_json(result.weights),
                "y": array_json(result.y),
                "ranking": array_json(None if result.weights is None else np.argsort(-result.weights)),
                "relation_states": array_json(result.relation_states),
                "llsm_weights": array_json(llsm),
            }
            append_csv(raw_path, output_row)
            with detail_path.open("a", encoding="utf-8") as handle:
                handle.write(json.dumps({
                    "instance_id": instance_id,
                    "alpha": result.alpha,
                    "relation_details": result.relation_details,
                }, separators=(",", ":")) + "\n")
            completed[key] = output_row
            completed_this_run += 1
            if completed_this_run % 100 == 0 or result.status != "OPTIMAL":
                print(
                    f"[{completed_this_run}/{total}] {instance_id} alpha={alpha_key} "
                    f"status={result.status} runtime={result.runtime:.3f}s "
                    f"NVR={result.nvr} GCI={result.gci_deviation}",
                    flush=True,
                )
            if result.has_solution:
                previous_y = result.y
                previous_states = result.relation_states
                previous_alpha = alpha
    summarize(results, raw_path)
    print(f"results -> {results}", flush=True)


if __name__ == "__main__":
    main()
