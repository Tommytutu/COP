from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _q90(x: pd.Series) -> float:
    return float(x.quantile(0.90))


def _q95(x: pd.Series) -> float:
    return float(x.quantile(0.95))


def _as_bool(series: pd.Series) -> pd.Series:
    """Read Boolean CSV columns robustly (including legacy string exports)."""
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def _runtime_summary(frame: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    rows = []
    for keys, group in frame.groupby(group_cols, dropna=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        solved = group[_as_bool(group["solved"])]
        times = group["runtime"].dropna()
        solved_times = solved["runtime"].dropna()
        row = dict(zip(group_cols, keys, strict=True))
        row.update({
            "instances": len(group), "solved_pct": 100 * len(solved) / len(group),
            "timeout_pct": 100 * (group["status"] == "TIME_LIMIT").mean(),
            "mean": times.mean(), "median": times.median(), "std": times.std(ddof=1),
            "p90": _q90(times) if len(times) else np.nan,
            "p95": _q95(times) if len(times) else np.nan,
            "maximum": times.max(),
            "solved_mean": solved_times.mean(),
        })
        rows.append(row)
    return pd.DataFrame(rows)


def build_all_reports(results: Path) -> None:
    raw, summary, figures = results / "raw", results / "summary", results / "figures"
    summary.mkdir(exist_ok=True, parents=True)
    figures.mkdir(exist_ok=True, parents=True)
    environment_path = results / "environment.json"
    stage_limit = 60.0
    if environment_path.exists():
        environment = json.loads(environment_path.read_text(encoding="utf-8"))
        stage_limit = float(environment.get("config", {}).get(
            "priority_time_limit_seconds", stage_limit
        ))

    dataset_path = raw / "datasets.csv"
    if dataset_path.exists():
        datasets = pd.read_csv(dataset_path)
        distribution = datasets.groupby(["n", "regime", "realized_cr_bin"], dropna=False).agg(
            instances=("instance_id", "count"), cr_mean=("cr", "mean"), cr_std=("cr", "std"),
            gci_mean=("gci", "mean"), gci_std=("gci", "std"),
        ).reset_index()
        distribution.to_csv(summary / "dataset_distribution.csv", index=False)

    priority_path = raw / "priority_results.csv"
    if priority_path.exists():
        priority = pd.read_csv(priority_path)
        solved = priority[_as_bool(priority["solved"])]
        quality = solved.groupby(["n", "regime", "method"]).agg(
            instances=("instance_id", "count"), nvr_mean=("nvr", "mean"), nvr_std=("nvr", "std"),
            kendall_mean=("kendall_tau_b", "mean"), kendall_std=("kendall_tau_b", "std"),
            best_choice_accuracy=("best_choice_accuracy", "mean"),
            lrmse_mean=("lrmse", "mean"), lrmse_std=("lrmse", "std"),
        ).reset_index()
        quality.to_csv(summary / "decision_quality.csv", index=False)
        solved.groupby(["regime", "method"]).agg(
            instances=("instance_id", "count"), nvr_mean=("nvr", "mean"), nvr_std=("nvr", "std"),
            kendall_mean=("kendall_tau_b", "mean"), kendall_std=("kendall_tau_b", "std"),
            best_choice_accuracy=("best_choice_accuracy", "mean"),
            lrmse_mean=("lrmse", "mean"), lrmse_std=("lrmse", "std"),
        ).reset_index().to_csv(summary / "decision_quality_overall.csv", index=False)
        comparison_methods = ["EM", "LLSM", "MNVEM", "MNVLLSM"]
        solved_comparison = solved[solved["method"].isin(comparison_methods)]
        common_ids = (
            solved_comparison.groupby("instance_id")["method"].nunique()
            .loc[lambda x: x == len(comparison_methods)].index
        )
        paired_quality = solved_comparison[solved_comparison["instance_id"].isin(common_ids)]
        paired_quality.groupby(["regime", "method"]).agg(
            instances=("instance_id", "count"), nvr_mean=("nvr", "mean"), nvr_std=("nvr", "std"),
            kendall_mean=("kendall_tau_b", "mean"), kendall_std=("kendall_tau_b", "std"),
            best_choice_accuracy=("best_choice_accuracy", "mean"),
            lrmse_mean=("lrmse", "mean"), lrmse_std=("lrmse", "std"),
        ).reset_index().to_csv(summary / "decision_quality_paired.csv", index=False)

        # A legacy run may have used a longer cap.  This certificate filter is
        # valid because every retained solve had already closed its gap before
        # the current per-stage limit; no post-limit incumbent is reused.
        mnv = solved_comparison[solved_comparison["method"].isin(["MNVEM", "MNVLLSM"])].copy()
        within_limit = mnv[
            (mnv["stage1_runtime"] <= stage_limit)
            & (mnv["stage2_runtime"] <= stage_limit)
        ]
        limited_ids = (
            within_limit.groupby("instance_id")["method"].nunique()
            .loc[lambda x: x == 2].index
        )
        limited_quality = solved_comparison[solved_comparison["instance_id"].isin(limited_ids)]
        limited_quality.groupby(["regime", "method"]).agg(
            instances=("instance_id", "count"), nvr_mean=("nvr", "mean"), nvr_std=("nvr", "std"),
            kendall_mean=("kendall_tau_b", "mean"), kendall_std=("kendall_tau_b", "std"),
            best_choice_accuracy=("best_choice_accuracy", "mean"),
            lrmse_mean=("lrmse", "mean"), lrmse_std=("lrmse", "std"),
        ).reset_index().to_csv(summary / "decision_quality_paired_60s.csv", index=False)
        _runtime_summary(priority, ["n", "method"]).to_csv(summary / "priority_runtime.csv", index=False)

        if dataset_path.exists():
            expected = pd.read_csv(dataset_path)[["instance_id", "n", "regime"]]
            methods = sorted(priority["method"].unique())
            coverage_rows = []
            for method in methods:
                observed = priority[priority["method"] == method][["instance_id", "n", "regime", "solved", "status"]]
                merged = expected.merge(observed, on=["instance_id", "n", "regime"], how="left")
                for (n, regime), group in merged.groupby(["n", "regime"]):
                    attempted = group["status"].notna()
                    certified = attempted & _as_bool(group["solved"])
                    coverage_rows.append({
                        "n": n, "regime": regime, "method": method,
                        "expected": len(group), "attempted": int(attempted.sum()),
                        "missing": int((~attempted).sum()), "certified": int(certified.sum()),
                        "attempted_pct": 100 * attempted.mean(),
                        "certified_pct_of_expected": 100 * certified.mean(),
                    })
            pd.DataFrame(coverage_rows).to_csv(summary / "execution_coverage.csv", index=False)

    stage1_path = raw / "stage1_cache.csv"
    if stage1_path.exists():
        stage1 = pd.read_csv(stage1_path)
        feasibility_rows = []
        for (n, regime), group in stage1.groupby(["n", "regime"]):
            certified = group[_as_bool(group["solved"])]
            feasibility_rows.append({
                "n": n, "regime": regime, "instances": len(group),
                "certified_instances": len(certified),
                "certified_pct": 100 * len(certified) / len(group),
                "p_nv_star_zero": float((certified["nv"] == 0).mean()) if len(certified) else np.nan,
                "nv_star_mean": certified["nv"].mean(),
            })
        feasibility = pd.DataFrame(feasibility_rows)
        feasibility.to_csv(summary / "representability.csv", index=False)

        certified_60 = stage1[
            _as_bool(stage1["solved"]) & (stage1["runtime"] <= stage_limit)
        ]
        representability_60 = certified_60.groupby("regime").agg(
            certified_instances=("instance_id", "count"),
            zero_nv_instances=("nv", lambda x: int((x == 0).sum())),
            p_nv_star_zero=("nv", lambda x: float((x == 0).mean())),
            nv_star_mean=("nv", "mean"),
        ).reset_index()
        representability_60["stage_limit_seconds"] = stage_limit
        representability_60.to_csv(summary / "representability_60s.csv", index=False)
        order = ["low", "moderate", "high", "cyclic"]
        plot_data = representability_60.set_index("regime").reindex(order).reset_index()
        fig, ax = plt.subplots(figsize=(6.5, 4.0))
        bars = ax.bar(plot_data["regime"], plot_data["p_nv_star_zero"],
                      color=["#4C78A8", "#F2CF5B", "#E45756", "#7A5195"])
        ax.set(ylabel=r"$P(NV^*=0)$", ylim=(0.0, 0.24))
        for bar, value in zip(bars, plot_data["p_nv_star_zero"], strict=True):
            ax.text(bar.get_x() + bar.get_width() / 2, value + 0.006,
                    f"{value:.3f}", ha="center", va="bottom", fontsize=9)
        fig.tight_layout()
        fig.savefig(figures / "representability_probability_60s.pdf")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(7, 4.5))
        for regime, group in feasibility.groupby("regime"):
            ax.plot(group["n"], group["p_nv_star_zero"], marker="o", label=regime)
        ax.set(xlabel="n", ylabel=r"$P(NV^*=0)$", ylim=(-0.03, 1.03))
        ax.legend(frameon=False, ncol=2)
        fig.tight_layout()
        fig.savefig(figures / "representability_probability.pdf")
        plt.close(fig)

    formulation_path = raw / "formulation_runtime.csv"
    if formulation_path.exists():
        formulation = pd.read_csv(formulation_path)
        _runtime_summary(formulation, ["n", "variant"]).to_csv(summary / "formulation_runtime.csv", index=False)
        paired = formulation.pivot(index="instance_id", columns="variant", values=["runtime", "solved"])
        if ("runtime", "basic") in paired and ("runtime", "strong") in paired:
            paired.columns = ["_".join(map(str, col)) for col in paired.columns]
            strong_runtime = pd.to_numeric(paired["runtime_strong"], errors="coerce").replace(0.0, np.nan)
            basic_runtime = pd.to_numeric(paired["runtime_basic"], errors="coerce")
            paired["speedup_basic_over_strong"] = basic_runtime / strong_runtime
            paired.reset_index().to_csv(summary / "formulation_paired.csv", index=False)

    evrim_path = raw / "evrim_results.csv"
    if evrim_path.exists():
        evrim = pd.read_csv(evrim_path)
        _runtime_summary(evrim, ["n", "regime", "variant"]).to_csv(summary / "evrim_runtime.csv", index=False)
        solved_e = evrim[_as_bool(evrim["solved"])]
        solved_e.groupby(["n", "regime", "variant"]).agg(
            instances=("instance_id", "count"), nrp_mean=("nrp", "mean"), nrp_std=("nrp", "std"),
            aoc_mean=("aoc", "mean"), aoc_std=("aoc", "std"),
            gci_mean=("gci", "mean"), nv_mean=("nv", "mean"), runtime_mean=("runtime", "mean"),
        ).reset_index().to_csv(summary / "evrim_quality.csv", index=False)
        variant_count = evrim["variant"].nunique()
        paired_ids = (
            solved_e.groupby("instance_id")["variant"].nunique()
            .loc[lambda x: x == variant_count].index
        )
        solved_e[solved_e["instance_id"].isin(paired_ids)].groupby(["regime", "variant"]).agg(
            instances=("instance_id", "count"), nrp_mean=("nrp", "mean"),
            aoc_mean=("aoc", "mean"), gci_mean=("gci", "mean"),
            nv_mean=("nv", "mean"), runtime_mean=("runtime", "mean"),
        ).reset_index().to_csv(summary / "evrim_paired.csv", index=False)

        certificate_variants = ["EVRIM-Direct", "EVRIM-OA"]
        certificate = evrim[evrim["variant"].isin(certificate_variants)].copy()
        certificate_solved = certificate[_as_bool(certificate["solved"])]
        common = (
            certificate_solved.groupby("instance_id")["variant"].nunique()
            .loc[lambda x: x == len(certificate_variants)].index
        )
        certificate = certificate_solved[certificate_solved["instance_id"].isin(common)]
        if len(certificate):
            wide = certificate.pivot(
                index=["instance_id", "n", "regime"], columns="variant",
                values=["nrp", "aoc", "gci", "nv", "runtime"],
            )
            wide.columns = ["_".join(map(str, column)) for column in wide.columns]
            wide["nrp_agree"] = wide["nrp_EVRIM-Direct"] == wide["nrp_EVRIM-OA"]
            wide["aoc_abs_difference"] = (
                wide["aoc_EVRIM-Direct"] - wide["aoc_EVRIM-OA"]
            ).abs()
            wide["speedup_direct_over_oa"] = (
                wide["runtime_EVRIM-Direct"]
                / wide["runtime_EVRIM-OA"].replace(0.0, np.nan)
            )
            wide.reset_index().to_csv(summary / "evrim_direct_oa_paired.csv", index=False)

    bnc_path = raw / "bnc_runtime.csv"
    if bnc_path.exists():
        bnc = pd.read_csv(bnc_path)
        _runtime_summary(bnc, ["n", "backend"]).to_csv(summary / "bnc_runtime.csv", index=False)
        paired = bnc.pivot(index="instance_id", columns="backend", values=["runtime", "solved"])
        if ("runtime", "direct") in paired and ("runtime", "oa_callback") in paired:
            paired.columns = ["_".join(map(str, col)) for col in paired.columns]
            paired["speedup_direct_over_callback"] = (
                pd.to_numeric(paired["runtime_direct"], errors="coerce")
                / pd.to_numeric(paired["runtime_oa_callback"], errors="coerce").replace(0.0, np.nan)
            )
            paired.reset_index().to_csv(summary / "bnc_paired.csv", index=False)

    sensitivity_path = raw / "epsilon_sensitivity.csv"
    if sensitivity_path.exists():
        sensitivity = pd.read_csv(sensitivity_path)
        rows = []
        for (n, regime), group in sensitivity.groupby(["n", "regime"]):
            per_instance = []
            for _instance, values in group.groupby("instance_id"):
                solved_values = values[_as_bool(values["solved"])]
                per_instance.append({
                    "all_solved": len(solved_values) == len(values),
                    "nv_stable": solved_values["nv_star"].nunique(dropna=False) <= 1,
                    "ranking_stable": solved_values["ranking"].nunique(dropna=False) <= 1,
                })
            checks = pd.DataFrame(per_instance)
            rows.append({
                "n": n, "regime": regime, "instances": len(checks),
                "all_solved_pct": 100 * checks["all_solved"].mean(),
                "nv_stable_pct": 100 * checks["nv_stable"].mean(),
                "ranking_stable_pct": 100 * checks["ranking_stable"].mean(),
                "runtime_mean": group["runtime"].mean(),
            })
        pd.DataFrame(rows).to_csv(summary / "epsilon_sensitivity.csv", index=False)

    example_sensitivity_path = raw / "example_epsilon_sensitivity.csv"
    if example_sensitivity_path.exists():
        pd.read_csv(example_sensitivity_path).to_csv(
            summary / "example_epsilon_sensitivity.csv", index=False
        )

    sanity_path = raw / "representation_sanity.csv"
    if sanity_path.exists():
        sanity = pd.read_csv(sanity_path)
        sanity.groupby("n").agg(
            instances=("replicate", "count"), feasible_pct=("feasible", lambda x: 100 * x.mean()),
            runtime_mean=("runtime", "mean"), runtime_max=("runtime", "max"),
        ).reset_index().to_csv(summary / "representation_sanity.csv", index=False)

    print(f"summaries -> {summary}; figures -> {figures}", flush=True)
