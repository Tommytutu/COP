from __future__ import annotations

from pathlib import Path
import json

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results_weighted_pareto_all_n_20260821"
SOURCE = RESULTS / "summary" / "alpha_mnvdm_decile_tables_wide.csv"
SUBMISSION = ROOT.parent / "latex_project"
REGIME_ORDER = ["low", "moderate", "high", "cyclic"]


def fmt(value: object) -> str:
    return f"{float(value):.4f}"


def pair_cells(row: pd.Series, keys: list[str]) -> list[str]:
    cells: list[str] = []
    for key in keys:
        cells.extend([fmt(row[f"{key}_nvr"]), fmt(row[f"{key}_gci"])])
    return cells


def panel(frame: pd.DataFrame, labels: list[str], keys: list[str]) -> str:
    if len(keys) != 6 or len(labels) != 6:
        raise RuntimeError("each panel must contain exactly six methods")
    header_top = ["Regime"] + [f"\\multicolumn{{2}}{{c}}{{{label}}}" for label in labels]
    header_bottom = [""] + [item for _ in keys for item in ("NVR", "GCI")]
    lines = [
        "\\setlength{\\tabcolsep}{2.3pt}%",
        "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}l*{12}{r}@{}}",
        "\\toprule",
        " & ".join(header_top) + r"\\",
        " & ".join(header_bottom) + r"\\",
        "\\midrule",
    ]
    indexed = frame.set_index("regime")
    for regime in REGIME_ORDER:
        row = indexed.loc[regime]
        cells = [regime.capitalize(), *pair_cells(row, keys)]
        lines.append(" & ".join(cells) + r"\\")
    lines.extend(["\\bottomrule", "\\end{tabular*}"])
    return "\n".join(lines)


def table_for_n(frame: pd.DataFrame, n: int, placement: str, label: str) -> str:
    current = frame[frame["n"] == n].copy()
    if set(current["regime"]) != set(REGIME_ORDER):
        raise RuntimeError(f"n={n} does not contain all four regimes")
    instances = sorted(current["instances"].astype(int).unique())
    if len(instances) != 1:
        raise RuntimeError(f"n={n} has inconsistent sample counts: {instances}")
    sample_count = instances[0]
    caption = (
        f"Mean $n={n}$ weighted-MNVDM results and conventional baselines over "
        f"{sample_count} PCMs per regime."
    )
    panel_a = panel(
        current,
        ["$\\alpha=1$", "$\\alpha=0.9$", "$\\alpha=0.8$", "$\\alpha=0.7$", "$\\alpha=0.6$", "$\\alpha=0.5$"],
        ["alpha_1p0", "alpha_0p9", "alpha_0p8", "alpha_0p7", "alpha_0p6", "alpha_0p5"],
    )
    panel_b = panel(
        current,
        ["$\\alpha=0.4$", "$\\alpha=0.3$", "$\\alpha=0.2$", "$\\alpha=0.1$", "LLSM", "EM"],
        ["alpha_0p4", "alpha_0p3", "alpha_0p2", "alpha_0p1", "alpha_0p0", "em"],
    )
    return "\n".join([
        f"\\begin{{table}}[{placement}]",
        "\\centering",
        "\\scriptsize",
        f"\\caption{{{caption}}}",
        f"\\label{{{label}}}",
        "\\textbf{Panel A: weighted MNVDM, $\\alpha=1$ to $0.5$}\\par\\smallskip",
        panel_a,
        "\\medskip",
        "\\textbf{Panel B: weighted MNVDM, $\\alpha=0.4$ to $0.1$, and conventional baselines}\\par\\smallskip",
        panel_b,
        "\\end{table}",
    ])


def main() -> None:
    frame = pd.read_csv(SOURCE)
    expected_sizes = set(range(3, 10))
    if set(frame["n"].astype(int)) != expected_sizes:
        raise RuntimeError("wide summary does not cover n=3,...,9")
    SUBMISSION.mkdir(parents=True, exist_ok=True)
    (SUBMISSION / "table_alpha_tradeoff_n7.tex").write_text(
        table_for_n(frame, 7, "t", "tab:alpha-tradeoff"), encoding="utf-8"
    )
    other_tables = [
        table_for_n(frame, n, "p", f"tab:alpha-tradeoff-n{n}")
        for n in (3, 4, 5, 6, 8, 9)
    ]
    (SUBMISSION / "table_alpha_tradeoff_other_n.tex").write_text(
        "\n\n".join(other_tables) + "\n", encoding="utf-8"
    )
    reported = pd.read_csv(RESULTS / "summary" / "alpha_mnvdm_reported_decile_points.csv")
    coverage = reported.groupby("instance_id")["alpha"].nunique()
    if len(coverage) != 424 or coverage.min() != 11 or coverage.max() != 11:
        raise RuntimeError("reported decile grid is incomplete")
    certification = reported.groupby(["n", "regime"], as_index=False).agg(
        runs=("alpha", "size"),
        instances=("instance_id", "nunique"),
        certified=("certified", "sum"),
        feasible=("python_feasible", "sum"),
        timeouts=("status", lambda values: (values == "TIME_LIMIT").sum()),
        mean_runtime=("runtime", "mean"),
        median_runtime=("runtime", "median"),
        p90_runtime=("runtime", lambda values: values.quantile(0.90)),
        p95_runtime=("runtime", lambda values: values.quantile(0.95)),
        max_runtime=("runtime", "max"),
    )
    certification["certified_rate"] = certification["certified"] / certification["runs"]
    certification.to_csv(
        RESULTS / "summary" / "alpha_mnvdm_decile_certification_by_n_regime.csv", index=False
    )
    monotonicity_violations: list[dict] = []
    # At alpha=0.9, a 1e-5 absolute objective tolerance permits up to 1e-4
    # movement in the GCI component because its coefficient is 0.1.
    tolerance = 1e-4
    for instance_id, group in reported.groupby("instance_id"):
        certified = group[
            group["certified"].astype(bool) & group["python_feasible"].astype(bool)
        ].sort_values("alpha", ascending=False)
        previous = None
        for _, point in certified.iterrows():
            if previous is not None:
                if float(point["nvr"]) < float(previous["nvr"]) - tolerance:
                    monotonicity_violations.append({
                        "instance_id": instance_id,
                        "criterion": "NVR",
                        "alpha_high": float(previous["alpha"]),
                        "alpha_low": float(point["alpha"]),
                    })
                if float(point["gci_deviation"]) > float(previous["gci_deviation"]) + tolerance:
                    monotonicity_violations.append({
                        "instance_id": instance_id,
                        "criterion": "GCI",
                        "alpha_high": float(previous["alpha"]),
                        "alpha_low": float(point["alpha"]),
                    })
            previous = point
    pd.DataFrame(monotonicity_violations).to_csv(
        RESULTS / "summary" / "alpha_mnvdm_decile_monotonicity_audit.csv", index=False
    )
    audit = {
        "instances": int(len(coverage)),
        "reported_rows": int(len(reported)),
        "alpha_values_per_instance": int(coverage.min()),
        "python_feasible": int(reported["python_feasible"].astype(bool).sum()),
        "certified": int(reported["certified"].astype(bool).sum()),
        "timeouts": int((reported["status"] == "TIME_LIMIT").sum()),
        "component_monotonicity_tolerance": tolerance,
        "certified_monotonicity_violations": len(monotonicity_violations),
    }
    (RESULTS / "summary" / "alpha_mnvdm_decile_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    print(SUBMISSION / "table_alpha_tradeoff_n7.tex")
    print(SUBMISSION / "table_alpha_tradeoff_other_n.tex")


if __name__ == "__main__":
    main()
