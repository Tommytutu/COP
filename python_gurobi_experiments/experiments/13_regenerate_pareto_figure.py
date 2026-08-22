from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results_weighted_pareto_all_n_20260821"
SOURCE = RESULTS / "summary" / "alpha_mnvdm_all_points.csv"
FIGURES = RESULTS / "figures"
LATEX = ROOT.parent / "latex_project"


def main() -> None:
    frame = pd.read_csv(SOURCE)
    # Figure 2 reports only the three computationally difficult sizes.
    frame = frame[frame["n"].isin([7, 8, 9])].copy()
    representatives = set(
        frame.sort_values(["n", "regime", "replicate"])
        .groupby(["n", "regime"], as_index=False)
        .first()["instance_id"]
    )

    figure, axis = plt.subplots(figsize=(8.2, 5.8))
    selected = frame[frame["instance_id"].isin(representatives)]
    for _instance_id, group in selected.groupby("instance_id"):
        valid = group[group["python_feasible"].astype(bool)].sort_values("nvr")
        if valid.empty:
            continue
        first = valid.iloc[0]
        label = f"{str(first['regime']).capitalize()}, n={int(first['n'])}"
        line = axis.plot(
            valid["nvr"],
            valid["gci_deviation"],
            marker="o",
            markersize=2.5,
            linewidth=1.0,
            label=label,
        )[0]
        uncertified = valid[~valid["certified"].astype(bool)]
        if not uncertified.empty:
            axis.scatter(
                uncertified["nvr"],
                uncertified["gci_deviation"],
                marker="x",
                s=30,
                linewidths=1.2,
                color=line.get_color(),
                zorder=4,
            )

    axis.set_xlabel("Normalized violation ratio (NVR)")
    axis.set_ylabel("GCI")
    axis.grid(alpha=0.25)
    axis.legend(
        fontsize=6.5,
        ncol=2,
        loc="upper left",
        frameon=True,
        framealpha=0.9,
    )
    figure.tight_layout()
    FIGURES.mkdir(exist_ok=True)
    pdf_path = FIGURES / "alpha_mnvdm_pareto_selected.pdf"
    png_path = FIGURES / "alpha_mnvdm_pareto_selected.png"
    figure.savefig(pdf_path)
    figure.savefig(png_path, dpi=220)
    LATEX.mkdir(exist_ok=True)
    figure.savefig(LATEX / "alpha_mnvdm_pareto_selected.pdf")
    plt.close(figure)
    print(pdf_path)
    print(LATEX / "alpha_mnvdm_pareto_selected.pdf")


if __name__ == "__main__":
    main()
