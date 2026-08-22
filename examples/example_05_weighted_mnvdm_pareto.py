"""Run alpha-MNVDM on one user-supplied PCM and print its Pareto points."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ahpcop import GurobiSettings, solve_weighted_mnvdm  # noqa: E402


A = np.array([
    [1, 2, 4, 7],
    [1 / 2, 1, 3, 5],
    [1 / 4, 1 / 3, 1, 2],
    [1 / 7, 1 / 5, 1 / 2, 1],
], dtype=float)

settings = GurobiSettings(time_limit=60, threads=12, output_flag=0)
previous_y = None
previous_states = None
for alpha in (1.0, 0.75, 0.5, 0.25, 0.0):
    # The public helper handles a single alpha.  The batch experiment additionally
    # passes the preceding solution as a continuation start for all 101 values.
    result = solve_weighted_mnvdm(A, alpha, settings=settings)
    print(
        f"alpha={alpha:.2f} status={result.status} "
        f"NV={result.nv} NVR={result.nvr:.6f} "
        f"GCI-deviation={result.gci_deviation:.6f} runtime={result.runtime:.3f}s"
    )
    print("weights:", np.round(result.weights, 6))
    previous_y = result.y
    previous_states = result.relation_states
