"""Imports and compact result printers shared by the examples."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC = PROJECT_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from cop_experiments.metrics import normalized_violation_rate, violation_score  # noqa: E402
from cop_experiments.pcm import consistency_ratio, geometric_consistency_index  # noqa: E402
from cop_experiments.priority import PriorityResult  # noqa: E402


def default_settings():
    from cop_experiments import GurobiSettings

    return GurobiSettings(
        epsilon=1e-4,
        y_bound=10.0,
        time_limit=60,
        mip_gap=1e-5,
        threads=1,
        seed=20260815,
        output_flag=0,
    )


def ranking(weights: np.ndarray) -> str:
    order = np.argsort(-weights) + 1
    return " > ".join(f"x{index}" for index in order)


def print_input_summary(matrix: np.ndarray) -> None:
    print("Input PCM:")
    print(np.array2string(matrix, precision=5, suppress_small=True))
    print(f"CR  = {consistency_ratio(matrix):.6f}")
    print(f"GCI = {geometric_consistency_index(matrix):.6f}")


def print_priority_result(matrix: np.ndarray, result: PriorityResult) -> None:
    print(f"\nMethod: {result.method}")
    print(f"Status: {result.status}; certified={result.solved}")
    print(f"Runtime: {result.runtime:.6f} s")
    if result.weights is None:
        print("No certified priority vector was returned.")
        return
    weights = result.weights
    print("Weights:", np.array2string(weights, precision=8, suppress_small=True))
    print("Ranking:", ranking(weights))
    print(f"NV from returned weights: {violation_score(matrix, weights):.6g}")
    print(f"NVR: {normalized_violation_rate(matrix, weights):.6g}")
    if result.nv_star is not None:
        print(f"Certified Stage-1 NV*: {result.nv_star:.6g}")
    if result.objective is not None:
        print(f"Stage-2 objective: {result.objective:.9g}")
    print(f"Stage-1 runtime: {result.stage1_runtime:.6f} s")
    print(f"Stage-2 runtime: {result.stage2_runtime:.6f} s")
    if result.gap is not None:
        print(f"Final MIP gap: {result.gap:.3e}")
