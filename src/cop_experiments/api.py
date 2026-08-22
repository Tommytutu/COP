"""Small public API for applying the manuscript models to one PCM.

The experiment pipeline is intentionally batch-oriented.  This module exposes
the same verified core solvers through a compact interface for readers who only
want to paste in one reciprocal pairwise-comparison matrix.
"""

from __future__ import annotations

import numpy as np

from .evrim import EVRIMResult, solve_evrim, solve_evrim_check_first
from .priority import (
    AlphaMNVDMResult,
    GurobiSettings,
    PriorityResult,
    classical_priorities,
    solve_alpha_mnvdm,
    solve_mnvem,
    solve_mnvllsm,
    solve_stage1,
)


def validate_pcm(matrix: np.ndarray, tolerance: float = 1e-8) -> np.ndarray:
    """Return a float PCM after checking positivity, diagonal, and reciprocity."""
    a = np.asarray(matrix, dtype=float)
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError("A PCM must be a square two-dimensional matrix.")
    if a.shape[0] < 3:
        raise ValueError("The manuscript models require at least three alternatives.")
    if not np.all(np.isfinite(a)) or np.any(a <= 0):
        raise ValueError("Every PCM entry must be finite and strictly positive.")
    if not np.allclose(np.diag(a), 1.0, atol=tolerance, rtol=0.0):
        raise ValueError("Every diagonal entry must equal 1.")
    if not np.allclose(a * a.T, 1.0, atol=tolerance, rtol=tolerance):
        maximum_error = float(np.max(np.abs(a * a.T - 1.0)))
        raise ValueError(
            "The matrix must be reciprocal: a[i,j] * a[j,i] = 1 "
            f"(maximum error {maximum_error:.3e})."
        )
    return a


def solve_mnvdm(
    matrix: np.ndarray,
    method: str = "LLSM",
    settings: GurobiSettings | None = None,
    formulation: str = "indicator",
) -> PriorityResult:
    """Solve lexicographic MNVDM for one PCM.

    Parameters
    ----------
    matrix:
        Positive reciprocal PCM.
    method:
        ``"LLSM"`` for MNVLLSM or ``"EM"`` for MNVEM.
    settings:
        Gurobi limits and numerical tolerances.  Defaults are suitable for a
        small illustrative matrix.
    formulation:
        Stage-1 formulation: ``"indicator"`` (recommended exact solver),
        ``"basic"``, or the strengthened ``"strong"`` big-M formulation.
    """
    a = validate_pcm(matrix)
    chosen = method.strip().upper()
    if chosen not in {"LLSM", "EM"}:
        raise ValueError("method must be 'LLSM' or 'EM'.")
    if formulation not in {"indicator", "strong", "basic"}:
        raise ValueError("formulation must be 'indicator', 'strong', or 'basic'.")
    actual_settings = settings or GurobiSettings()
    stage1 = solve_stage1(a, actual_settings, formulation)
    return (
        solve_mnvllsm(a, stage1, actual_settings)
        if chosen == "LLSM"
        else solve_mnvem(a, stage1, actual_settings)
    )


def solve_priority_methods(
    matrix: np.ndarray,
    settings: GurobiSettings | None = None,
) -> list[PriorityResult]:
    """Run EM, LLSM, MNVEM, and MNVLLSM on the same PCM."""
    a = validate_pcm(matrix)
    actual_settings = settings or GurobiSettings()
    stage1 = solve_stage1(a, actual_settings, "indicator")
    classical = classical_priorities(a)
    return classical + [
        solve_mnvem(a, stage1, actual_settings),
        solve_mnvllsm(a, stage1, actual_settings),
    ]


def solve_weighted_mnvdm(
    matrix: np.ndarray,
    alpha: float,
    settings: GurobiSettings | None = None,
    gci_normalizer: float = 1.0,
) -> AlphaMNVDMResult:
    """Solve the single-objective alpha-MNVDM model for one PCM."""
    a = validate_pcm(matrix)
    return solve_alpha_mnvdm(
        a,
        alpha=float(alpha),
        settings=settings or GurobiSettings(),
        gci_normalizer=float(gci_normalizer),
    )


def repair_with_evrim(
    matrix: np.ndarray,
    settings: GurobiSettings | None = None,
    threshold: float | None = None,
    value_protected: list[tuple[int, int]] | None = None,
    direction_protected: list[tuple[int, int]] | None = None,
    backend: str = "direct",
    check_first: bool = False,
) -> EVRIMResult:
    """Run EVRIM for one PCM; protected indices are zero-based Python pairs."""
    a = validate_pcm(matrix)
    actual_settings = settings or GurobiSettings()
    if check_first:
        return solve_evrim_check_first(
            a,
            actual_settings,
            threshold=threshold,
            value_protected=value_protected,
            direction_protected=direction_protected,
        )
    return solve_evrim(
        a, actual_settings, threshold=threshold,
        value_protected=value_protected,
        direction_protected=direction_protected, backend=backend,
    )
