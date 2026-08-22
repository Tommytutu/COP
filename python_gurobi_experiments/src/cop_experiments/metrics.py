from __future__ import annotations

from itertools import combinations

import numpy as np
from scipy.stats import kendalltau

from .pcm import upper_pairs


def order_relation_count(n: int) -> int:
    m = n * (n - 1) // 2
    # The additional item is the neutral comparison a_rr=1. Relations against
    # it encode POP, while relations among the m judgments encode POIP.
    return m * (m + 1) // 2


def _sign(value: float, tol: float) -> int:
    if value > tol:
        return 1
    if value < -tol:
        return -1
    return 0


def violation_score(a: np.ndarray, w: np.ndarray, tol: float = 1e-7) -> float:
    """Compute the manuscript's half-weighted NV on unique order relations."""
    pairs = upper_pairs(a.shape[0])
    items = pairs + [(0, 0)]
    xa = np.array([np.log(a[i, j]) for i, j in items])
    y = np.log(np.asarray(w, dtype=float))
    dw = np.array([y[i] - y[j] for i, j in items])
    nv = 0.0
    for p, q in combinations(range(len(items)), 2):
        sa = _sign(float(xa[p] - xa[q]), tol)
        sw = _sign(float(dw[p] - dw[q]), tol)
        if sa == 0:
            nv += 0.0 if sw == 0 else 0.5
        elif sw == 0:
            nv += 0.5
        elif sa != sw:
            nv += 1.0
    return nv


def normalized_violation_rate(a: np.ndarray, w: np.ndarray, tol: float = 1e-7) -> float:
    count = order_relation_count(a.shape[0])
    return violation_score(a, w, tol=tol) / count if count else 0.0


def kendall_recovery(w_hat: np.ndarray, w0: np.ndarray) -> float:
    tau = kendalltau(np.asarray(w_hat), np.asarray(w0), variant="b").statistic
    return float(0.0 if np.isnan(tau) else tau)


def best_choice_accuracy(w_hat: np.ndarray, w0: np.ndarray) -> int:
    return int(int(np.argmax(w_hat)) == int(np.argmax(w0)))


def log_ratio_rmse(w_hat: np.ndarray, w0: np.ndarray) -> float:
    yh = np.log(np.asarray(w_hat, dtype=float))
    y0 = np.log(np.asarray(w0, dtype=float))
    errors = [(yh[i] - yh[j]) - (y0[i] - y0[j]) for i, j in upper_pairs(len(y0))]
    return float(np.sqrt(np.mean(np.square(errors))))


def recovery_metrics(a: np.ndarray, w_hat: np.ndarray, w0: np.ndarray) -> dict[str, float | int]:
    nv = violation_score(a, w_hat)
    return {
        "nv": nv,
        "n_order_relations": order_relation_count(a.shape[0]),
        "nvr": normalized_violation_rate(a, w_hat),
        "kendall_tau_b": kendall_recovery(w_hat, w0),
        "best_choice_accuracy": best_choice_accuracy(w_hat, w0),
        "lrmse": log_ratio_rmse(w_hat, w0),
    }
