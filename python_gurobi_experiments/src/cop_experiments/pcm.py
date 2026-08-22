from __future__ import annotations

import json
from dataclasses import dataclass
from itertools import combinations
from typing import Iterable

import numpy as np


SAATY_SCALE = np.array(
    [1 / 9, 1 / 8, 1 / 7, 1 / 6, 1 / 5, 1 / 4, 1 / 3, 1 / 2,
     1, 2, 3, 4, 5, 6, 7, 8, 9],
    dtype=float,
)
LOG_SCALE = np.log(SAATY_SCALE)
RI = {3: 0.58, 4: 0.90, 5: 1.12, 6: 1.24, 7: 1.32, 8: 1.41, 9: 1.45}


@dataclass(frozen=True)
class SyntheticPCM:
    instance_id: str
    n: int
    regime: str
    replicate: int
    seed: int
    matrix: np.ndarray
    latent_weights: np.ndarray
    cr: float
    gci: float
    realized_cr_bin: str
    cycle: tuple[int, int, int] | None = None

    def as_record(self) -> dict:
        return {
            "instance_id": self.instance_id,
            "n": self.n,
            "regime": self.regime,
            "replicate": self.replicate,
            "seed": self.seed,
            "cr": self.cr,
            "gci": self.gci,
            "realized_cr_bin": self.realized_cr_bin,
            "cycle": json.dumps(self.cycle) if self.cycle else "",
            "matrix": json.dumps(self.matrix.tolist()),
            "latent_weights": json.dumps(self.latent_weights.tolist()),
        }


def upper_pairs(n: int) -> list[tuple[int, int]]:
    return list(combinations(range(n), 2))


def softmax(y: np.ndarray) -> np.ndarray:
    z = np.asarray(y, dtype=float)
    z = z - np.max(z)
    e = np.exp(z)
    return e / e.sum()


def llsm_weights(a: np.ndarray) -> np.ndarray:
    y = np.log(a).mean(axis=1)
    return softmax(y)


def em_weights(a: np.ndarray) -> np.ndarray:
    values, vectors = np.linalg.eig(a)
    idx = int(np.argmax(values.real))
    w = np.abs(vectors[:, idx].real)
    return w / w.sum()


def consistency_ratio(a: np.ndarray) -> float:
    n = a.shape[0]
    lam = float(np.max(np.linalg.eigvals(a).real))
    ci = max(0.0, (lam - n) / (n - 1))
    return ci / RI[n]


def geometric_consistency_index(a: np.ndarray, w: np.ndarray | None = None) -> float:
    n = a.shape[0]
    if w is None:
        w = llsm_weights(a)
    y = np.log(w)
    residuals = [np.log(a[i, j]) - y[i] + y[j] for i, j in upper_pairs(n)]
    return 2.0 * float(np.dot(residuals, residuals)) / ((n - 1) * (n - 2))


def cr_bin(cr: float) -> str:
    if cr <= 0.05 + 1e-12:
        return "CR<=0.05"
    if cr <= 0.10 + 1e-12:
        return "0.05<CR<=0.10"
    if cr <= 0.20 + 1e-12:
        return "0.10<CR<=0.20"
    return "CR>0.20"


def nearest_scale_log(x: float) -> float:
    return float(LOG_SCALE[int(np.argmin(np.abs(LOG_SCALE - x)))])


def matrix_from_upper_logs(n: int, logs: dict[tuple[int, int], float]) -> np.ndarray:
    a = np.ones((n, n), dtype=float)
    for (i, j), value in logs.items():
        a[i, j] = np.exp(value)
        a[j, i] = np.exp(-value)
    return a


def _set_directed(a: np.ndarray, i: int, j: int, value: float) -> None:
    a[i, j] = value
    a[j, i] = 1.0 / value


def _latent_vector(rng: np.random.Generator, n: int) -> tuple[np.ndarray, np.ndarray]:
    z = rng.normal(size=n)
    z -= z.mean()
    span = float(np.ptp(z))
    if span < 1e-8:
        z = np.linspace(-1.0, 1.0, n)
        span = 2.0
    target_span = rng.uniform(np.log(3.0), np.log(7.0))
    z *= target_span / span
    return z, softmax(z)


def _noisy_pcm(rng: np.random.Generator, z: np.ndarray, sigma: float) -> np.ndarray:
    n = z.size
    logs: dict[tuple[int, int], float] = {}
    for i, j in upper_pairs(n):
        noisy = z[i] - z[j] + rng.normal(0.0, sigma)
        noisy = float(np.clip(noisy, -np.log(9.0), np.log(9.0)))
        logs[i, j] = nearest_scale_log(noisy)
    return matrix_from_upper_logs(n, logs)


def _target_interval(regime: str) -> tuple[float, float]:
    return {
        "low": (0.0, 0.05),
        "moderate": (0.05, 0.10),
        "high": (0.10, 0.20),
    }[regime]


def _sigma_range(regime: str) -> tuple[float, float]:
    return {
        "low": (0.02, 0.22),
        "moderate": (0.22, 0.75),
        "high": (0.45, 1.25),
    }[regime]


def generate_pcm(
    n: int,
    regime: str,
    replicate: int,
    seed: int,
    max_attempts: int = 20_000,
) -> SyntheticPCM:
    """Generate a Saaty-scale reciprocal PCM and retain its latent weights.

    Low/moderate/high samples are rejection-classified by realized CR. Cyclic
    samples contain an explicit three-alternative directed cycle and are reported
    by their realized CR bin rather than forced into a CR interval.
    """
    stream_seed = int(seed + 10_000 * n + 1_000_000 * replicate + 10_000_000 *
                      {"low": 1, "moderate": 2, "high": 3, "cyclic": 4}[regime])
    rng = np.random.default_rng(stream_seed)
    cycle: tuple[int, int, int] | None = None

    if regime == "cyclic":
        z, w0 = _latent_vector(rng, n)
        a = _noisy_pcm(rng, z, rng.uniform(0.25, 0.65))
        ordered = np.argsort(-z)
        i, j, k = map(int, ordered[:3])
        strength = float(rng.choice(np.arange(3, 8)))
        _set_directed(a, i, j, strength)
        _set_directed(a, j, k, strength)
        _set_directed(a, k, i, strength)
        cycle = (i, j, k)
    else:
        low, high = _target_interval(regime)
        sig_lo, sig_hi = _sigma_range(regime)
        for _ in range(max_attempts):
            z, w0 = _latent_vector(rng, n)
            a = _noisy_pcm(rng, z, rng.uniform(sig_lo, sig_hi))
            value = consistency_ratio(a)
            lower_ok = value >= low - 1e-12 if low == 0 else value > low + 1e-12
            if lower_ok and value <= high + 1e-12:
                break
        else:
            raise RuntimeError(f"Could not generate {regime} n={n} after {max_attempts} attempts")

    value_cr = consistency_ratio(a)
    value_gci = geometric_consistency_index(a)
    return SyntheticPCM(
        instance_id=f"{regime}_n{n}_r{replicate:03d}",
        n=n,
        regime=regime,
        replicate=replicate,
        seed=stream_seed,
        matrix=a,
        latent_weights=w0,
        cr=value_cr,
        gci=value_gci,
        realized_cr_bin=cr_bin(value_cr),
        cycle=cycle,
    )


def load_record_matrix(value: str) -> np.ndarray:
    return np.asarray(json.loads(value), dtype=float)


def load_record_weights(value: str) -> np.ndarray:
    return np.asarray(json.loads(value), dtype=float)


def contains_directed_cycle(a: np.ndarray, triple: Iterable[int]) -> bool:
    i, j, k = triple
    return bool(a[i, j] > 1 and a[j, k] > 1 and a[k, i] > 1)
