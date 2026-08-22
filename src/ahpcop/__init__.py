"""AHPCOP: order-preserving pairwise preference models with Gurobi."""

__version__ = "1.0.0"
"""Order-preserving pairwise preference models and experiment utilities."""

from .api import (
    repair_with_evrim,
    solve_mnvdm,
    solve_priority_methods,
    solve_weighted_mnvdm,
    validate_pcm,
)
from .priority import GurobiSettings

__all__ = [
    "GurobiSettings",
    "repair_with_evrim",
    "solve_mnvdm",
    "solve_priority_methods",
    "solve_weighted_mnvdm",
    "validate_pcm",
]
