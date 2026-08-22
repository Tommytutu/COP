"""Repair a cyclic PCM with EVRIM and print the revised matrix."""

import numpy as np

from _setup import default_settings, print_input_summary, ranking
from ahpcop import repair_with_evrim


# x1 > x2, x2 > x3, and x3 > x1 is an explicit preference cycle.
A = np.array([
    [1, 3, 1 / 3],
    [1 / 3, 1, 3],
    [3, 1 / 3, 1],
], dtype=float)


def main() -> None:
    print_input_summary(A)
    result = repair_with_evrim(A, settings=default_settings())
    print(f"\nStatus: {result.status}; certified={result.solved}")
    print(f"NRP*: {result.nrp}; AOC: {result.aoc}; GCI: {result.gci}; NV: {result.nv}")
    print(f"Runtime: {result.runtime:.6f} s")
    if result.revised_matrix is not None and result.weights is not None:
        print("Revised PCM:")
        print(np.array2string(result.revised_matrix, precision=5, suppress_small=True))
        print("Weights:", np.array2string(result.weights, precision=8, suppress_small=True))
        print("Ranking:", ranking(result.weights))


if __name__ == "__main__":
    main()
