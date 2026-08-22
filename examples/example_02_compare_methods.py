"""Compare EM, LLSM, MNVEM, and MNVLLSM on one PCM."""

import numpy as np

from _setup import default_settings, print_input_summary, print_priority_result
from ahpcop import solve_priority_methods


A = np.array([
    [1, 2, 4, 9],
    [1 / 2, 1, 3, 7],
    [1 / 4, 1 / 3, 1, 5],
    [1 / 9, 1 / 7, 1 / 5, 1],
], dtype=float)


def main() -> None:
    print_input_summary(A)
    for result in solve_priority_methods(A, settings=default_settings()):
        print_priority_result(A, result)


if __name__ == "__main__":
    main()
