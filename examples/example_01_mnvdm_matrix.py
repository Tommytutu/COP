"""Paste one PCM below and solve it with MNVDM (MNVLLSM by default)."""

import numpy as np

from _setup import default_settings, print_input_summary, print_priority_result
from ahpcop import solve_mnvdm


# Replace this array with your own positive reciprocal matrix.
# This is Example 3 in the current manuscript.
A = np.array([
    [1, 2, 4, 9],
    [1 / 2, 1, 3, 7],
    [1 / 4, 1 / 3, 1, 5],
    [1 / 9, 1 / 7, 1 / 5, 1],
], dtype=float)



print_input_summary(A)
result = solve_mnvdm(A, method="LLSM", settings=default_settings())
print_priority_result(A, result)
