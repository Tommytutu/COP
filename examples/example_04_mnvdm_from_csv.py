"""Run MNVDM on a comma-separated reciprocal matrix file."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from _setup import default_settings, print_input_summary, print_priority_result
from cop_experiments import solve_mnvdm


def main() -> None:
    parser = argparse.ArgumentParser(description="Solve one CSV PCM with MNVDM")
    parser.add_argument("matrix", type=Path, help="CSV file with no header or row labels")
    parser.add_argument("--method", choices=["LLSM", "EM"], default="LLSM")
    args = parser.parse_args()
    matrix = np.loadtxt(args.matrix, delimiter=",")
    print_input_summary(matrix)
    result = solve_mnvdm(matrix, method=args.method, settings=default_settings())
    print_priority_result(matrix, result)


if __name__ == "__main__":
    main()
