# COP: order-preserving pairwise preference models

This repository contains the current Python/Gurobi implementation of MNVDM,
weighted MNVDM, and EVRIM, together with the experiment drivers and executable
single-matrix examples used for the Omega submission.

The original MATLAB/YALMIP files remain in the repository for archival use.
Their original documentation is available in
[`MATLAB_LEGACY_README.md`](MATLAB_LEGACY_README.md). New applications should
use the Python package in [`python_gurobi_experiments/`](python_gurobi_experiments/).

## Requirements

- Python 3.11 or later
- Gurobi Optimizer 12.0.0
- a valid Gurobi license

Install and test the project:

```powershell
cd python_gurobi_experiments
python -m pip install -e .
python -m pytest -q
```

The default numerical settings used by the public examples are a 60-second
time limit, MIP gap `1e-5`, and strict-order implementation tolerance `1e-4`.

## Concrete example: apply MNVLLSM to one PCM

The matrix below is Example 3 in the current manuscript:

```python
import numpy as np

from cop_experiments import GurobiSettings, solve_mnvdm

A = np.array([
    [1,     2,     4,     9],
    [1/2,   1,     3,     7],
    [1/4,   1/3,   1,     5],
    [1/9,   1/7,   1/5,   1],
], dtype=float)

settings = GurobiSettings(
    time_limit=60,
    mip_gap=1e-5,
    epsilon=1e-4,
    threads=1,
    output_flag=0,
)

result = solve_mnvdm(A, method="LLSM", settings=settings)

print("status:", result.status)
print("certified:", result.solved)
print("minimum NV:", result.nv_star)
print("weights:", result.weights)
print("runtime:", result.runtime)
```

The complete executable version is
[`examples/example_01_mnvdm_matrix.py`](python_gurobi_experiments/examples/example_01_mnvdm_matrix.py).
Its process is:

1. validate positivity, reciprocity, and the unit diagonal;
2. report the input CR and GCI;
3. solve MNVDM Stage 1 and certify the minimum order-violation value;
4. fix that value and solve the LLSM Stage 2;
5. report the priority vector, ranking, NV, NVR, objective values, runtime, and status.

The verified output on the current test machine is:

| Output | Value |
|---|---|
| Solver status | `OPTIMAL`, certified |
| Minimum NV | `0` |
| NVR | `0` |
| Priority vector | `(0.50197724, 0.31379938, 0.14331878, 0.04090460)` |
| Ranking | `x1 > x2 > x3 > x4` |
| Stage-2 LLSM objective | `0.132568674` |

Runtime depends on the processor, Gurobi version, license, and thread setting.

Run it with:

```powershell
cd python_gurobi_experiments
python examples/example_01_mnvdm_matrix.py
```

## Other examples

```powershell
# EM, LLSM, MNVEM, and MNVLLSM on the same PCM
python examples/example_02_compare_methods.py

# EVRIM revision of a cyclic PCM
python examples/example_03_evrim_repair.py

# Read a user PCM from CSV
python examples/example_04_mnvdm_from_csv.py examples/sample_matrix.csv

# Selected weighted-MNVDM trade-off points
python examples/example_05_weighted_mnvdm_pareto.py
```

For the cyclic PCM in Example 3, the verified EVRIM run changes one independent
judgment (`NRP = 1`), returns `AOC = 2.48490665`, `GCI = 0.21920261`, and
produces a revised PCM with `NV = 0`. The script prints the complete revised
matrix and its derived priority vector.

See [`examples/README.md`](python_gurobi_experiments/examples/README.md) for the
step-by-step purpose of every example and
[`README_CN.md`](python_gurobi_experiments/README_CN.md) for Chinese instructions
and the experiment-to-file mapping.

## Project structure

```text
python_gurobi_experiments/
  src/cop_experiments/   Public API and optimization models
  examples/              Single-matrix executable examples
  experiments/           One entry point for each manuscript experiment
  tests/                 Unit and solver integration tests
  run_all.py             Batch experiment dispatcher
  config_*.json          Reproducible solver configurations
```

Generated Gurobi logs and large result directories are intentionally excluded
from Git. Experiment scripts checkpoint raw CSV output locally and can resume
completed runs.

## Public Python API

- `validate_pcm`: validate a positive reciprocal PCM;
- `solve_mnvdm`: solve MNVLLSM or MNVEM;
- `solve_priority_methods`: compare EM, LLSM, MNVEM, and MNVLLSM;
- `solve_weighted_mnvdm`: solve one weighted MNVDM value of `alpha`;
- `repair_with_evrim`: revise a PCM with optional value or direction protection.

All protected judgment indices passed through Python are zero-based.
