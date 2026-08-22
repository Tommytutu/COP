# Single-matrix examples

Run these scripts from the repository root after installing the package and
activating a valid Gurobi license.

| Example | Purpose | Command |
|---|---|---|
| `example_01_mnvdm_matrix.py` | Paste a reciprocal PCM into Python and solve MNVLLSM | `python examples/example_01_mnvdm_matrix.py` |
| `example_02_compare_methods.py` | Compare EM, LLSM, MNVEM, and MNVLLSM | `python examples/example_02_compare_methods.py` |
| `example_03_evrim_repair.py` | Repair a cyclic PCM with EVRIM | `python examples/example_03_evrim_repair.py` |
| `example_04_mnvdm_from_csv.py` | Read a PCM from a CSV file | `python examples/example_04_mnvdm_from_csv.py examples/sample_matrix.csv` |
| `example_05_weighted_mnvdm_pareto.py` | Trace selected weighted-MNVDM trade-off points | `python examples/example_05_weighted_mnvdm_pareto.py` |

## Example 1: one PCM with MNVLLSM

Edit the NumPy array `A` in `example_01_mnvdm_matrix.py`. The script performs
the following operations:

1. checks that the matrix is square, positive, reciprocal, and has a unit diagonal;
2. reports CR and GCI for the input matrix;
3. solves MNVDM Stage 1 to minimize order violations;
4. fixes the certified Stage-1 optimum and solves the LLSM Stage 2;
5. reports the priority vector, ranking, NV, NVR, runtime, and solver status.

To use a file instead, place a numeric reciprocal matrix without headers in a
CSV file and run Example 4.

## Example 3: EVRIM revision

`example_03_evrim_repair.py` starts from a three-alternative preference cycle,
solves the two EVRIM objectives, and prints the revised PCM, NRP, AOC, GCI, NV,
priority vector, and ranking. Protected judgments can be supplied through the
public `repair_with_evrim` API as zero-based `(i, j)` index pairs.
