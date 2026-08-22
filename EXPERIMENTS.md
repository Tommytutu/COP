# AHPCOP Python/Gurobi experiment project

This project reproduces the revised experiments for the manuscript. It is a new
Python implementation informed by the public MATLAB/YALMIP code, but it replaces
the archived weighted objectives and fixed big-M constants with the revised
two-stage lexicographic formulations.

## Scope

- Matrix sizes are fixed at `n = 3,...,9`; no larger complete PCMs are used.
- Synthetic regimes: `low`, `moderate`, `high`, and explicit `cyclic` data.
- Low/moderate/high instances are retained by realized CR:
  - low: `CR <= 0.05`;
  - moderate: `0.05 < CR <= 0.10`;
  - high: `0.10 < CR <= 0.20`.
- Cyclic instances contain an explicit directed three-alternative preference cycle
  and are additionally classified by their realized CR bin.
- Prioritization methods: EM, LLSM, MNVEM, MNVLLSM, and the feasibility-only
  COP-LLSM sanity baseline.
- Recovery outcomes: NVR, Kendall tau-b, best-choice accuracy, and log-ratio RMSE.
- EVRIM experiments: high/cyclic PCMs, unprotected and protected variants, plus
  three actual house-buying cases.
- Runtime summaries report mean, median, standard deviation, P90, P95, maximum,
  solved percentage, and timeout percentage.

## Solver formulation

Stage 1 minimizes twice the half-weighted NV, so its objective is integer-valued.
The comparison domain contains the `m=n(n-1)/2` upper-triangular judgments and
one neutral item `a_rr=1`; its `m` neutral relations encode POP and its remaining
relations encode POIP. Thus `N_order_relations = choose(m+1,2)`.
Stage 2 fixes the certified Stage-1 optimum and minimizes the LLSM or EM criterion.
The primary solver uses explicit one-hot trichotomy and native Gurobi indicators.
The basic comparator uses a two-binary bounded big-M encoding with constants
derived for each relation from the declared `[-10,10]` log-weight domain. The
strengthened comparator uses exact coefficient aggregation, an explicit
three-state partition, native indicator constraints, and valid total-preorder
inequalities. It is integer-equivalent to the basic formulation; the paper does
not claim projection dominance between their continuous relaxations. EVRIM
similarly solves NRP first and AOC second, with a convex GCI constraint and
discrete Saaty-scale selection.

EVRIM keeps the POIP certificate `y` separate from the LLSM vector `u` used to
evaluate GCI. The `direct` backend adds the convex GCI quadratic constraint to
Gurobi. The check-first implementation used for Table 7 first certifies the
two-stage optimum without GCI and checks the resulting GCI. Only a failing
certified optimum starts a fresh outer-approximation solve; no MIP start is
passed. The fallback adds tangent cuts for GCI-infeasible solutions. A result is
reported as certified only when its MIP gap and independently recomputed GCI
meet the configured tolerances.

## Environment used for the official run

- Python: `D:\anaconda3\python.exe` (Python 3.12)
- Gurobi: 12.0.0, academic license
- Platform: Windows 11 (`10.0.26200`)
- CPU: Intel64 Family 6 Model 183, 32 logical processors
- RAM: 68,442,210,304 bytes (approximately 63.7 GiB)
- Threads: 12
- Master seed: 20260815
- MIP gap: 1e-6 for the archived lexicographic runs and 1e-5 for the
  weighted-path and Table 7 reruns
- Numerical strict-order tolerance: 1e-4
- Per-lexicographic-stage limit: 60 s for the model-aligned certificate runs

The exact machine and package details are written automatically to
`results_model_aligned_20260819/environment.json`.

## Commands

From this directory:

```powershell
D:\anaconda3\python.exe -m pip install -e .
D:\anaconda3\python.exe -m pytest -q
D:\anaconda3\python.exe run_all.py smoke
D:\anaconda3\python.exe run_all.py formulation --config config_model_aligned_20260819.json
D:\anaconda3\python.exe run_all.py evrim --config config_model_aligned_20260819.json
D:\anaconda3\python.exe run_all.py house --config config_model_aligned_20260819.json
D:\anaconda3\python.exe run_all.py summarize --config config_model_aligned_20260819.json
```

Individual phases can be resumed:

```powershell
D:\anaconda3\python.exe run_all.py generate
D:\anaconda3\python.exe run_all.py priority
D:\anaconda3\python.exe run_all.py formulation
D:\anaconda3\python.exe run_all.py evrim
D:\anaconda3\python.exe run_all.py bnc
D:\anaconda3\python.exe run_all.py house
D:\anaconda3\python.exe run_all.py summarize
```

Each manuscript experiment also has an independent entry point under
`experiments/`.  Ready-to-edit single-matrix examples are under `examples/`.
See `README_CN.md` for the experiment-to-file mapping and Chinese instructions.

The final additions are:

- `10_alpha_mnvdm_pareto.py`: weighted MNVDM paths;
- `11_generate_alpha_tables.py`: Table 4 and the size-specific supplement tables;
- `12_table7_evrim_check_first.py`: the two rerun OA variants in Table 7;
- `13_regenerate_pareto_figure.py`: Figure 2 for representative $n=7,8,9$
  paths from archived solver rows;
- `14_audit_citations.py`: DOI and bibliography audit.

Raw CSV files are checkpointed after every solved instance. Re-running a phase
skips completed `(instance_id, method/variant)` keys. Summary tables and figures
are generated under `results/summary` and `results/figures`.

## Model-aligned run coverage (2026-08-19)

The fixed design contains 560 prioritization instances (7 sizes x 4 regimes x
20 replications). Existing raw priority results are retained, while the current
paper applies a strict 60-second-per-stage certificate filter: 521 Stage-1
solutions were certified by 60 seconds, and 426 instances have common certified
EM/LLSM/MNVEM/MNVLLSM records with both MNVDM stages within that limit. No timed-out
or missing value is imputed.

The current formulation benchmark has 56 paired instances (112 solves). The
current EVRIM study has 42 fixed high/cyclic instances and three variants per
instance: Direct, B&C/OA, and B&C/OA with protected judgments (126 rows). All
three house cases are optimal under the 60-second-per-stage limit and each revised
PCM is subsequently re-prioritized by MNVLLSM. Raw rows, summaries, figures,
superseded partial runs, and the active environment are under
`results_model_aligned_20260819/`.

## Interpretation boundary

NVR directly measures the optimized violation criterion and therefore cannot by
itself establish decision quality. Kendall tau-b, best-choice accuracy, and LRMSE
compare each estimate with the retained generating vector and answer the separate
recovery question. Timing is environment-sensitive and should not be
compared byte-for-byte across machines.

## Weighted Pareto version (2026-08-20)

The separate weighted version solves
`alpha*NVR + (1-alpha)*D_GCI` at all 101 values `1.00, 0.99, ..., 0.00`
for every supported size `n=3,...,9`. It uses a 60-second total limit per
solve, `epsilon=1e-4`, Gurobi tolerances and MIP gap `1e-5`, and independent
post-solve residual checks. The production relation model uses valid
relation-specific big-M constants derived from the declared `y` bounds, with
`NumericFocus=3`, `ScaleFlag=2`, and `IntegralityFocus=1`; a certified point
that fails the independent residual audit aborts the run. The degenerate
`alpha=1` endpoint is started from the feasible `alpha=0.99` solution, while
each candidate endpoint state is verified in a big-M-free continuous subproblem.
Infeasible states receive no-good cuts; a 60-second endpoint timeout retains the
`alpha=0.99` vector only as an uncertified feasible incumbent. Saved tables are
sorted back to `1.00, 0.99, ..., 0.00`. The pilot contains one cyclic PCM for each of
`n=7,8,9` (303 solves):

```powershell
D:\anaconda3\python.exe experiments\10_alpha_mnvdm_pareto.py --mode pilot
D:\anaconda3\python.exe experiments\11_house_evrim_check_first.py
```

Use `--mode selected` for all 80 `n=7` PCMs plus one `n=8` and one `n=9`
PCM (8,282 solves). Completed pilot rows are skipped automatically.

Use `--mode full` for all 560 PCMs (56,560 solves). The check-first EVRIM
implementation passes no MIP start to its fresh B&C model. All new outputs are
kept under `results_weighted_pareto_selected_20260820`; prior results are untouched.
The completed selected run contains 8,282 unique rows, all independently
feasible, with 8,221 certified optima and 60 time-limit terminations.
