# COP Python/Gurobi 实验项目

本项目对应最终模型对齐版论文 `Omega_Submission_Model_Aligned_Final` 中的计算实验。
数学模型集中在 `src/cop_experiments`，论文中的每项实验分别放在
`experiments` 下的独立 Python 文件中；`examples` 用于处理读者自己的单个
PCM，不会启动整批 Monte Carlo 实验。

## 1. 代码流程

```text
输入 PCM / 合成 PCM
        |
        v
Stage 1: 最小化 NV
        |
        +--> NV*=0：PCM 存在零违反的 order-preserving representation
        |
        v
Stage 2: 固定 NV=NV*，最小化 LLSM 或 EM 偏差
        |
        v
MNVLLSM / MNVEM 权重、排序、NVR、运行时间

若需要修改判断：
PCM --> EVRIM Stage 1: 最小化 NRP
    --> EVRIM Stage 2: 固定 NRP=NRP*，最小化 AOC
    --> revised PCM --> 再运行 MNVDM
```

## 2. 安装

需要 Python 3.11 以上、Gurobi 12 和有效的 Gurobi license。正式发布后可以
直接从 PyPI 安装：

```powershell
python -m pip install cop-gurobi-experiments
```

从 GitHub 源码安装和测试：

```powershell
git clone https://github.com/Tommytutu/COP.git
cd COP
python -m pip install -e .
python -m pytest -q
```

## 3. 每个论文实验的独立文件

| 文件 | 对应实验 | 主要输出 |
|---|---|---|
| `experiments/01_generate_datasets.py` | Low/Moderate/High/Cyclic 数据 | `results/raw/datasets.csv` |
| `experiments/02_representation_sanity.py` | Proposition 1，n=3,...,9 sanity check | `representation_sanity.csv` |
| `experiments/03a_decision_quality.py` | EM/LLSM/MNVEM/MNVLLSM、COP-LLSM、NVR 与生成权重恢复 | `priority_results.csv`, `decision_quality*.csv` |
| `experiments/03b_order_representability.py` | 按 n 和 regime 计算 P(NV*=0) | `stage1_cache.csv`, `representability.csv` |
| `experiments/03c_mnvdm_runtime.py` | MNVDM runtime 与 certification rate | `priority_results.csv`, `priority_runtime.csv` |
| `experiments/04_basic_vs_strong.py` | basic 与 strong formulation | `formulation_runtime.csv` |
| `experiments/05_epsilon_sensitivity.py` | 4,683 个 weak orders 穷举与 epsilon sensitivity | `example_epsilon_sensitivity.csv` |
| `experiments/06_evrim_protected.py` | difficult PCM、EVRIM 与 EVRIM+T | `evrim_results.csv` |
| `experiments/07_evrim_bnc_benchmark.py` | direct MIQCP 与 OA/B&C | `bnc_runtime.csv` |
| `experiments/08_house_buying_evrim.py` | house-buying Cases A/B/C | `house_results.csv` |
| `experiments/09_build_tables_and_figures.py` | 汇总论文表格和图片 | `results/summary`, `results/figures` |
| `experiments/10_alpha_mnvdm_pareto.py` | 加权 MNVDM 的 alpha 路径 | `results_weighted_pareto_all_n_20260821` |
| `experiments/11_generate_alpha_tables.py` | 生成正文 Table 4 和补充表 | `latex_project/table_alpha_tradeoff_*.tex` |
| `experiments/12_table7_evrim_check_first.py` | 按新 Algorithm 1 重跑 Table 7 的两个 OA 变体 | `results_table7_check_first_20260821` |
| `experiments/13_regenerate_pareto_figure.py` | 从已保存结果重画仅含 $n=7,8,9$ 的 Figure 2 | `alpha_mnvdm_pareto_selected.pdf` |
| `experiments/14_audit_citations.py` | DOI、引文键和未引用条目检查 | `latex_project/citation_crossref_candidates.json` |

依次运行完整实验：

```powershell
D:\anaconda3\python.exe experiments\01_generate_datasets.py
D:\anaconda3\python.exe experiments\02_representation_sanity.py
D:\anaconda3\python.exe experiments\03a_decision_quality.py
D:\anaconda3\python.exe experiments\03b_order_representability.py
D:\anaconda3\python.exe experiments\03c_mnvdm_runtime.py
D:\anaconda3\python.exe experiments\04_basic_vs_strong.py
D:\anaconda3\python.exe experiments\05_epsilon_sensitivity.py
D:\anaconda3\python.exe experiments\06_evrim_protected.py
D:\anaconda3\python.exe experiments\07_evrim_bnc_benchmark.py
D:\anaconda3\python.exe experiments\08_house_buying_evrim.py
D:\anaconda3\python.exe experiments\09_build_tables_and_figures.py
```

快速检查入口可使用 `--limit`：

```powershell
D:\anaconda3\python.exe experiments\03a_decision_quality.py --limit 2
D:\anaconda3\python.exe experiments\03b_order_representability.py --limit 2
D:\anaconda3\python.exe experiments\03c_mnvdm_runtime.py --limit 2
D:\anaconda3\python.exe experiments\04_basic_vs_strong.py --limit 2
D:\anaconda3\python.exe experiments\06_evrim_protected.py --limit 1
```

`03a`、`03b` 和 `03c` 使用相同的 Stage-1/Stage-2 原始求解结果；首次运行其中
任何一个会完成所需求解，之后运行另外两个会直接复用 checkpoint，而不会重复
优化。所有脚本都会识别结果 CSV 中已经完成的键并跳过，便于中断后继续。若要建立完全独立
的新结果目录，复制 `config.json` 并添加例如
`"results_dir": "results_new"`，然后传入 `--config new_config.json`。

## 4. 我有一个矩阵，如何运行 MNVDM

最直接的文件是 `examples/example_01_mnvdm_matrix.py`。只需替换其中的 `A`：

```python
A = np.array([
    [1,   2,   4,   9],
    [1/2, 1,   3,   7],
    [1/4, 1/3, 1,   5],
    [1/9, 1/7, 1/5, 1],
], dtype=float)
```

然后运行：

```powershell
D:\anaconda3\python.exe examples\example_01_mnvdm_matrix.py
```

输出包括求解状态、权重、排序、返回权重的 NV/NVR、Stage-1 最优
`NV*`、Stage-2 偏差和两个阶段的运行时间。默认使用 MNVLLSM；将
`method="LLSM"` 改成 `method="EM"` 即运行 MNVEM。

该示例在当前环境中的实际 MNVLLSM 结果为：

| 项目 | 结果 |
|---|---|
| Status | `OPTIMAL`, certified=True |
| 权重 | `(0.50197724, 0.31379938, 0.14331878, 0.04090460)` |
| 排序 | `x1 > x2 > x3 > x4` |
| Stage-1 `NV*` | `0` |
| 返回权重的 NV / NVR | `0 / 0` |
| Stage-2 LLSM objective | `0.132568674` |

同一矩阵下，经典 EM 和 LLSM 的 NV 均为 1，而 MNVEM 与 MNVLLSM 的
NV 均为 0。运行时间依赖机器和 Gurobi 环境，不应要求与上述实跑时间逐位一致。

也可以从无表头的 CSV 读取：

```powershell
D:\anaconda3\python.exe examples\example_04_mnvdm_from_csv.py examples\sample_matrix.csv
D:\anaconda3\python.exe examples\example_04_mnvdm_from_csv.py examples\sample_matrix.csv --method EM
```

输入必须是正的 reciprocal PCM：对角线为 1，而且
`a[i,j] * a[j,i] = 1`。接口会先检查这些条件，错误矩阵不会静默求解。

## 5. 其他单矩阵示例

```powershell
# 同一矩阵比较 EM、LLSM、MNVEM 和 MNVLLSM
D:\anaconda3\python.exe examples\example_02_compare_methods.py

# 对一个显式循环矩阵运行 EVRIM
D:\anaconda3\python.exe examples\example_03_evrim_repair.py
```

核心单矩阵 API 为：

```python
from cop_experiments import GurobiSettings, solve_mnvdm

settings = GurobiSettings(time_limit=60, threads=1)
result = solve_mnvdm(A, method="LLSM", settings=settings)
print(result.weights, result.nv_star, result.objective)
```

最终证书实验使用 `config_model_aligned_20260819.json`，每个字典序阶段上限为
60 秒；结果在 `results_model_aligned_20260819`。`epsilon` 是有限精度下实现严格不等式的
数值容差；`y_bound` 是模型声明的 log-weight 计算域，不应在未说明的情况下随意
改变后与论文结果直接比较。

## 6. 最终结果口径

- 560 个固定矩阵不重新生成；统一按每阶段 60 秒筛选最优性证书。
- 521 个 Stage-1 结果在 60 秒内认证；四方法共同认证集为 426 个实例。
- basic/strong 使用 56 个配对实例，共 112 次求解。
- EVRIM 使用 42 个 high/cyclic 实例；论文 Table 7 保留 Direct 原值，并用 `12_table7_evrim_check_first.py` 重跑两个 OA 变体。
- check-first OA 先求不含 GCI 约束的两阶段最优解；仅当该认证解不满足 GCI 阈值时，才启动新的 OA 求解，且不传递 MIP start。
- EVRIM 的 POIP 证书变量 `y` 与 GCI 的 LLSM 变量 `u` 分开；house A/B/C 修订后都再次运行 MNVLLSM。
- 原始 CSV、汇总、图片、环境和被替代的中间运行均保存在 `results_model_aligned_20260819`。

## 7. 加权单目标与 Pareto 版本（2026-08-20）

本副本新增单目标模型 `min alpha*NVR + (1-alpha)*D_GCI`，其中所有
`n=3,...,9` 都统一使用 `alpha=1,0.99,...,0` 共 101 个取值。旧的
lexicographic 结果仅作为对照保留，不会被覆盖。数值设置为
`epsilon=1e-4`、Gurobi tolerances/MIPGap `1e-5`、每次求解 60 秒，并设置
`NumericFocus=3`、`ScaleFlag=2` 和 `IntegralityFocus=1`。POIP/COP 三分状态
使用由 `y` 显式边界推导的 relation-specific big-M，不使用统一常数，也不使用
indicator constraint。每个解还会由 Python 独立重算目标和约束残差；任何已认证
但未通过独立残差检查的点都会使批次立即停止。

先运行 `n=7,8,9` 各一个 cyclic PCM 的 303 次 pilot：

```powershell
D:\anaconda3\python.exe experiments\10_alpha_mnvdm_pareto.py --mode pilot
```

当前选定设计（`n=7` 全部 80 个 PCM，并保留 `n=8,9` 各一个 PCM，共
8,282 次）使用：

```powershell
D:\anaconda3\python.exe experiments\10_alpha_mnvdm_pareto.py --mode selected
```

运行全部 560 个 PCM（56,560 次）：

```powershell
D:\anaconda3\python.exe experiments\10_alpha_mnvdm_pareto.py --mode full
```

运行无热启动的 check-first House EVRIM：

```powershell
D:\anaconda3\python.exe experiments\11_house_evrim_check_first.py
```

结果写入 `results_weighted_pareto_selected_20260820`。原始逐次结果在 `raw/`，
关系级证据在 `details/`，Gurobi 日志在 `logs/`，Pareto 点、运行时间汇总和
图片在 `summary/` 与 `figures/`。程序按 `(instance_id, alpha)` 断点续跑。
为稳定处理纯 NVR 的退化端点，实际计算顺序为 `0.99, 1.00, 0.98, ..., 0`，
最终文件仍按 `1.00, 0.99, ..., 0` 排序，并检查 `alpha=1` 的 NVR 不得高于
`alpha=0.99`。每个 `alpha=1` 离散状态还会固定后用无 big-M 连续子问题验证；
不可行状态用 no-good cut 排除。若 60 秒内未证明端点最优，则只保留
`alpha=0.99` 的可行向量作为 `TIME_LIMIT` incumbent。

本次选定设计已完成 8,282 个唯一组合，全部通过独立可行性检查；其中
8,221 个具有最优性证书，60 个达到 60 秒上限，3 个纯 NVR 端点使用了明确
标记的可行回退值。

单矩阵示例：

```powershell
D:\anaconda3\python.exe examples\example_05_weighted_mnvdm_pareto.py
```
