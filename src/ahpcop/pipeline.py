from __future__ import annotations

import json
import os
import platform
import subprocess
import time
from pathlib import Path

import gurobipy as gp
import numpy as np
import pandas as pd
from scipy.optimize import linprog

from .evrim import HOUSE_PCM, default_protected_sets, solve_evrim
from .metrics import recovery_metrics
from .pcm import consistency_ratio, generate_pcm, load_record_matrix, load_record_weights, upper_pairs
from .priority import (
    GurobiSettings,
    Stage1Result,
    classical_priorities,
    cop_llsm_from_mnvllsm,
    solve_mnvem,
    solve_mnvllsm,
    solve_stage1,
)
from .reporting import build_all_reports
from .sensitivity import run_example_epsilon_sensitivity


def _json_array(value: np.ndarray | None) -> str:
    return "" if value is None else json.dumps(np.asarray(value).tolist())


def _append_row(path: Path, row: dict) -> None:
    frame = pd.DataFrame([row])
    frame.to_csv(path, mode="a", header=not path.exists(), index=False)


def _completed(path: Path, keys: list[str]) -> set[tuple]:
    if not path.exists() or path.stat().st_size == 0:
        return set()
    frame = pd.read_csv(path)
    return {tuple(row[k] for k in keys) for _, row in frame.iterrows()}


class ExperimentPipeline:
    def __init__(self, config_path: Path):
        self.root = config_path.resolve().parent
        self.config = json.loads(config_path.read_text(encoding="utf-8"))
        self.results = self.root / self.config.get("results_dir", "results")
        self.raw = self.results / "raw"
        self.summary = self.results / "summary"
        self.figures = self.results / "figures"
        self.logs = self.results / "logs"
        for directory in (self.raw, self.summary, self.figures, self.logs):
            directory.mkdir(parents=True, exist_ok=True)
        self.priority_settings = GurobiSettings(
            epsilon=self.config["epsilon"], y_bound=self.config["y_bound"],
            time_limit=self.config["priority_time_limit_seconds"], mip_gap=self.config["mip_gap"],
            threads=self.config["threads"], seed=self.config["gurobi_seed"], output_flag=self.config["output_flag"],
        )
        self.evrim_settings = GurobiSettings(
            epsilon=self.config["epsilon"], y_bound=self.config["y_bound"],
            time_limit=self.config["evrim_time_limit_seconds"], mip_gap=self.config["mip_gap"],
            threads=self.config["threads"], seed=self.config["gurobi_seed"], output_flag=self.config["output_flag"],
        )
        self.formulation_settings = GurobiSettings(
            epsilon=self.config["epsilon"], y_bound=self.config["y_bound"],
            time_limit=self.config["formulation_time_limit_seconds"], mip_gap=self.config["mip_gap"],
            threads=self.config["threads"], seed=self.config["gurobi_seed"], output_flag=self.config["output_flag"],
        )
        self._write_environment()

    @property
    def dataset_path(self) -> Path:
        return self.raw / "datasets.csv"

    def _write_environment(self) -> None:
        try:
            import psutil
            total_ram = int(psutil.virtual_memory().total)
            cpu_name = platform.processor() or subprocess.check_output(
                ["powershell", "-NoProfile", "-Command", "(Get-CimInstance Win32_Processor).Name"],
                text=True, timeout=10,
            ).strip()
        except (ImportError, OSError, subprocess.SubprocessError):
            total_ram = None
            cpu_name = platform.processor()
        environment = {
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "platform": platform.platform(),
            "python": platform.python_version(),
            "processor": cpu_name,
            "cpu_count": os.cpu_count(),
            "ram_bytes": total_ram,
            "gurobi": ".".join(map(str, gp.gurobi.version())),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "config": self.config,
        }
        (self.results / "environment.json").write_text(json.dumps(environment, indent=2), encoding="utf-8")

    def run_generate(self, samples_override: int | None = None) -> None:
        samples = samples_override or int(self.config["samples_per_cell"])
        records = []
        for n in self.config["n_values"]:
            for regime in self.config["regimes"]:
                for replicate in range(samples):
                    instance = generate_pcm(n, regime, replicate, self.config["master_seed"])
                    records.append(instance.as_record())
        pd.DataFrame(records).sort_values(["n", "regime", "replicate"]).to_csv(self.dataset_path, index=False)
        print(f"generated {len(records)} instances -> {self.dataset_path}", flush=True)

    def run_priority(self, limit: int | None = None) -> None:
        if not self.dataset_path.exists():
            self.run_generate()
        data = pd.read_csv(self.dataset_path)
        if limit is not None:
            data = data.head(limit)
        priority_path = self.raw / "priority_results.csv"
        stage1_path = self.raw / "stage1_cache.csv"
        done_priority = _completed(priority_path, ["instance_id", "method"])
        stage1_cache: dict[str, dict] = {}
        if stage1_path.exists():
            stage1_cache = {
                str(cached["instance_id"]): cached.to_dict()
                for _, cached in pd.read_csv(stage1_path).iterrows()
            }
        required_methods = {"EM", "LLSM", "MNVLLSM", "COP-LLSM", "MNVEM"}
        for row_number, row in data.iterrows():
            instance_id = row["instance_id"]
            priority_complete = all((instance_id, method) in done_priority for method in required_methods)
            if priority_complete:
                continue
            a = load_record_matrix(row["matrix"])
            w0 = load_record_weights(row["latent_weights"])
            print(f"priority {row_number + 1}/{len(data)} {instance_id}", flush=True)
            if instance_id in stage1_cache:
                cached = stage1_cache[instance_id]
                stage1 = Stage1Result(
                    status=str(cached["status"]),
                    solved=str(cached["solved"]).strip().lower() == "true",
                    nv=None if pd.isna(cached["nv"]) else float(cached["nv"]),
                    nv2=None if pd.isna(cached["nv2"]) else int(cached["nv2"]),
                    runtime=float(cached["runtime"]),
                    gap=None if pd.isna(cached["gap"]) else float(cached["gap"]),
                    node_count=None if pd.isna(cached["node_count"]) else float(cached["node_count"]),
                    signs=None, y=None, model_variant="indicator",
                )
            else:
                stage1 = solve_stage1(a, self.priority_settings, variant="indicator")
                cache_row = {
                    "instance_id": instance_id, "n": row["n"], "regime": row["regime"],
                    "status": stage1.status, "solved": stage1.solved, "nv": stage1.nv,
                    "nv2": stage1.nv2, "runtime": stage1.runtime, "gap": stage1.gap,
                    "node_count": stage1.node_count, "model_variant": stage1.model_variant,
                }
                _append_row(stage1_path, cache_row)
                stage1_cache[instance_id] = cache_row

            methods = [
                result for result in classical_priorities(a)
                if (instance_id, result.method) not in done_priority
            ]
            need_mnvllsm = any(
                (instance_id, method) not in done_priority for method in ("MNVLLSM", "COP-LLSM")
            )
            if need_mnvllsm:
                mnvllsm = solve_mnvllsm(a, stage1, self.priority_settings)
                methods.extend([mnvllsm, cop_llsm_from_mnvllsm(mnvllsm)])
            if (instance_id, "MNVEM") not in done_priority:
                methods.append(solve_mnvem(a, stage1, self.priority_settings))
            for result in methods:
                if (instance_id, result.method) in done_priority:
                    continue
                metrics = recovery_metrics(a, result.weights, w0) if result.weights is not None else {
                    "nv": None, "n_order_relations": None, "nvr": None,
                    "kendall_tau_b": None, "best_choice_accuracy": None, "lrmse": None,
                }
                _append_row(priority_path, {
                    "instance_id": instance_id, "n": row["n"], "regime": row["regime"],
                    "cr": row["cr"], "gci_input": row["gci"], "realized_cr_bin": row["realized_cr_bin"],
                    "method": result.method, "status": result.status, "solved": result.solved,
                    "runtime": result.runtime, "stage1_runtime": result.stage1_runtime,
                    "stage2_runtime": result.stage2_runtime, "nv_star": result.nv_star,
                    "objective": result.objective, "gap": result.gap,
                    "weights": _json_array(result.weights), **metrics,
                })

    def run_formulation(self, limit: int | None = None) -> None:
        """Compare the basic and strong Stage-1 formulations under equal settings."""
        if not self.dataset_path.exists():
            self.run_generate()
        data = pd.read_csv(self.dataset_path)
        data = data[
            data["replicate"] < int(self.config["formulation_samples_per_cell"])
        ].copy()
        if limit is not None:
            data = data.head(limit)
        path = self.raw / "formulation_runtime.csv"
        done = _completed(path, ["instance_id", "variant"])
        for row_number, row in data.reset_index(drop=True).iterrows():
            a = load_record_matrix(row["matrix"])
            instance_id = row["instance_id"]
            print(f"formulation {row_number + 1}/{len(data)} {instance_id}", flush=True)
            for variant in ("basic", "strong"):
                if (instance_id, variant) in done:
                    continue
                result = solve_stage1(a, self.formulation_settings, variant=variant)
                _append_row(path, {
                    "instance_id": instance_id,
                    "n": row["n"],
                    "regime": row["regime"],
                    "variant": variant,
                    "status": result.status,
                    "solved": result.solved,
                    "nv_star": result.nv,
                    "runtime": result.runtime,
                    "gap": result.gap,
                    "node_count": result.node_count,
                })

    def run_evrim(self, limit: int | None = None) -> None:
        if not self.dataset_path.exists():
            self.run_generate()
        data = pd.read_csv(self.dataset_path)
        subset = data[data["regime"].isin(["high", "cyclic"])].copy()
        subset = subset.groupby(["n", "regime"], group_keys=False).head(self.config["evrim_samples_per_cell"])
        if limit is not None:
            subset = subset.head(limit)
        path = self.raw / "evrim_results.csv"
        done = _completed(path, ["instance_id", "variant"])
        for row_number, row in subset.reset_index(drop=True).iterrows():
            a = load_record_matrix(row["matrix"])
            instance_id = row["instance_id"]
            print(f"evrim {row_number + 1}/{len(subset)} {instance_id}", flush=True)
            protected_value, protected_direction = default_protected_sets(a)
            cases = [
                ("EVRIM-Direct", "direct", [], []),
                ("EVRIM-OA", "oa_callback", [], []),
                ("EVRIM-OA+T", "oa_callback", protected_value, protected_direction),
            ]
            for variant, backend, value_set, direction_set in cases:
                if (instance_id, variant) in done:
                    continue
                result = solve_evrim(
                    a, self.evrim_settings, value_protected=value_set,
                    direction_protected=direction_set, variant=variant, backend=backend,
                )
                _append_row(path, {
                    "instance_id": instance_id, "n": row["n"], "regime": row["regime"],
                    "cr_input": row["cr"], "gci_input": row["gci"], "variant": variant,
                    "backend": backend,
                    "status": result.status, "solved": result.solved, "nrp": result.nrp,
                    "aoc": result.aoc, "gci": result.gci, "nv": result.nv,
                    "runtime": result.runtime, "stage1_runtime": result.stage1_runtime,
                    "stage2_runtime": result.stage2_runtime, "gap": result.gap,
                    "value_protected": json.dumps(result.value_protected),
                    "direction_protected": json.dumps(result.direction_protected),
                    "revised_matrix": _json_array(result.revised_matrix),
                    "weights": _json_array(result.weights),
                })

    def run_house(self) -> None:
        path = self.raw / "house_results.csv"
        done = _completed(path, ["case"])
        cases = [
            ("House-A-unprotected", []),
            ("House-B-protect-a37", [(2, 6)]),
            ("House-C-protect-a13-a37", [(0, 2), (2, 6)]),
        ]
        for case, value_protected in cases:
            if (case,) in done:
                continue
            result = solve_evrim(
                HOUSE_PCM, self.evrim_settings, value_protected=value_protected,
                variant=case, backend="oa_callback",
            )
            if result.solved and result.revised_matrix is not None:
                priority_stage1 = solve_stage1(
                    result.revised_matrix, self.priority_settings, variant="indicator"
                )
                priority = solve_mnvllsm(
                    result.revised_matrix, priority_stage1, self.priority_settings
                )
            else:
                priority_stage1 = None
                priority = None
            _append_row(path, {
                "case": result.variant, "backend": "oa_callback",
                "status": result.status, "solved": result.solved,
                "nrp": result.nrp, "aoc": result.aoc, "gci": result.gci, "nv": result.nv,
                "cr": None if result.revised_matrix is None else consistency_ratio(result.revised_matrix),
                "runtime": result.runtime, "stage1_runtime": result.stage1_runtime,
                "stage2_runtime": result.stage2_runtime, "gap": result.gap,
                "value_protected": json.dumps(result.value_protected),
                "direction_protected": json.dumps(result.direction_protected),
                "revised_matrix": _json_array(result.revised_matrix),
                "certificate_weights": _json_array(result.weights),
                "priority_status": None if priority is None else priority.status,
                "priority_solved": None if priority is None else priority.solved,
                "priority_nv_star": None if priority is None else priority.nv_star,
                "priority_objective": None if priority is None else priority.objective,
                "priority_runtime": None if priority is None else priority.runtime,
                "priority_weights": _json_array(None if priority is None else priority.weights),
            })
            print(f"house {result.variant}: {result.status}", flush=True)

    def run_bnc(self) -> None:
        """Equal-settings direct-MIQCP versus lazy-OA B&C benchmark."""
        if not self.dataset_path.exists():
            self.run_generate()
        data = pd.read_csv(self.dataset_path)
        subset = data[data["regime"] == "high"].groupby("n", group_keys=False).head(
            int(self.config.get("bnc_samples_per_n", 2))
        )
        path = self.raw / "bnc_runtime.csv"
        done = _completed(path, ["instance_id", "backend"])
        for _, row in subset.iterrows():
            a = load_record_matrix(row["matrix"])
            for backend in ("direct", "oa_callback"):
                if (row["instance_id"], backend) in done:
                    continue
                result = solve_evrim(
                    a, self.evrim_settings, variant=f"EVRIM-{backend}", backend=backend,
                )
                _append_row(path, {
                    "instance_id": row["instance_id"], "n": row["n"], "regime": row["regime"],
                    "backend": backend, "status": result.status, "solved": result.solved,
                    "nrp": result.nrp, "aoc": result.aoc, "gci": result.gci, "nv": result.nv,
                    "runtime": result.runtime, "stage1_runtime": result.stage1_runtime,
                    "stage2_runtime": result.stage2_runtime, "gap": result.gap,
                })
                print(f"bnc {row['instance_id']} {backend}: {result.status}", flush=True)

    def run_synthetic_sensitivity(self) -> None:
        """Optional synthetic tolerance batch; not used for the manuscript claim."""
        if not self.dataset_path.exists():
            self.run_generate()
        data = pd.read_csv(self.dataset_path)
        subset = data[
            data["n"].isin(self.config["sensitivity_n_values"])
            & (data["replicate"] == 0)
        ].copy()
        path = self.raw / "epsilon_sensitivity.csv"
        done = _completed(path, ["instance_id", "epsilon"])
        for _, row in subset.iterrows():
            a = load_record_matrix(row["matrix"])
            w0 = load_record_weights(row["latent_weights"])
            for epsilon in self.config["sensitivity_epsilons"]:
                key = (row["instance_id"], float(epsilon))
                if key in done:
                    continue
                settings = GurobiSettings(
                    epsilon=float(epsilon), y_bound=self.priority_settings.y_bound,
                    time_limit=self.priority_settings.time_limit, mip_gap=self.priority_settings.mip_gap,
                    threads=self.priority_settings.threads, seed=self.priority_settings.seed,
                    output_flag=self.priority_settings.output_flag,
                )
                stage1 = solve_stage1(a, settings, "indicator")
                result = solve_mnvllsm(a, stage1, settings)
                metrics = recovery_metrics(a, result.weights, w0) if result.weights is not None else {}
                ranking = None if result.weights is None else np.argsort(-result.weights).astype(int).tolist()
                _append_row(path, {
                    "instance_id": row["instance_id"], "n": row["n"], "regime": row["regime"],
                    "epsilon": float(epsilon), "status": result.status, "solved": result.solved,
                    "nv_star": result.nv_star, "runtime": result.runtime,
                    "objective": result.objective, "ranking": json.dumps(ranking),
                    "kendall_tau_b": metrics.get("kendall_tau_b"),
                    "lrmse": metrics.get("lrmse"), "nvr": metrics.get("nvr"),
                })
                print(f"sensitivity {row['instance_id']} epsilon={epsilon}: {result.status}", flush=True)

    def run_sensitivity(self) -> None:
        """Reproduce the manuscript's exhaustive 4x4 epsilon table."""
        path = self.raw / "example_epsilon_sensitivity.csv"
        frame = run_example_epsilon_sensitivity(
            [float(value) for value in self.config["sensitivity_epsilons"]],
            self.priority_settings,
            path,
        )
        print(frame.to_string(index=False), flush=True)
        print(f"example epsilon sensitivity -> {path}", flush=True)

    def run_sanity(self) -> None:
        """Feasibility sanity check for randomly constructed index-exchangeable PCMs."""
        path = self.raw / "representation_sanity.csv"
        if path.exists():
            path.unlink()
        epsilon = float(self.config["epsilon"])
        seed = int(self.config["master_seed"])
        for n in self.config["n_values"]:
            items = upper_pairs(n) + [(0, 0)]
            item_coefficients = []
            for i, j in items:
                coefficient = np.zeros(n)
                coefficient[i] += 1.0
                coefficient[j] -= 1.0
                item_coefficients.append(coefficient)
            for replicate in range(int(self.config["sanity_samples_per_n"])):
                rng = np.random.default_rng(seed + 100_000 * n + replicate)
                latent_y = rng.normal(size=n)
                latent_y -= latent_y.mean()
                differences = np.array([latent_y[i] - latent_y[j] for i, j in items])
                input_logs = differences + 0.05 * differences ** 3
                a_ub, b_ub, a_eq, b_eq = [], [], [np.ones(n)], [0.0]
                for p in range(len(items)):
                    for q in range(p + 1, len(items)):
                        coefficient = item_coefficients[p] - item_coefficients[q]
                        delta = input_logs[p] - input_logs[q]
                        if delta > 1e-12:
                            a_ub.append(-coefficient)
                            b_ub.append(-epsilon)
                        elif delta < -1e-12:
                            a_ub.append(coefficient)
                            b_ub.append(-epsilon)
                        else:
                            a_eq.append(coefficient)
                            b_eq.append(0.0)
                start = time.perf_counter()
                lp = linprog(
                    np.zeros(n), A_ub=np.asarray(a_ub) if a_ub else None,
                    b_ub=np.asarray(b_ub) if b_ub else None, A_eq=np.asarray(a_eq),
                    b_eq=np.asarray(b_eq), bounds=[(None, None)] * n, method="highs",
                )
                runtime = time.perf_counter() - start
                _append_row(path, {
                    "n": n, "replicate": replicate, "seed": seed + 100_000 * n + replicate,
                    "feasible": bool(lp.success), "status": int(lp.status),
                    "message": lp.message, "runtime": runtime,
                })
            print(f"sanity n={n}: completed", flush=True)

    def run_summarize(self) -> None:
        build_all_reports(self.results)

    def run_all(self) -> None:
        self.run_generate()
        self.run_sanity()
        self.run_priority()
        self.run_formulation()
        self.run_sensitivity()
        self.run_evrim()
        self.run_bnc()
        self.run_house()
        self.run_summarize()

    def run_smoke(self) -> None:
        smoke_config = dict(self.config)
        smoke_config.update({
            "n_values": [3, 4], "regimes": ["low", "moderate", "high", "cyclic"],
            "samples_per_cell": 1, "evrim_samples_per_cell": 1,
            "priority_time_limit_seconds": min(30, self.config["priority_time_limit_seconds"]),
            "evrim_time_limit_seconds": min(45, self.config["evrim_time_limit_seconds"]),
            "results_dir": "results_smoke_v2",
        })
        smoke_path = self.root / "smoke_config.json"
        smoke_path.write_text(json.dumps(smoke_config, indent=2), encoding="utf-8")
        smoke = ExperimentPipeline(smoke_path)
        smoke.run_generate()
        smoke.run_priority()
        smoke.run_formulation()
        smoke.run_evrim()
        smoke.run_summarize()
