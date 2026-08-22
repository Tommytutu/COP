from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from cop_experiments.pipeline import ExperimentPipeline  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the COP Python/Gurobi experiments")
    parser.add_argument(
        "command",
        choices=["generate", "sanity", "priority", "formulation", "sensitivity", "evrim", "bnc", "house", "summarize", "all", "smoke"],
    )
    parser.add_argument(
        "--config", default=str(ROOT / "config_model_aligned_20260819.json")
    )
    args = parser.parse_args()
    pipeline = ExperimentPipeline(Path(args.config))
    getattr(pipeline, f"run_{args.command}")()


if __name__ == "__main__":
    main()
