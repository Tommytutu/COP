"""Shared command-line setup for the independent experiment entry points."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC = PROJECT_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ahpcop.pipeline import ExperimentPipeline  # noqa: E402


def parser(description: str) -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(description=description)
    command.add_argument(
        "--config",
        type=Path,
        default=PROJECT_ROOT / "config_model_aligned_20260819.json",
        help="JSON configuration file (default: model-aligned 60-second configuration)",
    )
    return command


def pipeline(config: Path) -> ExperimentPipeline:
    return ExperimentPipeline(config.resolve())
