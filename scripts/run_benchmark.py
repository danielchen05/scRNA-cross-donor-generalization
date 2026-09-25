from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil

from _common import bootstrap_repo


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run one plot-free benchmark from a YAML DatasetConfig."
    )
    parser.add_argument("--config", required=True, type=Path, help="Benchmark YAML config")
    parser.add_argument(
        "--clean-output",
        action="store_true",
        help="Delete this dataset's existing result directory before running",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repo_root = bootstrap_repo()
    os.chdir(repo_root)

    from scrna_benchmark.config_io import load_dataset_config
    from scrna_benchmark.pipeline import run_dataset_benchmark

    config = load_dataset_config(args.config)
    if args.clean_output and config.dataset_output_dir.exists():
        print(f"[benchmark] removing {config.dataset_output_dir}")
        shutil.rmtree(config.dataset_output_dir)

    print(f"[benchmark] dataset={config.dataset_name}")
    print(f"[benchmark] adata={config.adata_path}")
    run_dataset_benchmark(config)
    print(f"[benchmark] complete -> {config.dataset_output_dir}")


if __name__ == "__main__":
    main()
