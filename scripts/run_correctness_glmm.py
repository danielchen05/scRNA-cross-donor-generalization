from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil

from _common import bootstrap_repo


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the prediction-level correctness GLMM for one benchmark config. "
            "This script saves numerical outputs only and creates no plots."
        )
    )
    parser.add_argument(
        "--benchmark-config",
        required=True,
        type=Path,
        help="Benchmark YAML config defining the frozen AnnData and representations",
    )
    parser.add_argument(
        "--site-col",
        default=None,
        help=(
            "Optional obs column copied into prediction metadata (for example 'batch'). "
            "It is not added to the GLMM fixed effects."
        ),
    )
    parser.add_argument(
        "--clean-output",
        action="store_true",
        help="Delete the existing mixed_models directory before fitting",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repo_root = bootstrap_repo()
    os.chdir(repo_root)

    from scrna_benchmark.config_io import load_dataset_config
    from scrna_benchmark.mixed_models import run_correctness_glmm_analysis
    from scrna_benchmark.pipeline import prepare_dataset

    config = load_dataset_config(args.benchmark_config)
    out_dir = config.dataset_output_dir / "mixed_models"
    if args.clean_output and out_dir.exists():
        print(f"[glmm] removing {out_dir}")
        shutil.rmtree(out_dir)

    print(f"[glmm] benchmark={config.dataset_name}")
    print(f"[glmm] adata={config.adata_path}")
    print(f"[glmm] classifier batch covariate={config.batch_col}")
    print(f"[glmm] metadata site column={args.site_col}")

    prepared = prepare_dataset(config)
    adata = prepared["adata"]

    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / "glmm_config.json").open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "benchmark_config": str(args.benchmark_config),
                "dataset_name": config.dataset_name,
                "adata_path": str(config.adata_path),
                "representations": config.representations,
                "celltype_col": config.celltype_col,
                "donor_col": config.donor_col,
                "site_col": args.site_col,
                "classifier_batch_col": config.batch_col,
                "test_size": config.test_size,
                "n_random_repeats": config.random_split_n_repeats,
                "n_folds": config.n_folds,
                "random_state": config.random_state,
            },
            handle,
            indent=2,
        )

    run_correctness_glmm_analysis(
        adata=adata,
        representations=config.representations,
        out_dir=out_dir,
        celltype_col=config.celltype_col,
        donor_col=config.donor_col,
        site_col=args.site_col,
        batch_col=config.batch_col,
        test_size=config.test_size,
        n_random_repeats=config.random_split_n_repeats,
        n_folds=config.n_folds,
        random_state=config.random_state,
        dataset_name=config.dataset_name,
        verbose=config.verbose,
    )
    print(f"[glmm] complete -> {out_dir}")


if __name__ == "__main__":
    main()
