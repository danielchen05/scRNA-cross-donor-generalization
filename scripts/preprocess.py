from __future__ import annotations

import argparse
import os
from pathlib import Path

from _common import bootstrap_repo


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create one frozen benchmark-ready AnnData object from a YAML config."
    )
    parser.add_argument("--config", required=True, type=Path, help="Preprocessing YAML config")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite an existing benchmark-ready .h5ad",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repo_root = bootstrap_repo()
    os.chdir(repo_root)

    from scrna_benchmark.preprocessing import load_preprocessing_config, run_preprocessing

    config = load_preprocessing_config(args.config)
    print(f"[preprocess] dataset={config.dataset_name}")
    print(f"[preprocess] input={config.input_path}")
    print(f"[preprocess] output={config.output_path}")
    adata = run_preprocessing(config, force=args.force)
    print(
        f"[preprocess] complete: {adata.n_obs:,} cells x {adata.n_vars:,} genes -> "
        f"{config.output_path}"
    )


if __name__ == "__main__":
    main()
