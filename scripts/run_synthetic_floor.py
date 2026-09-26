from __future__ import annotations

"""Run the PBMC synthetic-label floor negative control."""

import argparse
import json
import os
from pathlib import Path
import shutil

import anndata as ad

from _common import bootstrap_repo


DEFAULT_REPRESENTATIONS = {
    "hvg": "hvg",
    "pca": "X_pca",
    "harmony": "X_harmony",
    "scvi": "X_scVI",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run a synthetic-label floor control on the frozen PBMC benchmark "
            "object. Expression and donor structure stay real; only the outcome "
            "labels are replaced with donor-balanced random pseudo-classes."
        )
    )
    parser.add_argument(
        "--adata",
        type=Path,
        default=Path(
            "data/PBMC_Stephenson/stephenson_benchmark_ready.h5ad"
        ),
        help="Frozen benchmark-ready PBMC AnnData.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/synthetic_floor"),
        help="Directory for compact synthetic-floor numerical outputs.",
    )
    parser.add_argument(
        "--celltype-col",
        default="cell_type",
        help=(
            "Biological cell-type column. Used only to infer the default number "
            "of synthetic classes."
        ),
    )
    parser.add_argument(
        "--donor-col",
        default="patient_id",
        help="Donor identifier column.",
    )
    parser.add_argument(
        "--n-classes",
        type=int,
        default=None,
        help=(
            "Number of synthetic classes. Default: number of unique biological "
            "cell types in --celltype-col."
        ),
    )
    parser.add_argument(
        "--n-label-seeds",
        type=int,
        default=10,
        help="Number of independent synthetic-label randomizations.",
    )
    parser.add_argument(
        "--base-label-seed",
        type=int,
        default=1000,
        help="First synthetic-label seed.",
    )
    parser.add_argument(
        "--random-split-repeats",
        type=int,
        default=1,
        help=(
            "Random cell-level splits per synthetic-label seed. Default 1 keeps "
            "the negative control lightweight; label randomization provides the "
            "main replication."
        ),
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Base seed for evaluation splits and donor folds.",
    )
    parser.add_argument(
        "--n-folds",
        type=int,
        default=5,
        help="Number of donor-held-out folds.",
    )
    parser.add_argument(
        "--test-size",
        type=float,
        default=0.20,
        help="Random cell-level test fraction.",
    )
    parser.add_argument(
        "--clean-output",
        action="store_true",
        help="Delete the output directory before running.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repo_root = bootstrap_repo()
    os.chdir(repo_root)

    from scrna_benchmark.null_control import run_null_floor_benchmark
    from scrna_benchmark.representations import validate_representations

    adata_path = args.adata
    if not adata_path.is_absolute():
        adata_path = repo_root / adata_path

    output_dir = args.output_dir
    if not output_dir.is_absolute():
        output_dir = repo_root / output_dir

    if args.clean_output and output_dir.exists():
        print(f"[synthetic-floor] removing {output_dir}")
        shutil.rmtree(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    if not adata_path.exists():
        raise FileNotFoundError(adata_path)

    print(f"[synthetic-floor] loading {adata_path}")
    adata = ad.read_h5ad(adata_path)

    for col in [args.celltype_col, args.donor_col]:
        if col not in adata.obs:
            raise KeyError(
                f"{col!r} not found in adata.obs. "
                f"Available columns include: {adata.obs.columns.tolist()}"
            )

    validate_representations(
        adata,
        DEFAULT_REPRESENTATIONS,
    )

    n_classes = (
        args.n_classes
        if args.n_classes is not None
        else adata.obs[args.celltype_col].astype(str).nunique()
    )
    if n_classes < 2:
        raise ValueError("Need at least two synthetic classes.")

    if args.n_label_seeds < 1:
        raise ValueError("--n-label-seeds must be at least 1.")
    if args.random_split_repeats < 1:
        raise ValueError("--random-split-repeats must be at least 1.")

    label_seeds = [
        args.base_label_seed + i
        for i in range(args.n_label_seeds)
    ]
    random_split_seeds = [
        args.random_state + i
        for i in range(args.random_split_repeats)
    ]

    min_cells_per_donor = (
        adata.obs[args.donor_col]
        .astype(str)
        .value_counts()
        .min()
    )
    if min_cells_per_donor < n_classes:
        raise ValueError(
            "Cannot balance all synthetic classes within every donor: "
            f"minimum donor size={min_cells_per_donor}, n_classes={n_classes}."
        )

    print(
        "[synthetic-floor] "
        f"cells={adata.n_obs}, "
        f"donors={adata.obs[args.donor_col].astype(str).nunique()}, "
        f"synthetic_classes={n_classes}, "
        f"label_seeds={len(label_seeds)}, "
        f"random_splits_per_seed={len(random_split_seeds)}, "
        f"donor_folds={args.n_folds}"
    )

    runs, replicates, deltas, balance = run_null_floor_benchmark(
        adata,
        DEFAULT_REPRESENTATIONS,
        donor_col=args.donor_col,
        n_classes=n_classes,
        label_seeds=label_seeds,
        random_split_seeds=random_split_seeds,
        n_folds=args.n_folds,
        donor_fold_seed=args.random_state,
        test_size=args.test_size,
    )

    runs.to_csv(
        output_dir / "null_control_runs.csv",
        index=False,
    )
    replicates.to_csv(
        output_dir / "null_control_replicates.csv",
        index=False,
    )
    deltas.to_csv(
        output_dir / "null_control_deltas.csv",
        index=False,
    )
    balance.to_csv(
        output_dir / "null_label_balance.csv",
        index=False,
    )

    config = {
        "analysis": "synthetic_label_floor_control",
        "source_dataset": "PBMC Stephenson benchmark-ready object",
        "adata_path": str(adata_path.relative_to(repo_root)),
        "celltype_col_used_only_for_class_count": args.celltype_col,
        "donor_col": args.donor_col,
        "representations": DEFAULT_REPRESENTATIONS,
        "n_cells": int(adata.n_obs),
        "n_donors": int(
            adata.obs[args.donor_col].astype(str).nunique()
        ),
        "n_classes": int(n_classes),
        "chance_macro_f1": 1.0 / float(n_classes),
        "label_seeds": label_seeds,
        "random_split_seeds": random_split_seeds,
        "n_folds": int(args.n_folds),
        "donor_fold_seed": int(args.random_state),
        "test_size": float(args.test_size),
        "batch_covariate": None,
        "note": (
            "Synthetic labels are generated independently within every donor, "
            "with per-class counts differing by at most one cell. Frozen PBMC "
            "expression and representations are reused without retraining."
        ),
    }
    with open(
        output_dir / "config.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(config, handle, indent=2)

    # Hard sanity check: within every donor and label seed, class counts differ
    # by at most one cell.
    imbalance = (
        balance
        .groupby(["label_seed", args.donor_col])["n_cells"]
        .agg(lambda x: int(x.max() - x.min()))
    )
    if int(imbalance.max()) > 1:
        raise AssertionError(
            "Synthetic labels are not balanced within donor as intended."
        )

    summary = (
        replicates
        .groupby(["scheme", "representation"], as_index=False)
        .agg(
            macro_f1_mean=("macro_f1_mean", "mean"),
            macro_f1_sd_across_label_seeds=("macro_f1_mean", "std"),
        )
    )
    print("\n[synthetic-floor] summary across label randomizations")
    print(summary.to_string(index=False))
    print(
        "\n[synthetic-floor] mean paired delta by representation"
    )
    print(
        deltas
        .groupby("representation")["delta_macro_f1"]
        .agg(["mean", "std"])
        .to_string()
    )
    print(f"\n[synthetic-floor] complete -> {output_dir}")


if __name__ == "__main__":
    main()
