# src/scrna_benchmark/pipeline.py

"""End-to-end orchestration for one dataset benchmark notebooks."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import pandas as pd

from .ablation import run_donor_ablation, save_donor_ablation_outputs
from .config import DatasetConfig
from .cross_group import run_bidirectional_group_transfer, save_group_transfer_outputs
from .experiments import run_donor_cv_experiment, run_random_split_experiment
from .filtering import filter_celltypes_by_support, save_celltype_filtering_outputs
from .metadata import (
    save_metadata_summary,
    standardize_metadata,
    summarize_metadata,
)
from .representations import validate_representations


def save_config(config: DatasetConfig, out_file: str | Path) -> None:
    """Save DatasetConfig as JSON."""
    out_file = Path(out_file)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(config.to_dict(), f, indent=2)


def load_dataset(config: DatasetConfig):
    """Load the AnnData object for a benchmark dataset."""
    if not config.adata_path.exists():
        raise FileNotFoundError(f"AnnData file not found: {config.adata_path}")
    return ad.read_h5ad(config.adata_path)


def prepare_dataset(config: DatasetConfig):
    """Load, standardize metadata, filter cell types, and validate representations."""
    adata = load_dataset(config)

    adata = standardize_metadata(
        adata,
        required_cols=config.required_obs_columns(),
        rename_map=config.rename_obs,
        drop_missing=True,
        stringify=True,
    )

    summary_full = summarize_metadata(
        adata,
        celltype_col=config.celltype_col,
        donor_col=config.donor_col,
        batch_col=config.batch_col,
        group_col=config.group_col,
    )

    if config.min_cells is not None or config.min_donors is not None:
        adata_filtered, support_df, kept_celltypes, excluded_celltypes = (
            filter_celltypes_by_support(
                adata,
                celltype_col=config.celltype_col,
                donor_col=config.donor_col,
                min_cells=config.min_cells if config.min_cells is not None else 0,
                min_donors=config.min_donors if config.min_donors is not None else 0,
            )
        )
    else:
        adata_filtered = adata.copy()
        support_df = pd.DataFrame()
        kept_celltypes = sorted(adata.obs[config.celltype_col].astype(str).unique())
        excluded_celltypes = []

    summary_filtered = summarize_metadata(
        adata_filtered,
        celltype_col=config.celltype_col,
        donor_col=config.donor_col,
        batch_col=config.batch_col,
        group_col=config.group_col,
    )

    validate_representations(adata_filtered, config.representations)

    return {
        "adata_full": adata,
        "adata": adata_filtered,
        "support_df": support_df,
        "kept_celltypes": kept_celltypes,
        "excluded_celltypes": excluded_celltypes,
        "summary_full": summary_full,
        "summary_filtered": summary_filtered,
    }


def run_dataset_benchmark(config: DatasetConfig) -> dict[str, object]:
    """Run the configured benchmark for one dataset.

    This is the main function to call from each dataset-specific notebook.
    """
    out_dir = config.dataset_output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    save_config(config, out_dir / "config.json")

    prepared = prepare_dataset(config)
    adata = prepared["adata"]

    save_metadata_summary(prepared["summary_full"], out_dir / "metadata_full.csv")
    save_metadata_summary(prepared["summary_filtered"], out_dir / "metadata_filtered.csv")

    if not prepared["support_df"].empty:
        if config.save_filtered_adata:
            save_celltype_filtering_outputs(
                adata_full=prepared["adata_full"],
                adata_filtered=adata,
                support_df=prepared["support_df"],
                kept_celltypes=prepared["kept_celltypes"],
                excluded_celltypes=prepared["excluded_celltypes"],
                out_dir=config.filtering_output_dir,
                dataset_name=config.dataset_name,
            )
        else:
            config.filtering_output_dir.mkdir(parents=True, exist_ok=True)
            prepared["support_df"].to_csv(
                config.filtering_output_dir / "celltype_support.csv",
                index=False,
            )
            pd.Series(prepared["kept_celltypes"], name="cell_type").to_csv(
                config.filtering_output_dir / "kept_cell_types.csv",
                index=False,
            )
            pd.Series(prepared["excluded_celltypes"], name="cell_type").to_csv(
                config.filtering_output_dir / "excluded_cell_types.csv",
                index=False,
            )

    outputs: dict[str, object] = {"prepared": prepared}

    if config.run_random_split:
        if config.random_split_seeds is not None:
            random_seeds = config.random_split_seeds
        else:
            random_seeds = [
                config.random_state + i
                for i in range(config.random_split_n_repeats)
            ]

        random_tables = []

        for seed in random_seeds:
            seed_dir = config.random_output_dir / f"seed_{seed}"

            metrics_df = run_random_split_experiment(
                adata=adata,
                representations=config.representations,
                results_dir=seed_dir,
                scheme_label="random_split",
                celltype_col=config.celltype_col,
                batch_col=config.batch_col,
                test_size=config.test_size,
                random_state=seed,
            )

            metrics_df["repeat_seed"] = seed
            random_tables.append(metrics_df)

        random_metrics = pd.concat(random_tables, ignore_index=True)

        config.random_output_dir.mkdir(parents=True, exist_ok=True)
        random_metrics.to_csv(
            config.random_output_dir / "random_split_repeated_metrics.csv",
            index=False,
        )

        outputs["random_split_metrics"] = random_metrics

    if config.run_donor_cv:
        donor_metrics = run_donor_cv_experiment(
            adata=adata,
            representations=config.representations,
            results_dir=config.donor_cv_output_dir,
            scheme_label="donor_cv",
            celltype_col=config.celltype_col,
            donor_col=config.donor_col,
            batch_col=config.batch_col,
            n_folds=config.n_folds,
            random_state=config.random_state,
        )
        outputs["donor_cv_metrics"] = donor_metrics

    if config.run_group_transfer:
        if config.group_col is None:
            raise ValueError("run_group_transfer=True requires config.group_col.")
        if len(config.group_transfer_pairs) == 0:
            raise ValueError("run_group_transfer=True requires group_transfer_pairs.")

        transfer_tables = []
        transfer_results = {}

        for group_a, group_b in config.group_transfer_pairs:
            metrics_df, results = run_bidirectional_group_transfer(
                adata=adata,
                representations=config.representations,
                group_col=config.group_col,
                group_a=group_a,
                group_b=group_b,
                celltype_col=config.celltype_col,
                donor_col=config.donor_col,
                batch_col=config.batch_col,
                random_state=config.random_state,
                verbose=config.verbose,
            )
            pair_name = f"{group_a}_vs_{group_b}"
            metrics_df["pair"] = pair_name
            transfer_tables.append(metrics_df)
            transfer_results[pair_name] = results

            save_group_transfer_outputs(
                metrics_df=metrics_df,
                results=results,
                out_dir=config.group_transfer_output_dir / pair_name,
                prefix="group_transfer",
            )

        outputs["group_transfer_metrics"] = pd.concat(
            transfer_tables,
            ignore_index=True,
        )
        outputs["group_transfer_results"] = transfer_results

    if config.run_donor_ablation:
        if len(config.donor_ablation_k_values) == 0:
            raise ValueError(
                "run_donor_ablation=True requires donor_ablation_k_values."
            )

        raw_df, summary_df = run_donor_ablation(
            adata=adata,
            representations=config.representations,
            k_values=config.donor_ablation_k_values,
            n_repeats=config.donor_ablation_n_repeats,
            celltype_col=config.celltype_col,
            donor_col=config.donor_col,
            batch_col=config.batch_col,
            random_state=config.random_state,
            verbose=config.verbose,
        )

        save_donor_ablation_outputs(
            raw_df=raw_df,
            summary_df=summary_df,
            out_dir=config.donor_ablation_output_dir,
            prefix="donor_ablation",
        )

        outputs["donor_ablation_raw"] = raw_df
        outputs["donor_ablation_summary"] = summary_df

    return outputs
