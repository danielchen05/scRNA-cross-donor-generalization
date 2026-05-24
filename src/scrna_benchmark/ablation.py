# src/scrna_benchmark/ablation.py

from pathlib import Path

import numpy as np
import pandas as pd

from .models import run_donor_split_logreg
from .splits import sample_train_donors
from .filtering import subset_adata_by_celltypes


def run_donor_ablation(
    adata,
    representations,
    k_values,
    n_repeats,
    celltype_col="cell_type",
    donor_col="patient_id",
    batch_col=None,
    random_state=42,
    verbose=True,
):
    """
    Run donor ablation across representations.

    For each k, sample k training donors and test on all remaining donors.

    Parameters
    ----------
    adata : AnnData
        Input AnnData.
    representations : dict
        Mapping from representation label to AnnData representation key.
        Example: {"pca": "X_pca", "hvg": "hvg"}.
    k_values : list[int]
        Numbers of training donors to test.
    n_repeats : int
        Number of random repeats per k.
    celltype_col : str
        Cell-type column.
    donor_col : str
        Donor column.
    batch_col : str or None
        Optional batch/site covariate.
    random_state : int
        Base seed.
    verbose : bool
        Whether to print progress.

    Returns
    -------
    raw_df : pd.DataFrame
        One row per representation/k/repeat.
    summary_df : pd.DataFrame
        Mean/std summary per representation/k.
    """
    all_donors = np.array(sorted(adata.obs[donor_col].astype(str).unique()))

    results = []

    for rep_label, rep_key in representations.items():
        if verbose:
            print(f"\nRunning donor ablation for {rep_label} ({rep_key})")

        for k in k_values:
            if k >= len(all_donors):
                if verbose:
                    print(f"Skipping k={k}: must be smaller than n_donors={len(all_donors)}")
                continue

            for repeat in range(n_repeats):
                seed = random_state + 100 * k + repeat

                train_donors = sample_train_donors(
                    all_donors,
                    k=k,
                    random_state=seed,
                )
                test_donors = np.array([
                    donor for donor in all_donors
                    if donor not in set(train_donors)
                ])

                res = run_donor_split_logreg(
                    adata=adata,
                    rep_key=rep_key,
                    train_donors=train_donors,
                    test_donors=test_donors,
                    celltype_col=celltype_col,
                    donor_col=donor_col,
                    batch_col=batch_col,
                    random_state=random_state,
                )

                row = {
                    "representation": rep_label,
                    "rep_key": rep_key,
                    "k_train_donors": k,
                    "repeat": repeat + 1,
                    "macro_f1": res["macro_f1"],
                    "accuracy": res["accuracy"],
                    "n_train_cells": res["n_train_cells"],
                    "n_test_cells": res["n_test_cells"],
                    "n_classes_used": res["n_classes_used"],
                    "batch_covariate": batch_col if batch_col is not None else "None",
                }
                results.append(row)

                if verbose:
                    print(
                        f"  k={k:2d}, repeat={repeat + 1}: "
                        f"macro_f1={res['macro_f1']:.4f}, "
                        f"acc={res['accuracy']:.4f}"
                    )

    raw_df = pd.DataFrame(results)
    summary_df = summarize_donor_ablation(raw_df)

    return raw_df, summary_df


def summarize_donor_ablation(raw_df):
    """
    Summarize donor ablation results.
    """
    if raw_df.empty:
        return pd.DataFrame()

    summary_df = (
        raw_df
        .groupby(["representation", "rep_key", "k_train_donors"], as_index=False)
        .agg(
            macro_f1_mean=("macro_f1", "mean"),
            macro_f1_std=("macro_f1", "std"),
            accuracy_mean=("accuracy", "mean"),
            accuracy_std=("accuracy", "std"),
            mean_train_cells=("n_train_cells", "mean"),
            mean_test_cells=("n_test_cells", "mean"),
            mean_n_classes_used=("n_classes_used", "mean"),
        )
    )

    return summary_df


def run_subset_donor_ablation(
    adata,
    celltypes,
    representations,
    k_values,
    n_repeats,
    celltype_col="cell_type",
    donor_col="patient_id",
    batch_col=None,
    random_state=42,
    verbose=True,
):
    """
    Run donor ablation after restricting to selected cell types.

    Useful for T-cell subtype analyses or excluded-cell-type analyses.
    """
    adata_subset = subset_adata_by_celltypes(
        adata=adata,
        celltypes=celltypes,
        celltype_col=celltype_col,
    )

    raw_df, summary_df = run_donor_ablation(
        adata=adata_subset,
        representations=representations,
        k_values=k_values,
        n_repeats=n_repeats,
        celltype_col=celltype_col,
        donor_col=donor_col,
        batch_col=batch_col,
        random_state=random_state,
        verbose=verbose,
    )

    raw_df["analysis"] = "subset_celltypes"
    summary_df["analysis"] = "subset_celltypes"

    return raw_df, summary_df


def save_donor_ablation_outputs(
    raw_df,
    summary_df,
    out_dir,
    prefix="donor_ablation",
):
    """
    Save donor ablation raw and summary tables.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    raw_df.to_csv(out_dir / f"{prefix}_raw.csv", index=False)
    summary_df.to_csv(out_dir / f"{prefix}_summary.csv", index=False)