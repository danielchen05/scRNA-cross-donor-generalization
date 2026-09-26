from __future__ import annotations

"""Synthetic-label negative control for the cross-donor benchmark."""

from collections.abc import Iterable

import numpy as np
import pandas as pd

from .models import run_donor_split_logreg, run_random_split_logreg
from .splits import make_donor_folds


def make_donor_balanced_null_labels(
    adata,
    donor_col: str,
    n_classes: int,
    random_state: int,
    *,
    label_name: str = "_null_label",
    prefix: str = "null",
) -> pd.Series:
    """
    Generate labels that are random with respect to expression while being
    approximately balanced within every donor.

    Within each donor, pseudo-class counts differ by at most one cell. This
    intentionally removes donor-level class-composition signal while preserving
    the real expression matrix, donor structure, and frozen representations.

    Parameters
    ----------
    adata
        AnnData object.
    donor_col
        Column in ``adata.obs`` identifying donors.
    n_classes
        Number of synthetic classes.
    random_state
        Seed controlling the random label assignment.
    label_name
        Name of the returned Series.
    prefix
        Prefix used for pseudo-class labels.

    Returns
    -------
    pandas.Series
        Synthetic labels aligned exactly to ``adata.obs_names``.
    """
    if donor_col not in adata.obs:
        raise KeyError(f"{donor_col!r} not found in adata.obs")
    if n_classes < 2:
        raise ValueError("n_classes must be at least 2.")

    donors = adata.obs[donor_col].astype(str).to_numpy()
    unique_donors = np.array(sorted(np.unique(donors)))

    donor_counts = pd.Series(donors).value_counts()
    too_small = donor_counts[donor_counts < n_classes]
    if not too_small.empty:
        raise ValueError(
            "Every donor must contain at least n_classes cells so that every "
            "synthetic class is represented within every donor. "
            f"Too-small donors: {too_small.to_dict()}"
        )

    rng = np.random.default_rng(random_state)
    labels = np.empty(adata.n_obs, dtype=object)
    class_codes = np.arange(n_classes)

    for donor in unique_donors:
        positions = np.flatnonzero(donors == donor)
        n_cells = len(positions)
        quotient, remainder = divmod(n_cells, n_classes)

        # Start with exactly equal counts.
        codes = np.repeat(class_codes, quotient)

        # Randomly choose which classes receive the at-most-one extra cell.
        if remainder:
            extras = rng.choice(
                class_codes,
                size=remainder,
                replace=False,
            )
            codes = np.concatenate([codes, extras])

        # Break any connection between row order/expression and pseudo-label.
        rng.shuffle(codes)

        labels[positions] = [
            f"{prefix}_{int(code) + 1:02d}"
            for code in codes
        ]

    return pd.Series(
        labels,
        index=adata.obs_names,
        name=label_name,
        dtype="object",
    )


def summarize_null_label_balance(
    adata,
    donor_col: str,
    label_col: str,
) -> pd.DataFrame:
    """Return donor-by-pseudo-class counts in tidy format."""
    if donor_col not in adata.obs:
        raise KeyError(f"{donor_col!r} not found in adata.obs")
    if label_col not in adata.obs:
        raise KeyError(f"{label_col!r} not found in adata.obs")

    balance = (
        pd.crosstab(
            adata.obs[donor_col].astype(str),
            adata.obs[label_col].astype(str),
        )
        .rename_axis(index=donor_col, columns="synthetic_label")
        .reset_index()
        .melt(
            id_vars=donor_col,
            var_name="synthetic_label",
            value_name="n_cells",
        )
    )
    return balance


def run_null_floor_benchmark(
    adata,
    representations: dict[str, str],
    *,
    donor_col: str,
    n_classes: int,
    label_seeds: Iterable[int],
    random_split_seeds: Iterable[int] = (42,),
    n_folds: int = 5,
    donor_fold_seed: int = 42,
    test_size: float = 0.20,
    label_col: str = "_null_label",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Run the synthetic-label floor control without saving prediction-level files.

    For each independent synthetic-label seed:
      1. create donor-balanced random labels;
      2. evaluate one or more random cell-level splits;
      3. evaluate the same 5-fold donor-held-out scheme used in the benchmark.

    Returns
    -------
    runs_df
        One row per representation x evaluation split/fold x label seed.
    replicates_df
        One row per representation x scheme x label seed, averaging the
        evaluation units within that seed.
    deltas_df
        Paired random-minus-donor macro-F1 values for each representation and
        label seed.
    balance_df
        Donor-by-pseudo-class counts for all label seeds.
    """
    label_seeds = [int(x) for x in label_seeds]
    random_split_seeds = [int(x) for x in random_split_seeds]

    if not label_seeds:
        raise ValueError("label_seeds must be non-empty.")
    if not random_split_seeds:
        raise ValueError("random_split_seeds must be non-empty.")

    chance_macro_f1 = 1.0 / float(n_classes)
    donor_folds = make_donor_folds(
        adata,
        donor_col=donor_col,
        n_folds=n_folds,
        random_state=donor_fold_seed,
    )

    run_rows: list[dict[str, object]] = []
    balance_tables: list[pd.DataFrame] = []

    had_existing_label_col = label_col in adata.obs
    existing_label_col = (
        adata.obs[label_col].copy()
        if had_existing_label_col
        else None
    )

    try:
        for label_seed in label_seeds:
            null_labels = make_donor_balanced_null_labels(
                adata,
                donor_col=donor_col,
                n_classes=n_classes,
                random_state=label_seed,
                label_name=label_col,
            )
            adata.obs[label_col] = null_labels

            balance = summarize_null_label_balance(
                adata,
                donor_col=donor_col,
                label_col=label_col,
            )
            balance.insert(0, "label_seed", label_seed)
            balance_tables.append(balance)

            # Random cell-level split(s).
            for rep_label, rep_key in representations.items():
                for split_seed in random_split_seeds:
                    result = run_random_split_logreg(
                        adata=adata,
                        rep_key=rep_key,
                        celltype_col=label_col,
                        batch_col=None,
                        test_size=test_size,
                        random_state=split_seed,
                    )
                    run_rows.append(
                        {
                            "label_seed": label_seed,
                            "scheme": "random_split",
                            "representation": rep_label,
                            "evaluation_type": "split",
                            "evaluation_id": split_seed,
                            "macro_f1": result["macro_f1"],
                            "accuracy": result["accuracy"],
                            "n_test_cells": len(result["y_test"]),
                            "n_classes": result["n_classes_used"],
                            "chance_macro_f1": chance_macro_f1,
                        }
                    )

            # Donor-held-out CV.
            for rep_label, rep_key in representations.items():
                for fold_idx, test_donors in enumerate(donor_folds):
                    train_donors = np.concatenate(
                        [
                            donor_folds[j]
                            for j in range(len(donor_folds))
                            if j != fold_idx
                        ]
                    )
                    result = run_donor_split_logreg(
                        adata=adata,
                        rep_key=rep_key,
                        train_donors=train_donors,
                        test_donors=test_donors,
                        celltype_col=label_col,
                        donor_col=donor_col,
                        batch_col=None,
                        random_state=donor_fold_seed,
                    )
                    run_rows.append(
                        {
                            "label_seed": label_seed,
                            "scheme": "donor_cv",
                            "representation": rep_label,
                            "evaluation_type": "fold",
                            "evaluation_id": fold_idx + 1,
                            "macro_f1": result["macro_f1"],
                            "accuracy": result["accuracy"],
                            "n_test_cells": len(result["y_test"]),
                            "n_classes": result["n_classes_used"],
                            "chance_macro_f1": chance_macro_f1,
                        }
                    )
    finally:
        if had_existing_label_col:
            adata.obs[label_col] = existing_label_col
        elif label_col in adata.obs:
            del adata.obs[label_col]

    runs_df = pd.DataFrame(run_rows)
    balance_df = pd.concat(balance_tables, ignore_index=True)

    replicates_df = (
        runs_df
        .groupby(
            ["label_seed", "scheme", "representation"],
            as_index=False,
        )
        .agg(
            macro_f1_mean=("macro_f1", "mean"),
            macro_f1_std=("macro_f1", "std"),
            accuracy_mean=("accuracy", "mean"),
            accuracy_std=("accuracy", "std"),
            n_evaluations=("macro_f1", "size"),
            n_classes=("n_classes", "max"),
            chance_macro_f1=("chance_macro_f1", "first"),
        )
    )

    paired = (
        replicates_df
        .pivot(
            index=["label_seed", "representation"],
            columns="scheme",
            values="macro_f1_mean",
        )
        .reset_index()
    )

    required_schemes = {"random_split", "donor_cv"}
    missing_schemes = required_schemes.difference(paired.columns)
    if missing_schemes:
        raise RuntimeError(
            f"Could not form paired delta table; missing schemes: {missing_schemes}"
        )

    deltas_df = paired.rename(
        columns={
            "random_split": "random_split_macro_f1",
            "donor_cv": "donor_cv_macro_f1",
        }
    )
    deltas_df["delta_macro_f1"] = (
        deltas_df["random_split_macro_f1"]
        - deltas_df["donor_cv_macro_f1"]
    )
    deltas_df["n_classes"] = n_classes
    deltas_df["chance_macro_f1"] = chance_macro_f1

    return runs_df, replicates_df, deltas_df, balance_df
