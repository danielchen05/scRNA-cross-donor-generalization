# src/scrna_benchmark/splits.py

import numpy as np
import pandas as pd


def make_donor_folds(adata, donor_col, n_folds=5, random_state=42):
    """
    Split unique donors into shuffled donor folds.
    """
    donors = np.array(sorted(adata.obs[donor_col].astype(str).unique()))

    if n_folds > len(donors):
        raise ValueError(
            f"n_folds={n_folds} is larger than number of donors={len(donors)}"
        )

    rng = np.random.default_rng(random_state)
    shuffled = donors.copy()
    rng.shuffle(shuffled)

    return np.array_split(shuffled, n_folds)


def restrict_to_train_labels(y_train, y_test, X_test, obs_test):
    """
    Remove test cells whose labels are absent from the training set.
    """
    train_labels = set(pd.Series(y_train).astype(str).unique())
    keep_mask = pd.Series(y_test).astype(str).isin(train_labels).to_numpy()

    return (
        np.asarray(y_test)[keep_mask],
        X_test[keep_mask],
        obs_test.iloc[keep_mask].copy(),
    )


def sample_train_donors(all_donors, k, random_state=42):
    """
    Randomly sample k training donors without replacement.
    """
    all_donors = np.asarray(all_donors)

    if k >= len(all_donors):
        raise ValueError("k must be smaller than the number of donors.")

    rng = np.random.default_rng(random_state)
    return rng.choice(all_donors, size=k, replace=False)