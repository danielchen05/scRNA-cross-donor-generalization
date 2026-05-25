# src/scrna_benchmark/covariates.py

from __future__ import annotations

import numpy as np
import pandas as pd


def make_covariate_matrix(
    obs: pd.DataFrame,
    covariate_cols: list[str] | str | None,
    drop_first: bool = False,
    prefix_sep: str = "=",
) -> tuple[np.ndarray | None, list[str]]:
    """One-hot encode one or more categorical covariates.

    Returns None and an empty feature list if covariate_cols is None or empty.
    """
    if covariate_cols is None:
        return None, []

    if isinstance(covariate_cols, str):
        covariate_cols = [covariate_cols]

    if len(covariate_cols) == 0:
        return None, []

    missing = [col for col in covariate_cols if col not in obs.columns]
    if missing:
        raise KeyError(f"Missing covariate columns: {missing}")

    cov_df = pd.get_dummies(
        obs[covariate_cols].astype(str),
        columns=covariate_cols,
        prefix=covariate_cols,
        prefix_sep=prefix_sep,
        drop_first=drop_first,
    )

    return cov_df.to_numpy(dtype=float), cov_df.columns.tolist()


def align_covariate_frames(
    train_obs: pd.DataFrame,
    test_obs: pd.DataFrame,
    covariate_cols: list[str] | str | None,
    drop_first: bool = False,
    prefix_sep: str = "=",
) -> tuple[np.ndarray | None, np.ndarray | None, list[str]]:
    """One-hot encode train/test covariates and align their columns."""
    if covariate_cols is None or covariate_cols == []:
        return None, None, []

    if isinstance(covariate_cols, str):
        covariate_cols = [covariate_cols]

    missing_train = [col for col in covariate_cols if col not in train_obs.columns]
    missing_test = [col for col in covariate_cols if col not in test_obs.columns]
    if missing_train or missing_test:
        raise KeyError(
            f"Missing covariates. train missing={missing_train}, test missing={missing_test}"
        )

    train_df = pd.get_dummies(
        train_obs[covariate_cols].astype(str),
        columns=covariate_cols,
        prefix=covariate_cols,
        prefix_sep=prefix_sep,
        drop_first=drop_first,
    )
    test_df = pd.get_dummies(
        test_obs[covariate_cols].astype(str),
        columns=covariate_cols,
        prefix=covariate_cols,
        prefix_sep=prefix_sep,
        drop_first=drop_first,
    )

    train_df, test_df = train_df.align(test_df, join="outer", axis=1, fill_value=0)

    return (
        train_df.to_numpy(dtype=float),
        test_df.to_numpy(dtype=float),
        train_df.columns.tolist(),
    )


def append_covariates(
    X: np.ndarray,
    covariate_matrix: np.ndarray | None,
) -> np.ndarray:
    """Append a covariate matrix to a representation matrix."""
    if covariate_matrix is None:
        return X

    if X.shape[0] != covariate_matrix.shape[0]:
        raise ValueError(
            f"Row mismatch: X has {X.shape[0]} rows but covariates have "
            f"{covariate_matrix.shape[0]} rows."
        )

    return np.hstack([X, covariate_matrix])
