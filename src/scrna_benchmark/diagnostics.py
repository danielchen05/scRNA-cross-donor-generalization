"""Support diagnostics for PBMC directional-transfer analyses.

These helpers quantify cell-type support and test whether directional
site-transfer differences persist after matching cell-type support or
donor counts.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from .evaluation import summarize_classification
from .representations import get_representation
from .splits import restrict_to_train_labels

# -----------------------------------------------------------------------------
# PBMC site-transfer support diagnostics
# -----------------------------------------------------------------------------


def group_celltype_support(
    obs: pd.DataFrame,
    *,
    group_col: str,
    celltype_col: str,
    donor_col: str,
) -> pd.DataFrame:
    """Cell and donor support for every group x cell-type combination."""
    work = obs[[group_col, celltype_col, donor_col]].dropna().copy()
    for col in [group_col, celltype_col, donor_col]:
        work[col] = work[col].astype(str)

    cells = (
        work.groupby([group_col, celltype_col], observed=True)
        .size()
        .rename("n_cells")
    )
    donors = (
        work.groupby([group_col, celltype_col], observed=True)[donor_col]
        .nunique()
        .rename("n_donors")
    )
    return pd.concat([cells, donors], axis=1).reset_index()


def _fit_group_transfer(
    adata,
    *,
    rep_key: str,
    group_col: str,
    train_group: str,
    test_group: str,
    celltype_col: str,
    random_state: int,
    train_obs_names: Iterable[str] | None = None,
    allowed_donors: Mapping[str, Sequence[str]] | None = None,
) -> dict:
    """Fit one directional group-transfer model using the primary model recipe."""
    obs = adata.obs.copy()
    X_rep = get_representation(adata, rep_key)
    groups = obs[group_col].astype(str)
    train_mask = groups.eq(str(train_group))
    test_mask = groups.eq(str(test_group))

    if allowed_donors is not None:
        # The caller supplies group-specific donor lists.
        donor_col = allowed_donors.get("donor_col")
        if donor_col is None:
            raise ValueError("allowed_donors must contain key 'donor_col'.")
        donor_series = obs[str(donor_col)].astype(str)
        train_mask &= donor_series.isin(
            pd.Series(allowed_donors[str(train_group)]).astype(str)
        )
        test_mask &= donor_series.isin(
            pd.Series(allowed_donors[str(test_group)]).astype(str)
        )

    if train_obs_names is not None:
        train_obs_names = set(map(str, train_obs_names))
        train_mask &= pd.Index(obs.index.astype(str)).isin(train_obs_names)

    X_train = X_rep[train_mask.to_numpy()]
    X_test = X_rep[test_mask.to_numpy()]
    y_train = obs.loc[train_mask, celltype_col].astype(str).to_numpy()
    y_test = obs.loc[test_mask, celltype_col].astype(str).to_numpy()
    obs_test = obs.loc[test_mask].copy()

    train_counts = pd.Series(y_train).value_counts()
    keep_train_labels = train_counts[train_counts >= 2].index
    keep_train = pd.Series(y_train).isin(keep_train_labels).to_numpy()
    X_train = X_train[keep_train]
    y_train = y_train[keep_train]

    y_test, X_test, obs_test = restrict_to_train_labels(
        y_train=y_train,
        y_test=y_test,
        X_test=X_test,
        obs_test=obs_test,
    )
    if len(y_test) == 0:
        raise ValueError("No valid test cells remain after restricting to training labels.")

    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)
    clf = LogisticRegression(max_iter=5000, random_state=random_state)
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    labels = np.unique(np.concatenate([y_train, y_test]))
    summary = summarize_classification(y_test, y_pred, labels=labels)
    return {
        "macro_f1": float(summary["macro_f1"]),
        "accuracy": float(summary["accuracy"]),
        "per_class_f1": summary["per_class_f1"],
        "n_train_cells": int(len(y_train)),
        "n_test_cells": int(len(y_test)),
        "n_classes": int(len(labels)),
    }


def _sample_balanced_training_names(
    obs: pd.DataFrame,
    *,
    group_col: str,
    celltype_col: str,
    group_a: str,
    group_b: str,
    random_state: int,
) -> tuple[list[str], list[str], pd.DataFrame]:
    """Match training cell counts by cell type across two groups."""
    work = obs[[group_col, celltype_col]].copy()
    work[group_col] = work[group_col].astype(str)
    work[celltype_col] = work[celltype_col].astype(str)
    rng = np.random.default_rng(random_state)

    labels = sorted(
        set(work.loc[work[group_col] == str(group_a), celltype_col])
        & set(work.loc[work[group_col] == str(group_b), celltype_col])
    )
    names_a: list[str] = []
    names_b: list[str] = []
    rows = []
    for label in labels:
        idx_a = work.index[
            (work[group_col] == str(group_a)) & (work[celltype_col] == label)
        ].astype(str).to_numpy()
        idx_b = work.index[
            (work[group_col] == str(group_b)) & (work[celltype_col] == label)
        ].astype(str).to_numpy()
        n = min(len(idx_a), len(idx_b))
        if n < 2:
            continue
        chosen_a = rng.choice(idx_a, size=n, replace=False)
        chosen_b = rng.choice(idx_b, size=n, replace=False)
        names_a.extend(chosen_a.tolist())
        names_b.extend(chosen_b.tolist())
        rows.append(
            {
                "cell_type": label,
                f"original_{group_a}": len(idx_a),
                f"original_{group_b}": len(idx_b),
                "matched_per_group": n,
            }
        )
    return names_a, names_b, pd.DataFrame(rows)


def run_support_matched_group_transfer(
    adata,
    representations: Mapping[str, str],
    *,
    group_col: str,
    group_a: str,
    group_b: str,
    celltype_col: str,
    n_repeats: int = 10,
    seed: int = 5000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Repeat bidirectional transfer with matched training support by class."""
    rows = []
    support_parts = []
    for repeat in range(n_repeats):
        state = seed + repeat
        names_a, names_b, support = _sample_balanced_training_names(
            adata.obs,
            group_col=group_col,
            celltype_col=celltype_col,
            group_a=group_a,
            group_b=group_b,
            random_state=state,
        )
        support.insert(0, "repeat", repeat)
        support_parts.append(support)

        for rep_label, rep_key in representations.items():
            for train_group, test_group, names in [
                (group_a, group_b, names_a),
                (group_b, group_a, names_b),
            ]:
                result = _fit_group_transfer(
                    adata,
                    rep_key=rep_key,
                    group_col=group_col,
                    train_group=train_group,
                    test_group=test_group,
                    celltype_col=celltype_col,
                    random_state=state,
                    train_obs_names=names,
                )
                rows.append(
                    {
                        "analysis": "support_matched",
                        "repeat": repeat,
                        "representation": rep_label,
                        "rep_key": rep_key,
                        "train_group": str(train_group),
                        "test_group": str(test_group),
                        "direction": f"{train_group}_to_{test_group}",
                        **{k: v for k, v in result.items() if k != "per_class_f1"},
                    }
                )
    return pd.DataFrame(rows), pd.concat(support_parts, ignore_index=True)


def run_donor_count_matched_group_transfer(
    adata,
    representations: Mapping[str, str],
    *,
    group_col: str,
    group_a: str,
    group_b: str,
    celltype_col: str,
    donor_col: str,
    n_repeats: int = 10,
    seed: int = 6000,
) -> pd.DataFrame:
    """Repeat bidirectional transfer after matching the number of donors/site."""
    obs = adata.obs.copy()
    groups = obs[group_col].astype(str)
    donors = obs[donor_col].astype(str)
    donors_a = np.array(sorted(donors[groups == str(group_a)].unique()))
    donors_b = np.array(sorted(donors[groups == str(group_b)].unique()))
    n = min(len(donors_a), len(donors_b))
    rng = np.random.default_rng(seed)

    rows = []
    for repeat in range(n_repeats):
        chosen_a = np.sort(rng.choice(donors_a, n, replace=False))
        chosen_b = np.sort(rng.choice(donors_b, n, replace=False))
        allowed = {
            "donor_col": donor_col,
            str(group_a): chosen_a,
            str(group_b): chosen_b,
        }
        for rep_label, rep_key in representations.items():
            for train_group, test_group in [
                (group_a, group_b),
                (group_b, group_a),
            ]:
                result = _fit_group_transfer(
                    adata,
                    rep_key=rep_key,
                    group_col=group_col,
                    train_group=train_group,
                    test_group=test_group,
                    celltype_col=celltype_col,
                    random_state=seed + repeat,
                    allowed_donors=allowed,
                )
                rows.append(
                    {
                        "analysis": "donor_count_matched",
                        "repeat": repeat,
                        "representation": rep_label,
                        "rep_key": rep_key,
                        "train_group": str(train_group),
                        "test_group": str(test_group),
                        "direction": f"{train_group}_to_{test_group}",
                        "n_donors_per_group": n,
                        "sampled_group_a_donors": ";".join(chosen_a.tolist()),
                        "sampled_group_b_donors": ";".join(chosen_b.tolist()),
                        **{k: v for k, v in result.items() if k != "per_class_f1"},
                    }
                )
    return pd.DataFrame(rows)


def summarize_transfer_sensitivity(raw: pd.DataFrame) -> pd.DataFrame:
    """Summarize repeated transfer sensitivity runs."""
    group_cols = [
        c
        for c in ["analysis", "representation", "train_group", "test_group", "direction"]
        if c in raw.columns
    ]
    return (
        raw.groupby(group_cols, as_index=False)
        .agg(
            macro_f1_mean=("macro_f1", "mean"),
            macro_f1_std=("macro_f1", "std"),
            accuracy_mean=("accuracy", "mean"),
            accuracy_std=("accuracy", "std"),
            n_runs=("macro_f1", "size"),
        )
    )