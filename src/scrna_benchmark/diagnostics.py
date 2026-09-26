"""Focused diagnostic helpers for the cross-donor generalization project.

These functions are intentionally analysis-oriented rather than part of the
primary benchmark API. They support:

* Blood Atlas metadata-alignment and donor-composition audits
* reduced-donor and label-permutation sensitivity analyses
* PBMC directional-transfer support diagnostics

They reuse the project's existing logistic-regression runners so the
sensitivity analyses stay aligned with the primary benchmark implementation.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from .evaluation import summarize_classification
from .models import run_donor_split_logreg, run_random_split_logreg
from .representations import get_representation
from .splits import make_donor_folds, restrict_to_train_labels


# -----------------------------------------------------------------------------
# Metadata alignment
# -----------------------------------------------------------------------------


def audit_metadata_alignment(
    obs_names: Sequence[str],
    metadata: pd.DataFrame,
    id_col: str,
) -> pd.DataFrame:
    """Audit whether external metadata IDs match AnnData observation names.

    The function deliberately separates *set identity* from *row-order
    identity*. A metadata table can contain exactly the correct barcodes while
    still being unsafe to assign by row order.
    """
    obs = pd.Index(pd.Series(obs_names, dtype="string").astype(str))
    if id_col not in metadata.columns:
        raise KeyError(f"Metadata ID column not found: {id_col!r}")

    meta_ids = pd.Index(metadata[id_col].astype(str))
    same_length = len(obs) == len(meta_ids)
    order_match_fraction = np.nan
    if same_length and len(obs) > 0:
        order_match_fraction = float(np.mean(obs.to_numpy() == meta_ids.to_numpy()))

    obs_set = set(obs)
    meta_set = set(meta_ids)
    intersection = obs_set & meta_set
    union = obs_set | meta_set

    row = {
        "n_obs": len(obs),
        "n_metadata_rows": len(meta_ids),
        "obs_names_unique": bool(obs.is_unique),
        "metadata_ids_unique": bool(meta_ids.is_unique),
        "same_length": bool(same_length),
        "set_equal": bool(obs_set == meta_set),
        "n_intersection": len(intersection),
        "intersection_fraction_obs": (
            len(intersection) / len(obs_set) if obs_set else np.nan
        ),
        "intersection_fraction_metadata": (
            len(intersection) / len(meta_set) if meta_set else np.nan
        ),
        "jaccard_id_overlap": (
            len(intersection) / len(union) if union else np.nan
        ),
        "order_match_fraction": order_match_fraction,
        "safe_for_row_order_assignment": bool(
            same_length
            and obs.is_unique
            and meta_ids.is_unique
            and order_match_fraction == 1.0
        ),
        "safe_for_id_reindex": bool(
            obs.is_unique and meta_ids.is_unique and obs_set == meta_set
        ),
    }
    return pd.DataFrame([row])


def align_metadata_to_obs(
    obs_names: Sequence[str],
    metadata: pd.DataFrame,
    id_col: str,
    *,
    require_exact_set: bool = True,
) -> pd.DataFrame:
    """Reindex metadata to AnnData observation order using explicit IDs.

    This is the safe replacement for positional metadata assignment.
    """
    obs = pd.Index(pd.Series(obs_names, dtype="string").astype(str))
    if id_col not in metadata.columns:
        raise KeyError(f"Metadata ID column not found: {id_col!r}")
    if not obs.is_unique:
        raise ValueError("AnnData observation names are not unique.")

    meta = metadata.copy()
    meta[id_col] = meta[id_col].astype(str)
    if meta[id_col].duplicated().any():
        examples = meta.loc[meta[id_col].duplicated(), id_col].head().tolist()
        raise ValueError(f"Metadata IDs are not unique. Examples: {examples}")

    meta = meta.set_index(id_col, drop=False)
    if require_exact_set and set(meta.index) != set(obs):
        missing = obs.difference(meta.index)
        extra = meta.index.difference(obs)
        raise ValueError(
            "Metadata and AnnData ID sets differ. "
            f"Missing metadata for {len(missing)} observations; "
            f"metadata has {len(extra)} extra IDs."
        )

    aligned = meta.reindex(obs)
    if aligned[id_col].isna().any():
        n_missing = int(aligned[id_col].isna().sum())
        raise ValueError(f"Alignment left {n_missing} observations without metadata.")

    aligned.index = obs
    return aligned


# -----------------------------------------------------------------------------
# Donor x cell-type composition
# -----------------------------------------------------------------------------


def donor_celltype_counts(
    obs: pd.DataFrame,
    donor_col: str,
    celltype_col: str,
    *,
    labels: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Return a donor x cell-type count matrix."""
    work = obs[[donor_col, celltype_col]].dropna().copy()
    work[donor_col] = work[donor_col].astype(str)
    work[celltype_col] = work[celltype_col].astype(str)
    table = pd.crosstab(work[donor_col], work[celltype_col])
    if labels is not None:
        table = table.reindex(columns=list(labels), fill_value=0)
    return table.sort_index()


def donor_celltype_proportions(counts: pd.DataFrame) -> pd.DataFrame:
    """Row-normalize a donor x cell-type count matrix."""
    denom = counts.sum(axis=1).replace(0, np.nan)
    return counts.div(denom, axis=0).fillna(0.0)


def _entropy_base2(p: np.ndarray) -> float:
    p = np.asarray(p, dtype=float)
    p = p[p > 0]
    if p.size == 0:
        return np.nan
    return float(-(p * np.log2(p)).sum())


def _js_divergence_base2(p: np.ndarray, q: np.ndarray) -> float:
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    if p.sum() <= 0 or q.sum() <= 0:
        return np.nan
    p = p / p.sum()
    q = q / q.sum()
    m = 0.5 * (p + q)

    def kl(a: np.ndarray, b: np.ndarray) -> float:
        mask = a > 0
        return float(np.sum(a[mask] * np.log2(a[mask] / b[mask])))

    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


def summarize_donor_composition(
    counts: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Summarize donor-specific composition heterogeneity.

    Returns
    -------
    donor_metrics
        One row per donor with total cells, entropy, normalized entropy and
        Jensen-Shannon divergence from the cohort-wide composition.
    celltype_metrics
        One row per cell type with cross-donor proportion summaries.
    """
    props = donor_celltype_proportions(counts)
    global_p = counts.sum(axis=0).to_numpy(dtype=float)
    global_p = global_p / global_p.sum()
    k = max(1, counts.shape[1])
    max_entropy = np.log2(k) if k > 1 else 1.0

    donor_rows = []
    for donor, row in props.iterrows():
        p = row.to_numpy(dtype=float)
        ent = _entropy_base2(p)
        donor_rows.append(
            {
                "donor": str(donor),
                "n_cells": int(counts.loc[donor].sum()),
                "entropy_bits": ent,
                "entropy_normalized": (
                    ent / max_entropy if np.isfinite(ent) and max_entropy > 0 else np.nan
                ),
                "js_divergence_from_global": _js_divergence_base2(p, global_p),
            }
        )
    donor_metrics = pd.DataFrame(donor_rows)

    cell_rows = []
    for cell_type in props.columns:
        s = props[cell_type]
        q1 = s.quantile(0.25)
        q3 = s.quantile(0.75)
        cell_rows.append(
            {
                "cell_type": str(cell_type),
                "mean_proportion": float(s.mean()),
                "sd_proportion": float(s.std(ddof=1)),
                "median_proportion": float(s.median()),
                "iqr_proportion": float(q3 - q1),
                "min_proportion": float(s.min()),
                "max_proportion": float(s.max()),
                "fraction_donors_present": float((counts[cell_type] > 0).mean()),
            }
        )
    celltype_metrics = pd.DataFrame(cell_rows)
    return donor_metrics, celltype_metrics


def composition_summary_row(
    counts: pd.DataFrame,
    *,
    dataset: str | None = None,
) -> dict:
    """Return compact dataset-level composition heterogeneity metrics."""
    donor_metrics, celltype_metrics = summarize_donor_composition(counts)
    return {
        "dataset": dataset,
        "n_donors": int(counts.shape[0]),
        "n_cell_types": int(counts.shape[1]),
        "n_cells": int(counts.to_numpy().sum()),
        "mean_js_divergence": float(
            donor_metrics["js_divergence_from_global"].mean()
        ),
        "median_js_divergence": float(
            donor_metrics["js_divergence_from_global"].median()
        ),
        "mean_normalized_entropy": float(
            donor_metrics["entropy_normalized"].mean()
        ),
        "mean_celltype_sd": float(celltype_metrics["sd_proportion"].mean()),
    }


# -----------------------------------------------------------------------------
# Core random-vs-donor sensitivity helper
# -----------------------------------------------------------------------------


def evaluate_random_vs_donor(
    adata,
    representations: Mapping[str, str],
    *,
    celltype_col: str,
    donor_col: str,
    random_state: int = 42,
    n_folds: int = 5,
    test_size: float = 0.2,
    batch_col: str | None = None,
) -> pd.DataFrame:
    """Evaluate random split and donor-held-out CV on one AnnData object.

    The donor score is the mean across folds. This compact paired output is
    convenient for repeated cohort/permutation sensitivity experiments.
    """
    folds = make_donor_folds(
        adata,
        donor_col=donor_col,
        n_folds=n_folds,
        random_state=random_state,
    )
    all_donors = set(adata.obs[donor_col].astype(str).unique())

    rows = []
    for rep_label, rep_key in representations.items():
        random_result = run_random_split_logreg(
            adata=adata,
            rep_key=rep_key,
            celltype_col=celltype_col,
            batch_col=batch_col,
            test_size=test_size,
            random_state=random_state,
        )

        donor_scores = []
        donor_accuracies = []
        for fold_id, test_donors in enumerate(folds):
            test_donors = [str(x) for x in test_donors]
            train_donors = sorted(all_donors - set(test_donors))
            result = run_donor_split_logreg(
                adata=adata,
                rep_key=rep_key,
                train_donors=train_donors,
                test_donors=test_donors,
                celltype_col=celltype_col,
                donor_col=donor_col,
                batch_col=batch_col,
                random_state=random_state + fold_id,
            )
            donor_scores.append(float(result["macro_f1"]))
            donor_accuracies.append(float(result["accuracy"]))

        donor_scores_arr = np.asarray(donor_scores, dtype=float)
        donor_acc_arr = np.asarray(donor_accuracies, dtype=float)
        rows.append(
            {
                "representation": str(rep_label),
                "rep_key": str(rep_key),
                "random_state": int(random_state),
                "random_macro_f1": float(random_result["macro_f1"]),
                "random_accuracy": float(random_result["accuracy"]),
                "donor_macro_f1_mean": float(donor_scores_arr.mean()),
                "donor_macro_f1_std": float(donor_scores_arr.std(ddof=1)),
                "donor_accuracy_mean": float(donor_acc_arr.mean()),
                "donor_accuracy_std": float(donor_acc_arr.std(ddof=1)),
                "delta_macro_f1_random_minus_donor": float(
                    random_result["macro_f1"] - donor_scores_arr.mean()
                ),
                "n_cells": int(adata.n_obs),
                "n_donors": int(adata.obs[donor_col].astype(str).nunique()),
                "n_cell_types": int(adata.obs[celltype_col].astype(str).nunique()),
                "n_folds": int(n_folds),
            }
        )
    return pd.DataFrame(rows)


def run_reduced_donor_cohorts(
    adata,
    representations: Mapping[str, str],
    *,
    celltype_col: str,
    donor_col: str,
    n_donors: int = 20,
    n_cohorts: int = 10,
    n_folds: int = 5,
    seed: int = 2026,
) -> pd.DataFrame:
    """Repeatedly sample reduced donor cohorts and rerun both schemes."""
    donors = np.array(sorted(adata.obs[donor_col].astype(str).unique()))
    if n_donors > len(donors):
        raise ValueError(f"Requested {n_donors} donors but only {len(donors)} exist.")

    rng = np.random.default_rng(seed)
    parts = []
    for cohort_id in range(n_cohorts):
        chosen = np.sort(rng.choice(donors, size=n_donors, replace=False))
        mask = adata.obs[donor_col].astype(str).isin(chosen).to_numpy()
        sub = adata[mask].copy()
        cohort_seed = int(seed + cohort_id)
        result = evaluate_random_vs_donor(
            sub,
            representations,
            celltype_col=celltype_col,
            donor_col=donor_col,
            random_state=cohort_seed,
            n_folds=min(n_folds, n_donors),
        )
        result.insert(0, "cohort_id", cohort_id)
        result["sampled_donors"] = ";".join(chosen.tolist())
        parts.append(result)
    return pd.concat(parts, ignore_index=True)


# -----------------------------------------------------------------------------
# Label-permutation and label-coarsening controls
# -----------------------------------------------------------------------------


def permute_labels(
    obs: pd.DataFrame,
    *,
    label_col: str,
    donor_col: str,
    mode: str,
    random_state: int,
) -> pd.Series:
    """Return globally or within-donor permuted labels aligned to obs.index."""
    labels = obs[label_col].astype(str).to_numpy(copy=True)
    rng = np.random.default_rng(random_state)

    if mode == "global":
        out = rng.permutation(labels)
    elif mode == "within_donor":
        out = labels.copy()
        donors = obs[donor_col].astype(str).to_numpy()
        for donor in pd.unique(donors):
            idx = np.flatnonzero(donors == donor)
            out[idx] = rng.permutation(out[idx])
    else:
        raise ValueError("mode must be 'global' or 'within_donor'.")

    return pd.Series(out, index=obs.index, name=f"{label_col}_{mode}_permuted")


def run_label_permutation_sensitivity(
    adata,
    representations: Mapping[str, str],
    *,
    celltype_col: str,
    donor_col: str,
    modes: Sequence[str] = ("global", "within_donor"),
    n_permutations: int = 10,
    n_folds: int = 5,
    seed: int = 3000,
) -> pd.DataFrame:
    """Evaluate random vs donor-held-out performance after label permutation."""
    scratch = "__diagnostic_permuted_label__"
    if scratch in adata.obs.columns:
        raise ValueError(f"Scratch column already exists: {scratch}")

    parts = []
    try:
        for mode in modes:
            for perm_id in range(n_permutations):
                random_state = int(seed + 1000 * list(modes).index(mode) + perm_id)
                adata.obs[scratch] = permute_labels(
                    adata.obs,
                    label_col=celltype_col,
                    donor_col=donor_col,
                    mode=mode,
                    random_state=random_state,
                )
                result = evaluate_random_vs_donor(
                    adata,
                    representations,
                    celltype_col=scratch,
                    donor_col=donor_col,
                    random_state=random_state,
                    n_folds=n_folds,
                )
                result.insert(0, "permutation_id", perm_id)
                result.insert(0, "permutation_mode", mode)
                parts.append(result)
    finally:
        if scratch in adata.obs.columns:
            del adata.obs[scratch]

    return pd.concat(parts, ignore_index=True)


def run_coarse_label_sensitivity(
    adata,
    representations: Mapping[str, str],
    *,
    celltype_col: str,
    donor_col: str,
    mapping: Mapping[str, str],
    n_repeats: int = 5,
    n_folds: int = 5,
    seed: int = 4000,
) -> pd.DataFrame:
    """Evaluate a coarsened label hierarchy without modifying the source file."""
    scratch = "__diagnostic_coarse_label__"
    if scratch in adata.obs.columns:
        raise ValueError(f"Scratch column already exists: {scratch}")

    original = adata.obs[celltype_col].astype(str)
    adata.obs[scratch] = original.map(lambda x: mapping.get(x, x)).astype(str)
    try:
        parts = []
        for repeat in range(n_repeats):
            part = evaluate_random_vs_donor(
                adata,
                representations,
                celltype_col=scratch,
                donor_col=donor_col,
                random_state=seed + repeat,
                n_folds=n_folds,
            )
            part.insert(0, "repeat", repeat)
            parts.append(part)
        return pd.concat(parts, ignore_index=True)
    finally:
        del adata.obs[scratch]


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


# -----------------------------------------------------------------------------
# Marker sanity helper
# -----------------------------------------------------------------------------


def marker_group_means(
    adata,
    markers: Sequence[str],
    *,
    group_col: str,
    use_raw: bool = True,
) -> tuple[pd.DataFrame, list[str]]:
    """Mean marker expression by group, using raw when available.

    This intentionally relies on var_names as gene symbols. It returns missing
    markers explicitly rather than silently dropping them.
    """
    source = adata.raw if use_raw and adata.raw is not None else adata
    var_names = pd.Index(source.var_names.astype(str))
    present = [m for m in markers if m in var_names]
    missing = [m for m in markers if m not in var_names]
    if not present:
        raise ValueError("None of the requested marker genes are present in var_names.")

    X = source[:, present].X
    groups = adata.obs[group_col].astype(str)
    rows = []
    for group in sorted(groups.unique()):
        mask = groups.eq(group).to_numpy()
        sub = X[mask]
        means = np.asarray(sub.mean(axis=0)).ravel()
        rows.append(pd.Series(means, index=present, name=group))
    return pd.DataFrame(rows), missing


def zscore_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Column-wise z scores, guarding constant columns."""
    centered = df - df.mean(axis=0)
    sd = df.std(axis=0, ddof=0).replace(0, np.nan)
    return centered.div(sd, axis=1).fillna(0.0)


# -----------------------------------------------------------------------------
# Small I/O helpers
# -----------------------------------------------------------------------------


def write_csv(df: pd.DataFrame, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return path
