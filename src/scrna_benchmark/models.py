# src/scrna_benchmark/models.py

import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from .representations import get_representation
from .evaluation import summarize_classification
from .splits import restrict_to_train_labels


def get_batch_matrix(obs, batch_col):
    """
    One-hot encode a batch/site covariate.
    """
    batch_df = pd.get_dummies(
        obs[batch_col].astype(str),
        prefix=batch_col,
        drop_first=False,
    )
    return batch_df.to_numpy(), batch_df.columns.tolist()


def _fit_predict_logreg(X_train, X_test, y_train, random_state=42):
    """
    Scale features, fit logistic regression, and predict.

    This helper is currently unused by the main runners because the main
    runners need to scale representation features separately from optional
    batch covariates.
    """
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    clf = LogisticRegression(
        max_iter=5000,
        random_state=random_state,
        n_jobs=-1,
    )
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)

    return clf, scaler, y_pred


def run_random_split_logreg(
    adata,
    rep_key,
    celltype_col,
    batch_col=None,
    test_size=0.2,
    random_state=42,
):
    """
    Run random cell-level split logistic regression.

    Returns
    -------
    dict
        Model, scaler, predictions, metrics, and train/test metadata.
    """
    X_rep = get_representation(adata, rep_key)
    y = adata.obs[celltype_col].astype(str)

    counts = y.value_counts()
    keep_labels = counts[counts >= 2].index
    keep_mask = y.isin(keep_labels)

    X_rep = X_rep[keep_mask.to_numpy()]
    y = y[keep_mask].to_numpy()
    obs_sub = adata.obs.loc[keep_mask].copy()

    if batch_col is not None:
        X_batch, batch_feature_names = get_batch_matrix(obs_sub, batch_col)
    else:
        X_batch = None
        batch_feature_names = []

    idx = np.arange(len(y))

    idx_train, idx_test = train_test_split(
        idx,
        test_size=test_size,
        random_state=random_state,
        stratify=y,
    )

    X_train_rep = X_rep[idx_train]
    X_test_rep = X_rep[idx_test]

    y_train = y[idx_train]
    y_test = y[idx_test]

    obs_train = obs_sub.iloc[idx_train].copy()
    obs_test = obs_sub.iloc[idx_test].copy()

    if X_batch is not None:
        X_train_batch = X_batch[idx_train]
        X_test_batch = X_batch[idx_test]
    else:
        X_train_batch = None
        X_test_batch = None

    scaler = StandardScaler()
    X_train_rep = scaler.fit_transform(X_train_rep)
    X_test_rep = scaler.transform(X_test_rep)

    if X_train_batch is not None:
        X_train = np.hstack([X_train_rep, X_train_batch])
        X_test = np.hstack([X_test_rep, X_test_batch])
    else:
        X_train = X_train_rep
        X_test = X_test_rep

    clf = LogisticRegression(
        max_iter=5000,
        random_state=random_state,
        n_jobs=-1,
    )
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)

    labels = np.unique(y)
    summary = summarize_classification(y_test, y_pred, labels=labels)

    return {
        "model": clf,
        "scaler": scaler,
        "obs_train": obs_train,
        "obs_test": obs_test,
        "y_train": y_train,
        "y_test": y_test,
        "y_pred": y_pred,
        "labels": labels,
        "macro_f1": summary["macro_f1"],
        "accuracy": summary["accuracy"],
        "cm": summary["cm"],
        "per_class_f1": summary["per_class_f1"],
        "n_cells_used": len(y),
        "n_classes_used": len(labels),
        "batch_feature_names": batch_feature_names,
    }


def run_donor_split_logreg(
    adata,
    rep_key,
    train_donors,
    test_donors,
    celltype_col,
    donor_col,
    batch_col=None,
    random_state=42,
):
    """
    Train on one donor set and test on a disjoint donor set.

    Returns
    -------
    dict
        Model, scaler, predictions, metrics, donor IDs, and train/test metadata.
    """
    obs = adata.obs.copy()
    X_rep = get_representation(adata, rep_key)

    train_mask = obs[donor_col].astype(str).isin(pd.Series(train_donors).astype(str))
    test_mask = obs[donor_col].astype(str).isin(pd.Series(test_donors).astype(str))

    X_train_rep = X_rep[train_mask.to_numpy()]
    X_test_rep = X_rep[test_mask.to_numpy()]

    y_train = obs.loc[train_mask, celltype_col].astype(str).to_numpy()
    y_test = obs.loc[test_mask, celltype_col].astype(str).to_numpy()

    obs_train = obs.loc[train_mask].copy()
    obs_test = obs.loc[test_mask].copy()

    train_counts = pd.Series(y_train).value_counts()
    keep_train_labels = train_counts[train_counts >= 2].index

    keep_train_mask = pd.Series(y_train).isin(keep_train_labels).to_numpy()

    X_train_rep = X_train_rep[keep_train_mask]
    y_train = y_train[keep_train_mask]
    obs_train = obs_train.iloc[keep_train_mask].copy()

    y_test, X_test_rep, obs_test = restrict_to_train_labels(
        y_train=y_train,
        y_test=y_test,
        X_test=X_test_rep,
        obs_test=obs_test,
    )

    if len(y_test) == 0:
        raise ValueError("No valid test cells remain after restricting to training labels.")

    if batch_col is not None:
        train_batch_df = pd.get_dummies(
            obs_train[batch_col].astype(str),
            prefix=batch_col,
            drop_first=False,
        )
        test_batch_df = pd.get_dummies(
            obs_test[batch_col].astype(str),
            prefix=batch_col,
            drop_first=False,
        )

        train_batch_df, test_batch_df = train_batch_df.align(
            test_batch_df,
            join="outer",
            axis=1,
            fill_value=0,
        )

        X_train_batch = train_batch_df.to_numpy()
        X_test_batch = test_batch_df.to_numpy()
        batch_feature_names = train_batch_df.columns.tolist()
    else:
        X_train_batch = None
        X_test_batch = None
        batch_feature_names = []

    scaler = StandardScaler()
    X_train_rep = scaler.fit_transform(X_train_rep)
    X_test_rep = scaler.transform(X_test_rep)

    if X_train_batch is not None:
        X_train = np.hstack([X_train_rep, X_train_batch])
        X_test = np.hstack([X_test_rep, X_test_batch])
    else:
        X_train = X_train_rep
        X_test = X_test_rep

    clf = LogisticRegression(
        max_iter=5000,
        random_state=random_state,
        n_jobs=-1,
    )
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)

    labels = np.unique(np.concatenate([y_train, y_test]))
    summary = summarize_classification(y_test, y_pred, labels=labels)

    return {
        "model": clf,
        "scaler": scaler,
        "obs_train": obs_train,
        "obs_test": obs_test,
        "train_donors": list(train_donors),
        "test_donors": list(test_donors),
        "y_train": y_train,
        "y_test": y_test,
        "y_pred": y_pred,
        "labels": labels,
        "macro_f1": summary["macro_f1"],
        "accuracy": summary["accuracy"],
        "cm": summary["cm"],
        "per_class_f1": summary["per_class_f1"],
        "n_train_cells": len(y_train),
        "n_test_cells": len(y_test),
        "n_classes_used": len(labels),
        "batch_feature_names": batch_feature_names,
    }