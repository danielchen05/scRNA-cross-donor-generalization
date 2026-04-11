import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    f1_score,
    confusion_matrix,
    accuracy_score,
    classification_report,
)


CELLTYPE_COL = "cell_type"
DONOR_COL = "patient_id"
BATCH_COL = "Site"
RANDOM_STATE = 42

REPRESENTATIONS = {
    "hvg": "hvg",
    "pca": "X_pca",
    "harmony": "X_harmony",
    "scvi": "X_scVI",
}


# get feature matrix for a chosen representation
def get_representation(adata, rep_name):
    if rep_name == "hvg":
        X = adata.X
    else:
        X = adata.obsm[rep_name]

    if hasattr(X, "toarray"):
        X = X.toarray()
    else:
        X = np.asarray(X)

    return X


# one-hot encode batch/site covariate
def get_batch_matrix(obs, batch_col):
    batch_df = pd.get_dummies(
        obs[batch_col].astype(str),
        prefix=batch_col,
        drop_first=False,
    )
    return batch_df.to_numpy(), batch_df.columns.tolist()


# remove test labels absent from training set
def restrict_to_train_labels(y_train, y_test, X_test, obs_test):
    train_labels = set(pd.Series(y_train).astype(str).unique())
    keep_mask = pd.Series(y_test).astype(str).isin(train_labels).to_numpy()

    y_test = np.asarray(y_test)[keep_mask]
    X_test = X_test[keep_mask]
    obs_test = obs_test.loc[keep_mask].copy()

    return y_test, X_test, obs_test


# summarize multiclass classification results
def summarize_classification(y_true, y_pred, labels):
    cm = confusion_matrix(y_true, y_pred, labels=labels)
    macro_f1 = f1_score(y_true, y_pred, average="macro")
    acc = accuracy_score(y_true, y_pred)

    report = classification_report(
        y_true,
        y_pred,
        labels=labels,
        output_dict=True,
        zero_division=0,
    )

    per_class_f1 = pd.DataFrame({
        "cell_type": labels,
        "f1": [report[label]["f1-score"] for label in labels],
        "precision": [report[label]["precision"] for label in labels],
        "recall": [report[label]["recall"] for label in labels],
        "support": [report[label]["support"] for label in labels],
    })

    return {
        "macro_f1": macro_f1,
        "accuracy": acc,
        "cm": cm,
        "per_class_f1": per_class_f1,
    }


# train on one donor set, test on another
def run_donor_split_logreg(
    adata,
    rep_name,
    train_donors,
    test_donors,
    celltype_col=CELLTYPE_COL,
    donor_col=DONOR_COL,
    batch_col=None,
    random_state=RANDOM_STATE,
):
    obs = adata.obs.copy()
    X_rep = get_representation(adata, rep_name)

    train_mask = obs[donor_col].astype(str).isin(pd.Series(train_donors).astype(str))
    test_mask = obs[donor_col].astype(str).isin(pd.Series(test_donors).astype(str))

    X_train_rep = X_rep[train_mask.to_numpy()]
    X_test_rep = X_rep[test_mask.to_numpy()]

    y_train = obs.loc[train_mask, celltype_col].astype(str).to_numpy()
    y_test = obs.loc[test_mask, celltype_col].astype(str).to_numpy()

    obs_train = obs.loc[train_mask].copy()
    obs_test = obs.loc[test_mask].copy()

    # remove rare labels from training
    train_counts = pd.Series(y_train).value_counts()
    keep_train_labels = train_counts[train_counts >= 2].index

    keep_train_mask = pd.Series(y_train).isin(keep_train_labels).to_numpy()
    X_train_rep = X_train_rep[keep_train_mask]
    y_train = y_train[keep_train_mask]
    obs_train = obs_train.loc[keep_train_mask].copy()

    # keep only test labels seen in training
    y_test, X_test_rep, obs_test = restrict_to_train_labels(
        y_train=y_train,
        y_test=y_test,
        X_test=X_test_rep,
        obs_test=obs_test,
    )

    if len(y_test) == 0:
        raise ValueError("No valid test cells remain after restricting to training labels.")

    # optional batch features
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

    # scale representation only
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
    summary = summarize_classification(y_test, y_pred, labels)

    return {
        "model": clf,
        "scaler": scaler,
        "train_donors": list(train_donors),
        "test_donors": list(test_donors),
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


# split donors into folds
def make_donor_folds(
    adata,
    donor_col=DONOR_COL,
    n_folds=5,
    random_state=RANDOM_STATE,
):
    donors = np.array(sorted(adata.obs[donor_col].astype(str).unique()))
    rng = np.random.default_rng(random_state)
    shuffled = donors.copy()
    rng.shuffle(shuffled)
    return np.array_split(shuffled, n_folds)


# run donor-held-out cross-validation for one representation
def run_donor_cv_experiment(
    adata,
    rep_name,
    donor_col=DONOR_COL,
    batch_col=None,
    n_folds=5,
    random_state=RANDOM_STATE,
):
    donor_folds = make_donor_folds(
        adata,
        donor_col=donor_col,
        n_folds=n_folds,
        random_state=random_state,
    )

    fold_metrics = []
    all_pred_rows = []
    per_class_tables = []

    cm_sum = None
    final_labels = None

    for fold_idx in range(len(donor_folds)):
        test_donors = donor_folds[fold_idx]
        train_donors = np.concatenate(
            [donor_folds[j] for j in range(len(donor_folds)) if j != fold_idx]
        )

        res = run_donor_split_logreg(
            adata=adata,
            rep_name=rep_name,
            train_donors=train_donors,
            test_donors=test_donors,
            donor_col=donor_col,
            batch_col=batch_col,
            random_state=random_state,
        )

        pred_df = pd.DataFrame({
            "fold": fold_idx + 1,
            "y_true": res["y_test"],
            "y_pred": res["y_pred"],
        })
        all_pred_rows.append(pred_df)

        fold_per_class = res["per_class_f1"].copy()
        fold_per_class["fold"] = fold_idx + 1
        per_class_tables.append(fold_per_class)

        fold_metrics.append({
            "fold": fold_idx + 1,
            "representation": rep_name,
            "macro_f1": res["macro_f1"],
            "accuracy": res["accuracy"],
            "n_train_cells": res["n_train_cells"],
            "n_test_cells": res["n_test_cells"],
            "n_classes_used": res["n_classes_used"],
            "n_train_donors": len(res["train_donors"]),
            "n_test_donors": len(res["test_donors"]),
            "batch_covariate": batch_col if batch_col is not None else "None",
            "n_batch_features": len(res["batch_feature_names"]),
        })

        if cm_sum is None:
            cm_sum = res["cm"].copy()
            final_labels = res["labels"]
        else:
            cm_sum = cm_sum + res["cm"]

    fold_metrics_df = pd.DataFrame(fold_metrics)
    all_preds_df = pd.concat(all_pred_rows, ignore_index=True)
    per_class_all_df = pd.concat(per_class_tables, ignore_index=True)

    per_class_mean_df = (
        per_class_all_df
        .groupby("cell_type", as_index=False)[["f1", "precision", "recall", "support"]]
        .mean()
        .sort_values("f1", ascending=False)
        .reset_index(drop=True)
    )

    summary_row = {
        "representation": rep_name,
        "macro_f1_mean": fold_metrics_df["macro_f1"].mean(),
        "macro_f1_std": fold_metrics_df["macro_f1"].std(ddof=1),
        "accuracy_mean": fold_metrics_df["accuracy"].mean(),
        "accuracy_std": fold_metrics_df["accuracy"].std(ddof=1),
        "n_folds": len(fold_metrics_df),
        "mean_train_cells": fold_metrics_df["n_train_cells"].mean(),
        "mean_test_cells": fold_metrics_df["n_test_cells"].mean(),
        "mean_n_classes_used": fold_metrics_df["n_classes_used"].mean(),
        "batch_covariate": batch_col if batch_col is not None else "None",
        "mean_n_batch_features": fold_metrics_df["n_batch_features"].mean(),
    }

    return {
        "summary_row": summary_row,
        "fold_metrics_df": fold_metrics_df,
        "all_predictions_df": all_preds_df,
        "per_class_all_df": per_class_all_df,
        "per_class_mean_df": per_class_mean_df,
        "cm_sum": cm_sum,
        "labels": final_labels,
    }


# sample k training donors for donor ablation
def sample_train_donors(all_donors, k, random_state=42):
    rng = np.random.default_rng(random_state)
    return rng.choice(all_donors, size=k, replace=False)


# run donor ablation across k values and repeats
def run_donor_ablation_experiment(
    adata,
    rep_name,
    donor_col=DONOR_COL,
    k_values=(3, 5, 8, 10, 15, 20),
    n_repeats=5,
    random_state=RANDOM_STATE,
):
    all_donors = np.array(sorted(adata.obs[donor_col].astype(str).unique()))
    results = []

    for k in k_values:
        for repeat in range(n_repeats):
            seed = random_state + 100 * k + repeat

            train_donors = sample_train_donors(
                all_donors,
                k,
                random_state=seed,
            )
            test_donors = np.array([d for d in all_donors if d not in train_donors])

            res = run_donor_split_logreg(
                adata=adata,
                rep_name=rep_name,
                train_donors=train_donors,
                test_donors=test_donors,
                donor_col=donor_col,
                batch_col=None,
                random_state=random_state,
            )

            results.append({
                "representation": rep_name,
                "k_train_donors": k,
                "repeat": repeat + 1,
                "macro_f1": res["macro_f1"],
                "accuracy": res["accuracy"],
                "n_train_cells": res["n_train_cells"],
                "n_test_cells": res["n_test_cells"],
                "n_classes_used": res["n_classes_used"],
            })

    results_df = pd.DataFrame(results)

    summary_df = (
        results_df
        .groupby(["representation", "k_train_donors"], as_index=False)
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

    return {
        "results_df": results_df,
        "summary_df": summary_df,
    }


# run cross-site train/test in one direction
def run_cross_site_experiment(
    adata,
    rep_name,
    train_site_name,
    test_site_name,
    train_donors,
    test_donors,
):
    out = run_donor_split_logreg(
        adata=adata,
        rep_name=rep_name,
        train_donors=train_donors,
        test_donors=test_donors,
        donor_col=DONOR_COL,
        batch_col=None,
        random_state=RANDOM_STATE,
    )

    return {
        "representation": rep_name,
        "train_site": train_site_name,
        "test_site": test_site_name,
        "direction": f"{train_site_name} → {test_site_name}",
        "macro_f1": out["macro_f1"],
        "accuracy": out["accuracy"],
        "n_train_donors": len(train_donors),
        "n_test_donors": len(test_donors),
        "n_train_cells": out["n_train_cells"],
        "n_test_cells": out["n_test_cells"],
        "n_classes_used": out["n_classes_used"],
        "per_class_f1": out["per_class_f1"],
        "cm": out["cm"],
        "labels": out["labels"],
        "y_test": out["y_test"],
        "y_pred": out["y_pred"],
    }