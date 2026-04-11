import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    f1_score,
    confusion_matrix,
    accuracy_score,
    classification_report,
)


CELLTYPE_COL = "cell_type"
RANDOM_STATE = 42
TEST_SIZE = 0.2


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


# run random cell-level split logistic regression
def run_random_split_logreg(
    adata,
    rep_name,
    celltype_col=CELLTYPE_COL,
    batch_col=None,
    test_size=TEST_SIZE,
    random_state=RANDOM_STATE,
):
    X_rep = get_representation(adata, rep_name)
    y = adata.obs[celltype_col].astype(str)

    # remove labels with fewer than 2 cells
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

    if X_batch is not None:
        X_train_rep, X_test_rep, X_train_batch, X_test_batch, y_train, y_test = train_test_split(
            X_rep,
            X_batch,
            y,
            test_size=test_size,
            random_state=random_state,
            stratify=y,
        )
    else:
        X_train_rep, X_test_rep, y_train, y_test = train_test_split(
            X_rep,
            y,
            test_size=test_size,
            random_state=random_state,
            stratify=y,
        )

    # scale representation only
    scaler = StandardScaler()
    X_train_rep = scaler.fit_transform(X_train_rep)
    X_test_rep = scaler.transform(X_test_rep)

    if X_batch is not None:
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
    summary = summarize_classification(y_test, y_pred, labels)

    return {
        "model": clf,
        "scaler": scaler,
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


# run random split for multiple representations
def run_random_split_experiment(
    adata,
    representations,
    celltype_col=CELLTYPE_COL,
    batch_col=None,
    test_size=TEST_SIZE,
    random_state=RANDOM_STATE,
):
    metrics_rows = []
    results = {}

    for rep_label, rep_key in representations.items():
        res = run_random_split_logreg(
            adata=adata,
            rep_name=rep_key,
            celltype_col=celltype_col,
            batch_col=batch_col,
            test_size=test_size,
            random_state=random_state,
        )

        results[rep_label] = res

        metrics_rows.append({
            "representation": rep_label,
            "batch_covariate": batch_col if batch_col is not None else "None",
            "macro_f1": res["macro_f1"],
            "accuracy": res["accuracy"],
            "n_cells_used": res["n_cells_used"],
            "n_classes_used": res["n_classes_used"],
            "n_batch_features": len(res["batch_feature_names"]),
        })

    metrics_df = pd.DataFrame(metrics_rows).sort_values(
        "macro_f1",
        ascending=False,
    ).reset_index(drop=True)

    return {
        "metrics_df": metrics_df,
        "results": results,
    }