# src/scrna_benchmark/experiments.py

from pathlib import Path

import numpy as np
import pandas as pd

from .evaluation import save_confusion_matrix, save_prediction_outputs
from .models import run_random_split_logreg, run_donor_split_logreg
from .splits import make_donor_folds


def run_random_split_experiment(
    adata,
    representations,
    results_dir,
    scheme_label,
    celltype_col,
    batch_col=None,
    test_size=0.2,
    random_state=42,
):
    """
    Run random split evaluation for all representations.
    """
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    rows = []

    for rep_label, rep_key in representations.items():
        print(f"Running {scheme_label}: {rep_label} ({rep_key})")

        res = run_random_split_logreg(
            adata=adata,
            rep_key=rep_key,
            celltype_col=celltype_col,
            batch_col=batch_col,
            test_size=test_size,
            random_state=random_state,
        )

        prefix = f"{scheme_label}_{rep_label}"
        save_prediction_outputs(res, results_dir, prefix)

        save_confusion_matrix(
            res["cm"],
            res["labels"],
            results_dir / f"{prefix}_confusion_matrix_normalized.png",
            title=f"{scheme_label} - {rep_label}",
            normalize=True,
        )

        rows.append({
            "scheme": scheme_label,
            "representation": rep_label,
            "macro_f1": res["macro_f1"],
            "accuracy": res["accuracy"],
            "n_cells_used": res["n_cells_used"],
            "n_classes_used": res["n_classes_used"],
            "batch_covariate": batch_col if batch_col is not None else "None",
            "n_batch_features": len(res["batch_feature_names"]),
        })

    metrics = pd.DataFrame(rows)
    metrics.to_csv(results_dir / "metrics.csv", index=False)

    return metrics


def run_donor_cv_experiment(
    adata,
    representations,
    results_dir,
    scheme_label,
    celltype_col,
    donor_col,
    batch_col=None,
    n_folds=5,
    random_state=42,
):
    """
    Run donor-held-out CV for all representations.
    """
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    donor_folds = make_donor_folds(
        adata,
        donor_col=donor_col,
        n_folds=n_folds,
        random_state=random_state,
    )

    summary_rows = []

    for rep_label, rep_key in representations.items():
        print(f"Running {scheme_label}: {rep_label} ({rep_key})")

        fold_rows = []
        all_pred_rows = []
        per_class_tables = []

        for fold_idx, test_donors in enumerate(donor_folds):
            train_donors = np.concatenate([
                donor_folds[j]
                for j in range(len(donor_folds))
                if j != fold_idx
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

            fold = fold_idx + 1
            prefix = f"{scheme_label}_{rep_label}_fold{fold}"

            save_prediction_outputs(res, results_dir, prefix)

            fold_rows.append({
                "fold": fold,
                "representation": rep_label,
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

            all_pred_rows.append(pd.DataFrame({
                "fold": fold,
                "y_true": res["y_test"],
                "y_pred": res["y_pred"],
            }))

            per_class = res["per_class_f1"].copy()
            per_class["fold"] = fold
            per_class_tables.append(per_class)

        fold_metrics = pd.DataFrame(fold_rows)
        fold_metrics.to_csv(
            results_dir / f"{scheme_label}_{rep_label}_fold_metrics.csv",
            index=False,
        )

        pd.concat(all_pred_rows, ignore_index=True).to_csv(
            results_dir / f"{scheme_label}_{rep_label}_all_predictions.csv",
            index=False,
        )

        per_class_all = pd.concat(per_class_tables, ignore_index=True)
        per_class_all.to_csv(
            results_dir / f"{scheme_label}_{rep_label}_all_per_class_f1.csv",
            index=False,
        )

        label_col = celltype_col

        if label_col not in per_class_all.columns:
            if "cell_type" in per_class_all.columns:
                label_col = "cell_type"
            else:
                raise KeyError(
                    f"Could not find cell-type label column in per-class metrics. "
                    f"Tried '{celltype_col}' and 'cell_type'. "
                    f"Available columns: {per_class_all.columns.tolist()}"
                )

        per_class_mean = (
            per_class_all
            .groupby(label_col, as_index=False)[["f1", "precision", "recall", "support"]]
            .mean()
            .sort_values("f1", ascending=False)
        )

        if label_col != "cell_type":
            per_class_mean = per_class_mean.rename(columns={label_col: "cell_type"})

        per_class_mean.to_csv(
            results_dir / f"{scheme_label}_{rep_label}_mean_per_class_f1.csv",
            index=False,
        )

        summary_rows.append({
            "scheme": scheme_label,
            "representation": rep_label,
            "macro_f1_mean": fold_metrics["macro_f1"].mean(),
            "macro_f1_std": fold_metrics["macro_f1"].std(ddof=1),
            "accuracy_mean": fold_metrics["accuracy"].mean(),
            "accuracy_std": fold_metrics["accuracy"].std(ddof=1),
            "n_folds": len(fold_metrics),
            "mean_train_cells": fold_metrics["n_train_cells"].mean(),
            "mean_test_cells": fold_metrics["n_test_cells"].mean(),
            "mean_n_classes_used": fold_metrics["n_classes_used"].mean(),
            "batch_covariate": batch_col if batch_col is not None else "None",
            "mean_n_batch_features": fold_metrics["n_batch_features"].mean(),
        })

    metrics = pd.DataFrame(summary_rows)
    metrics.to_csv(results_dir / "metrics.csv", index=False)

    return metrics