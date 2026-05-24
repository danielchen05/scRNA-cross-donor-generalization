# src/scrna_benchmark/evaluation.py

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
    ConfusionMatrixDisplay,
)


def summarize_classification(y_true, y_pred, labels=None):
    """
    Compute macro F1, accuracy, confusion matrix, and per-class metrics.
    """
    if labels is None:
        labels = np.unique(np.concatenate([y_true, y_pred]))

    cm = confusion_matrix(y_true, y_pred, labels=labels)
    macro_f1 = f1_score(y_true, y_pred, average="macro")
    accuracy = accuracy_score(y_true, y_pred)

    report = classification_report(
        y_true,
        y_pred,
        labels=labels,
        output_dict=True,
        zero_division=0,
    )

    per_class = pd.DataFrame({
        "cell_type": labels,
        "f1": [report[label]["f1-score"] for label in labels],
        "precision": [report[label]["precision"] for label in labels],
        "recall": [report[label]["recall"] for label in labels],
        "support": [report[label]["support"] for label in labels],
    })

    return {
        "macro_f1": macro_f1,
        "accuracy": accuracy,
        "cm": cm,
        "labels": labels,
        "per_class_f1": per_class,
    }


def save_confusion_matrix(cm, labels, out_file, title=None, normalize=True):
    """
    Save a confusion matrix figure.
    """
    out_file = Path(out_file)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    if normalize:
        row_sums = cm.sum(axis=1, keepdims=True)
        cm_plot = np.divide(
            cm.astype(float),
            row_sums,
            out=np.zeros_like(cm, dtype=float),
            where=row_sums != 0,
        )
    else:
        cm_plot = cm

    fig, ax = plt.subplots(figsize=(12, 10))
    disp = ConfusionMatrixDisplay(
        confusion_matrix=cm_plot,
        display_labels=labels,
    )
    disp.plot(
        ax=ax,
        xticks_rotation=90,
        cmap="Blues",
        colorbar=True,
        values_format=None,
    )

    if disp.text_ is not None:
        for text in disp.text_.ravel():
            if text is not None:
                text.set_visible(False)

    if title is not None:
        ax.set_title(title)

    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close()


def save_prediction_outputs(result, results_dir, prefix):
    """
    Save predictions, per-class metrics, and confusion matrix table.
    """
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    pd.DataFrame({
        "y_true": result["y_test"],
        "y_pred": result["y_pred"],
    }).to_csv(results_dir / f"{prefix}_predictions.csv", index=False)

    result["per_class_f1"].to_csv(
        results_dir / f"{prefix}_per_class_f1.csv",
        index=False,
    )

    cm_df = pd.DataFrame(
        result["cm"],
        index=result["labels"],
        columns=result["labels"],
    )
    cm_df.to_csv(results_dir / f"{prefix}_confusion_matrix.csv")