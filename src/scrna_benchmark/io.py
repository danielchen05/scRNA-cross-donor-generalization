# src/scrna_benchmark/io.py

from pathlib import Path

import pandas as pd


def ensure_dir(path):
    """
    Create directory if it does not exist.
    """
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def save_table(df, path, index=False):
    """
    Save a DataFrame to CSV.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=index)


def read_table(path, **kwargs):
    """
    Read a CSV table.
    """
    return pd.read_csv(path, **kwargs)


def save_metrics_table(df, out_dir, filename="metrics.csv"):
    """
    Save metrics table.
    """
    out_dir = ensure_dir(out_dir)
    df.to_csv(out_dir / filename, index=False)


def save_predictions_table(y_true, y_pred, out_dir, filename="predictions.csv", extra_cols=None):
    """
    Save prediction table.
    """
    out_dir = ensure_dir(out_dir)

    pred_df = pd.DataFrame({
        "y_true": y_true,
        "y_pred": y_pred,
    })

    if extra_cols is not None:
        extra_df = pd.DataFrame(extra_cols)
        pred_df = pd.concat([extra_df.reset_index(drop=True), pred_df], axis=1)

    pred_df.to_csv(out_dir / filename, index=False)

    return pred_df


def save_per_class_table(per_class_df, out_dir, filename="per_class_f1.csv"):
    """
    Save per-class metrics table.
    """
    out_dir = ensure_dir(out_dir)
    per_class_df.to_csv(out_dir / filename, index=False)


def save_confusion_matrix_table(cm, labels, out_dir, filename="confusion_matrix.csv"):
    """
    Save confusion matrix as CSV.
    """
    out_dir = ensure_dir(out_dir)
    cm_df = pd.DataFrame(cm, index=labels, columns=labels)
    cm_df.to_csv(out_dir / filename)
    return cm_df


def save_experiment_result(
    result,
    out_dir,
    prefix,
):
    """
    Save standard result dict from model-running functions.
    """
    out_dir = ensure_dir(out_dir)

    save_predictions_table(
        y_true=result["y_test"],
        y_pred=result["y_pred"],
        out_dir=out_dir,
        filename=f"{prefix}_predictions.csv",
    )

    save_per_class_table(
        per_class_df=result["per_class_f1"],
        out_dir=out_dir,
        filename=f"{prefix}_per_class_f1.csv",
    )

    save_confusion_matrix_table(
        cm=result["cm"],
        labels=result["labels"],
        out_dir=out_dir,
        filename=f"{prefix}_confusion_matrix.csv",
    )


def load_metrics_from_dirs(
    result_dirs,
    filename="metrics.csv",
    source_col="source",
):
    """
    Load metrics.csv from multiple result directories and concatenate them.
    """
    dfs = []

    for result_dir in result_dirs:
        result_dir = Path(result_dir)
        path = result_dir / filename

        if not path.exists():
            raise FileNotFoundError(f"Missing metrics file: {path}")

        df = pd.read_csv(path)
        df[source_col] = result_dir.name
        dfs.append(df)

    return pd.concat(dfs, ignore_index=True)