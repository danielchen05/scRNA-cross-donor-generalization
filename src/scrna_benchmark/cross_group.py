from pathlib import Path

import pandas as pd

from .models import run_donor_split_logreg
from .evaluation import save_prediction_outputs


def run_group_transfer_experiment(
    adata,
    representations,
    group_col,
    train_group,
    test_group,
    celltype_col="cell_type",
    donor_col="patient_id",
    batch_col=None,
    random_state=42,
    verbose=True,
):
    """Train on one group/site and evaluate on another."""
    obs = adata.obs.copy()

    if group_col not in obs.columns:
        raise KeyError(f"{group_col} not found in adata.obs")

    train_mask = obs[group_col].astype(str) == str(train_group)
    test_mask = obs[group_col].astype(str) == str(test_group)

    train_donors = sorted(obs.loc[train_mask, donor_col].astype(str).unique())
    test_donors = sorted(obs.loc[test_mask, donor_col].astype(str).unique())

    if len(train_donors) == 0:
        raise ValueError(f"No train donors found for {group_col}={train_group}")
    if len(test_donors) == 0:
        raise ValueError(f"No test donors found for {group_col}={test_group}")

    rows = []
    results = {}
    direction = f"{train_group}_to_{test_group}"

    for rep_label, rep_key in representations.items():
        if verbose:
            print(f"Running group transfer {direction}: {rep_label} ({rep_key})")

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

        results[rep_label] = res
        rows.append(
            {
                "direction": direction,
                "group_col": group_col,
                "train_group": train_group,
                "test_group": test_group,
                "representation": rep_label,
                "rep_key": rep_key,
                "macro_f1": res["macro_f1"],
                "accuracy": res["accuracy"],
                "n_train_cells": res["n_train_cells"],
                "n_test_cells": res["n_test_cells"],
                "n_classes_used": res["n_classes_used"],
                "n_train_donors": len(res["train_donors"]),
                "n_test_donors": len(res["test_donors"]),
                "batch_covariate": batch_col if batch_col is not None else "None",
                "n_batch_features": len(res["batch_feature_names"]),
            }
        )

    return pd.DataFrame(rows), results


def run_bidirectional_group_transfer(
    adata,
    representations,
    group_col,
    group_a,
    group_b,
    celltype_col="cell_type",
    donor_col="patient_id",
    batch_col=None,
    random_state=42,
    verbose=True,
):
    """Run group A -> group B and group B -> group A transfer."""
    metrics_ab, results_ab = run_group_transfer_experiment(
        adata=adata,
        representations=representations,
        group_col=group_col,
        train_group=group_a,
        test_group=group_b,
        celltype_col=celltype_col,
        donor_col=donor_col,
        batch_col=batch_col,
        random_state=random_state,
        verbose=verbose,
    )
    metrics_ba, results_ba = run_group_transfer_experiment(
        adata=adata,
        representations=representations,
        group_col=group_col,
        train_group=group_b,
        test_group=group_a,
        celltype_col=celltype_col,
        donor_col=donor_col,
        batch_col=batch_col,
        random_state=random_state,
        verbose=verbose,
    )

    metrics_df = pd.concat([metrics_ab, metrics_ba], ignore_index=True)
    results = {
        f"{group_a}_to_{group_b}": results_ab,
        f"{group_b}_to_{group_a}": results_ba,
    }
    return metrics_df, results


def save_group_transfer_outputs(
    metrics_df,
    results,
    out_dir,
    prefix="group_transfer",
):
    """Save group-transfer numerical outputs only; plotting is downstream."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    metrics_df.to_csv(out_dir / f"{prefix}_summary.csv", index=False)

    for direction, direction_results in results.items():
        for rep_label, res in direction_results.items():
            out_prefix = f"{prefix}_{direction}_{rep_label}"
            save_prediction_outputs(
                result=res,
                results_dir=out_dir,
                prefix=out_prefix,
            )
