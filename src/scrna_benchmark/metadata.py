# src/scrna_benchmark/metadata.py

"""Metadata validation and standardization helpers."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def require_obs_columns(adata, columns: list[str]) -> None:
    """Raise a clear error if required columns are missing from adata.obs."""
    missing = [col for col in columns if col not in adata.obs.columns]
    if missing:
        raise KeyError(f"Missing required adata.obs columns: {missing}")


def rename_obs_columns(adata, rename_map: dict[str, str] | None = None, copy: bool = True):
    """Rename adata.obs columns using a dataset-specific mapping."""
    if copy:
        adata = adata.copy()

    if not rename_map:
        return adata

    missing = [old for old in rename_map if old not in adata.obs.columns]
    if missing:
        raise KeyError(f"Cannot rename missing adata.obs columns: {missing}")

    adata.obs = adata.obs.rename(columns=rename_map)
    return adata


def coerce_obs_to_str(adata, columns: list[str], copy: bool = True):
    """Convert selected adata.obs columns to strings.

    This avoids subtle mismatches when donor IDs are mixed int/string/categorical.
    """
    if copy:
        adata = adata.copy()

    require_obs_columns(adata, columns)
    for col in columns:
        adata.obs[col] = adata.obs[col].astype(str)

    return adata


def drop_missing_obs_values(adata, columns: list[str], copy: bool = True):
    """Drop cells with missing values in required metadata columns."""
    if copy:
        adata = adata.copy()

    require_obs_columns(adata, columns)

    obs = adata.obs[columns].copy()
    mask = obs.notna().all(axis=1)
    for col in columns:
        values = obs[col].astype(str).str.lower()
        mask &= ~values.isin(["nan", "none", "na", "unknown", ""])

    return adata[mask.to_numpy()].copy()


def standardize_metadata(
    adata,
    required_cols: list[str],
    rename_map: dict[str, str] | None = None,
    drop_missing: bool = True,
    stringify: bool = True,
):
    """Apply the standard metadata preparation used before benchmarking."""
    adata = rename_obs_columns(adata, rename_map=rename_map, copy=True)
    require_obs_columns(adata, required_cols)

    if drop_missing:
        adata = drop_missing_obs_values(adata, required_cols, copy=False)

    if stringify:
        adata = coerce_obs_to_str(adata, required_cols, copy=False)

    return adata


def summarize_metadata(
    adata,
    celltype_col: str,
    donor_col: str,
    batch_col: str | None = None,
    group_col: str | None = None,
) -> dict[str, int | str | None]:
    """Return a small metadata summary for sanity checks."""
    require_obs_columns(adata, [celltype_col, donor_col])

    summary = {
        "n_cells": int(adata.n_obs),
        "n_features": int(adata.n_vars),
        "celltype_col": celltype_col,
        "donor_col": donor_col,
        "n_celltypes": int(adata.obs[celltype_col].nunique()),
        "n_donors": int(adata.obs[donor_col].nunique()),
        "batch_col": batch_col,
        "n_batches": None,
        "group_col": group_col,
        "n_groups": None,
    }

    if batch_col is not None:
        require_obs_columns(adata, [batch_col])
        summary["n_batches"] = int(adata.obs[batch_col].nunique())

    if group_col is not None:
        require_obs_columns(adata, [group_col])
        summary["n_groups"] = int(adata.obs[group_col].nunique())

    return summary


def donor_celltype_coverage(adata, celltype_col: str, donor_col: str) -> pd.DataFrame:
    """Return cell counts for each donor x cell-type combination."""
    require_obs_columns(adata, [celltype_col, donor_col])
    return pd.crosstab(
        adata.obs[donor_col].astype(str),
        adata.obs[celltype_col].astype(str),
    )


def save_metadata_summary(summary: dict, out_file: str | Path) -> None:
    """Save a metadata summary dict as a two-column CSV."""
    out_file = Path(out_file)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    pd.Series(summary).rename_axis("field").reset_index(name="value").to_csv(
        out_file,
        index=False,
    )