# src/scrna_benchmark/filtering.py

from pathlib import Path

import pandas as pd


def summarize_celltype_support(
    adata,
    celltype_col="cell_type",
    donor_col="patient_id",
    min_cells=None,
    min_donors=None,
):
    """
    Summarize cell-type support by total cell count and donor coverage.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object.
    celltype_col : str
        Column in adata.obs containing cell-type labels.
    donor_col : str
        Column in adata.obs containing donor/sample IDs.
    min_cells : int or None
        Optional cell-count threshold.
    min_donors : int or None
        Optional donor-coverage threshold.

    Returns
    -------
    pd.DataFrame
        Table with cell_type, n_cells, n_donors, and optional keep flags.
    """
    obs = adata.obs.copy()

    if celltype_col not in obs.columns:
        raise KeyError(f"{celltype_col} not found in adata.obs")

    if donor_col not in obs.columns:
        raise KeyError(f"{donor_col} not found in adata.obs")

    cell_counts = (
        obs[celltype_col]
        .astype(str)
        .value_counts()
        .rename("n_cells")
    )

    ct_donor = pd.crosstab(
        obs[celltype_col].astype(str),
        obs[donor_col].astype(str),
    )
    donor_counts = (ct_donor > 0).sum(axis=1).rename("n_donors")

    support_df = (
        pd.concat([cell_counts, donor_counts], axis=1)
        .reset_index()
        .rename(columns={"index": "cell_type"})
        .sort_values(["n_cells", "n_donors"], ascending=[False, False])
        .reset_index(drop=True)
    )

    if min_cells is not None:
        support_df["keep_by_cell_count"] = support_df["n_cells"] >= min_cells

    if min_donors is not None:
        support_df["keep_by_donor_coverage"] = support_df["n_donors"] >= min_donors

    if min_cells is not None and min_donors is not None:
        support_df["keep"] = (
            support_df["keep_by_cell_count"]
            & support_df["keep_by_donor_coverage"]
        )

    return support_df


def filter_celltypes_by_support(
    adata,
    celltype_col="cell_type",
    donor_col="patient_id",
    min_cells=200,
    min_donors=5,
):
    """
    Filter AnnData to cell types satisfying minimum cell count and donor coverage.

    Returns
    -------
    adata_filtered : AnnData
        Filtered AnnData object.
    support_df : pd.DataFrame
        Cell-type support table.
    kept_celltypes : list[str]
        Cell types retained.
    excluded_celltypes : list[str]
        Cell types excluded.
    """
    support_df = summarize_celltype_support(
        adata=adata,
        celltype_col=celltype_col,
        donor_col=donor_col,
        min_cells=min_cells,
        min_donors=min_donors,
    )

    kept_celltypes = support_df.loc[support_df["keep"], "cell_type"].tolist()
    excluded_celltypes = support_df.loc[~support_df["keep"], "cell_type"].tolist()

    adata_filtered = adata[
        adata.obs[celltype_col].astype(str).isin(kept_celltypes)
    ].copy()

    return adata_filtered, support_df, kept_celltypes, excluded_celltypes


def save_celltype_filtering_outputs(
    adata_full,
    adata_filtered,
    support_df,
    kept_celltypes,
    excluded_celltypes,
    out_dir,
    dataset_name=None,
):
    """
    Save full/filtered AnnData objects and cell-type filtering tables.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    prefix = f"{dataset_name}_" if dataset_name else ""

    adata_full.write(out_dir / f"{prefix}adata_full_celltypes.h5ad")
    adata_filtered.write(out_dir / f"{prefix}adata_filtered_celltypes.h5ad")

    support_df.to_csv(out_dir / f"{prefix}celltype_support.csv", index=False)

    pd.Series(kept_celltypes, name="cell_type").to_csv(
        out_dir / f"{prefix}kept_cell_types.csv",
        index=False,
    )

    pd.Series(excluded_celltypes, name="cell_type").to_csv(
        out_dir / f"{prefix}excluded_cell_types.csv",
        index=False,
    )


def subset_adata_by_celltypes(
    adata,
    celltypes,
    celltype_col="cell_type",
):
    """
    Return a copy of AnnData restricted to selected cell types.
    """
    return adata[
        adata.obs[celltype_col].astype(str).isin(celltypes)
    ].copy()