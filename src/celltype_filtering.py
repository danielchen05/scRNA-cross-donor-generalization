import scanpy as sc
import pandas as pd


# count number of cells per cell type
def compute_celltype_counts(adata, celltype_col="cell_type"):
    return adata.obs[celltype_col].value_counts().sort_values(ascending=False)


# count number of donors each cell type appears in
def compute_donor_counts(adata, celltype_col="cell_type", donor_col="patient_id"):
    ct_donor = pd.crosstab(adata.obs[celltype_col], adata.obs[donor_col])
    return (ct_donor > 0).sum(axis=1)


# determine which cell types to keep based on thresholds
def get_kept_celltypes(
    adata,
    celltype_col="cell_type",
    donor_col="patient_id",
    cell_cutoff=200,
    donor_cutoff=5,
):
    counts = compute_celltype_counts(adata, celltype_col)
    donor_counts = compute_donor_counts(adata, celltype_col, donor_col)

    # align order
    donor_counts = donor_counts.loc[counts.index]

    keep_types = counts[
        (counts >= cell_cutoff) & (donor_counts >= donor_cutoff)
    ].index

    return keep_types, counts, donor_counts


# subset adata to selected cell types
def filter_adata_by_celltypes(adata, keep_types, celltype_col="cell_type"):
    return adata[adata.obs[celltype_col].isin(keep_types)].copy()


# full pipeline: load, filter, save outputs
def save_celltype_filtered_outputs(
    input_path="data/adata_processed.h5ad",
    full_output_path="data/adata_full_celltypes.h5ad",
    filtered_output_path="data/adata_filtered_celltypes.h5ad",
    kept_types_output_path="data/kept_cell_types.txt",
    celltype_col="cell_type",
    donor_col="patient_id",
    cell_cutoff=200,
    donor_cutoff=5,
):
    adata = sc.read_h5ad(input_path)

    keep_types, counts, donor_counts = get_kept_celltypes(
        adata,
        celltype_col,
        donor_col,
        cell_cutoff,
        donor_cutoff,
    )

    adata_filtered = filter_adata_by_celltypes(
        adata,
        keep_types,
        celltype_col,
    )

    # save outputs
    adata.write(full_output_path)
    adata_filtered.write(filtered_output_path)
    pd.Series(keep_types, name=celltype_col).to_csv(kept_types_output_path, index=False)

    return {
        "n_kept_celltypes": len(keep_types),
        "n_total_celltypes": len(counts),
        "full_shape": adata.shape,
        "filtered_shape": adata_filtered.shape,
    }