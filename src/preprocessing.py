import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import harmonypy as hm


# -------------------------
# STEP 1: Load data
# -------------------------
def load_data(path=None, use_pertpy=True):
    if use_pertpy:
        import pertpy as pt
        adata = pt.data.stephenson_2021_subsampled()
    else:
        adata = sc.read_h5ad(path)
    return adata


# -------------------------
# STEP 2: Filter to healthy donors + downsample
# -------------------------
def filter_and_downsample(adata, n_cells_per_donor=500, seed=42):
    adata_ctrl = adata[adata.obs["Status"] == "Healthy"].copy()

    np.random.seed(seed)
    cells_to_keep = []

    for donor, idx in adata_ctrl.obs.groupby("patient_id", observed=True).groups.items():
        idx = np.array(list(idx))
        if len(idx) > n_cells_per_donor:
            cells_to_keep.extend(np.random.choice(idx, n_cells_per_donor, replace=False))
        else:
            cells_to_keep.extend(idx)

    return adata_ctrl[cells_to_keep].copy()


# -------------------------
# STEP 3: HVG selection
# -------------------------
def compute_hvg(adata, n_top_genes=2000):
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        batch_key="patient_id"
    )

    adata.raw = adata.copy()
    adata = adata[:, adata.var["highly_variable"]].copy()

    return adata


# -------------------------
# STEP 4: PCA
# -------------------------
def compute_pca(adata, n_comps=50, keep_dims=15):
    sc.tl.pca(adata, n_comps=n_comps, random_state=42)
    adata.obsm["X_pca"] = adata.obsm["X_pca"][:, :keep_dims]
    return adata


# -------------------------
# STEP 5: Harmony
# -------------------------
def compute_harmony(adata, keep_dims=15):
    # recompute PCA for Harmony input
    sc.tl.pca(adata, n_comps=50, random_state=42)

    ho = hm.run_harmony(
        adata.obsm["X_pca"],
        adata.obs,
        "patient_id"
    )

    adata.obsm["X_harmony"] = ho.Z_corr[:, :keep_dims]
    return adata


# -------------------------
# STEP 6: Verify representations
# -------------------------
def verify_adata(adata):
    summary = {
        "X_shape": adata.X.shape,
        "raw_shape": adata.raw.X.shape if adata.raw is not None else None,
        "pca_shape": adata.obsm["X_pca"].shape,
        "harmony_shape": adata.obsm["X_harmony"].shape,
        "scvi_shape": adata.obsm["X_scVI"].shape if "X_scVI" in adata.obsm else None,
        "n_cell_types": adata.obs["cell_type"].nunique(),
        "n_donors": adata.obs["patient_id"].nunique(),
        "n_cells": adata.n_obs,
    }
    return summary


# -------------------------
# STEP 7: Full pipeline
# -------------------------
def run_preprocessing_pipeline(
    input_path=None,
    output_path="data/adata_processed.h5ad",
    use_pertpy=True
):
    adata = load_data(input_path, use_pertpy)
    adata = filter_and_downsample(adata)
    adata = compute_hvg(adata)
    adata = compute_pca(adata)
    adata = compute_harmony(adata)

    adata.write(output_path)

    return adata