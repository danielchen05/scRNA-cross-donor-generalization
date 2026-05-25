# src/scrna_benchmark/embedding.py

"""Utilities for building and registering AnnData representations.

This file is intentionally separate from representations.py:
- representations.py retrieves already-existing matrices.
- embedding.py creates or registers new matrices in adata.X / adata.obsm.
"""

from __future__ import annotations

import numpy as np


def compute_hvg_subset(
    adata,
    n_top_genes: int = 2000,
    batch_key: str | None = None,
    flavor: str = "seurat_v3",
    layer: str | None = None,
    copy: bool = True,
):
    """Return an AnnData object subset to highly variable genes.

    Your current pipeline treats rep_key='hvg' as adata.X, so this function is
    useful for creating a filtered AnnData where X is the HVG feature matrix.
    """
    try:
        import scanpy as sc
    except ImportError as exc:
        raise ImportError("compute_hvg_subset requires scanpy.") from exc

    adata_hvg = adata.copy() if copy else adata

    sc.pp.highly_variable_genes(
        adata_hvg,
        n_top_genes=n_top_genes,
        batch_key=batch_key,
        flavor=flavor,
        layer=layer,
        subset=False,
    )

    if "highly_variable" not in adata_hvg.var.columns:
        raise RuntimeError("scanpy did not create adata.var['highly_variable'].")

    return adata_hvg[:, adata_hvg.var["highly_variable"].to_numpy()].copy()


def compute_pca(
    adata,
    n_comps: int = 50,
    key_added: str = "X_pca",
    zero_center: bool = True,
    random_state: int = 42,
):
    """Compute PCA and store it in adata.obsm[key_added]."""
    try:
        import scanpy as sc
    except ImportError as exc:
        raise ImportError("compute_pca requires scanpy.") from exc

    sc.pp.pca(
        adata,
        n_comps=n_comps,
        zero_center=zero_center,
        random_state=random_state,
    )

    if key_added != "X_pca":
        adata.obsm[key_added] = adata.obsm["X_pca"].copy()

    return adata


def compute_harmony(
    adata,
    batch_col: str,
    basis: str = "X_pca",
    key_added: str = "X_pca_harmony",
    **kwargs,
):
    """Run Harmony batch correction on a PCA-like representation."""
    if basis not in adata.obsm:
        raise KeyError(f"{basis} not found in adata.obsm. Run PCA first.")
    if batch_col not in adata.obs.columns:
        raise KeyError(f"{batch_col} not found in adata.obs.")

    try:
        import scanpy.external as sce
    except ImportError as exc:
        raise ImportError("compute_harmony requires scanpy.external/harmonypy.") from exc

    sce.pp.harmony_integrate(
        adata,
        key=batch_col,
        basis=basis,
        adjusted_basis=key_added,
        **kwargs,
    )
    return adata


def compute_scanorama(
    adata,
    batch_col: str,
    basis: str = "X_pca",
    key_added: str = "X_scanorama",
    **kwargs,
):
    """Run Scanorama integration and store the corrected representation.

    This uses scanpy.external.pp.scanorama_integrate if available.
    """
    if basis not in adata.obsm:
        raise KeyError(f"{basis} not found in adata.obsm. Run PCA first.")
    if batch_col not in adata.obs.columns:
        raise KeyError(f"{batch_col} not found in adata.obs.")

    try:
        import scanpy.external as sce
    except ImportError as exc:
        raise ImportError("compute_scanorama requires scanpy.external and scanorama.") from exc

    sce.pp.scanorama_integrate(
        adata,
        key=batch_col,
        basis=basis,
        adjusted_basis=key_added,
        **kwargs,
    )
    return adata


def compute_scvi_latent(
    adata,
    batch_col: str | None = None,
    layer: str | None = None,
    key_added: str = "X_scVI",
    n_latent: int = 30,
    max_epochs: int | None = None,
    random_state: int = 42,
    **model_kwargs,
):
    """Train scVI and store the latent embedding in adata.obsm[key_added]."""
    try:
        import scvi
    except ImportError as exc:
        raise ImportError("compute_scvi_latent requires scvi-tools.") from exc

    scvi.settings.seed = random_state

    setup_kwargs = {}
    if batch_col is not None:
        if batch_col not in adata.obs.columns:
            raise KeyError(f"{batch_col} not found in adata.obs.")
        setup_kwargs["batch_key"] = batch_col
    if layer is not None:
        setup_kwargs["layer"] = layer

    scvi.model.SCVI.setup_anndata(adata, **setup_kwargs)
    model = scvi.model.SCVI(adata, n_latent=n_latent, **model_kwargs)
    model.train(max_epochs=max_epochs)
    adata.obsm[key_added] = model.get_latent_representation()

    return adata, model


def register_obsm_representation(adata, key: str, X) -> None:
    """Register an externally computed matrix as adata.obsm[key]."""
    X = np.asarray(X)
    if X.shape[0] != adata.n_obs:
        raise ValueError(
            f"Representation has {X.shape[0]} rows, but adata has {adata.n_obs} cells."
        )
    adata.obsm[key] = X


def available_representations(adata) -> list[str]:
    """Return representation keys currently available to the benchmark."""
    return ["hvg"] + sorted(list(adata.obsm.keys()))
