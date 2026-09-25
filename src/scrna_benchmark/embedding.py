from __future__ import annotations

"""Utilities for building and registering AnnData representations."""

import numpy as np


def compute_hvg_subset(
    adata,
    n_top_genes: int = 2000,
    batch_key: str | None = None,
    layer: str | None = None,
    copy: bool = True,
):
    try:
        import scanpy as sc
    except ImportError as exc:
        raise ImportError("compute_hvg_subset requires scanpy.") from exc

    adata_hvg = adata.copy() if copy else adata
    sc.pp.highly_variable_genes(
        adata_hvg,
        n_top_genes=n_top_genes,
        batch_key=batch_key,
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
    try:
        import scanpy as sc
    except ImportError as exc:
        raise ImportError("compute_pca requires scanpy.") from exc

    sc.tl.pca(
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
    key_added: str = "X_harmony",
    n_comps: int | None = 15,
    **kwargs,
):
    import harmonypy as hm

    if batch_col not in adata.obs.columns:
        raise KeyError(f"{batch_col} not found in adata.obs.")
    if basis not in adata.obsm:
        raise KeyError(f"{basis} not found in adata.obsm.")

    X = np.asarray(adata.obsm[basis])
    ho = hm.run_harmony(X, adata.obs, batch_col, **kwargs)
    Z = ho.Z_corr
    if Z.shape[0] != adata.n_obs and Z.shape[1] == adata.n_obs:
        Z = Z.T
    if Z.shape[0] != adata.n_obs:
        raise ValueError(
            f"Harmony output has shape {Z.shape}, expected first dimension {adata.n_obs}."
        )
    if n_comps is not None:
        Z = Z[:, :n_comps]
    adata.obsm[key_added] = Z.copy()
    return adata


def compute_scvi_latent(
    adata,
    batch_col: str | None = None,
    layer: str | None = None,
    key_added: str = "X_scVI",
    n_latent: int = 30,
    max_epochs: int | None = None,
    random_state: int = 42,
    train_kwargs: dict | None = None,
    **model_kwargs,
):
    """Train scVI and store the latent embedding in ``adata.obsm[key_added]``.

    ``model_kwargs`` are passed to ``scvi.model.SCVI`` while ``train_kwargs`` are
    passed to ``model.train``. This keeps training controls such as
    ``early_stopping=True`` out of the model constructor.
    """
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

    train_kwargs = dict(train_kwargs or {})
    if max_epochs is not None:
        train_kwargs.setdefault("max_epochs", max_epochs)
    model.train(**train_kwargs)
    adata.obsm[key_added] = model.get_latent_representation()
    return adata, model


def register_obsm_representation(adata, key: str, X) -> None:
    X = np.asarray(X)
    if X.shape[0] != adata.n_obs:
        raise ValueError(
            f"Representation has {X.shape[0]} rows, but adata has {adata.n_obs} cells."
        )
    adata.obsm[key] = X


def available_representations(adata) -> list[str]:
    return ["hvg"] + sorted(list(adata.obsm.keys()))
