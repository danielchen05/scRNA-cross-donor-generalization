# src/scrna_benchmark/representations.py

import numpy as np


def get_representation(adata, rep_key):
    """
    Return a dense feature matrix for a given representation.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    rep_key : str
        'hvg' uses adata.X. Other values are keys in adata.obsm.

    Returns
    -------
    np.ndarray
        Feature matrix with shape (n_cells, n_features).
    """
    if rep_key == "hvg":
        X = adata.X
    else:
        if rep_key not in adata.obsm:
            raise KeyError(f"{rep_key} not found in adata.obsm")
        X = adata.obsm[rep_key]

    if hasattr(X, "toarray"):
        X = X.toarray()
    else:
        X = np.asarray(X)

    return X


def validate_representations(adata, representations):
    """
    Check that all requested representations exist.
    """
    missing = []

    for rep_label, rep_key in representations.items():
        if rep_key == "hvg":
            if adata.X is None:
                missing.append(rep_label)
        elif rep_key not in adata.obsm:
            missing.append(rep_label)

    if missing:
        raise KeyError(f"Missing representations: {missing}")

    return True