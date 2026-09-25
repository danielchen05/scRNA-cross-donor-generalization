from __future__ import annotations

"""Deterministic preprocessing for benchmark-ready scRNA-seq AnnData objects.

The preprocessing stage owns all choices that define the frozen cell population and
representations used by downstream benchmarking. It intentionally performs no plotting.
"""

from dataclasses import asdict, dataclass, field
from importlib import metadata as importlib_metadata
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from .embedding import compute_harmony, compute_pca, compute_scvi_latent
from .filtering import filter_celltypes_by_support, summarize_celltype_support
from .representations import validate_representations


@dataclass
class FilterRule:
    column: str
    op: str = "eq"
    value: Any = None


@dataclass
class SupportConfig:
    min_cells: int = 200
    min_donors: int = 5
    # Most datasets define cell-type eligibility before downsampling. PBMC's
    # legacy workflow applies the same support rule after its donor cap.
    stage: str = "before_downsampling"


@dataclass
class DownsampleConfig:
    group_cols: list[str] = field(default_factory=lambda: ["donor_id"])
    max_cells_per_group: int = 500
    random_state: int = 42
    rng: str = "random_state"  # random_state (legacy) or default_rng
    # None preserves the existing dataset-specific behavior: RandomState samples
    # every group; default_rng samples only oversized groups. PBMC explicitly
    # sets False to reproduce its original donor-cap loop.
    sample_full_groups: bool | None = None


@dataclass
class GeneFilterConfig:
    exclude_var_if_true: list[str] = field(default_factory=list)
    min_cells: int | None = None


@dataclass
class CountsConfig:
    source: str = "raw"  # raw, X, or layer:<name>
    target_sum: float = 1e4
    require_integer_like: bool = True
    # True: rebuild normalized/log1p X from raw counts (kidney/pancreas/lung).
    # False: retain the input X as the feature matrix (PBMC legacy object).
    prepare_from_counts: bool = True


@dataclass
class HVGConfig:
    n_top_genes: int = 1000
    batch_key: str | None = "donor_id"


@dataclass
class PCAConfig:
    n_comps_compute: int = 50
    n_comps_keep: int = 15
    random_state: int = 42


@dataclass
class HarmonyConfig:
    enabled: bool = True
    batch_key: str = "donor_id"
    basis: str = "X_pca"
    key_added: str = "X_harmony"
    n_comps: int = 15
    kwargs: dict[str, Any] = field(default_factory=dict)


@dataclass
class SCVIConfig:
    enabled: bool = True
    batch_key: str = "donor_id"
    layer: str = "counts"
    key_added: str = "X_scVI"
    n_latent: int = 15
    max_epochs: int = 200
    early_stopping: bool = True
    random_state: int = 42
    model_kwargs: dict[str, Any] = field(default_factory=dict)
    train_kwargs: dict[str, Any] = field(default_factory=dict)


@dataclass
class ExpectedConfig:
    n_cells: int | None = None
    n_donors: int | None = None
    n_celltypes: int | None = None
    n_hvg: int | None = None


@dataclass
class PreprocessingConfig:
    dataset_name: str
    input_path: str | Path
    output_path: str | Path
    qc_output_dir: str | Path
    celltype_col: str = "cell_type"
    donor_col: str = "donor_id"
    filters: list[FilterRule] = field(default_factory=list)
    metadata_from_obs: dict[str, str] = field(default_factory=dict)
    support: SupportConfig = field(default_factory=SupportConfig)
    downsample: DownsampleConfig = field(default_factory=DownsampleConfig)
    gene_filter: GeneFilterConfig = field(default_factory=GeneFilterConfig)
    counts: CountsConfig = field(default_factory=CountsConfig)
    hvg: HVGConfig = field(default_factory=HVGConfig)
    pca: PCAConfig = field(default_factory=PCAConfig)
    harmony: HarmonyConfig = field(default_factory=HarmonyConfig)
    scvi: SCVIConfig = field(default_factory=SCVIConfig)
    expected: ExpectedConfig = field(default_factory=ExpectedConfig)
    representations: dict[str, str] = field(
        default_factory=lambda: {
            "hvg": "hvg",
            "pca": "X_pca",
            "harmony": "X_harmony",
            "scvi": "X_scVI",
        }
    )

    def __post_init__(self) -> None:
        self.input_path = Path(self.input_path)
        self.output_path = Path(self.output_path)
        self.qc_output_dir = Path(self.qc_output_dir)

    def to_dict(self) -> dict[str, Any]:
        out = asdict(self)
        out["input_path"] = str(self.input_path)
        out["output_path"] = str(self.output_path)
        out["qc_output_dir"] = str(self.qc_output_dir)
        return out


def _coerce_filter_rules(items: list[dict[str, Any]] | None) -> list[FilterRule]:
    return [FilterRule(**item) for item in (items or [])]


def load_preprocessing_config(path: str | Path) -> PreprocessingConfig:
    path = Path(path)
    with path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)

    return PreprocessingConfig(
        dataset_name=raw["dataset_name"],
        input_path=raw["input_path"],
        output_path=raw["output_path"],
        qc_output_dir=raw["qc_output_dir"],
        celltype_col=raw.get("celltype_col", "cell_type"),
        donor_col=raw.get("donor_col", "donor_id"),
        filters=_coerce_filter_rules(raw.get("filters")),
        metadata_from_obs=raw.get("metadata_from_obs", {}),
        support=SupportConfig(**raw.get("support", {})),
        downsample=DownsampleConfig(**raw.get("downsample", {})),
        gene_filter=GeneFilterConfig(**raw.get("gene_filter", {})),
        counts=CountsConfig(**raw.get("counts", {})),
        hvg=HVGConfig(**raw.get("hvg", {})),
        pca=PCAConfig(**raw.get("pca", {})),
        harmony=HarmonyConfig(**raw.get("harmony", {})),
        scvi=SCVIConfig(**raw.get("scvi", {})),
        expected=ExpectedConfig(**raw.get("expected", {})),
        representations=raw.get(
            "representations",
            {"hvg": "hvg", "pca": "X_pca", "harmony": "X_harmony", "scvi": "X_scVI"},
        ),
    )


def _require_obs_columns(adata, columns: list[str]) -> None:
    missing = [c for c in columns if c not in adata.obs.columns]
    if missing:
        raise KeyError(f"Missing required adata.obs columns: {missing}")


def _apply_filters(adata, rules: list[FilterRule]):
    out = adata
    for rule in rules:
        _require_obs_columns(out, [rule.column])
        series = out.obs[rule.column]
        if rule.op == "eq":
            mask = series.astype(str) == str(rule.value)
        elif rule.op == "ne":
            mask = series.astype(str) != str(rule.value)
        elif rule.op == "in":
            allowed = {str(x) for x in rule.value}
            mask = series.astype(str).isin(allowed)
        elif rule.op == "not_in":
            excluded = {str(x) for x in rule.value}
            mask = ~series.astype(str).isin(excluded)
        elif rule.op == "notna":
            mask = series.notna()
        elif rule.op == "isna":
            mask = series.isna()
        else:
            raise ValueError(f"Unsupported filter op: {rule.op}")
        out = out[mask.to_numpy()].copy()
    return out


def _assign_metadata(adata, mapping: dict[str, str]) -> None:
    if not mapping:
        return
    _require_obs_columns(adata, list(mapping.values()))
    for target, source in mapping.items():
        adata.obs[target] = adata.obs[source].astype(str)


def _downsample(
    adata,
    group_cols: list[str],
    max_cells: int,
    random_state: int,
    rng: str = "random_state",
    sample_full_groups: bool | None = None,
):
    """Deterministically downsample within metadata groups.

    ``random_state`` preserves the legacy kidney/pancreas NumPy RandomState
    behavior. ``default_rng`` reproduces the lung preparation notebook.
    ``sample_full_groups=False`` reproduces PBMC, where ``choice`` was called
    only when a donor exceeded the cap.
    """
    _require_obs_columns(adata, group_cols)
    if max_cells <= 0:
        raise ValueError("max_cells_per_group must be positive")

    if rng == "random_state":
        generator = np.random.RandomState(random_state)
    elif rng == "default_rng":
        generator = np.random.default_rng(random_state)
    else:
        raise ValueError("downsample.rng must be 'random_state' or 'default_rng'")

    if sample_full_groups is None:
        # Backward-compatible defaults used by the already-refactored datasets.
        sample_full_groups = rng == "random_state"

    keep: list[str] = []
    grouped = adata.obs.groupby(group_cols, observed=True, sort=True).groups
    for _, idx in grouped.items():
        idx = np.asarray(list(idx), dtype=object)
        n = min(len(idx), max_cells)
        if sample_full_groups or len(idx) > max_cells:
            idx = generator.choice(idx, n, replace=False)
        keep.extend(idx.tolist())
    return adata[keep].copy()


def _apply_gene_filters(adata, config: GeneFilterConfig, sc):
    """Apply optional gene-level eligibility filters before count recovery/HVG selection."""
    out = adata
    for column in config.exclude_var_if_true:
        if column in out.var.columns:
            mask = ~out.var[column].astype(bool).to_numpy()
            out = out[:, mask].copy()
    if config.min_cells is not None:
        if config.min_cells <= 0:
            raise ValueError("gene_filter.min_cells must be positive when provided")
        sc.pp.filter_genes(out, min_cells=config.min_cells)
    return out


def _extract_counts(adata, source: str):
    if source == "raw":
        if adata.raw is None:
            raise ValueError("counts.source='raw' requested, but adata.raw is None")
        raw_names = pd.Index(adata.raw.var_names)
        missing = pd.Index(adata.var_names).difference(raw_names)
        if len(missing):
            raise KeyError(f"{len(missing)} current genes are absent from adata.raw")
        return adata.raw[:, adata.var_names].X.copy()
    if source == "X":
        return adata.X.copy()
    if source.startswith("layer:"):
        layer = source.split(":", 1)[1]
        if layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' not found in adata.layers")
        return adata.layers[layer].copy()
    raise ValueError("counts.source must be one of: raw, X, layer:<name>")


def _sample_values(matrix, n: int = 10000) -> np.ndarray:
    if hasattr(matrix, "data") and not isinstance(matrix, np.ndarray):
        vals = np.asarray(matrix.data)
    else:
        vals = np.asarray(matrix).ravel()
    return vals[: min(n, len(vals))]


def _is_integer_like(matrix, atol: float = 1e-8) -> bool:
    vals = _sample_values(matrix)
    if len(vals) == 0:
        return True
    return bool(np.allclose(vals, np.round(vals), atol=atol))


def _software_versions() -> dict[str, str | None]:
    packages = [
        "numpy",
        "pandas",
        "scipy",
        "anndata",
        "scanpy",
        "harmonypy",
        "scvi-tools",
        "scikit-learn",
        "pyyaml",
    ]
    versions: dict[str, str | None] = {}
    for package in packages:
        try:
            versions[package] = importlib_metadata.version(package)
        except importlib_metadata.PackageNotFoundError:
            versions[package] = None
    return versions


def _summary(adata, celltype_col: str, donor_col: str) -> dict[str, int]:
    return {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_donors": int(adata.obs[donor_col].astype(str).nunique()),
        "n_celltypes": int(adata.obs[celltype_col].astype(str).nunique()),
    }


def _validate_expected(adata, config: PreprocessingConfig) -> dict[str, Any]:
    observed = _summary(adata, config.celltype_col, config.donor_col)
    observed["n_hvg"] = int(adata.n_vars)
    expected = asdict(config.expected)
    checks = {}
    for key, exp in expected.items():
        if exp is None:
            continue
        obs = observed[key]
        checks[key] = {"expected": exp, "observed": obs, "ok": bool(obs == exp)}
        if obs != exp:
            raise ValueError(f"Expected {key}={exp}, observed {obs} for {config.dataset_name}")
    return checks


def _save_qc_artifacts(
    *,
    adata,
    config: PreprocessingConfig,
    support_before: pd.DataFrame,
    support_after: pd.DataFrame,
    stage_summaries: dict[str, dict[str, int]],
    expected_checks: dict[str, Any],
) -> None:
    out = config.qc_output_dir
    out.mkdir(parents=True, exist_ok=True)

    support_before.to_csv(out / "celltype_support_before_filter.csv", index=False)
    support_after.to_csv(out / "celltype_support_final.csv", index=False)

    adata.obs[config.donor_col].astype(str).value_counts().rename_axis(config.donor_col).rename(
        "n_cells"
    ).reset_index().to_csv(out / "donor_counts_final.csv", index=False)

    adata.obs[config.celltype_col].astype(str).value_counts().rename_axis(config.celltype_col).rename(
        "n_cells"
    ).reset_index().to_csv(out / "celltype_counts_final.csv", index=False)

    pd.crosstab(
        adata.obs[config.donor_col].astype(str),
        adata.obs[config.celltype_col].astype(str),
    ).to_csv(out / "donor_by_celltype_final.csv")

    pd.Series(adata.obs_names.astype(str), name="obs_name").to_csv(
        out / "selected_cells.csv", index=False
    )
    pd.Series(adata.var_names.astype(str), name="gene").to_csv(out / "hvg_genes.csv", index=False)

    with (out / "preprocessing_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "dataset_name": config.dataset_name,
                "stages": stage_summaries,
                "expected_checks": expected_checks,
                "representations": {
                    key: (
                        list(adata.X.shape)
                        if rep == "hvg"
                        else list(adata.obsm[rep].shape)
                    )
                    for key, rep in config.representations.items()
                },
            },
            handle,
            indent=2,
        )

    with (out / "software_versions.json").open("w", encoding="utf-8") as handle:
        json.dump(_software_versions(), handle, indent=2)

    with (out / "resolved_config.yaml").open("w", encoding="utf-8") as handle:
        yaml.safe_dump(config.to_dict(), handle, sort_keys=False)


def run_preprocessing(config: PreprocessingConfig, *, force: bool = False):
    """Create one frozen benchmark-ready AnnData object and QC/provenance tables."""
    try:
        import anndata as ad
        import scanpy as sc
    except ImportError as exc:
        raise ImportError("run_preprocessing requires anndata and scanpy") from exc

    if config.output_path.exists() and not force:
        raise FileExistsError(
            f"Output already exists: {config.output_path}. Use force=True / --force to overwrite."
        )
    if not config.input_path.exists():
        raise FileNotFoundError(config.input_path)

    adata = ad.read_h5ad(config.input_path)
    _require_obs_columns(adata, [config.celltype_col, config.donor_col])
    stage_summaries: dict[str, dict[str, int]] = {
        "loaded": _summary(adata, config.celltype_col, config.donor_col)
    }

    adata = _apply_filters(adata, config.filters)
    _assign_metadata(adata, config.metadata_from_obs)
    for col in [config.celltype_col, config.donor_col, *config.metadata_from_obs.keys()]:
        if col in adata.obs:
            adata.obs[col] = adata.obs[col].astype(str)
    stage_summaries["after_dataset_filters"] = _summary(
        adata, config.celltype_col, config.donor_col
    )

    support_stage = config.support.stage
    if support_stage not in {"before_downsampling", "after_downsampling"}:
        raise ValueError(
            "support.stage must be 'before_downsampling' or 'after_downsampling'"
        )

    support_before: pd.DataFrame
    supported_celltypes: set[str]

    if support_stage == "before_downsampling":
        support_before = summarize_celltype_support(
            adata,
            celltype_col=config.celltype_col,
            donor_col=config.donor_col,
            min_cells=config.support.min_cells,
            min_donors=config.support.min_donors,
        )
        adata, _, _, _ = filter_celltypes_by_support(
            adata,
            celltype_col=config.celltype_col,
            donor_col=config.donor_col,
            min_cells=config.support.min_cells,
            min_donors=config.support.min_donors,
        )
        supported_celltypes = set(adata.obs[config.celltype_col].astype(str).unique())
        stage_summaries["after_celltype_support_filter"] = _summary(
            adata, config.celltype_col, config.donor_col
        )

        adata = _downsample(
            adata,
            group_cols=config.downsample.group_cols,
            max_cells=config.downsample.max_cells_per_group,
            random_state=config.downsample.random_state,
            rng=config.downsample.rng,
            sample_full_groups=config.downsample.sample_full_groups,
        )
        stage_summaries["after_downsampling"] = _summary(
            adata, config.celltype_col, config.donor_col
        )

        downsampled_celltypes = set(adata.obs[config.celltype_col].astype(str).unique())
        unexpected_celltypes = sorted(downsampled_celltypes - supported_celltypes)
        if unexpected_celltypes:
            raise ValueError(
                "Downsampled object contains cell types that were not present "
                f"after the support filter: {unexpected_celltypes}"
            )
        lost_celltypes = sorted(supported_celltypes - downsampled_celltypes)
        if lost_celltypes:
            raise ValueError(
                "Downsampling completely removed cell types that passed the "
                f"support filter: {lost_celltypes}"
            )
    else:
        # PBMC legacy order: apply the donor cap first, then define the supported
        # benchmark label set from the capped object.
        adata = _downsample(
            adata,
            group_cols=config.downsample.group_cols,
            max_cells=config.downsample.max_cells_per_group,
            random_state=config.downsample.random_state,
            rng=config.downsample.rng,
            sample_full_groups=config.downsample.sample_full_groups,
        )
        stage_summaries["after_downsampling"] = _summary(
            adata, config.celltype_col, config.donor_col
        )

        support_before = summarize_celltype_support(
            adata,
            celltype_col=config.celltype_col,
            donor_col=config.donor_col,
            min_cells=config.support.min_cells,
            min_donors=config.support.min_donors,
        )
        adata, _, _, _ = filter_celltypes_by_support(
            adata,
            celltype_col=config.celltype_col,
            donor_col=config.donor_col,
            min_cells=config.support.min_cells,
            min_donors=config.support.min_donors,
        )
        supported_celltypes = set(adata.obs[config.celltype_col].astype(str).unique())
        stage_summaries["after_celltype_support_filter"] = _summary(
            adata, config.celltype_col, config.donor_col
        )

    adata = _apply_gene_filters(adata, config.gene_filter, sc)
    stage_summaries["after_gene_filtering"] = _summary(
        adata, config.celltype_col, config.donor_col
    )

    source_counts = None
    if config.counts.prepare_from_counts:
        source_counts = _extract_counts(adata, config.counts.source)
        if config.counts.require_integer_like and not _is_integer_like(source_counts):
            raise ValueError(
                f"Counts source '{config.counts.source}' is not integer-like for {config.dataset_name}"
            )

        # Kidney/pancreas/lung rebuild the feature matrix from raw counts.
        adata.X = source_counts.copy()
        sc.pp.normalize_total(adata, target_sum=config.counts.target_sum)
        sc.pp.log1p(adata)
    elif config.scvi.enabled:
        raise ValueError(
            "counts.prepare_from_counts=False cannot be combined with scvi.enabled=True; "
            "training scVI requires a raw-count layer."
        )

    # Preserve the full feature matrix before HVG restriction. For PBMC this is
    # the normalized/log-transformed X already supplied by the source object.
    adata.raw = adata.copy()

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=config.hvg.n_top_genes,
        batch_key=config.hvg.batch_key,
        subset=False,
    )
    if "highly_variable" not in adata.var:
        raise RuntimeError("scanpy did not create adata.var['highly_variable']")
    hvg_mask = adata.var["highly_variable"].to_numpy(dtype=bool)
    if int(hvg_mask.sum()) != config.hvg.n_top_genes:
        raise ValueError(
            f"Expected {config.hvg.n_top_genes} HVGs, got {int(hvg_mask.sum())}"
        )

    if source_counts is not None:
        source_counts_hvg = source_counts[:, hvg_mask].copy()
    else:
        source_counts_hvg = None

    adata = adata[:, hvg_mask].copy()
    if source_counts_hvg is not None:
        adata.layers[config.scvi.layer] = source_counts_hvg

    compute_pca(
        adata,
        n_comps=config.pca.n_comps_compute,
        random_state=config.pca.random_state,
    )
    adata.obsm["X_pca"] = np.asarray(adata.obsm["X_pca"])[:, : config.pca.n_comps_keep].copy()

    if config.harmony.enabled:
        compute_harmony(
            adata,
            batch_col=config.harmony.batch_key,
            basis=config.harmony.basis,
            key_added=config.harmony.key_added,
            n_comps=config.harmony.n_comps,
            **config.harmony.kwargs,
        )

    if config.scvi.enabled:
        train_kwargs = dict(config.scvi.train_kwargs)
        train_kwargs.setdefault("early_stopping", config.scvi.early_stopping)
        adata, _ = compute_scvi_latent(
            adata,
            batch_col=config.scvi.batch_key,
            layer=config.scvi.layer,
            key_added=config.scvi.key_added,
            n_latent=config.scvi.n_latent,
            max_epochs=config.scvi.max_epochs,
            random_state=config.scvi.random_state,
            train_kwargs=train_kwargs,
            **config.scvi.model_kwargs,
        )
    elif config.scvi.key_added not in adata.obsm:
        raise KeyError(
            f"scvi.enabled=False but existing embedding '{config.scvi.key_added}' "
            "was not found in the input AnnData."
        )

    validate_representations(adata, config.representations)
    support_after = summarize_celltype_support(
        adata,
        celltype_col=config.celltype_col,
        donor_col=config.donor_col,
        min_cells=config.support.min_cells,
        min_donors=config.support.min_donors,
    )

    # Verify that later feature/representation steps did not alter the frozen label set.
    final_celltypes = set(adata.obs[config.celltype_col].astype(str).unique())
    unexpected_celltypes = sorted(final_celltypes - supported_celltypes)
    if unexpected_celltypes:
        raise ValueError(
            "Final object contains cell types that were not present after "
            f"the support filter: {unexpected_celltypes}"
        )
    lost_celltypes = sorted(supported_celltypes - final_celltypes)
    if lost_celltypes:
        raise ValueError(
            "Final object is missing cell types that passed the support "
            f"filter: {lost_celltypes}"
        )

    expected_checks = _validate_expected(adata, config)
    stage_summaries["benchmark_ready"] = _summary(adata, config.celltype_col, config.donor_col)

    config.output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write(config.output_path)
    _save_qc_artifacts(
        adata=adata,
        config=config,
        support_before=support_before,
        support_after=support_after,
        stage_summaries=stage_summaries,
        expected_checks=expected_checks,
    )
    return adata

