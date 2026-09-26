#!/usr/bin/env python
"""Run compute-heavy Blood Atlas sensitivity analyses.

Examples
--------
Core probes on the frozen benchmark-ready object::

    python scripts/run_blood_atlas_sensitivity.py --all-core

Generate the donor-only, composition-preserving Stage-0 object::

    python scripts/run_blood_atlas_sensitivity.py \
        --prepare-composition-stage0 \
        --full-h5ad /path/to/all_pbmcs_rna.h5ad \
        --metadata-csv /path/to/all_pbmcs_metadata.csv \
        --metadata-id-col "Unnamed: 0"

Then run the normal preprocessing entrypoint::

    python scripts/preprocess.py \
        --config configs/preprocessing/blood_atlas_composition_preserving.yaml

Finally compare the primary and composition-preserving checkpoints::

    python scripts/run_blood_atlas_sensitivity.py \
        --composition-eval --composition-permutation
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import scanpy as sc


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC = REPO_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from scrna_benchmark.ablation import run_donor_ablation  # noqa: E402
from scrna_benchmark.diagnostics import (  # noqa: E402
    align_metadata_to_obs,
    audit_metadata_alignment,
    evaluate_random_vs_donor,
    run_coarse_label_sensitivity,
    run_label_permutation_sensitivity,
    run_reduced_donor_cohorts,
)


PRIMARY_PATH = REPO_ROOT / "data/blood_atlas/blood_atlas_benchmark_ready.h5ad"
ALT_STAGE0_PATH = REPO_ROOT / "data/blood_atlas/blood_atlas_donor_only_subsampled.h5ad"
ALT_READY_PATH = (
    REPO_ROOT / "data/blood_atlas/blood_atlas_composition_preserving_benchmark_ready.h5ad"
)
OUT_ROOT = REPO_ROOT / "results/blood_atlas_diagnostics"

CELLTYPE_COL = "cell_type"
DONOR_COL = "donor_id"
REPRESENTATIONS = {
    "hvg": "hvg",
    "pca": "X_pca",
    "harmony": "X_harmony",
}


def _write(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    print(f"wrote {path.relative_to(REPO_ROOT)}")


def _available_representations(adata) -> dict[str, str]:
    reps: dict[str, str] = {}
    for label, key in REPRESENTATIONS.items():
        if key == "hvg" or key in adata.obsm:
            reps[label] = key
    if not reps:
        raise ValueError("No configured Blood Atlas representations are available.")
    return reps


def _coarse_mapping(labels: list[str]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for label in labels:
        low = label.lower()
        if (
            "cd4" in low
            or "cd8" in low
            or "mait" in low
            or "gd t" in low
            or "γδ" in low
            or "gamma delta" in low
        ):
            mapping[label] = "T cells"
        else:
            mapping[label] = label
    return mapping


def run_donor20(adata, args) -> pd.DataFrame:
    reps = _available_representations(adata)
    raw = run_reduced_donor_cohorts(
        adata,
        reps,
        celltype_col=CELLTYPE_COL,
        donor_col=DONOR_COL,
        n_donors=args.n_donors,
        n_cohorts=args.n_cohorts,
        n_folds=args.n_folds,
        seed=args.seed,
    )
    _write(raw, OUT_ROOT / "donor20/donor20_random_vs_donor_raw.csv")
    summary = (
        raw.groupby("representation", as_index=False)
        .agg(
            random_macro_f1_mean=("random_macro_f1", "mean"),
            random_macro_f1_std=("random_macro_f1", "std"),
            donor_macro_f1_mean=("donor_macro_f1_mean", "mean"),
            donor_macro_f1_between_cohort_std=("donor_macro_f1_mean", "std"),
            delta_macro_f1_mean=("delta_macro_f1_random_minus_donor", "mean"),
            delta_macro_f1_std=("delta_macro_f1_random_minus_donor", "std"),
            n_cohorts=("cohort_id", "nunique"),
        )
    )
    _write(summary, OUT_ROOT / "donor20/donor20_random_vs_donor_summary.csv")
    return raw


def run_donor20_ablation(adata, args, donor20_raw: pd.DataFrame | None = None) -> None:
    if donor20_raw is None:
        path = OUT_ROOT / "donor20/donor20_random_vs_donor_raw.csv"
        if not path.exists():
            donor20_raw = run_donor20(adata, args)
        else:
            donor20_raw = pd.read_csv(path)

    reps = _available_representations(adata)
    cohort_info = (
        donor20_raw[["cohort_id", "sampled_donors"]]
        .drop_duplicates()
        .sort_values("cohort_id")
        .head(args.max_ablation_cohorts)
    )
    parts = []
    for _, row in cohort_info.iterrows():
        cohort_id = int(row["cohort_id"])
        donors = str(row["sampled_donors"]).split(";")
        mask = adata.obs[DONOR_COL].astype(str).isin(donors).to_numpy()
        sub = adata[mask].copy()
        raw, _ = run_donor_ablation(
            adata=sub,
            representations=reps,
            k_values=args.ablation_k,
            n_repeats=args.ablation_repeats,
            celltype_col=CELLTYPE_COL,
            donor_col=DONOR_COL,
            random_state=args.seed + 100 * cohort_id,
            verbose=False,
        )
        raw.insert(0, "cohort_id", cohort_id)
        parts.append(raw)

    if not parts:
        return
    all_raw = pd.concat(parts, ignore_index=True)
    _write(all_raw, OUT_ROOT / "donor20_ablation/donor20_ablation_raw.csv")
    summary = (
        all_raw.groupby(["representation", "rep_key", "k_train_donors"], as_index=False)
        .agg(
            macro_f1_mean=("macro_f1", "mean"),
            macro_f1_std=("macro_f1", "std"),
            accuracy_mean=("accuracy", "mean"),
            accuracy_std=("accuracy", "std"),
            n_runs=("macro_f1", "size"),
            n_cohorts=("cohort_id", "nunique"),
        )
    )
    _write(summary, OUT_ROOT / "donor20_ablation/donor20_ablation_summary.csv")


def run_permutations(adata, args) -> None:
    raw = run_label_permutation_sensitivity(
        adata,
        _available_representations(adata),
        celltype_col=CELLTYPE_COL,
        donor_col=DONOR_COL,
        modes=("global", "within_donor"),
        n_permutations=args.n_permutations,
        n_folds=args.n_folds,
        seed=args.seed + 10000,
    )
    _write(raw, OUT_ROOT / "label_permutation/label_permutation_raw.csv")
    summary = (
        raw.groupby(["permutation_mode", "representation"], as_index=False)
        .agg(
            random_macro_f1_mean=("random_macro_f1", "mean"),
            random_macro_f1_std=("random_macro_f1", "std"),
            donor_macro_f1_mean=("donor_macro_f1_mean", "mean"),
            donor_macro_f1_std=("donor_macro_f1_mean", "std"),
            delta_macro_f1_mean=("delta_macro_f1_random_minus_donor", "mean"),
            delta_macro_f1_std=("delta_macro_f1_random_minus_donor", "std"),
            n_permutations=("permutation_id", "nunique"),
        )
    )
    _write(summary, OUT_ROOT / "label_permutation/label_permutation_summary.csv")


def run_coarse(adata, args) -> None:
    labels = sorted(adata.obs[CELLTYPE_COL].astype(str).unique())
    mapping = _coarse_mapping(labels)
    print("Coarse-label mapping:")
    for old, new in mapping.items():
        print(f"  {old} -> {new}")
    raw = run_coarse_label_sensitivity(
        adata,
        _available_representations(adata),
        celltype_col=CELLTYPE_COL,
        donor_col=DONOR_COL,
        mapping=mapping,
        n_repeats=args.coarse_repeats,
        n_folds=args.n_folds,
        seed=args.seed + 20000,
    )
    _write(raw, OUT_ROOT / "coarse_labels/coarse_4class_raw.csv")
    summary = (
        raw.groupby("representation", as_index=False)
        .agg(
            random_macro_f1_mean=("random_macro_f1", "mean"),
            random_macro_f1_std=("random_macro_f1", "std"),
            donor_macro_f1_mean=("donor_macro_f1_mean", "mean"),
            donor_macro_f1_std=("donor_macro_f1_mean", "std"),
            delta_macro_f1_mean=("delta_macro_f1_random_minus_donor", "mean"),
            delta_macro_f1_std=("delta_macro_f1_random_minus_donor", "std"),
        )
    )
    _write(summary, OUT_ROOT / "coarse_labels/coarse_4class_summary.csv")


def prepare_composition_preserving_stage0(args, primary) -> None:
    if args.full_h5ad is None or args.metadata_csv is None:
        raise ValueError(
            "--prepare-composition-stage0 requires --full-h5ad and --metadata-csv."
        )
    full_h5ad = Path(args.full_h5ad).expanduser()
    metadata_csv = Path(args.metadata_csv).expanduser()
    meta = pd.read_csv(metadata_csv)
    full = sc.read_h5ad(full_h5ad, backed="r")
    try:
        audit = audit_metadata_alignment(full.obs_names, meta, args.metadata_id_col)
        _write(audit, OUT_ROOT / "alignment/full_source_metadata_alignment_audit.csv")
        if not bool(audit.iloc[0]["safe_for_id_reindex"]):
            raise ValueError(
                "Full Blood Atlas expression and metadata do not have identical, unique "
                "cell-ID sets. Do not construct a composition-preserving checkpoint until "
                "the source alignment is resolved. See the saved audit table."
            )

        aligned = align_metadata_to_obs(
            full.obs_names,
            meta,
            args.metadata_id_col,
            require_exact_set=True,
        )
        rename_map = {
            "Donor_id": "donor_id",
            "Cluster_names": "cell_type",
            "Batch": "batch",
            "File_name": "sample_id",
            "Age": "age",
            "Sex": "sex",
        }
        for source, target in rename_map.items():
            if source in aligned.columns and target not in aligned.columns:
                aligned[target] = aligned[source]

        for required in [DONOR_COL, CELLTYPE_COL]:
            if required not in aligned.columns:
                raise KeyError(
                    f"Aligned metadata lacks {required!r}. Available: {aligned.columns.tolist()}"
                )

        final_labels = set(primary.obs[CELLTYPE_COL].astype(str).unique())
        aligned[DONOR_COL] = aligned[DONOR_COL].astype(str)
        aligned[CELLTYPE_COL] = aligned[CELLTYPE_COL].astype(str)
        eligible = aligned[aligned[CELLTYPE_COL].isin(final_labels)].copy()

        rng = np.random.default_rng(args.stage0_seed)
        selected: list[str] = []
        for _, group in eligible.groupby(DONOR_COL, sort=True):
            names = group.index.astype(str).to_numpy()
            n = min(len(names), args.max_cells_per_donor)
            if len(names) > n:
                names = rng.choice(names, size=n, replace=False)
            selected.extend(names.tolist())

        selected_set = set(selected)
        obs_index = pd.Index(full.obs_names.astype(str))
        keep_mask = obs_index.isin(selected_set)
        if int(keep_mask.sum()) != len(selected_set):
            raise RuntimeError("Not all selected cell IDs were found in the full AnnData.")

        print(
            f"Creating donor-only Stage-0 object with {int(keep_mask.sum()):,} cells "
            f"across {eligible[DONOR_COL].nunique()} donors."
        )
        subset = full[keep_mask, :].to_memory()
        aligned_subset = aligned.reindex(subset.obs_names.astype(str)).copy()
        aligned_subset.index = subset.obs_names
        subset.obs = aligned_subset

        ALT_STAGE0_PATH.parent.mkdir(parents=True, exist_ok=True)
        subset.write_h5ad(ALT_STAGE0_PATH)
        print(f"wrote {ALT_STAGE0_PATH.relative_to(REPO_ROOT)}")

        counts = (
            subset.obs.groupby([DONOR_COL, CELLTYPE_COL], observed=True)
            .size()
            .rename("n_cells")
            .reset_index()
        )
        _write(counts, OUT_ROOT / "composition_preserving/stage0_group_counts.csv")
    finally:
        try:
            full.file.close()
        except Exception:
            pass


def run_composition_eval(primary, args) -> None:
    if not ALT_READY_PATH.exists():
        raise FileNotFoundError(
            f"Missing {ALT_READY_PATH}. First generate Stage-0 and run the provided "
            "composition-preserving preprocessing config."
        )
    alt = sc.read_h5ad(ALT_READY_PATH)
    parts = []
    for cohort_name, adata in [("primary_capped", primary), ("composition_preserving", alt)]:
        for repeat in range(args.composition_repeats):
            result = evaluate_random_vs_donor(
                adata,
                _available_representations(adata),
                celltype_col=CELLTYPE_COL,
                donor_col=DONOR_COL,
                random_state=args.seed + 30000 + repeat,
                n_folds=args.n_folds,
            )
            result.insert(0, "repeat", repeat)
            result.insert(0, "cohort_design", cohort_name)
            parts.append(result)
    raw = pd.concat(parts, ignore_index=True)
    _write(raw, OUT_ROOT / "composition_preserving/random_vs_donor_raw.csv")
    summary = (
        raw.groupby(["cohort_design", "representation"], as_index=False)
        .agg(
            random_macro_f1_mean=("random_macro_f1", "mean"),
            donor_macro_f1_mean=("donor_macro_f1_mean", "mean"),
            delta_macro_f1_mean=("delta_macro_f1_random_minus_donor", "mean"),
            delta_macro_f1_std=("delta_macro_f1_random_minus_donor", "std"),
        )
    )
    _write(summary, OUT_ROOT / "composition_preserving/random_vs_donor_summary.csv")


def run_composition_permutation(args) -> None:
    if not ALT_READY_PATH.exists():
        raise FileNotFoundError(f"Missing {ALT_READY_PATH}")
    alt = sc.read_h5ad(ALT_READY_PATH)
    raw = run_label_permutation_sensitivity(
        alt,
        _available_representations(alt),
        celltype_col=CELLTYPE_COL,
        donor_col=DONOR_COL,
        modes=("within_donor",),
        n_permutations=args.n_permutations,
        n_folds=args.n_folds,
        seed=args.seed + 40000,
    )
    raw.insert(0, "cohort_design", "composition_preserving")
    _write(
        raw,
        OUT_ROOT / "composition_preserving/within_donor_permutation_raw.csv",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", type=Path, default=PRIMARY_PATH)
    parser.add_argument("--all-core", action="store_true")
    parser.add_argument("--donor20", action="store_true")
    parser.add_argument("--donor20-ablation", action="store_true")
    parser.add_argument("--permutations", action="store_true")
    parser.add_argument("--coarse", action="store_true")
    parser.add_argument("--prepare-composition-stage0", action="store_true")
    parser.add_argument("--composition-eval", action="store_true")
    parser.add_argument("--composition-permutation", action="store_true")

    parser.add_argument("--n-donors", type=int, default=20)
    parser.add_argument("--n-cohorts", type=int, default=10)
    parser.add_argument("--n-folds", type=int, default=5)
    parser.add_argument("--n-permutations", type=int, default=10)
    parser.add_argument("--coarse-repeats", type=int, default=5)
    parser.add_argument("--composition-repeats", type=int, default=5)
    parser.add_argument("--seed", type=int, default=2026)

    parser.add_argument("--ablation-k", type=int, nargs="+", default=[5, 10, 15])
    parser.add_argument("--ablation-repeats", type=int, default=5)
    parser.add_argument("--max-ablation-cohorts", type=int, default=5)

    parser.add_argument("--full-h5ad", type=Path)
    parser.add_argument("--metadata-csv", type=Path)
    parser.add_argument("--metadata-id-col", default="Unnamed: 0")
    parser.add_argument("--max-cells-per-donor", type=int, default=650)
    parser.add_argument("--stage0-seed", type=int, default=0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.all_core:
        args.donor20 = True
        args.donor20_ablation = True
        args.permutations = True
        args.coarse = True

    primary = sc.read_h5ad(args.benchmark)
    donor20_raw = None
    if args.donor20:
        donor20_raw = run_donor20(primary, args)
    if args.donor20_ablation:
        run_donor20_ablation(primary, args, donor20_raw)
    if args.permutations:
        run_permutations(primary, args)
    if args.coarse:
        run_coarse(primary, args)
    if args.prepare_composition_stage0:
        prepare_composition_preserving_stage0(args, primary)
    if args.composition_eval:
        run_composition_eval(primary, args)
    if args.composition_permutation:
        run_composition_permutation(args)

    if not any(
        [
            args.donor20,
            args.donor20_ablation,
            args.permutations,
            args.coarse,
            args.prepare_composition_stage0,
            args.composition_eval,
            args.composition_permutation,
        ]
    ):
        print("No action requested. Use --all-core or one of the analysis flags.")


if __name__ == "__main__":
    main()
