# Data

This directory contains AnnData (`.h5ad`) files for the cross-donor cell-type classification benchmarks. The complete `data/` directory is distributed through [Zenodo](https://doi.org/10.5281/zenodo.19717558); expression files are excluded from GitHub. Restore this folder at the repository root, preserving the paths below.

## Files

Each cohort includes a **preprocessing input** and a **benchmark-ready object**. Inputs are source objects or intermediate checkpoints, not necessarily raw-count datasets. Benchmark-ready objects contain the final filtered cell population and expression representations used for classification.

All paths below are relative to `data/`.

| Cohort | Preprocessing input | Benchmark-ready object |
| --- | --- | --- |
| PBMC | `PBMC_Stephenson/stephenson_2021_subsampled.h5ad` | `PBMC_Stephenson/stephenson_benchmark_ready.h5ad` |
| Lung | `lung/lung.h5ad` | `lung/lung_benchmark_ready.h5ad` |
| Kidney | `kidney/kidney.h5ad` | `kidney/kidney_benchmark_ready.h5ad` |
| Pancreas | `pancreas/pancreas.h5ad` | `pancreas/pancreas_benchmark_ready.h5ad` |
| Blood Atlas | `blood_atlas/blood_atlas_subsampled.h5ad` | `blood_atlas/blood_atlas_benchmark_ready.h5ad` |

## Benchmark-ready contents

| Cohort | Cells | Donors | Cell types |
| --- | ---: | ---: | ---: |
| PBMC | 10,107 | 23 | 13 |
| Lung | 83,966 | 104 | 31 |
| Kidney | 8,746 | 19 | 9 |
| Pancreas | 8,500 | 17 | 9 |
| Blood Atlas | 108,682 | 166 | 7 |

Objects contain expression for 1,000 highly variable genes in `.X`, cell-type labels (`cell_type`) and donor identifiers in `.obs`, and embeddings in `.obsm` (`X_pca`, `X_harmony`, and `X_scVI` where available). The donor column is `patient_id` for PBMC and `donor_id` for the other cohorts.

- **PBMC:** retains the supplied processed expression and precomputed scVI embedding.
- **Blood Atlas:** uses a subsampled Stage-0 checkpoint with supplied transformed expression. Gene-level raw counts are unavailable, so no scVI embedding is included.
- **Lung, kidney, and pancreas:** expression is normalized to 10,000 counts per cell and log-transformed during preprocessing.

## Use

Use the benchmark-ready files directly to reproduce classification results; preprocessing can be skipped. To rebuild them from the input files, run from the repository root:

```bash
python scripts/run_all_preprocessing.py
```

Settings are in [configs/preprocessing/](../configs/preprocessing/). Existing benchmark-ready files are protected from overwriting unless `--force` is supplied. See the [main README](../README.md) for environment setup and analysis commands, and the [manuscript](https://doi.org/10.64898/2026.10.26.761949) for dataset provenance and original-study references.