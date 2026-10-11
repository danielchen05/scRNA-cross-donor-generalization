# Data

Expression data are excluded from GitHub. All datasets are available with the repository snapshot on Zenodo.

## Expected files

Paths below are relative to `data/`.

| Cohort | Preprocessing input | Benchmark-ready output |
| --- | --- | --- |
| PBMC | `PBMC_Stephenson/stephenson_2021_subsampled.h5ad` | `PBMC_Stephenson/stephenson_benchmark_ready.h5ad` |
| Lung | `lung/lung.h5ad` | `lung/lung_benchmark_ready.h5ad` |
| Kidney | `kidney/kidney.h5ad` | `kidney/kidney_benchmark_ready.h5ad` |
| Pancreas | `pancreas/pancreas.h5ad` | `pancreas/pancreas_benchmark_ready.h5ad` |
| Blood Atlas | `blood_atlas/blood_atlas_subsampled.h5ad` | `blood_atlas/blood_atlas_benchmark_ready.h5ad` |

## Preprocessing

Run `scripts/preprocess.py` with the corresponding configuration in
[configs/preprocessing](../configs/preprocessing).  

PBMC uses an existing subsampled object with a precomputed scVI embedding.
Blood Atlas uses the checkpoint produced by the
[Stage-0 notebook](../notebooks/preprocessing/blood_atlas_stage0_subsample.ipynb).

If benchmark-ready files are already available, preprocessing can be skipped.