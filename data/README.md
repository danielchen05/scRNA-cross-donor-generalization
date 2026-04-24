## Data Access and Structure

This project uses the Stephenson et al. scRNA-seq dataset, accessed via the `pertpy` package.  
To ensure reproducibility and platform compatibility, we provide both programmatic access and pre-saved data formats.

### Data Loading Options

The dataset can be loaded in two ways:

(1) Using `pertpy`:  
```python
import pertpy as pt
adata = pt.data.stephenson_2021_subsampled()
```

(2) Using pre-saved `.h5ad` files:  
Preprocessed data objects are provided with source code submission and Zenodo snapshot to avoid dependency and compatibility issues with the `pertpy` package.

---

### Provided Data Files

All `.h5ad` files will be included in the source code submission (but are not hosted on GitHub due to size constraints).

- `stephenson_2021_subsampled.h5ad`  
  Raw downsampled dataset equivalent to the output from `pertpy`.

- `adata_processed.h5ad` / `adata_full_celltypes.h5ad`  
  Output from `01_preprocessing.ipynb`.  
  Contains normalized expression, highly variable genes (HVGs), and learned representations (e.g., PCA, scVI).  
  Used for downstream model training and evaluation.

- `adata_filtered_celltypes.h5ad`  
  Output from `02_celltype.ipynb`.  
  Contains a subset of well-represented cell types after filtering.  
  This is the primary dataset used for most analyses and experiments.

- `kept_cell_types.txt`  
  List of cell types retained after filtering.

---

### Usage Notes

- Users may choose either data loading method depending on their environment.  
- Pre-saved `.h5ad` files are recommended for faster execution and consistent reproducibility.  
- All downstream analyses, models, and figures in this project are generated from the processed `.h5ad` files.

---

### Notes on Data Availability

- Due to file size limitations, `.h5ad` files are not included in the GitHub repository.  
- These files will be provided alongside the final source code submission.