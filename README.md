# scRNA-cross-donor-generalization
Final Research Project Code Repository for JHU EN.580.448 Computational Genomics: Data Analysis  
Manuscript: https://docs.google.com/document/d/1ljOWM2J1JYLTrstWga6W-QoZ3kqSHeyz_d5JUgY7Kwo/edit?usp=sharing  
Presentation:  
Proposal: https://docs.google.com/document/d/1RDp3KSHCICNNnTg9bc4NFLe74MFOHukmjl1v7DhUTKw/edit?usp=sharing  

---

## Project Overview

This project investigates how evaluation strategy affects the measured performance of scRNA-seq cell type classification models.

We compare:
- **Random cell-level splits (Scheme A)** — a commonly used but potentially biased evaluation
- **Donor-held-out evaluation (Scheme B)** — a more realistic cross-donor generalization setting

Across multiple transcriptomic representations:
- HVG (highly variable genes)
- PCA
- Harmony
- scVI

The goal is to quantify how data leakage and dataset structure impact reported model performance.

---

## Repository Structure
scRNA-cross-donor-generalization/
│
├── data/ # Processed AnnData objects and metadata
├── notebooks/ # Main analysis pipeline (run in order)
├── src/ # Modular Python functions used in notebooks
├── results/ # All outputs (figures, tables, metrics)
├── requirements.txt # Python dependencies
└── README.md # This file

### Key Notes

- The **main pipeline starts from `01_preprocessing.ipynb`**
- `00_prelim_viz.ipynb` is **exploratory only** and not required
- All results in the manuscript are generated from the notebooks

---

## How to Run the Pipeline

### 1. Setup environment

```bash
pip install -r requirements.txt
```

### 2. Run notebooks in order

Navigate to the `notebooks/` directory and run:  

`01_preprocessing.ipynb`  
`02_celltype.ipynb`  
`03_random_split.ipynb`  
`04_donor_held_out.ipynb`  
`05_visualizations.ipynb`  
`06_other_explorations.ipynb`  

### 3. Outputs

All outputs are saved to `results/`.  
Including:
- Figures
- Model metrics (CSV)
- Confusion matrices
- Per-class F1 scores

## Data

Data is loaded via:

`import pertpy as pt`  
`adata = pt.data.stephenson_2021_subsampled()`

Alternatively, preprocessed `.h5ad` files are provided in:

`data/`

### Key files:

- `adata_processed.h5ad`: main processed dataset  
- `adata_filtered_celltypes.h5ad`: filtered dataset used in analysis  

---

## Results

- Main figures are located in: `results/figures/`  
- Summary tables: `results/tables/`  

Detailed breakdowns of evaluation schemes are in:

- `schemeA_*`: random split experiments  
- `schemeB_*`: donor-held-out experiments  

Additional documentation is provided within subdirectories:
- `data/README.md`: dataset details and preprocessing
- `results/README.md`: detailed description of output structure and evaluation schemes

---

## Reproducibility Notes

- All results in the manuscript can be reproduced by running the notebooks sequentially  
- Precomputed outputs are included to avoid long runtimes  
- No additional scripts are required beyond the notebooks  

---

## Authors

- Xingyi (Daniel) Chen  
- Jiabei (Carol) Li  