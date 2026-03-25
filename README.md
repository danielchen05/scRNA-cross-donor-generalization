# scRNA-cross-donor-generalization
Final Research Project Code Repository for JHU EN.580.448 Computational Genomics: Data Analysis  
Proposal Link: https://docs.google.com/document/d/1RDp3KSHCICNNnTg9bc4NFLe74MFOHukmjl1v7DhUTKw/edit?usp=sharing  

## Repository Structure (Initial Setup)

```
scRNA-cross-donor-generalization/
│
├── notebooks/          # Jupyter notebooks for EDA, modeling, and evaluation
├── src/                # Reusable Python scripts (data loading, models, metrics). ALL core logic should appear here again.
├── results/
│   ├── figures/        # Final plots for reports/presentations
│   └── tables/         # Evaluation outputs if needed (CSV)
├── data/               # Data access instructions (raw data will be stored in this folder, but will not be tracked by Github)
│
├── requirements.txt    # Python dependencies
├── README.md           # Project overview
└── .gitignore
```

## Carol: Preprocessing & Representations

- [ ] Load dataset
- [ ] Filter healthy donors, downsample (refer to 00_prelim, status=healthy, downsample 1 donor)
- [ ] Normalize + log transform
- [ ] Subset adata to HVGs (adata.X should only contain HVGs)
- [ ] Store full data in adata.raw before subsetting
- [ ] Compute PCA
- [ ] Run Harmony
- [ ] Extract scVI embeddings from processed data
- [ ] Sanity-check plots (UMAP, heatmap), generate additional plots if necessary -> Figure 1 (results saved to results/figures using relative paths)
- [ ] Export processed data (.h5ad, saved in /data. Will not be tracked by Github but should save locally)

Note: Jupyter notebooks go in notebooks/ (ex: 01_preprocessing), also copy the core logic as .py files and save in src/. Import all packages at the beginning of notebooks.

**Suggested Output format (for Daniel):**
- adata.X -> HVG expression
- adata.raw -> full expression before HVG subsetting
- adata.obsm["X_pca"]
- adata.obsm["X_harmony"]
- adata.obsm["X_scvi"]

## Daniel: Modeling & Evaluation
- [ ] Load processed adata object (re-run notebooks and generate local copy of processed adata)
- [ ] Implement logistic regression
- [ ] Scheme A (random split)
- [ ] Scheme B (donor-held-out)
- [ ] Compute macro F1
- [ ] Confusion matrices
- [ ] Donor ablation analysis
- [ ] Generate result plots
- [ ] Save metrics to results/metrics.csv
- [ ] Save ablation results to results/ablation.csv
- [ ] Final results table

## Current Figure Plans (Subject to Change)
Figure Ideas:  
Figure 1: Dataset overview: UMAP by cell type, UMAP by donor, donor/cell type composition bar plot, additional?
Figure 2: Main benchmark: cell-type annotation performance under random split/donor held-out split across representations: bar plot   comparing macro F1 across representations under cell-level and donor-held-out evaluation.  
Figure 3: Per cell type performance: heatmap of F1 by cell type and method: per-class F1 by method and evaluation scheme   
Figure 4: ? donor ablation plots  
Figure 5: confusion matrix for best donor-held-out method (dependent on figure 2, may be skip)  