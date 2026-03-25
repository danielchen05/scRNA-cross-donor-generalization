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
├── data/               # Data access instructions (no raw data stored)
│
├── requirements.txt    # Python dependencies
├── README.md           # Project overview
└── .gitignore
```

## Current Figure Plans (Subject to Change)
Figure Ideas:  
Figure 1: Dataset overview: UMAP by cell type, UMAP by donor, donor/cell type composition bar plot  
Figure 2: Main benchmark: cell-type annotation performance under random split/donor held-out split across representations: bar plot   comparing macro F1 across representations under cell-level and donor-held-out evaluation.  
Figure 3: Per cell type performance: heatmap of F1 by cell type and method: per-class F1 by method and evaluation scheme   
Figure 4: ? donor ablation plots  
Figure 5: confusion matrix for best donor-held-out method (dependent on figure 2, may be skip)  