# Results Directory

This directory contains all outputs generated from the evaluation pipelines described in:

- `03-Random-Split.ipynb` (Scheme A: random cell-level split) 
- `04-donor-held-out.ipynb` (Scheme B: donor-held-out evaluation)
- `06-other-explorations.ipynb` (Additional Data Explorations)

These results support all figures, tables, and analyses in the manuscript.

---

## Directory Structure

### Scheme A: Random Cell-Level Split

- **schemeA_case1/**  
  Filtered dataset, logistic regression (no batch covariate)

- **schemeA_case2/**  
  Filtered dataset, logistic regression + batch covariate (`Site`)

- **schemeA_case3/**  
  Full dataset (all cell types), logistic regression (no batch covariate)

- **schemeA_case4/**  
  Full dataset, logistic regression + batch covariate (`Site`)

---

### Scheme B: Donor-Held-Out Evaluation

- **schemeB_case1/**  
  Filtered dataset, 5-fold donor cross-validation (no batch covariate)

- **schemeB_case2/**  
  Filtered dataset, 5-fold donor cross-validation + batch covariate (`Site`)

- **schemeB_full_comparison/**  
  Donor-held-out evaluation on the **full dataset**, used to assess the impact of rare cell types

---

### Supplementary Analyses

- **supp_excluded_celltypes/**  
  Donor-held-out evaluation on excluded (filtered-out) cell types only.  
  Used to assess how sparsity and low representation impact classification performance.

- **supp_glmm_s3/**  
  Mixed-effects modeling analysis (GLMM) evaluating correctness while accounting for donor-level variation as a random effect.  
  Used to validate robustness of evaluation strategy beyond standard metrics.

### Aggregated Outputs

- **tables/**
  - `schemeA_summary.csv`: Summary across all Scheme A cases
  - `schemeB_summary.csv`: Summary across Scheme B filtered experiments
  - `schemeB_filtered_vs_full.csv`: Comparison of filtered vs full datasets
  - `donor_ablation_summary.csv`: Donor ablation results (performance vs # training donors)

- **figures/**
  Final figures generated from all analyses, including:
  - Main figures (`figure*`)
  - Supplementary figures (`supplementary*`)
  - Dataset overview visualizations

  These correspond directly to figures referenced in the manuscript.

---

## Outputs Within Each Case Folder

Each `schemeX_caseY/` directory contains:

- `metrics.csv`  
  Summary performance metrics (macro F1, accuracy)

- `*_predictions.csv`  
  Model predictions (`y_true`, `y_pred`)

- `*_per_class_f1.csv`  
  Per-cell-type performance (F1, precision, recall, support)

- `*_confusion_matrix.csv`  
  Raw confusion matrix

- `*_confusion_matrix_normalized.png`  
  Row-normalized confusion matrix visualization

For Scheme B (cross-validation):
- Fold-level outputs (`fold1`, ..., `fold5`)
- Aggregated confusion matrices across folds
- Mean per-class performance across folds

For supplementary GLMM analysis (`supp_glmm_s3/`):

- `glmm_input_predictions.csv`  
  Combined prediction-level dataset across all splits and representations

- `glmm_predicted_probabilities.csv`  
  Estimated probability of correct classification under each evaluation scheme

These outputs are used to quantify evaluation bias while accounting for donor-level random effects.

---

## Donor Ablation Analysis

Performed under Scheme B using the filtered dataset:

- Varies number of training donors (`k = 3–20`)
- Repeated sampling for robustness
- Results saved in:
  - `donor_ablation_raw.csv`
  - `donor_ablation_summary.csv`

This analysis quantifies how training donor diversity affects generalization.

---

## Reproducibility

To regenerate these results:

1. Run preprocessing (`01_preprocessing.ipynb`)
2. Run cell type filtering (`02_celltype_filtering.ipynb`)
3. Execute:
   - `03-Random-Split.ipynb`
   - `04-donor-held-out.ipynb`
   - `05-visualizations.ipynb`
   - `06-other-explorations.ipynb`
4. Outputs will be written automatically to this `results/` directory

Refer to the main repository `README.md` for full pipeline setup.

---

## Naming Conventions

- `schemeA`: Random cell-level split (potential donor leakage)
- `schemeB`: Donor-held-out evaluation (no donor overlap)
- `filtered`: Dataset restricted to well-represented cell types
- `full`: All cell types included
- `supp_*`: Supplementary analyses not part of the main evaluation pipeline