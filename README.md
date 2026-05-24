# scRNA-cross-donor-generalization

**Random Cell-Level Splits Introduce Systematic Bias in scRNA-seq Cell Type Annotation**  

---

## Overview

This project investigates how evaluation strategy affects the measured performance of scRNA-seq cell type classification models.

We compare:
- **Scheme A (Random cell-level split)** – commonly used but biased
- **Scheme B (Donor-held-out)** – more realistic cross-donor evaluation

Representations: HVG, PCA, Harmony, scVI

The goal is to quantify how data leakage and dataset structure impact reported model performance.

---

## Project history

This repository began as a Spring 2026 Computational Genomics final project. 
The original course-project materials, including the proposal, presentation, report draft, notebooks, and results, are preserved in:

`archive/`

The active top-level repository is being actively developed as a cleaned, extended, multi-cohort benchmark.