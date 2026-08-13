# DRACC cell-type labeling

This folder transfers MSK subcelltype labels to the final DRACC cells in a
shared PCA space and records holdout validation of the label-transfer rule.

## Contents

- [`01_prepare_pca.R`](01_prepare_pca.R) aligns direct MSK labels to DRACC
  cells, selects 2,000 variable genes, and computes 30 principal components.
- [`02_label_dracc_cells_knn.R`](02_label_dracc_cells_knn.R) validates and
  applies 15-nearest-neighbor label transfer with a 70% vote threshold.
- [`01_prepare_pca.rds`](outputs/01_prepare_pca.rds) is the compact PCA and
  metadata checkpoint needed to rerun labeling without the private inputs.
- [`cell_type_annotations.csv`](outputs/02_label_dracc_cells_knn/cell_type_annotations.csv)
  contains the final label, source, and kNN vote fraction for every DRACC cell.
- [`knn_validation_summary.csv`](outputs/02_label_dracc_cells_knn/knn_validation_summary.csv)
  reports overall and stage-level holdout performance.

## Main Findings

- Direct barcode overlap supplies MSK labels for 113,891 of 121,444 DRACC
  cells; kNN assigns 7,178 additional cells and leaves 375 unassigned.
- In a stage-and-label-stratified 20% holdout, 98.10% of cells pass the vote
  threshold and accepted predictions are 99.53% accurate.
- Stage-level accepted-prediction accuracy ranges from 98.71% to 100%.

## Scope and Data Availability

The combined DRACC Seurat object and canonical MSK metadata are not distributed
in this repository. The retained PCA checkpoint makes the labeling and
validation step reproducible from public files, while rebuilding the PCA still
requires those source objects. Direct MSK labels always take precedence; kNN
is applied only to cells without a direct barcode match.
