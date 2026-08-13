# Genotype effects on cell-type composition

This folder tests whether each knockout changes cell-type composition within
post-endocrine DRACC samples while accounting for LARRY clone structure.

## Contents

- [`01_genotype_effects_on_cell_type_composition.R`](01_genotype_effects_on_cell_type_composition.R)
  analyzes PP1, PP2, S5_1, S5_2, S6_1, and S6_2 separately. It excludes
  unassigned cells, requires at least 200 cells per modeled cell type, and fits
  `is_cell_type ~ geneBC_type + (1 | larryBC)` with WT as the reference. The
  Bayesian binomial mixed models use a normal fixed-effect prior; response-scale
  KO-versus-WT contrasts use Dunnett adjustment.
- [`genotype_effects_on_cell_type_composition.csv`](outputs/01_genotype_effects_on_cell_type_composition/genotype_effects_on_cell_type_composition.csv)
  contains the complete public contrast table with odds ratios, confidence
  intervals, adjusted p-values, and supporting cell counts.

## Main Findings

- Thirty-two sample-by-cell-type models produce 960 KO-versus-WT contrasts.
- Of these contrasts, 224 have a Dunnett-adjusted p-value below 0.05.

## Scope and Data Availability

The analysis uses the combined QC Seurat object and final
[cell-type annotations](../dracc_cell_type_labeling_knn/outputs/02_label_dracc_cells_knn/cell_type_annotations.csv).
The annotations and complete contrast table are included. The combined Seurat
object and fitted model objects are not distributed.
