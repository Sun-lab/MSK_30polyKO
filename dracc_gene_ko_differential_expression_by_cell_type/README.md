# Gene knockout differential expression by cell type

This folder tests KO-versus-WT expression changes within developmental stage
and DRACC cell type while accounting for LARRY clone structure.

## Contents

- [`01_gene_ko_differential_expression_by_cell_type.R`](01_gene_ko_differential_expression_by_cell_type.R)
  groups samples as S0 = ES, S1 = DE + DE_XM, S3 = PP1, S4 = PP2,
  S5 = S5_1 + S5_2, and S6 = S6_1 + S6_2. It excludes unassigned cells and
  requires at least 30 KO and 30 WT cells per stage-group, cell-type contrast.
  NEBULA negative-binomial gamma mixed models use LARRY clone as the subject,
  total RNA count as the offset, and a sample covariate for S1, S5, and S6;
  p-values are Bonferroni-adjusted across genes within each contrast.
- [`02_collect_significant_differential_expression.R`](02_collect_significant_differential_expression.R)
  collects converged results with Bonferroni-adjusted p-values below 0.05.
- [`02_collect_significant_differential_expression.csv`](outputs/02_collect_significant_differential_expression.csv)
  is the compact public result table.

## Main Findings

- The eligibility rule produces 311 KO-versus-WT contrasts across S0, S1, S3,
  S4, S5, and S6.
- The compact table contains 94,322 significant, converged gene-level results.

## Scope and Data Availability

The analysis uses the combined QC Seurat object and final
[cell-type annotations](../dracc_cell_type_labeling_knn/outputs/02_label_dracc_cells_knn/cell_type_annotations.csv).
The annotations and compact result table are included. The combined Seurat
object, 311 complete gene-level tables, and fitted models are not distributed;
the complete local outputs total about 1.8 GB.
