# MSK 30-polyKO DRACC reanalysis

This repository provides an auditable R workflow for the revised
DRACC-processed MSK 30-polyKO single-cell dataset, from input validation and
quality control through cell-type labeling and knockout-effect analysis.

## Repository Contents

| Folder | Purpose | Main public content |
|---|---|---|
| [`msk_audit`](msk_audit) | Reconstruct and validate the original MSK reference analysis. | Workflow reconstruction, source Rmds, and rendered reports. |
| [`dracc_audit`](dracc_audit) | Inspect the processed DRACC inputs and validate clone and genotype feature calls. | Source Rmds and rendered input-audit reports. |
| [`msk_dracc_consistency`](msk_dracc_consistency) | Compare MSK and DRACC cell selection, feature identities, and gene space. | Source analyses, rendered report, and compact comparison outputs. |
| [`dracc_qc`](dracc_qc) | Build the nine-sample, raw-count DRACC dataset used downstream. | QC scripts and the combined cell-count summary. |
| [`dracc_cell_type_labeling_knn`](dracc_cell_type_labeling_knn) | Transfer MSK cell-type labels to DRACC cells and validate the transfer. | Labeling scripts, PCA checkpoint, final annotations, and validation summary. |
| [`dracc_genotype_effects_on_cell_type_composition`](dracc_genotype_effects_on_cell_type_composition) | Test KO effects on cell-type abundance. | Clone-aware composition model and complete contrast table. |
| [`dracc_gene_ko_differential_expression_by_cell_type`](dracc_gene_ko_differential_expression_by_cell_type) | Test KO expression effects within stage and cell type. | Clone-aware DE model, result collector, and compact significant-result table. |

## Analysis Workflow

The main DRACC workflow is:

1. Inspect the inputs in [`dracc_audit`](dracc_audit).
2. Validate the transition from the MSK reference in
   [`msk_dracc_consistency`](msk_dracc_consistency).
3. Construct the filtered dataset in [`dracc_qc`](dracc_qc).
4. Assign and validate cell types in
   [`dracc_cell_type_labeling_knn`](dracc_cell_type_labeling_knn).
5. Run the two biological analyses:
   [`dracc_genotype_effects_on_cell_type_composition`](dracc_genotype_effects_on_cell_type_composition)
   and
   [`dracc_gene_ko_differential_expression_by_cell_type`](dracc_gene_ko_differential_expression_by_cell_type).

[`msk_audit`](msk_audit) documents the original reference analysis used by the
consistency and cell-labeling steps.

## Main Findings

- [QC summary](dracc_qc/outputs/03_combine_qc_stages/dracc_qc_summary.csv):
  121,444 cells across nine samples, including 116,237 KO and 5,207 WT cells.
- [Cell-type annotations](dracc_cell_type_labeling_knn/outputs/02_label_dracc_cells_knn/cell_type_annotations.csv):
  113,891 direct labels, 7,178 kNN labels, and 375 unassigned cells; accepted
  holdout predictions are 99.53% accurate.
- [Composition results](dracc_genotype_effects_on_cell_type_composition/outputs/01_genotype_effects_on_cell_type_composition/genotype_effects_on_cell_type_composition.csv):
  630 KO-versus-WT contrasts, including 178 with Dunnett-adjusted p-values
  below 0.05.
- [Differential-expression results](dracc_gene_ko_differential_expression_by_cell_type/outputs/02_collect_significant_differential_expression.csv):
  94,353 converged, Bonferroni-significant gene-level results from 311 eligible
  contrasts.

## Scope and Data Availability

The repository includes analysis code, rendered audit reports, compact result
tables, final cell annotations, and the PCA checkpoint used for label transfer.
Raw H5AD files, the combined Seurat object, internal source communications,
full differential-expression tables, fitted models, and cluster logs are not
distributed.

The scripts retain absolute input paths from the analysis environment for
provenance. To rerun them elsewhere, replace the input and output paths defined
near the top of each script. R package requirements are declared by the scripts;
the DRACC H5AD steps also require Python `anndata` through `reticulate`.
