# MSK reference audit

This folder documents and validates the original MSK analysis that supplies
the reference cells and labels for the DRACC reanalysis. It is supporting
provenance work, not a replacement for the original MSK pipeline.

## Contents

- [`01_workflow.md`](01_workflow.md) reconstructs the available workflow from
  Cell Ranger processing through the final reference metadata and marks each
  step as documented, reconstructed, or unavailable.
- [`02_validate_workflow.Rmd`](02_validate_workflow.Rmd) and its
  [rendered report](02_validate_workflow.html) verify the supplied QC fields
  and reconstruct the final cell set.
- [`03_msk_eda.Rmd`](03_msk_eda.Rmd) and its
  [rendered report](03_msk_eda.html) summarize the integrated Seurat object,
  key metadata fields, and subcelltype UMAP.

## Main Findings

- The supplied and canonical final metadata objects are identical.
- The reconstructed transcriptome-QC and barcode-selection rules reproduce all
  nine sample-specific cell sets and all 132,087 final cells exactly.
- The integrated reference contains 132,087 cells, 35,365 RNA features, and 19
  subcelltype labels.

## Scope and Data Availability

The audit relies on internal MSK source material and reference objects that are
not distributed in this repository. The reports therefore preserve the exact
analysis and findings but are not independently runnable from public files
alone. See [`01_workflow.md`](01_workflow.md) for the available evidence,
unavailable steps, QC thresholds, barcode rules, and stage mapping.
