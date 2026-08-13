# DRACC input audit

This folder inspects the processed DRACC inputs and verifies the clone and
genotype feature calls used by the downstream QC workflow.

## Contents

- [`01_list_processed_dracc_files.Rmd`](01_list_processed_dracc_files.Rmd) and
  its [rendered report](01_list_processed_dracc_files.html) inventory the main
  DRACC release plus the corrected DE_XM release.
- [`02_inspect_larry_feature_h5ad.Rmd`](02_inspect_larry_feature_h5ad.Rmd) and
  its [rendered report](02_inspect_larry_feature_h5ad.html) inspect the LARRY
  feature matrices and metadata across nine stages.
- [`03_inspect_gene_feature_h5ad.Rmd`](03_inspect_gene_feature_h5ad.Rmd) and
  its [rendered report](03_inspect_gene_feature_h5ad.html) inspect the geneBC
  feature matrices and metadata.
- [`04_check_clone_genotype_cell_meta_consistency.Rmd`](04_check_clone_genotype_cell_meta_consistency.Rmd)
  and its [rendered report](04_check_clone_genotype_cell_meta_consistency.html)
  compare decoded feature calls with the extracted count metadata.

## Main Findings

- All nine stage-level LARRY and gene-feature H5AD files were readable with
  the metadata fields required downstream.
- The corrected DE_XM release is used instead of DE_XM from the main release.
- Every nonempty clone and genotype call in count metadata matches the feature
  name decoded from its corresponding feature H5AD.

## Scope and Data Availability

The processed DRACC H5AD files and extracted count metadata are not distributed
in this repository. The rendered reports preserve the inspected schemas,
examples, and validation result, but the audit is not independently runnable
from public files alone. Tables are embedded in the reports; this folder does
not write a separate outputs directory.
