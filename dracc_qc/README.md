# DRACC quality control

This folder constructs the nine-stage, raw-count DRACC Seurat object used by
cell-type labeling and downstream knockout analyses.

## Contents

- [`01_extract_raw_count_cell_meta.R`](01_extract_raw_count_cell_meta.R)
  extracts count-file cell metadata for each stage.
- [`02_build_qc_stage.R`](02_build_qc_stage.R) applies the cumulative filters
  and constructs one raw-count Seurat object per stage.
- [`03_combine_qc_stages.R`](03_combine_qc_stages.R) verifies matching feature
  spaces, merges the nine stage objects, and produces the final QC summary.
- [`dracc_qc_summary.csv`](outputs/03_combine_qc_stages/dracc_qc_summary.csv)
  reports stage-level and combined cell counts.

## Main Findings

- Of 19,564,447 raw droplets, 260,228 are called cells, 212,386 are singlets,
  and 199,057 pass transcriptome QC.
- The final dataset contains 121,444 cells: 116,237 KO and 5,207 WT cells
  across nine stages.
- Every stage retains the same 38,606 Ensembl features; gene symbols are stored
  as feature metadata rather than replacing stable feature IDs.

## Scope and Data Availability

Filtering is cumulative: called cell, singlet, at least 200 RNA counts, at
least 100 detected genes, at most 20% mitochondrial reads, then the documented
geneBC/LARRY eligibility rules. The raw H5AD files, internal MSK WT-LARRY
mapping, stage objects, and combined Seurat object are not distributed. The
compact QC summary is included; large objects remain local workflow handoffs.
