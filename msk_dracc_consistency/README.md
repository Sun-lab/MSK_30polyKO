# MSK–DRACC consistency

This folder measures how closely DRACC processing reproduces the original MSK
cell selection, feature identities, and gene space.

## Contents

- [`01_compare_cell_selection.Rmd`](01_compare_cell_selection.Rmd) and its
  [rendered report](01_compare_cell_selection.html) compare cumulative cell
  selection, doublet calls, transcriptome QC, and final feature identities.
- [`02_compare_gene_space.R`](02_compare_gene_space.R) compares RNA features
  from the integrated MSK object with gene symbols in each DRACC count file.
- [`outputs/`](outputs) contains the comparison tables and diagnostic figures.

## Main Findings

- The two workflows share 252,999 initially called cells: 97.14% of MSK and
  97.22% of DRACC called cells.
- They share 105,945 final geneBC/LARRY-eligible cells: 95.82% of the MSK set
  and 87.24% of the DRACC set.
- Among shared final cells, genotype labels agree for 99.96% and LARRY clone
  IDs agree for 99.40%.
- The integrated MSK object has 35,365 RNA features; 35,348 are present in each
  DRACC gene space.

## Scope and Data Availability

Cell comparisons use the 16-base 10x barcode within an explicit stage mapping;
MSK-specific guide suffixes and WT clone wrappers are removed only for the
final identity comparison. The source MSK and DRACC objects are not distributed,
so the included reports and summary outputs document the comparison but are not
independently reproducible from public files alone.
