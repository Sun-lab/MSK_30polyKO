# MSK Cell-Selection Workflow

This document records the available MSK workflow from FASTQ processing to the
132,087-cell reference metadata set (`df.meta.rds`).

Status is `documented` when supplied source material states the step,
`reconstructed` when the rule was recovered from supplied metadata, and
`unavailable` when the required source was not provided.

## Workflow

1. **FASTQ and library inputs (`unavailable`)**

   The FASTQ manifests and per-sample Cell Ranger library CSV files were not
   supplied. The Cell Ranger command receives a library CSV as `--libraries`
   and a feature reference named `ref_feature_larryBC.csv`, indicating gene
   expression plus geneBC/LARRY feature-barcode libraries.

2. **Per-sample Cell Ranger count (`documented`)**

   MSK used Cell Ranger 9.0.1 with
   `refdata-gex-GRCh38-2024-A` and:

   ```text
   cellranger count \
     --id=<sample> \
     --transcriptome=<refdata-gex-GRCh38-2024-A> \
     --libraries=<library CSV> \
     --feature-ref=<ref_feature_larryBC.csv> \
     --create-bam=false \
     --nosecondary \
     --output-dir=<output directory> \
     --disable-ui
   ```

   The Cell Ranger aggregation report identifies DE_NZ as Single Cell 3' v4
   and the other eight inputs as Single Cell 3' v3.

3. **Cell Ranger aggregation (`documented`)**

   The nine Cell Ranger outputs were combined in a run named `aggrall`. The
   supplied report records Cell Ranger 9.0.1, GRCh38, and 260,440 estimated
   cells before downstream MSK QC.

   The exact `cellranger aggr` command and aggregation CSV are `unavailable`.

4. **Transcriptome metadata and doublet calling (`documented`)**

   MSK read the Cell Ranger `filtered_feature_bc_matrix`, selected the
   `Gene Expression` assay, and created a Seurat object with `min.cells = 3`
   and `min.features = 0`. Mitochondrial and ribosomal percentages used
   `^MT-` and `^RP[SL]`, respectively.

   scDblFinder 1.20.2 was run on cells with `nCount_RNA >= 200`. Its score and
   class were joined back to the full Seurat metadata by cell barcode.

5. **Transcriptome QC (`documented`)**

   A cell has `Qualified == "Yes"` only when all four conditions hold:

   ```text
   nCount_RNA >= 200
   scDblFinder.class == "singlet"
   nFeature_RNA >= 100
   percent_mito <= 20
   ```

6. **geneBC and LARRY assignment (`documented`)**

   Each modality was classified from its largest and second-largest feature
   UMI counts:

   ```text
   Type 1: maximum UMI >= 5 and second/maximum <= 0.5
   Type 2: maximum UMI >= 5 and second/maximum > 0.5
   Type 3: maximum UMI < 5
   ```

   The supplied metadata contains the resulting top barcode, top UMI,
   category, and barcode type. The original R source and second-largest UMI
   values are `unavailable`, so the ratio calculation cannot be rerun.

7. **Final cell selection (`reconstructed`)**

   Transcriptome-QC and barcode metadata are joined by the 10x cell barcode.
   A retained cell must:

   ```text
   Qualified == "Yes"
   larryBCCategory is present and is not the literal string "NA"
   geneBCCategory is present and is not the literal string "NA"
   ```

   It must then satisfy either:

   ```text
   KO cell:
     geneBCCategory == "Type 1: Max >=5 & 2nd/Max <= 0.5"

   WT cell:
     geneBCCategory == "Type 3: Max < 5"
     larryBCCategory == "Type 1: Max >=5 & 2nd/Max <= 0.5"
     larryBC_type == "WT"
   ```

   The geneBC rule has precedence: a Type-1 geneBC cell is classified as a KO
   cell even when its LARRY barcode is labeled WT.

8. **Cell identifiers and final sample labels (`reconstructed`)**

   The `-1` suffix is removed from the 10x barcode and the final sample label
   is appended. Barcode overlap gives this mapping:

   | Source stage | Final `orig.ident` |
   |---|---|
   | DE_NZ | DE_GemX |
   | DE_XM | DE |
   | ES | ES |
   | PP1 | PP1 |
   | PP2 | PP2 |
   | S5_1 | S5_1 |
   | S5_2 | S5_2 |
   | S6_1 | S6_1 |
   | S6_2 | S6_2 |

   The DE mapping is based only on exact cell-barcode overlap and should not
   be interpreted as a technology annotation.

9. **Endpoint (`documented`)**

   Nan's supplied `30KO_scRNAseq/df.meta.rds` is byte-identical to the
   canonical MSK `df.meta.rds`. The reconstructed rule
   reproduces all nine sample-specific barcode sets and the combined 132,087
   cells exactly; see `02_validate_workflow.html`.

## Internal Source Material

These inputs were used for the audit but are not distributed in this
repository.

- `../external_sources/msk/0.alignment.larryBC.sbatch`
- `../external_sources/msk/1.QC.TX.DE_NZ.html`
- `../external_sources/msk/1.BC.assign.DE_NZ.html`
- `../external_sources/msk/30KO_scRNAseq/ReadMe.md`
- `../external_sources/msk/30KO_scRNAseq/1.aggrall.web_summary.html`
- `../external_sources/msk/30KO_scRNAseq/30KO_QC.pptx`
- `../external_sources/msk/30KO_scRNAseq/2.demulti_BC/`
- `../external_sources/msk/30KO_scRNAseq/3.TX_QC/`
- `../external_sources/msk/30KO_scRNAseq/df.meta.rds`
