library(anndata)
library(Matrix)
library(Seurat)

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript 02_build_qc_stage.R <stage_id>")
}

stage_id = args[[1]]
stage_order = c(
  "ES", "DE", "DE_XM", "PP1", "PP2",
  "S5_1", "S5_2", "S6_1", "S6_2"
)
if (!stage_id %in% stage_order) {
  stop("Unknown DRACC stage: ", stage_id)
}

output_dir = "/fh/fast/sun_w/yub/github/MSK_30polyKO/dracc_qc/outputs/02_build_qc_stage"
stage_rds_dir = file.path(output_dir, "stage_rds")
stage_summary_dir = file.path(output_dir, "stage_summary")
dir.create(stage_rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stage_summary_dir, recursive = TRUE, showWarnings = FALSE)

dracc_dir = "/fh/fast/sun_w/MorPhiC/data/MorPhiC_internal_releases/MorPhiC_Release_June_2026/MSK_30polyKO_revised_DRACC_processed_v3/processed"
corrected_de_xm_dir = "/fh/fast/sun_w/MorPhiC/data/Hong_shared_globus_20260706/MSK_30KO_DE_XM_corrected_2026-07-06_preRelease"
corrected_de_xm_path = file.path(corrected_de_xm_dir, "counts.h5ad")
msk_source_dir = "/fh/fast/sun_w/yub/github/MSK_30polyKO/external_sources/msk/30KO_scRNAseq"
msk_bc_dir = file.path(
  msk_source_dir,
  "2.demulti_BC"
)
msk_source = c(
  ES = "ES",
  DE = "DE_XM",
  DE_XM = "DE_NZ",
  PP1 = "PP1",
  PP2 = "PP2",
  S5_1 = "S5_1",
  S5_2 = "S5_2",
  S6_1 = "S6_1",
  S6_2 = "S6_2"
)
feature_suffix = c(
  ES = "es",
  DE = "de",
  DE_XM = "de_gemx",
  PP1 = "pp1",
  PP2 = "pp2",
  S5_1 = "s5_1",
  S5_2 = "s5_2",
  S6_1 = "s6_1",
  S6_2 = "s6_2"
)

# 1. Load the stage and its verified WT-LARRY mapping
h5ad_path = if (stage_id == "DE_XM") {
  corrected_de_xm_path
} else {
  file.path(dracc_dir, paste0("30_KO_", stage_id), "counts.h5ad")
}
msk_bc_path = file.path(
  msk_bc_dir,
  paste0("df.all.larryBC.geneBC.", msk_source[[stage_id]], ".rds")
)
stopifnot(file.exists(h5ad_path), file.exists(msk_bc_path))

msk_bc = readRDS(msk_bc_path)
stopifnot(all(c("larryBC", "larryBC_type") %in% colnames(msk_bc)))
wt_larry_labels = sort(unique(msk_bc$larryBC[msk_bc$larryBC_type == "WT"]))
stopifnot(
  length(wt_larry_labels) == 12L,
  all(grepl("^WTBC[0-9]+_BC[0-9]+$", wt_larry_labels))
)
wt_larry_barcodes = sort(sub("^WTBC[0-9]+_", "", wt_larry_labels))
stopifnot(
  !anyDuplicated(wt_larry_barcodes),
  !any(msk_bc$larryBC %in% wt_larry_barcodes)
)
rm(msk_bc, wt_larry_labels)

adata = anndata::read_h5ad(h5ad_path)
obs = as.data.frame(adata$obs, stringsAsFactors = FALSE)
var = as.data.frame(adata$var, stringsAsFactors = FALSE)
raw_cell_ids = rownames(obs)
ensembl_ids = rownames(var)
gene_symbols = as.character(var$gene_symbols)

guide_prefix = paste0(
  "CRISPR_Guide_Capture_grna_",
  feature_suffix[[stage_id]]
)
larry_prefix = paste0("Custom_larry_", feature_suffix[[stage_id]])
guide_fields = c(
  "__num_umis",
  "__feature_call",
  "__feature1_count",
  "__feature2_count"
)
larry_fields = c(
  "__feature_call",
  "__feature1_count",
  "__feature2_count"
)
required_obs = c(
  "is_cell", "singlet", "total_counts", "n_genes", "mt_pct",
  paste0(guide_prefix, guide_fields),
  paste0(larry_prefix, larry_fields)
)
stopifnot(
  all(required_obs %in% colnames(obs)),
  "gene_symbols" %in% colnames(var),
  !anyNA(raw_cell_ids),
  !anyDuplicated(raw_cell_ids),
  !anyNA(ensembl_ids),
  !anyDuplicated(ensembl_ids),
  !anyNA(gene_symbols),
  is.logical(obs$is_cell),
  is.logical(obs$singlet),
  !anyNA(obs$is_cell),
  !anyNA(obs$singlet),
  is.numeric(obs$total_counts),
  is.numeric(obs$n_genes),
  is.numeric(obs$mt_pct),
  all(is.finite(obs$total_counts)),
  all(is.finite(obs$n_genes)),
  all(is.finite(obs$mt_pct))
)

# 2. Apply the agreed cumulative transcriptome and feature rules
pass_called = obs$is_cell
pass_singlet = pass_called & obs$singlet
pass_transcriptome_qc = (
  pass_singlet &
    obs$total_counts >= 200 &
    obs$n_genes >= 100 &
    obs$mt_pct <= 20
)

guide_call = as.character(obs[[paste0(guide_prefix, "__feature_call")]])
guide_total = obs[[paste0(guide_prefix, "__num_umis")]]
guide_top = obs[[paste0(guide_prefix, "__feature1_count")]]
guide_second = obs[[paste0(guide_prefix, "__feature2_count")]]
larry_call = as.character(obs[[paste0(larry_prefix, "__feature_call")]])
larry_top = obs[[paste0(larry_prefix, "__feature1_count")]]
larry_second = obs[[paste0(larry_prefix, "__feature2_count")]]

stopifnot(
  !anyNA(guide_call),
  !anyNA(guide_total),
  !anyNA(guide_top),
  !anyNA(guide_second),
  !anyNA(larry_call),
  !anyNA(larry_top),
  !anyNA(larry_second),
  all(wt_larry_barcodes %in% levels(
    obs[[paste0(larry_prefix, "__feature_call")]]
  ))
)

guide_type_1 = (
  guide_top >= 5 &
    guide_second <= 0.5 * guide_top
)
guide_type_3 = guide_total > 0 & guide_top < 5
larry_type_1 = (
  larry_top >= 5 &
    larry_second <= 0.5 * larry_top
)
stopifnot(
  !any(guide_type_1 & guide_call == ""),
  !any(larry_type_1 & larry_call == "")
)

larry_is_wt = larry_call %in% wt_larry_barcodes
is_ko = guide_type_1 & larry_type_1 & !larry_is_wt
is_wt = guide_type_3 & larry_type_1 & larry_is_wt
pass_final = pass_transcriptome_qc & (is_ko | is_wt)

geneBC_type = rep(NA_character_, nrow(obs))
geneBC_type[is_ko] = guide_call[is_ko]
geneBC_type[is_wt] = "WT"
stopifnot(
  sum(pass_final) > 0L,
  any(pass_transcriptome_qc & is_wt),
  !anyNA(geneBC_type[pass_final]),
  !any(larry_call[pass_final] == "")
)

# 3. Build the stage Seurat object from raw counts
counts = adata$X[pass_final, , drop = FALSE]
counts = as(Matrix::t(counts), "CsparseMatrix")
rownames(counts) = ensembl_ids

sample_id = paste0("30_KO_", stage_id)
cell_names = paste(sample_id, raw_cell_ids[pass_final], sep = ":")
colnames(counts) = cell_names
meta = data.frame(
  stage_id = stage_id,
  dracc_barcode = substr(raw_cell_ids[pass_final], 1, 16),
  geneBC_type = geneBC_type[pass_final],
  larryBC = larry_call[pass_final],
  mt_pct = obs$mt_pct[pass_final],
  row.names = cell_names,
  stringsAsFactors = FALSE
)
stopifnot(
  identical(colnames(counts), rownames(meta)),
  !anyDuplicated(meta$dracc_barcode)
)

obj = Seurat::CreateSeuratObject(
  counts = counts,
  meta.data = meta,
  assay = "RNA",
  project = sample_id
)
obj[["RNA"]] = SeuratObject::AddMetaData(
  obj[["RNA"]],
  metadata = data.frame(
    gene_symbol = gene_symbols,
    row.names = ensembl_ids,
    stringsAsFactors = FALSE
  )
)
stopifnot(
  nrow(obj) == length(ensembl_ids),
  ncol(obj) == sum(pass_final),
  identical(colnames(obj), rownames(obj@meta.data)),
  all(obj$nCount_RNA == obs$total_counts[pass_final]),
  all(obj$nFeature_RNA == obs$n_genes[pass_final])
)

# 4. Save the stage object and cumulative counts
stage_summary = data.frame(
  stage_id = stage_id,
  n_raw_cells = nrow(obs),
  n_initial_called = sum(pass_called),
  n_singlets = sum(pass_singlet),
  n_transcriptome_qc = sum(pass_transcriptome_qc),
  n_ko = sum(pass_transcriptome_qc & is_ko),
  n_wt = sum(pass_transcriptome_qc & is_wt),
  n_final = sum(pass_final),
  n_features = nrow(obj),
  stringsAsFactors = FALSE
)

saveRDS(obj, file.path(stage_rds_dir, paste0(stage_id, ".rds")))
write.csv(
  stage_summary,
  file.path(stage_summary_dir, paste0(stage_id, ".csv")),
  row.names = FALSE
)
print(stage_summary)
