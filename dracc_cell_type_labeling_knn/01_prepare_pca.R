library(Seurat)

output_path = "/fh/fast/sun_w/yub/github/MSK_30polyKO/dracc_cell_type_labeling_knn/outputs/01_prepare_pca.rds"
dracc_path = "/fh/fast/sun_w/yub/github/MSK_30polyKO/dracc_qc/outputs/03_combine_qc_stages/dracc_qc_seurat.rds"
msk_meta_path = "/fh/fast/sun_w/yub/data/msk_30ko/raw/df.meta.rds"
n_variable_features = 2000L
n_pcs = 30L

stage_map = c(
  ES = "ES",
  DE = "DE",
  DE_XM = "DE_GemX",
  PP1 = "PP1",
  PP2 = "PP2",
  S5_1 = "S5_1",
  S5_2 = "S5_2",
  S6_1 = "S6_1",
  S6_2 = "S6_2"
)

# 1. Align direct MSK labels to the final DRACC cells
obj = readRDS(dracc_path)
msk_meta = as.data.frame(readRDS(msk_meta_path), stringsAsFactors = FALSE)
stopifnot(
  identical(colnames(obj), rownames(obj@meta.data)),
  all(c("stage_id", "dracc_barcode") %in% colnames(obj@meta.data)),
  all(c("cellid", "subcelltype") %in% colnames(msk_meta))
)

msk_labels = unique(data.frame(
  cellid = as.character(msk_meta$cellid),
  subcelltype = as.character(msk_meta$subcelltype),
  stringsAsFactors = FALSE
))
msk_labels = msk_labels[
  !is.na(msk_labels$cellid) &
    !is.na(msk_labels$subcelltype) &
    nzchar(msk_labels$subcelltype),
  ,
  drop = FALSE
]
if (anyDuplicated(msk_labels$cellid) > 0L) {
  stop("MSK cell IDs are not unique after label selection.")
}

msk_stage = unname(stage_map[as.character(obj$stage_id)])
if (anyNA(msk_stage)) {
  stop(
    "Missing MSK stage mapping for: ",
    paste(sort(unique(obj$stage_id[is.na(msk_stage)])), collapse = ", ")
  )
}
label_lookup = setNames(msk_labels$subcelltype, msk_labels$cellid)
msk_cellid = paste0(as.character(obj$dracc_barcode), "-", msk_stage)
direct_subcelltype = unname(label_lookup[msk_cellid])
if (sum(!is.na(direct_subcelltype)) < 15L) {
  stop("Too few direct MSK labels for kNN annotation.")
}

# 2. Compute the reusable PCA matrix
DefaultAssay(obj) = "RNA"
obj = NormalizeData(obj, verbose = FALSE)
obj = FindVariableFeatures(
  obj,
  selection.method = "vst",
  nfeatures = n_variable_features,
  verbose = FALSE
)
obj = ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj = RunPCA(
  obj,
  features = VariableFeatures(obj),
  npcs = n_pcs,
  verbose = FALSE
)
pca = Embeddings(obj, "pca")[, seq_len(n_pcs), drop = FALSE]
stopifnot(identical(rownames(pca), colnames(obj)))

# 3. Save only what the kNN step needs
pca_input = list(
  pca = pca,
  cell_meta = data.frame(
    cell_id = colnames(obj),
    stage_id = as.character(obj$stage_id),
    direct_subcelltype = direct_subcelltype,
    stringsAsFactors = FALSE
  )
)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(pca_input, output_path)

cat("Cells: ", nrow(pca), "\n", sep = "")
cat("Direct labels: ", sum(!is.na(direct_subcelltype)), "\n", sep = "")
cat("Unlabeled cells: ", sum(is.na(direct_subcelltype)), "\n", sep = "")
cat("Wrote: ", output_path, "\n", sep = "")
