library(Seurat)

stage_output_dir = "/fh/fast/sun_w/yub/github/MSK_30polyKO/dracc_qc/outputs/02_build_qc_stage"
stage_rds_dir = file.path(stage_output_dir, "stage_rds")
stage_summary_dir = file.path(stage_output_dir, "stage_summary")
output_dir = "/fh/fast/sun_w/yub/github/MSK_30polyKO/dracc_qc/outputs/03_combine_qc_stages"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stage_order = c(
  "ES", "DE", "DE_XM", "PP1", "PP2",
  "S5_1", "S5_2", "S6_1", "S6_2"
)
stage_rds_paths = file.path(stage_rds_dir, paste0(stage_order, ".rds"))
stage_summary_paths = file.path(
  stage_summary_dir,
  paste0(stage_order, ".csv")
)
missing_paths = c(
  stage_rds_paths[!file.exists(stage_rds_paths)],
  stage_summary_paths[!file.exists(stage_summary_paths)]
)
if (length(missing_paths) > 0L) {
  stop("Missing stage outputs:\n", paste(missing_paths, collapse = "\n"))
}

# 1. Load and validate the stage outputs
objects = lapply(stage_rds_paths, readRDS)
names(objects) = stage_order
stage_summary = do.call(
  rbind,
  lapply(
    stage_summary_paths,
    read.csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
)
stage_summary = stage_summary[
  match(stage_order, stage_summary$stage_id),
  c(
    "stage_id", "n_raw_cells", "n_initial_called", "n_singlets",
    "n_transcriptome_qc", "n_ko", "n_wt", "n_final", "n_features"
  ),
  drop = FALSE
]

feature_meta = objects[[1]][["RNA"]][[]]
stopifnot(
  identical(stage_summary$stage_id, stage_order),
  "gene_symbol" %in% colnames(feature_meta),
  all(vapply(
    objects[-1],
    \(x) identical(rownames(x), rownames(objects[[1]])),
    logical(1)
  )),
  all(vapply(
    objects[-1],
    \(x) identical(x[["RNA"]][[]]$gene_symbol, feature_meta$gene_symbol),
    logical(1)
  ))
)

all_cells = unlist(lapply(objects, colnames), use.names = FALSE)
if (anyDuplicated(all_cells) > 0L) {
  stop("Duplicated cell IDs across stage objects.")
}

# 2. Merge raw-count objects without normalization
combined = merge(
  objects[[1]],
  y = objects[-1],
  project = "DRACC_QC",
  merge.data = FALSE
)
feature_meta = feature_meta[rownames(combined), , drop = FALSE]
combined[["RNA"]] = SeuratObject::AddMetaData(
  combined[["RNA"]],
  metadata = feature_meta
)

required_meta = c(
  "orig.ident", "nCount_RNA", "nFeature_RNA", "stage_id",
  "dracc_barcode", "geneBC_type", "larryBC", "mt_pct"
)
stopifnot(
  ncol(combined) == sum(stage_summary$n_final),
  identical(colnames(combined), rownames(combined@meta.data)),
  all(required_meta %in% colnames(combined@meta.data)),
  !anyNA(combined$geneBC_type),
  !anyNA(combined$larryBC)
)

# 3. Add the overall cumulative summary and save
overall_summary = data.frame(
  stage_id = "All",
  n_raw_cells = sum(stage_summary$n_raw_cells),
  n_initial_called = sum(stage_summary$n_initial_called),
  n_singlets = sum(stage_summary$n_singlets),
  n_transcriptome_qc = sum(stage_summary$n_transcriptome_qc),
  n_ko = sum(stage_summary$n_ko),
  n_wt = sum(stage_summary$n_wt),
  n_final = sum(stage_summary$n_final),
  n_features = nrow(combined),
  stringsAsFactors = FALSE
)
final_summary = rbind(stage_summary, overall_summary)

object_path = file.path(output_dir, "dracc_qc_seurat.rds")
summary_path = file.path(output_dir, "dracc_qc_summary.csv")
saveRDS(combined, object_path)
write.csv(final_summary, summary_path, row.names = FALSE)

cat("Final cells: ", ncol(combined), "\n", sep = "")
cat("Final features: ", nrow(combined), "\n", sep = "")
cat("Wrote: ", object_path, "\n", sep = "")
