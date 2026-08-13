library(hdf5r)
library(SeuratObject)

# Read gene symbols without loading a DRACC count matrix.
read_dracc_genes <- function(path) {
  h5 = hdf5r::H5File$new(path, mode = "r")
  on.exit(h5$close_all())

  gene_group = h5[["var"]][["gene_symbols"]]
  categories = as.character(gene_group[["categories"]][])
  codes = as.integer(gene_group[["codes"]][])
  keep = !is.na(codes) & codes >= 0
  unique(categories[codes[keep] + 1])
}

# 1. Define inputs
output_dir = "/fh/fast/sun_w/yub/github/MSK_30polyKO/msk_dracc_consistency/outputs/02_compare_gene_space"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

msk_seurat_path = "/fh/fast/sun_w/yub/data/msk_30ko/raw/all.seu.merged.integrated.rds"
dracc_counts_dir = "/fh/fast/sun_w/MorPhiC/data/MorPhiC_internal_releases/MorPhiC_Release_June_2026/MSK_30polyKO_revised_DRACC_processed_v3/processed"
corrected_de_xm_dir = "/fh/fast/sun_w/MorPhiC/data/Hong_shared_globus_20260706/MSK_30KO_DE_XM_corrected_2026-07-06_preRelease"
corrected_de_xm_path = file.path(corrected_de_xm_dir, "counts.h5ad")

stage_map = data.frame(
  stage = c(
    "ES", "DE", "DE_GemX", "PP1", "PP2",
    "S5_1", "S5_2", "S6_1", "S6_2"
  ),
  dracc_source = c(
    "ES", "DE", "DE_XM", "PP1", "PP2",
    "S5_1", "S5_2", "S6_1", "S6_2"
  ),
  stringsAsFactors = FALSE
)
dracc_counts_paths = ifelse(
  stage_map$dracc_source == "DE_XM",
  corrected_de_xm_path,
  file.path(
    dracc_counts_dir,
    paste0("30_KO_", stage_map$dracc_source),
    "counts.h5ad"
  )
)
stopifnot(
  file.exists(msk_seurat_path),
  all(file.exists(dracc_counts_paths))
)

# 2. Read the gene lists
msk_seurat = readRDS(msk_seurat_path)
msk_genes = unique(rownames(msk_seurat[["RNA"]]))
rm(msk_seurat)
invisible(gc())

dracc_genes = lapply(dracc_counts_paths, read_dracc_genes)

# 3. Compare gene symbols by stage
gene_space = do.call(
  rbind,
  lapply(seq_len(nrow(stage_map)), \(i) {
    shared_genes = intersect(msk_genes, dracc_genes[[i]])

    data.frame(
      stage = stage_map$stage[[i]],
      n_msk = length(msk_genes),
      n_dracc = length(dracc_genes[[i]]),
      n_shared = length(shared_genes),
      shared_over_msk = 100 * length(shared_genes) / length(msk_genes),
      shared_over_dracc = 100 * length(shared_genes) /
        length(dracc_genes[[i]]),
      stringsAsFactors = FALSE
    )
  })
)

# 4. Save output
write.csv(
  gene_space,
  file.path(output_dir, "gene_space.csv"),
  row.names = FALSE
)
