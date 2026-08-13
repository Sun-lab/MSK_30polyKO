library(anndata)
library(tibble)

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript 01_extract_raw_count_cell_meta.R <stage_id>")
}

stage_id = args[[1]]
stage_order = c("ES", "DE", "DE_XM", "PP1", "PP2", "S5_1", "S5_2", "S6_1", "S6_2")
if (!stage_id %in% stage_order) {
  stop("Unknown DRACC stage: ", stage_id)
}

output_dir = "/fh/fast/sun_w/yub/github/MSK_30polyKO/dracc_qc/outputs/01_extract_raw_count_cell_meta"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dracc_dir = "/fh/fast/sun_w/MorPhiC/data/MorPhiC_internal_releases/MorPhiC_Release_June_2026/MSK_30polyKO_revised_DRACC_processed_v3/processed"
corrected_de_xm_dir = "/fh/fast/sun_w/MorPhiC/data/Hong_shared_globus_20260706/MSK_30KO_DE_XM_corrected_2026-07-06_preRelease"

# Keep the DE_XM prerelease override aligned with the singlet workflow.
counts_path_for_stage <- function(stage_id) {
  if (stage_id == "DE_XM") {
    return(file.path(corrected_de_xm_dir, "counts.h5ad"))
  }
  file.path(dracc_dir, paste0("30_KO_", stage_id), "counts.h5ad")
}

# 1. Load raw count cell metadata
sample_id = paste0("30_KO_", stage_id)
h5ad_path = counts_path_for_stage(stage_id)
if (!file.exists(h5ad_path)) {
  stop("Missing DRACC counts file: ", h5ad_path)
}

h5ad = anndata::read_h5ad(h5ad_path, backed = "r")
obs = as.data.frame(h5ad$obs, stringsAsFactors = FALSE)

raw_cell_ids = rownames(obs)
if (anyDuplicated(raw_cell_ids) > 0) {
  stop("Duplicated raw DRACC cell IDs in stage: ", stage_id)
}

# 2. Build tibble
reserved_cols = c("sample_id", "stage_id", "dracc_cell_id", "dracc_barcode")
obs = obs[, setdiff(colnames(obs), reserved_cols), drop = FALSE]

cell_meta = tibble::as_tibble(
  data.frame(
    sample_id = sample_id,
    stage_id = stage_id,
    dracc_cell_id = raw_cell_ids,
    dracc_barcode = substr(raw_cell_ids, 1, 16),
    obs,
    check.names = FALSE
  )
)

# 3. Save output
saveRDS(cell_meta, file.path(output_dir, paste0(stage_id, ".rds")))
