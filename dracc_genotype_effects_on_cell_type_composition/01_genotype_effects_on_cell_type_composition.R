suppressPackageStartupMessages({
  library(SeuratObject)
  library(blme)
  library(emmeans)
})

object_path = "/fh/fast/sun_w/yub/github/MSK_30polyKO/dracc_qc/outputs/03_combine_qc_stages/dracc_qc_seurat.rds"
annotation_path = "/fh/fast/sun_w/yub/github/MSK_30polyKO/dracc_cell_type_labeling_knn/outputs/02_label_dracc_cells_knn/cell_type_annotations.csv"
output_dir = "/fh/fast/sun_w/yub/github/MSK_30polyKO/dracc_genotype_effects_on_cell_type_composition/outputs/01_genotype_effects_on_cell_type_composition"
model_dir = file.path(output_dir, "models")

samples_to_run = c("PP1", "PP2", "S5_1", "S5_2", "S6_1", "S6_2")
min_cell_type_cells = 200L

# 1. Join the final cell-type annotations to the QC metadata
obj = readRDS(object_path)
annotations = read.csv(annotation_path, stringsAsFactors = FALSE)
meta = obj[[]]

stopifnot(
  inherits(obj, "Seurat"),
  identical(rownames(meta), colnames(obj)),
  all(c("stage_id", "geneBC_type", "larryBC") %in% colnames(meta)),
  all(c("cell_id", "subcelltype") %in% colnames(annotations)),
  nrow(annotations) == ncol(obj),
  anyDuplicated(annotations$cell_id) == 0L
)
annotation_idx = match(colnames(obj), annotations$cell_id)
stopifnot(!anyNA(annotation_idx))
annotations = annotations[annotation_idx, , drop = FALSE]
stopifnot(identical(annotations$cell_id, colnames(obj)))

analysis_meta = data.frame(
  sample_id = as.character(meta$stage_id),
  cell_type = as.character(annotations$subcelltype),
  genotype = as.character(meta$geneBC_type),
  clone_id = as.character(meta$larryBC),
  stringsAsFactors = FALSE
)
rm(obj, meta, annotations)
analysis_meta = analysis_meta[
  analysis_meta$cell_type != "unassigned",
  ,
  drop = FALSE
]
stopifnot(
  !anyNA(analysis_meta),
  all(nzchar(analysis_meta$cell_type)),
  all(nzchar(analysis_meta$genotype)),
  all(nzchar(analysis_meta$clone_id)),
  "WT" %in% analysis_meta$genotype,
  all(samples_to_run %in% analysis_meta$sample_id)
)

# 2. Test genotype effects within each sample and cell type
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
result_list = list()
result_index = 1L

for (sample_id in samples_to_run) {
  sample_meta = analysis_meta[
    analysis_meta$sample_id == sample_id,
    ,
    drop = FALSE
  ]
  cell_type_counts = table(sample_meta$cell_type)
  cell_types = names(cell_type_counts[cell_type_counts >= min_cell_type_cells])
  genotype_counts = table(sample_meta$genotype)
  sample_models = list()

  for (cell_type in cell_types) {
    message("Fitting ", sample_id, ": ", cell_type)
    model_data = data.frame(
      is_cell_type = as.integer(sample_meta$cell_type == cell_type),
      geneBC_type = relevel(factor(sample_meta$genotype), ref = "WT"),
      larryBC = factor(sample_meta$clone_id)
    )
    stopifnot(
      nlevels(model_data$geneBC_type) >= 2L,
      nlevels(model_data$larryBC) >= 2L
    )

    fixed_design = model.matrix(~ geneBC_type, data = model_data)
    fit = bglmer(
      is_cell_type ~ geneBC_type + (1 | larryBC),
      data = model_data,
      family = binomial,
      fixef.prior = normal(cov = diag(2.5^2, ncol(fixed_design))),
      control = lme4::glmerControl(
        optCtrl = list(maxfun = 100000)
      )
    )
    marginal_means = emmeans(fit, "geneBC_type", type = "response")
    result = as.data.frame(summary(
      contrast(marginal_means, "trt.vs.ctrl", ref = "WT"),
      infer = TRUE
    ))
    ko_genotypes = setdiff(levels(model_data$geneBC_type), "WT")
    stopifnot(nrow(result) == length(ko_genotypes))

    result$sample_id = sample_id
    result$cell_type = cell_type
    result$genotype = ko_genotypes
    result$n_cell_type = unname(cell_type_counts[[cell_type]])
    result$n_genotype = unname(genotype_counts[ko_genotypes])
    result$n_wt = unname(genotype_counts[["WT"]])
    first_cols = c(
      "sample_id", "cell_type", "genotype", "contrast",
      "n_cell_type", "n_genotype", "n_wt"
    )
    result = result[c(first_cols, setdiff(colnames(result), first_cols))]

    sample_models[[cell_type]] = fit
    result_list[[result_index]] = result
    result_index = result_index + 1L
  }

  saveRDS(
    sample_models,
    file.path(
      model_dir,
      paste0(sample_id, "_cell_type_composition_models.rds")
    )
  )
}

# 3. Save Dunnett-adjusted KO-versus-WT results
if (length(result_list) == 0L) {
  stop("No cell types met the minimum cell-count threshold.")
}
write.csv(
  do.call(rbind, result_list),
  file.path(output_dir, "genotype_effects_on_cell_type_composition.csv"),
  row.names = FALSE
)
