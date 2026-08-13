suppressPackageStartupMessages({
  library(SeuratObject)
  library(nebula)
})

object_path = "/fh/fast/sun_w/yub/github/MSK_30polyKO/dracc_qc/outputs/03_combine_qc_stages/dracc_qc_seurat.rds"
annotation_path = "/fh/fast/sun_w/yub/github/MSK_30polyKO/dracc_cell_type_labeling_knn/outputs/02_label_dracc_cells_knn/cell_type_annotations.csv"
output_dir = "/fh/fast/sun_w/yub/github/MSK_30polyKO/dracc_gene_ko_differential_expression_by_cell_type/outputs/01_gene_ko_differential_expression_by_cell_type"

# Native DRACC sample IDs grouped as in the example analysis.
stage_samples = list(
  S0 = "ES",
  S1 = c("DE", "DE_XM"),
  S3 = "PP1",
  S4 = "PP2",
  S5 = c("S5_1", "S5_2"),
  S6 = c("S6_1", "S6_2")
)
min_ko_cells = 30L
min_wt_cells = 30L
nebula_cores = 4L

args = commandArgs(trailingOnly = TRUE)
if (length(args) > 1L || (length(args) == 1L && !args %in% names(stage_samples))) {
  stop("Optional argument must be one of: ", paste(names(stage_samples), collapse = ", "))
}
stages_to_run = if (length(args) == 1L) args else names(stage_samples)

safe_name <- function(x) {
  gsub("(^_+|_+$)", "", gsub("[^A-Za-z0-9_.-]+", "_", x))
}

# 1. Join the final cell-type annotations to the QC object
obj = readRDS(object_path)
annotations = read.csv(annotation_path, stringsAsFactors = FALSE)
meta = obj[[]]

stopifnot(
  inherits(obj, "Seurat"),
  identical(rownames(meta), colnames(obj)),
  "RNA" %in% Assays(obj),
  all(
    c("stage_id", "geneBC_type", "larryBC", "nCount_RNA") %in%
      colnames(meta)
  ),
  all(c("cell_id", "subcelltype") %in% colnames(annotations)),
  nrow(annotations) == ncol(obj),
  anyDuplicated(annotations$cell_id) == 0L
)
annotation_idx = match(colnames(obj), annotations$cell_id)
stopifnot(!anyNA(annotation_idx))
annotations = annotations[annotation_idx, , drop = FALSE]
stopifnot(identical(annotations$cell_id, colnames(obj)))

obj$subcelltype = annotations$subcelltype
obj = subset(obj, cells = colnames(obj)[obj$subcelltype != "unassigned"])
meta = obj[[]]
stopifnot(
  !anyNA(meta[c("stage_id", "subcelltype", "geneBC_type", "larryBC")]),
  all(nzchar(as.character(meta$subcelltype))),
  all(nzchar(as.character(meta$geneBC_type))),
  all(nzchar(as.character(meta$larryBC))),
  "WT" %in% meta$geneBC_type,
  all(unlist(stage_samples, use.names = FALSE) %in% meta$stage_id),
  all(is.finite(meta$nCount_RNA)),
  all(meta$nCount_RNA > 0)
)

# 2. Test each eligible KO-versus-WT contrast within stage and cell type
for (stage_name in stages_to_run) {
  sample_ids = stage_samples[[stage_name]]
  stage_meta = meta[meta$stage_id %in% sample_ids, , drop = FALSE]

  for (cell_type in sort(unique(stage_meta$subcelltype))) {
    cell_type_meta = stage_meta[
      stage_meta$subcelltype == cell_type,
      ,
      drop = FALSE
    ]
    ko_genotypes = sort(setdiff(unique(cell_type_meta$geneBC_type), "WT"))

    for (genotype in ko_genotypes) {
      ko_cells = rownames(cell_type_meta)[
        cell_type_meta$geneBC_type == genotype
      ]
      wt_cells = rownames(cell_type_meta)[
        cell_type_meta$geneBC_type == "WT"
      ]
      if (length(ko_cells) < min_ko_cells || length(wt_cells) < min_wt_cells) {
        next
      }

      obj_sub = subset(obj, cells = c(ko_cells, wt_cells))
      obj_sub$geneBC_type = factor(
        obj_sub$geneBC_type,
        levels = c("WT", genotype)
      )
      predictor_cols = "geneBC_type"
      use_sample_covariate = length(sample_ids) > 1L
      if (use_sample_covariate) {
        obj_sub$stage_id = factor(obj_sub$stage_id)
        predictor_cols = c("stage_id", predictor_cols)
      }

      nebula_input = scToNeb(
        obj = obj_sub,
        assay = "RNA",
        id = "larryBC",
        pred = predictor_cols,
        offset = "nCount_RNA"
      )
      design = if (use_sample_covariate) {
        model.matrix(~ stage_id + geneBC_type, data = nebula_input$pred)
      } else {
        model.matrix(~ geneBC_type, data = nebula_input$pred)
      }
      if (qr(design)$rank < ncol(design)) {
        stop(
          "Non-estimable sample/genotype design for ",
          stage_name, ", ", cell_type, ", ", genotype, " versus WT."
        )
      }

      grouped = group_cell(
        count = nebula_input$count,
        id = nebula_input$id,
        pred = design,
        offset = nebula_input$offset
      )
      fit = nebula(
        grouped$count,
        grouped$id,
        pred = grouped$pred,
        method = "HL",
        offset = grouped$offset,
        ncore = nebula_cores
      )

      genotype_term = grep("^geneBC_type", colnames(design), value = TRUE)
      stopifnot(length(genotype_term) == 1L)
      results = as.data.frame(fit$summary)
      log_fc_col = paste0("logFC_", genotype_term)
      p_value_col = paste0("p_", genotype_term)
      stopifnot(all(c(log_fc_col, p_value_col) %in% colnames(results)))

      results$stage = stage_name
      results$cell_type = cell_type
      results$genotype = genotype
      results$contrast = paste0(genotype, "_vs_WT")
      results$n_ko_cells = length(ko_cells)
      results$n_wt_cells = length(wt_cells)
      results$genotype_log_fc = results[[log_fc_col]]
      results$genotype_fold_change = exp(results$genotype_log_fc)
      results$genotype_p_value = results[[p_value_col]]
      results$genotype_p_value_bonferroni = p.adjust(
        results$genotype_p_value,
        method = "bonferroni"
      )
      if (length(fit$convergence) == nrow(results)) {
        results$convergence = fit$convergence
      }
      if (length(fit$algorithm) == nrow(results)) {
        results$algorithm = fit$algorithm
      }
      first_cols = c(
        "stage", "cell_type", "genotype", "contrast",
        "n_ko_cells", "n_wt_cells",
        "genotype_log_fc", "genotype_fold_change",
        "genotype_p_value", "genotype_p_value_bonferroni"
      )
      results = results[c(first_cols, setdiff(colnames(results), first_cols))]

      contrast_dir = file.path(
        output_dir,
        stage_name,
        safe_name(cell_type)
      )
      dir.create(contrast_dir, recursive = TRUE, showWarnings = FALSE)
      file_stem = paste(
        stage_name,
        safe_name(cell_type),
        safe_name(genotype),
        "vs_WT",
        sep = "_"
      )
      saveRDS(
        fit,
        file.path(
          contrast_dir,
          paste0(file_stem, "_differential_expression_model.rds")
        )
      )
      write.csv(
        results,
        file.path(
          contrast_dir,
          paste0(file_stem, "_differential_expression.csv")
        ),
        row.names = FALSE
      )
    }
  }
}
