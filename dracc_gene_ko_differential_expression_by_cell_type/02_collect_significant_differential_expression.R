input_dir = "/fh/fast/sun_w/yub/github/MSK_30polyKO/dracc_gene_ko_differential_expression_by_cell_type/outputs/01_gene_ko_differential_expression_by_cell_type"
output_path = "/fh/fast/sun_w/yub/github/MSK_30polyKO/dracc_gene_ko_differential_expression_by_cell_type/outputs/02_collect_significant_differential_expression.csv"
significance_threshold = 0.05

# 1. Collect converged, Bonferroni-significant genotype effects
input_files = list.files(
  input_dir,
  pattern = "_vs_WT_differential_expression\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
stopifnot(length(input_files) > 0L)

output_columns = c(
  "stage", "cell_type", "genotype", "contrast",
  "n_ko_cells", "n_wt_cells", "gene_id", "gene",
  "genotype_log_fc", "genotype_fold_change",
  "genotype_p_value", "genotype_p_value_bonferroni",
  "convergence", "algorithm"
)

significant_results = lapply(input_files, \(input_file) {
  results = read.csv(input_file, check.names = FALSE)
  stopifnot(all(output_columns %in% colnames(results)))

  keep = !is.na(results$convergence) &
    results$convergence == 1L &
    !is.na(results$genotype_p_value_bonferroni) &
    results$genotype_p_value_bonferroni < significance_threshold

  results[keep, output_columns, drop = FALSE]
})
significant_results = do.call(rbind, significant_results)
significant_results = significant_results[
  order(
    significant_results$stage,
    significant_results$cell_type,
    significant_results$genotype,
    significant_results$genotype_p_value_bonferroni,
    significant_results$gene
  ),
  ,
  drop = FALSE
]
row.names(significant_results) = NULL

# 2. Save the compact public result table
write.csv(significant_results, output_path, row.names = FALSE)
message("Saved ", nrow(significant_results), " results to ", output_path)
