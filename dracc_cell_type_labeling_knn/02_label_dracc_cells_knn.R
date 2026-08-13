library(RANN)

pca_path = "/fh/fast/sun_w/yub/github/MSK_30polyKO/dracc_cell_type_labeling_knn/outputs/01_prepare_pca.rds"
output_dir = "/fh/fast/sun_w/yub/github/MSK_30polyKO/dracc_cell_type_labeling_knn/outputs/02_label_dracc_cells_knn"

# ponytail: one fixed configuration; reconsider only if validation is poor.
k = 15L
vote_threshold = 0.7
holdout_fraction = 0.2
set.seed(1)

# Predict labels by majority vote among global PCA neighbors.
predict_knn <- function(pca, query_idx, reference_idx, reference_labels) {
  if (length(query_idx) == 0L) {
    return(data.frame(
      row_idx = integer(),
      predicted_label = character(),
      vote_fraction = double(),
      stringsAsFactors = FALSE
    ))
  }
  if (length(reference_idx) < k) {
    stop("Fewer directly labeled reference cells than k.")
  }

  neighbor_idx = RANN::nn2(
    data = pca[reference_idx, , drop = FALSE],
    query = pca[query_idx, , drop = FALSE],
    k = k
  )$nn.idx
  neighbor_labels = matrix(
    reference_labels[reference_idx[neighbor_idx]],
    nrow = nrow(neighbor_idx),
    ncol = ncol(neighbor_idx)
  )
  votes = lapply(seq_len(nrow(neighbor_labels)), \(i) {
    counts = sort(table(neighbor_labels[i, ]), decreasing = TRUE)
    list(
      label = names(counts)[[1]],
      fraction = as.integer(counts[[1]]) / k
    )
  })

  data.frame(
    row_idx = query_idx,
    predicted_label = vapply(votes, `[[`, character(1), "label"),
    vote_fraction = vapply(votes, `[[`, numeric(1), "fraction"),
    stringsAsFactors = FALSE
  )
}

# Summarize accepted-label coverage and accuracy for one validation subset.
summarize_validation <- function(validation, stage_id) {
  accepted = validation$vote_fraction >= vote_threshold
  data.frame(
    stage_id = stage_id,
    n_holdout = nrow(validation),
    n_accepted = sum(accepted),
    coverage = mean(accepted),
    accuracy_accepted = if (any(accepted)) {
      mean(validation$predicted_label[accepted] == validation$true_label[accepted])
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

# 1. Load the compact PCA checkpoint
pca_input = readRDS(pca_path)
pca = pca_input$pca
cell_meta = pca_input$cell_meta
stopifnot(
  is.matrix(pca),
  identical(rownames(pca), cell_meta$cell_id),
  all(c("cell_id", "stage_id", "direct_subcelltype") %in% colnames(cell_meta))
)

direct_label = cell_meta$direct_subcelltype
reference_idx = which(!is.na(direct_label))
target_idx = which(is.na(direct_label))

# 2. Validate once, then label cells without direct overlap
validation_groups = split(
  reference_idx,
  paste(cell_meta$stage_id[reference_idx], direct_label[reference_idx], sep = "\r")
)
holdout_idx = unlist(lapply(validation_groups, \(idx) {
  if (length(idx) < 2L) {
    return(integer())
  }
  sample(idx, min(length(idx) - 1L, max(1L, floor(length(idx) * holdout_fraction))))
}), use.names = FALSE)
training_idx = setdiff(reference_idx, holdout_idx)

validation = predict_knn(pca, holdout_idx, training_idx, direct_label)
validation$true_label = direct_label[validation$row_idx]
validation$stage_id = cell_meta$stage_id[validation$row_idx]
validation_summary = rbind(
  summarize_validation(validation, "All"),
  do.call(rbind, lapply(
    split(validation, validation$stage_id),
    \(x) summarize_validation(x, unique(x$stage_id))
  ))
)

target_predictions = predict_knn(pca, target_idx, reference_idx, direct_label)
accepted = target_predictions$vote_fraction >= vote_threshold
accepted_idx = target_predictions$row_idx[accepted]

final_label = direct_label
label_source = ifelse(!is.na(direct_label), "MSK_direct", "unassigned")
knn_vote_fraction = rep(NA_real_, nrow(cell_meta))
knn_vote_fraction[target_predictions$row_idx] = target_predictions$vote_fraction
final_label[accepted_idx] = target_predictions$predicted_label[accepted]
label_source[accepted_idx] = "kNN"
final_label[is.na(final_label)] = "unassigned"

assignments = data.frame(
  cell_id = cell_meta$cell_id,
  subcelltype = final_label,
  subcelltype_source = label_source,
  subcelltype_knn_vote_fraction = knn_vote_fraction,
  stringsAsFactors = FALSE
)

# 3. Save the minimal annotation and validation tables
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  assignments,
  file.path(output_dir, "cell_type_annotations.csv"),
  row.names = FALSE
)
write.csv(
  validation_summary,
  file.path(output_dir, "knn_validation_summary.csv"),
  row.names = FALSE
)

print(validation_summary)
print(table(assignments$subcelltype_source))
