# Build a reduced MaAsLin results table for cluster metadata only.

input_file <- "Output/MAASLIN3/001_disc_cluster/001_cluster_gc/all_results.tsv"
output_file <- "Output/MAASLIN3/001_disc_cluster/001_cluster_gc/all_results_single_table.tsv"

all_results <- read.delim(
  file = input_file,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!all(c("feature", "metadata", "value", "coef", "model") %in% names(all_results))) {
  stop("Missing one or more required columns: feature, metadata, value, coef, model")
}

pval_col <- if ("p.val_individual" %in% names(all_results)) {
  "p.val_individual"
} else if ("pval_individual" %in% names(all_results)) {
  "pval_individual"
} else {
  stop("Missing p-value column: expected p.val_individual or pval_individual")
}

cluster_only <- all_results[tolower(all_results$metadata) == "cluster", ]

final_table <- data.frame(
  feature_model = paste(cluster_only$feature, cluster_only$model, sep = "_"),
  value = cluster_only$value,
  coef = cluster_only$coef,
  p.val_individual = cluster_only[[pval_col]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Remove rows with NA in coef or p.val_individual
final_table <- final_table[!is.na(final_table$coef) & !is.na(final_table$p.val_individual), ]

write.table(
  final_table,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Saved: ", output_file)
