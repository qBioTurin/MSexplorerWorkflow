# Load shared helpers and prepare the output directory.
source("Settings/utilities.R")
output_folder <- "Output/DAS_ALPHA_MAASLIN3_clean/"
createFolder(output_folder)

# Load compositional bacterial abundances and retain samples with a valid treatment group.
bact_baselines_ds_abund <- readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
baselines_decB_table <- as.data.frame(abundances(bact_baselines_ds_abund, transform = "compositional"))
metaData <- read.csv("InputData/Metadata_MS.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
gc_patients <- metaData$id[metaData$gc_treatment %in% c("positive", "negative")]
selected_samples <- colnames(baselines_decB_table) %in% gc_patients
baselines_decB_table <- baselines_decB_table[, selected_samples, drop = FALSE]

# Calculate the Shannon diversity contribution using a shared sample denominator.
shannon_index <- function(x, denominator) {
  x <- x[!is.na(x) & x > 0]
  p <- x / denominator
  -sum(p * log(p))
}

# Load the four anatomical taxon sets for the OO1 subset.
lesion <- readRDS("Output/BOI/001/Bacteria_lesion_001.rds")
spinal <- readRDS("Output/BOI/001/Bacteria_spinal_cord_001.rds")
gado <- readRDS("Output/BOI/001/Bacteria_gadolinium_001.rds")
sub <- readRDS("Output/BOI/001/Bacteria_subtentorial_lesions_001.rds")

# Calculate each anatomical group's contribution to sample-level Shannon diversity.
createTab <- function(lesion, spinal, gado, sub) {
  group_taxa <- list(
    Lesion = rownames(lesion@tax_table),
    spinal_Cord = rownames(spinal@tax_table),
    Gadolinium = rownames(gado@tax_table),
    Subtentorial = rownames(sub@tax_table)
  )
  group_tables <- lapply(group_taxa, function(taxa) {
    baselines_decB_table[rownames(baselines_decB_table) %in% taxa, , drop = FALSE]
  })

  shannon_rows <- as.data.frame(matrix(
    0,
    nrow = length(group_tables),
    ncol = ncol(baselines_decB_table),
    dimnames = list(names(group_tables), colnames(baselines_decB_table))
  ))

  for (i in seq_len(ncol(shannon_rows))) {
    denominator <- sum(vapply(group_tables, function(group) sum(group[, i], na.rm = TRUE), numeric(1)))
    for (j in seq_along(group_tables)) {
      shannon_rows[j, i] <- shannon_index(group_tables[[j]][, i], denominator)
    }
  }
  shannon_rows
}

# Add the original analysis labels and arrange samples as rows.
tabMod <- function(tab, Alpha, Method, Subset) {
  tab$Alpha <- Alpha
  tab$Method <- Method
  tab$Discriminant <- rownames(tab)
  tab$Subset <- Subset
  tab <- t(tab)
  tab <- cbind(id = rownames(tab), tab)
  as.data.frame(tab)
}

# Export the Shannon table for the OO1 subset only.
shannon_table <- createTab(lesion, spinal, gado, sub)
alpha_table <- tabMod(shannon_table, "Shannon", "Both", "OO1")
write.table(alpha_table, file = paste0(output_folder, "merged_MSHD_alpha.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)