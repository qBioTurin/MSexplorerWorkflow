source("Settings/utilities.R")
output_folder <- "Output/DAS_ALPHA_MAASLIN3/"
createFolder(output_folder)
# abundance
bact_baselines_ds_abund <- readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
baselines_decB_table <- as.data.frame(abundances(bact_baselines_ds_abund, transform = "compositional"))
# setPatientIds
metaData <- read.csv("InputData/merged_data_cluster.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
metadava_revised <- metaData[!is.na(metaData$gc_treatment), ] %>%
  filter(gc_treatment == "positive" | gc_treatment == "negative") %>%
  select(id, gc_treatment)
GC_patient <- metadava_revised[metadava_revised$gc_treatment == "positive", "id"]
NO_GC_patient <- metadava_revised[metadava_revised$gc_treatment == "negative", "id"]

NT_HV_patient <- union(GC_patient, NO_GC_patient)
# setTaxid
shannon_index <- function(x, denom) {
  x <- x[!is.na(x) & x > 0]
  p <- x / denom
  -sum(p * log(p))
}

normalize_abundance <- function(abundance_table) {
  for (i in 1:ncol(abundance_table)) {
    col_sum <- sum(abundance_table[, i], na.rm = TRUE)
    if (col_sum != 0) { # Avoid division by zero
      abundance_table[, i] <- abundance_table[, i] / col_sum
    }
  }
  return(abundance_table) # Return the normalized table
}

createTab <- function(lesion, spinal, gado, sub, filtered_baselines_decB_table) {
  lesion <- rownames(lesion@tax_table)
  spinal_cord <- rownames(spinal@tax_table)
  gadolinium <- rownames(gado@tax_table)
  subtentorial <- rownames(sub@tax_table)

  filtered_baselines_decB_table <- baselines_decB_table[, colnames(baselines_decB_table) %in% NT_HV_patient, drop = FALSE]
  lesion_table <- filtered_baselines_decB_table[rownames(baselines_decB_table) %in% lesion, ]
  spinal_table <- filtered_baselines_decB_table[rownames(baselines_decB_table) %in% spinal_cord, ]
  gado_table <- filtered_baselines_decB_table[rownames(baselines_decB_table) %in% gadolinium, ]
  sub_table <- filtered_baselines_decB_table[rownames(baselines_decB_table) %in% subtentorial, ]


  shannon_rows <- as.data.frame(matrix(0, nrow = 4, ncol = ncol(filtered_baselines_decB_table)))
  colnames(shannon_rows) <- colnames(filtered_baselines_decB_table)
  rownames(shannon_rows) <- c("Lesion", "spinal_Cord", "Gadolinium", "Subtentorial")
  for (i in 1:ncol(shannon_rows)) {
    denom <- sum(lesion_table[, i]) + sum(spinal_table[, i]) + sum(gado_table[, i]) + sum(sub_table[, i])
    shannon_rows[nrow(shannon_rows) - 3, i] <- shannon_index(lesion_table[, i], denom)
    shannon_rows[nrow(shannon_rows) - 2, i] <- shannon_index(spinal_table[, i], denom)
    shannon_rows[nrow(shannon_rows) - 1, i] <- shannon_index(gado_table[, i], denom)
    shannon_rows[nrow(shannon_rows), i] <- shannon_index(sub_table[, i], denom)
  }
  return(list(shannon = shannon_rows ))
}

lesion_001 <- readRDS("Output/MAASLIN3_model/001_disc/001_lesion_burden_onlygc_plus_gc/lesion_burden_MAASLIN3_SR_features_001_T.rds")
bm_001 <- readRDS("Output/MAASLIN3_model/001_disc/001_spinal_cord_lesion_onlygc_plus_gc/spinal_cord_lesion_MAASLIN3_SR_features_001_T.rds")
gado_001 <- readRDS("Output/MAASLIN3_model/001_disc/001_gadolinium_contrast_onlygc_plus_gc/gadolinium_contrast_MAASLIN3_SR_features_001_T.rds")
sub_001 <- readRDS("Output/MAASLIN3_model/001_disc/001_subtentorial_lesions_onlygc_plus_gc/subtentorial_lesions_MAASLIN3_SR_features_001_T.rds")
tax <- as.data.frame(lesion_001@tax_table)
otu <- as.data.frame(lesion_001@otu_table)
metaData <- as.data.frame(as.matrix(sample_data(lesion_001)))

metaData <- sample_data(lesion_001) %>% as.data.frame()
print(lesion_001)
tab_list_001 <- createTab(lesion = lesion_001, spinal = bm_001, gado = gado_001, sub = sub_001, filtered_baselines_decB_table)

tabMod <- function(tab, Alpha, Method, Subset) {
  tab$Alpha <- Alpha
  tab$Method <- Method
  tab$Discriminant <- rownames(tab)
  tab$Subset <- Subset
  tab <- t(tab)
  tab <- cbind(id = rownames(tab), tab)
  tab <- as.data.frame(tab)
  return(tab)
}

constructTabSingle <- function(tab001, Alpha) {
  tab_001 <- tabMod(tab001, Alpha, "Both", "001")
  merged_sh <- t(tab_001)
  return(merged_sh)
}


Shannon <- constructTabSingle(tab_list_001$shannon, "Shannon")
Shannon <- t(Shannon)
merged_sh<-Shannon
write.table(merged_sh, file = paste0(output_folder, "merged_shannon_union.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

mshd_rds_candidates <- c(
  "Output/MAASLIN3_model/Bacteria_MSHD_disc_001/category_MAASLIN3_SR_features_K_T.rds",
  "Image/Rebuttal2/MSHD_Bacteria_heatmap.rds"
)
mshd_rds_path <- mshd_rds_candidates[file.exists(mshd_rds_candidates)][1]
if (is.na(mshd_rds_path)) {
  stop("Nessun file MSHD .rds trovato tra i path candidati.")
}
tab2 <- readRDS(mshd_rds_path)

createTabone <- function(tab1, dataset, output_mod2, filtered_baselines_decB_table) {
  taxa_names1 <- rownames(tab1@tax_table)


  filtered_baselines_decB_table <- baselines_decB_table
  filtered_baselines_decB_table1 <- filtered_baselines_decB_table[rownames(baselines_decB_table) %in% taxa_names1, ]


  shannon_rows <- as.data.frame(matrix(0, nrow = 1, ncol = ncol(filtered_baselines_decB_table)))
  colnames(shannon_rows) <- colnames(filtered_baselines_decB_table)
  rownames(shannon_rows) <- c("HD")


  for (i in 1:ncol(shannon_rows)) {
    denom <- sum(filtered_baselines_decB_table1[, i], na.rm = TRUE)
    shannon_rows[nrow(shannon_rows), i] <- shannon_index(filtered_baselines_decB_table1[, i], denom)
  }
  return(list(shannon = shannon_rows))
}

HD_001 <- createTabone(tab2, "001", output_folder, filtered_baselines_decB_table)

MSSh <- constructTabSingle(HD_001$shannon, "Shannon")
MSSh <- t(MSSh)
merged_sh <- as.data.frame(MSSh)
merged_sh <- t(merged_sh)
write.table(merged_sh, file = paste0(output_folder, "merged_MSHD_alpha.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

