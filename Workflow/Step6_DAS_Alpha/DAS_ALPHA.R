source("Settings/utilities.R")
output_folder <- "Output/DAS_ALPHA_ANCOMBC2/"
createFolder(output_folder)
# abundance
bact_baselines_ds_abund <- readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
baselines_decB_table <- as.data.frame(abundances(bact_baselines_ds_abund, transform = "compositional"))
# setPatientIds
metaData <- read.csv("InputData/Metadata_MS.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
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

# simpson_index <- function(x) {
#   x <- x[!is.na(x) & x > 0]
#   p <- x / sum(x)
#   1 - sum(p^2)
# }

# EH_index <- function(x, denom) {
#   shannon_index(x, denom) / log(length(x[x > 0]))
# }

# observed_index <- function(x) {
#   sum(x > 0, na.rm = TRUE)
# }

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

  # simpson_rows <- as.data.frame(matrix(0, nrow = 4, ncol = ncol(filtered_baselines_decB_table)))
  # colnames(simpson_rows) <- colnames(filtered_baselines_decB_table)
  # rownames(simpson_rows) <- c("Lesion", "Spinal_Cord", "Gadolinium", "Subtentorial")

  # EH_rows <- as.data.frame(matrix(0, nrow = 4, ncol = ncol(filtered_baselines_decB_table)))
  # colnames(EH_rows) <- colnames(filtered_baselines_decB_table)
  # rownames(EH_rows) <- c("Lesion", "Spinal_Cord", "Gadolinium", "Subtentorial")

  # observed_rows <- as.data.frame(matrix(0, nrow = 4, ncol = ncol(filtered_baselines_decB_table)))
  # colnames(observed_rows) <- colnames(filtered_baselines_decB_table)
  # rownames(observed_rows) <- c("Lesion", "Spinal_Cord", "Gadolinium", "Subtentorial")
  for (i in 1:ncol(shannon_rows)) {
    denom <- sum(lesion_table[, i]) + sum(spinal_table[, i]) + sum(gado_table[, i]) + sum(sub_table[, i])
    shannon_rows[nrow(shannon_rows) - 3, i] <- shannon_index(lesion_table[, i], denom)
    shannon_rows[nrow(shannon_rows) - 2, i] <- shannon_index(spinal_table[, i], denom)
    shannon_rows[nrow(shannon_rows) - 1, i] <- shannon_index(gado_table[, i], denom)
    shannon_rows[nrow(shannon_rows), i] <- shannon_index(sub_table[, i], denom)

    # simpson_rows[nrow(simpson_rows) - 3, i] <- simpson_index(lesion_table[, i])
    # simpson_rows[nrow(simpson_rows) - 2, i] <- simpson_index(spinal_table[, i])
    # simpson_rows[nrow(simpson_rows) - 1, i] <- simpson_index(gado_table[, i])
    # simpson_rows[nrow(simpson_rows), i] <- simpson_index(sub_table[, i])

    # EH_rows[nrow(EH_rows) - 3, i] <- EH_index(lesion_table[, i], denom)
    # EH_rows[nrow(EH_rows) - 2, i] <- EH_index(spinal_table[, i], denom)
    # EH_rows[nrow(EH_rows) - 1, i] <- EH_index(gado_table[, i], denom)
    # EH_rows[nrow(EH_rows), i] <- EH_index(sub_table[, i], denom)

    # observed_rows[1, i] <- observed_index(lesion_table[, i])
    # observed_rows[2, i] <- observed_index(spinal_table[, i])
    # observed_rows[3, i] <- observed_index(gado_table[, i])
    # observed_rows[4, i] <- observed_index(sub_table[, i])
  }
  return(list(shannon = shannon_rows ))#, simpson = simpson_rows, EH = EH_rows, observed = observed_rows))
}


print(lesion_05)
tax<-as.data.frame(lesion_05@tax_table)
otu<-as.data.frame(lesion_05@otu_table)
metaData <- as.data.frame(as.matrix(sample_data(lesion_05)))

# lesion_05 <- readRDS("Output/BOI/05/Bacteria_lesion_05.rds")
# bm_05 <- readRDS("Output/BOI/05/Bacteria_spinal_cord_05.rds")
# gado_05 <- readRDS("Output/BOI/05/Bacteria_gadolinium_05.rds")
# sub_05 <- readRDS("Output/BOI/05/Bacteria_subtentorial_lesions_05.rds")

# lesion_01 <- readRDS("Output/BOI/01/Bacteria_lesion_01.rds")
# bm_01 <- readRDS("Output/BOI/01/Bacteria_spinal_cord_01.rds")
# gado_01 <- readRDS("Output/BOI/01/Bacteria_gadolinium_01.rds")
# sub_01 <- readRDS("Output/BOI/01/Bacteria_subtentorial_lesions_01.rds")

# lesion_001 <- readRDS("Output/BOI/001/Bacteria_lesion_001.rds")
# bm_001 <- readRDS("Output/BOI/001/Bacteria_spinal_cord_001.rds")
# gado_001 <- readRDS("Output/BOI/001/Bacteria_gadolinium_001.rds")
# sub_001 <- readRDS("Output/BOI/001/Bacteria_subtentorial_lesions_001.rds")

# lesion_05 <- readRDS("Output/BOI2/05_unioned/Bacteria_subtentorial_05.rds")
# bm_05 <- readRDS("Output/BOI2/05_unioned/Bacteria_subtentorial_05.rds")
# gado_05 <- readRDS("Output/BOI2/05_unioned/Bacteria_subtentorial_05.rds")
# sub_05 <- readRDS("Output/BOI2/05_unioned/Bacteria_subtentorial_05.rds")
# lesion_01 <- readRDS("Output/BOI2/01_unioned/Bacteria_lesion_01.rds")
# bm_01 <- readRDS("Output/BOI2/01_unioned/Bacteria_spinal_cord_01.rds")
# gado_01 <- readRDS("Output/BOI2/01_unioned/Bacteria_gadolinium_01.rds")
# sub_01 <- readRDS("Output/BOI2/01_unioned/Bacteria_subtentorial_01.rds")
# lesion_001 <- readRDS("Output/BOI2/001_unioned/Bacteria_lesion_001.rds")
# bm_001 <- readRDS("Output/BOI2/001_unioned/Bacteria_spinal_cord_001.rds")
# gado_001 <- readRDS("Output/BOI2/001_unioned/Bacteria_gadolinium_001.rds")
# sub_001 <- readRDS("Output/BOI2/001_unioned/Bacteria_subtentorial_001.rds")
lesion_01 <- readRDS("Output/BOI_AncomBC2/01/Bacteria_lesion_burden_01.rds")
bm_01 <- readRDS("Output/BOI_AncomBC2/01/Bacteria_spinal_cord_lesion_01.rds")
gado_01 <- readRDS("Output/BOI_AncomBC2/01/Bacteria_gadolinium_contrast_01.rds")
sub_01 <- readRDS("Output/BOI_AncomBC2/01/Bacteria_subtentorial_lesions_01.rds")
lesion_001 <- readRDS("Output/BOI_AncomBC2/001/Bacteria_lesion_burden_001.rds")
bm_001 <- readRDS("Output/BOI_AncomBC2/001/Bacteria_spinal_cord_lesion_001.rds")
gado_001 <- readRDS("Output/BOI_AncomBC2/001/Bacteria_gadolinium_contrast_001.rds")
sub_001 <- readRDS("Output/BOI_AncomBC2/001/Bacteria_subtentorial_lesions_001.rds")



tab_list_05 <- createTab(lesion = lesion_05,
 spinal = bm_05, gado = gado_05, sub = sub_05, filtered_baselines_decB_table)
tab_list_01 <- createTab(lesion = lesion_01, spinal = bm_01, gado = gado_01, sub = sub_01, filtered_baselines_decB_table)
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

constructTab <- function(tab05, tab01, tab001, Alpha) {
  tab_05 <- tabMod(tab05, Alpha, "Both", "O5")
  tab_01 <- tabMod(tab01, Alpha, "Both", "O1")
  tab_001 <- tabMod(tab001, Alpha, "Both", "OO1")

  merged_sh <- left_join(tab_05, tab_01, by = "id") %>%
    left_join(tab_001, by = "id")

  merged_sh <- t(merged_sh)

  return(merged_sh)
}
constructTab2 <- function(tab05, tab001, Alpha) {
  tab_05 <- tabMod(tab05, Alpha, "Both", "05")
  tab_001 <- tabMod(tab001, Alpha, "Both", "001")

  merged_sh <- left_join(tab_05, tab_001, by = "id")

  merged_sh <- t(merged_sh)

  return(merged_sh)
}


Shannon <- constructTab(tab_list_05$shannon, tab_list_01$shannon, tab_list_001$shannon, "Shannon")
Shannon <- t(Shannon)
# Simpson <- constructTab(tab_list_05$simpson, tab_list_01$simpson, tab_list_001$simpson, "Simpson")
# Simpson <- t(Simpson)
# EH <- constructTab(tab_list_05$EH, tab_list_01$EH, tab_list_001$EH, "EH")
# EH <- t(EH)
# Observed <- constructTab(tab_list_05$observed, tab_list_01$observed, tab_list_001$observed, "Observed")
# Observed <- t(Observed)
merged_sh<-Shannon
# merged_sh <- left_join(as.data.frame(Shannon), as.data.frame(Simpson), by = "id") %>%
#   left_join(as.data.frame(EH), by = "id") %>%
#   left_join(as.data.frame(Observed), by = "id")
write.table(merged_sh, file = paste0(output_folder, "merged_shannon_union.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

tab1 <- readRDS("Output/merge_DAS/MSHD/Bacteria_MsHd_05_merged.rds")
tab2 <- readRDS("Output/merge_DAS/MSHD/Bacteria_MsHd_merged.rds")

createTabone <- function(tab1, dataset, output_mod2, filtered_baselines_decB_table) {
  taxa_names1 <- rownames(tab1@tax_table)


  filtered_baselines_decB_table <- baselines_decB_table
  filtered_baselines_decB_table1 <- filtered_baselines_decB_table[rownames(baselines_decB_table) %in% taxa_names1, ]


  shannon_rows <- as.data.frame(matrix(0, nrow = 1, ncol = ncol(filtered_baselines_decB_table)))
  colnames(shannon_rows) <- colnames(filtered_baselines_decB_table)
  rownames(shannon_rows) <- c("HD")

  simpson_rows <- as.data.frame(matrix(0, nrow = 1, ncol = ncol(filtered_baselines_decB_table)))
  colnames(simpson_rows) <- colnames(filtered_baselines_decB_table)
  rownames(simpson_rows) <- c("HD")

  EH_rows <- as.data.frame(matrix(0, nrow = 1, ncol = ncol(filtered_baselines_decB_table)))
  colnames(EH_rows) <- colnames(filtered_baselines_decB_table)
  rownames(EH_rows) <- c("HD")

  observed_rows <- as.data.frame(matrix(0, nrow = 1, ncol = ncol(filtered_baselines_decB_table)))
  colnames(observed_rows) <- colnames(filtered_baselines_decB_table)
  rownames(observed_rows) <- c("HD")


  for (i in 1:ncol(EH_rows)) {
    shannon_rows[nrow(shannon_rows), i] <- shannon_index(filtered_baselines_decB_table1[, i])

    simpson_rows[nrow(simpson_rows), i] <- simpson_index(filtered_baselines_decB_table1[, i])

    EH_rows[nrow(EH_rows), i] <- EH_index(filtered_baselines_decB_table1[, i])

    observed_rows[nrow(observed_rows), i] <- observed_index(filtered_baselines_decB_table1[, i])
  }
  return(list(shannon = shannon_rows, simpson = simpson_rows, EH = EH_rows, observed = observed_rows))
}

HD_05 <- createTabone(tab1, "05", output_folder, filtered_baselines_decB_table)
HD_001 <- createTabone(tab2, "001", output_folder, filtered_baselines_decB_table)

MSSh <- constructTab2(HD_05$shannon, HD_001$shannon, "Shannon")
MSSh <- t(MSSh)
MSSi <- constructTab2(HD_05$simpson, HD_001$simpson, "Simpson")
MSSi <- t(MSSi)
MSEh <- constructTab2(HD_05$EH, HD_001$EH, "EH")
MSEh <- t(MSEh)
MSOb <- constructTab2(HD_05$observed, HD_001$observed, "Observed")
MSOb <- t(MSOb)
merged_sh <- left_join(as.data.frame(MSSh), as.data.frame(MSSi), by = "id") %>%
  left_join(as.data.frame(MSEh), by = "id") %>%
  left_join(as.data.frame(MSOb), by = "id")
merged_sh <- t(merged_sh)
write.table(merged_sh, file = paste0(output_folder, "merged_MSHD_alpha.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

t1<-read_tsv("Output/DAS_ALPHA_ANCOMBC2/merged_shannon_union.tsv",col_names=FALSE,row_names=TRUE)
rownames
t2<-t(t1)
