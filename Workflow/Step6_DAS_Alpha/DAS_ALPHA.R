source("Settings/utilities.R")
output_folder = "Output/DAS_ALPHA/"
createFolder(output_folder)
#abundance
bact_baselines_ds_abund = readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
baselines_decB_table = as.data.frame(abundances(bact_baselines_ds_abund, transform = "compositional"))
#setPatientIds
metaData <- read.csv("InputData/metadataMS.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
metadava_revised <- metaData[!is.na(metaData$gc_treatment),]%>%
  filter(gc_treatment == "positive" | gc_treatment == "negative")%>%
  select(id, gc_treatment)
GC_patient <- metadava_revised[metadava_revised$gc_treatment == "positive", "id"]
NO_GC_patient <- metadava_revised[metadava_revised$gc_treatment == "negative", "id"]

NT_HV_patient <- union(GC_patient, NO_GC_patient)
#setTaxid
shannon_index <- function(x) {
  x <- x[!is.na(x) & x > 0]
  p <- x / sum(x)
  -sum(p * log(p))
}

simpson_index <- function(x) {
  x <- x[!is.na(x) & x > 0]
  p <- x / sum(x)
  1 - sum(p^2)
}

EH_index <- function(x) {
  shannon_index(x) / log(length(x[x > 0]))
}

normalize_abundance <- function(abundance_table) {
  for (i in 1:ncol(abundance_table)) {
    col_sum <- sum(abundance_table[, i], na.rm = TRUE)
    if (col_sum != 0) {  # Avoid division by zero
      abundance_table[, i] <- abundance_table[, i] / col_sum
    }
  }
  return(abundance_table)  # Return the normalized table
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

  simpson_rows <- as.data.frame(matrix(0, nrow = 4, ncol = ncol(filtered_baselines_decB_table)))
  colnames(simpson_rows) <- colnames(filtered_baselines_decB_table)
  rownames(simpson_rows) <- c("Lesion", "Spinal_Cord", "Gadolinium", "Subtentorial")

  EH_rows <- as.data.frame(matrix(0, nrow = 4, ncol = ncol(filtered_baselines_decB_table)))
  colnames(EH_rows) <- colnames(filtered_baselines_decB_table)
  rownames(EH_rows) <- c("Lesion", "Spinal_Cord", "Gadolinium", "Subtentorial")

  for (i in 1:ncol(EH_rows)) {
    shannon_rows[nrow(shannon_rows) - 3, i] <- shannon_index(lesion_table[, i])
    shannon_rows[nrow(shannon_rows) - 2, i] <- shannon_index(spinal_table[, i])
    shannon_rows[nrow(shannon_rows) - 1, i] <- shannon_index(gado_table[, i])
    shannon_rows[nrow(shannon_rows), i] <- shannon_index(sub_table[, i])

    simpson_rows[nrow(simpson_rows) - 3, i] <- simpson_index(lesion_table[, i])
    simpson_rows[nrow(simpson_rows) - 2, i] <- simpson_index(spinal_table[, i])
    simpson_rows[nrow(simpson_rows) - 1, i] <- simpson_index(gado_table[, i])
    simpson_rows[nrow(simpson_rows), i] <- simpson_index(sub_table[, i])

    EH_rows[nrow(EH_rows) - 3, i] <- EH_index(lesion_table[, i])
    EH_rows[nrow(EH_rows) - 2, i] <- EH_index(spinal_table[, i])
    EH_rows[nrow(EH_rows) - 1, i] <- EH_index(gado_table[, i])
    EH_rows[nrow(EH_rows), i] <- EH_index(sub_table[, i])
  }
  return(list(shannon=shannon_rows, simpson= simpson_rows, EH= EH_rows))
}

lesion_05<- readRDS("Output/merge_DAS/05/Bacteria_lesion_burden_05_merged.rds")
bm_05<- readRDS("Output/merge_DAS/05/Bacteria_spinal_cord_lesion_05_merged.rds")
gado_05<- readRDS("Output/merge_DAS/05/Bacteria_gadolinium_contrast_05_merged.rds")
sub_05<- readRDS("Output/merge_DAS/05/Bacteria_subtentorial_lesions_05_merged.rds")

lesion_01<- readRDS("Output/merge_DAS/01/Bacteria_lesion_burden_01_merged.rds")
bm_01<- readRDS("Output/merge_DAS/01/Bacteria_spinal_cord_lesion_01_merged.rds")
gado_01<- readRDS("Output/merge_DAS/01/Bacteria_gadolinium_contrast_01_merged.rds")
sub_01<- readRDS("Output/merge_DAS/01/Bacteria_subtentorial_lesions_01_merged.rds")

lesion_001<- readRDS("Output/merge_DAS/001/Bacteria_lesion_burden_001_merged.rds")
bm_001<- readRDS("Output/merge_DAS/001/Bacteria_spinal_cord_lesion_001_merged.rds")
gado_001<- readRDS("Output/merge_DAS/001/Bacteria_gadolinium_contrast_001_merged.rds")
sub_001<- readRDS("Output/merge_DAS/001/Bacteria_subtentorial_lesions_001_merged.rds")


tab_list_05=createTab(lesion_05,bm_05,gado_05,sub_05,filtered_baselines_decB_table)
tab_list_01=createTab(lesion_01,bm_01,gado_01,sub_01,filtered_baselines_decB_table)
tab_list_001=createTab(lesion_001,bm_001,gado_001,sub_001, filtered_baselines_decB_table)

tabMod<-function(tab,Alpha,Method,Subset){
  tab$Alpha <- Alpha
  tab$Method <- Method
  tab$Discriminant <- rownames(tab)
  tab$Subset <- Subset
  tab=t(tab)
  tab <- cbind(id = rownames(tab), tab)
  tab <- as.data.frame(tab) 
  return(tab)
}

constructTab <- function(tab05, tab01, tab001, Alpha){
  tab_05 <- tabMod(tab05, Alpha, "Alpha", "05")
  tab_01 <- tabMod(tab01, Alpha, "Alpha", "01")
  tab_001 <- tabMod(tab001, Alpha, "Alpha", "001")

  merged_sh <- left_join(tab_05, tab_01, by = "id") %>%
    left_join(tab_001, by = "id")

  merged_sh <- t(merged_sh)

  write.table(
    merged_sh,
    file = paste0(output_folder, "merged_", Alpha, "_alpha.tsv"),
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}


constructTab(tab_list_05$shannon,tab_list_01$shannon,tab_list_001$shannon,"Shannon")
constructTab(tab_list_05$simpson,tab_list_01$simpson,tab_list_001$simpson,"Simpson")
constructTab(tab_list_05$EH,tab_list_01$EH,tab_list_001$EH,"EH")
merged_sh <- left_join(sh_05, sh_01, by = "id") %>%
  left_join(sh_001, by = "id")
merged_sh<-t(merged_sh)
write.table(merged_sh, file = paste0(output_folder, "merged_shannon_alpha.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE)

sh_05=tab_list_05$shannon
sh_05=tabMod(sh_05,"Shannon", "Alpha", "05")
sh_01=tab_list_01$shannon
sh_01=tabMod(sh_01,"Shannon", "Alpha", "01")
sh_001=tab_list_001$shannon
sh_001=tabMod(sh_001,"Shannon", "Alpha", "001")

merged_sh <- left_join(sh_05, sh_01, by = "id") %>%
  left_join(sh_001, by = "id")
merged_sh<-t(merged_sh)
write.table(merged_sh, file = paste0(output_folder, "merged_shannon_alpha.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE)

