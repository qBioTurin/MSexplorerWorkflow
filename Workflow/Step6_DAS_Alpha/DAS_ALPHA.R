source("utilities.R")
output_folder = "Output/DAS_ALPHA/"
output_mod="Output/modDAS_ALPHA/"
createFolder(output_folder)
createFolder(output_mod)

#setPath
base_path <- "Output/MERGED_DAS/"

output_folder39 = "Output/Alpha/39/"
output_folder116 = "Output/Alpha/116/"
output_folder167 = "Output/Alpha/167/"
output_folder264 = "Output/Alpha/264/"

output_folder39LF="Output/DAS_ONLY_LEFSE/Alpha39/"
output_folder116LF="Output/DAS_ONLY_LEFSE/Alpha116/"
output_folder264LF="Output/DAS_ONLY_LEFSE/Alpha264/"

createFolder(output_folder39LF)
createFolder(output_folder116LF)
createFolder(output_folder264LF)

createFolder(output_folder39)
createFolder(output_folder116)
createFolder(output_folder167)
createFolder(output_folder264)


lesion_39<- readRDS("Output/MERGED_DAS/39/Bacteria_Lesion39_merged.rds")
bm_39<- readRDS("Output/MERGED_DAS/39/Bacteria_BM_Lesion39_merged.rds")
gado_39<- readRDS("Output/MERGED_DAS/39/Bacteria_Gadolinium39_merged.rds")
sub_39<- readRDS("Output/MERGED_DAS/39/Bacteria_Subtentorial39_merged.rds")


lesion_116<- readRDS("Output/MERGED_DAS/116/Bacteria_Lesion116_merged.rds")
bm_116<- readRDS("Output/MERGED_DAS/116/Bacteria_BM_Lesion116_merged.rds")
gado_116<- readRDS("Output/MERGED_DAS/116/Bacteria_Gadolinium116_merged.rds")
sub_116<- readRDS("Output/MERGED_DAS/116/Bacteria_Subtentorial116_merged.rds")

lesion_167<- readRDS("Output/MERGED_DAS/167/Bacteria_Lesion167_merged.rds")
bm_167<- readRDS("Output/MERGED_DAS/167/Bacteria_BM_Lesion167_merged.rds")
gado_167<- readRDS("Output/MERGED_DAS/167/Bacteria_Gadolinium167_merged.rds")
sub_167<- readRDS("Output/MERGED_DAS/167/Bacteria_Subtentorial167_merged.rds")


lesion_264<- readRDS("Output/MERGED_DAS/264/Bacteria_LesionAll_merged.rds")
bm_264<- readRDS("Output/MERGED_DAS/264/Bacteria_BM_LesionAll_merged.rds")
gado_264<- readRDS("Output/MERGED_DAS/264/Bacteria_GadoliniumAll_merged.rds")
sub_264<- readRDS("Output/MERGED_DAS/264/Bacteria_SubtentorialAll_merged.rds")


#lef_lesion_39<-readRDS("Output/DAS_ONLY_LEFSE/39RDS/Bacteria_Lesion39_merged.rds")
#lef_bm_39<-readRDS("Output/DAS_ONLY_LEFSE/39RDS/Bacteria_BM_Lesion39_merged.rds")
#lef_gado_39<-readRDS("Output/DAS_ONLY_LEFSE/39RDS/Bacteria_Gadolinium39_merged.rds")
#lef_sub_39<-readRDS("Output/DAS_ONLY_LEFSE/39RDS/Bacteria_Subtentorial39_merged.rds")

#lef_lesion_116<-readRDS("Output/DAS_ONLY_LEFSE/116RDS/Bacteria_Lesion116_merged.rds")
#lef_bm_116<-readRDS("Output/DAS_ONLY_LEFSE/116RDS/Bacteria_BM_Lesion116_merged.rds")
#lef_gado_116<-readRDS("Output/DAS_ONLY_LEFSE/116RDS/Bacteria_Gadolinium116_merged.rds")
#lef_sub_116<-readRDS("Output/DAS_ONLY_LEFSE/116RDS/Bacteria_Subtentorial116_merged.rds")

#lef_lesion_264<-readRDS("Output/DAS_ONLY_LEFSE/264RDS/Bacteria_LesionAll_merged.rds")
#lef_bm_264<-readRDS("Output/DAS_ONLY_LEFSE/264RDS/Bacteria_BM_LesionAll_merged.rds")
#lef_gado_264<-readRDS("Output/DAS_ONLY_LEFSE/264RDS/Bacteria_GadoliniumAll_merged.rds")
#lef_sub_264<-readRDS("Output/DAS_ONLY_LEFSE/264RDS/Bacteria_SubtentorialAll_merged.rds")

createTab(lesion_39,bm_39,gado_39,sub_39,"39",output_folder39,filtered_baselines_decB_table)
createTab(lesion_116,bm_116,gado_116,sub_116,"116",output_folder116, filtered_baselines_decB_table)
createTab(lesion_264,bm_264,gado_264,sub_264,"264",output_folder264, filtered_baselines_decB_table)
createTab(lesion_167,bm_167,gado_167,sub_167,"167",output_folder167, filtered_baselines_decB_table)

createTab(lef_lesion_39,lef_bm_39,lef_gado_39,lef_sub_39,"39",output_folder39LF,filtered_baselines_decB_table)
createTab(lef_lesion_116,lef_bm_116,lef_gado_116,lef_sub_116,"116",output_folder116LF, filtered_baselines_decB_table)
createTab(lef_lesion_264,lef_bm_264,lef_gado_264,lef_sub_264,"264",output_folder264LF, filtered_baselines_decB_table)

#abundance
bact_baselines_ds_abund = readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontamOLD.rds")
print(bact_baselines_ds_abund)
baselines_decB_table = as.data.frame(abundances(bact_baselines_ds_abund, transform = "compositional"))
#setPatientIds
metaData <- read.csv("input/20241205_MetadataHSvsT0_modified.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
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

createTab <- function(lesion, spinal, gado, sub, dataset, output_mod2, filtered_baselines_decB_table) {
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

  colnames(shannon_rows) <- paste0(colnames(shannon_rows), "_39")
  colnames(simpson_rows) <- paste0(colnames(simpson_rows), "_39")
  colnames(EH_rows) <- paste0(colnames(EH_rows), "_39")

  write.csv(shannon_rows, file = paste0(output_mod2, "shannon_alpha", "_", dataset, ".csv"), row.names = TRUE)
  write.csv(simpson_rows, file = paste0(output_mod2, "simpson_alpha", "_", dataset, ".csv"), row.names = TRUE)
  write.csv(EH_rows, file = paste0(output_mod2, "EH_alpha", "_", dataset, ".csv"), row.names = TRUE)
}


