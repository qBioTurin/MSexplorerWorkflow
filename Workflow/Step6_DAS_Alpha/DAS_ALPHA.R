source("utilities.R")
output_folder = "Output/DAS_ALPHA/"
output_mod="Output/modDAS_ALPHA/"
createFolder(output_folder)
createFolder(output_mod)
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
    if (col_sum != 0) {  # Per evitare divisioni per zero
      abundance_table[, i] <- abundance_table[, i] / col_sum
    }
  }
  return(abundance_table)  # Restituisci la tabella normalizzata
}



DAS_alpha<-function(file_path,domain,file_name,patient,output_folder){

  data<- readRDS(gsub(" ","",paste(file_path,domain,"/",file_name,"_merged.rds")))
  taxa_names <- rownames(data@tax_table)
  filtered_baselines_decB_table <- baselines_decB_table[, colnames(baselines_decB_table) %in% patient, drop = FALSE]
  filtered_baselines_decB_table<- filtered_baselines_decB_table[rownames(baselines_decB_table) %in% taxa_names, ]
  new_rows <- as.data.frame(matrix(0, nrow = 3, ncol = ncol(filtered_baselines_decB_table)))
  rownames(new_rows)<-c("Shannon","Simpson","EH")
  colnames(new_rows) <- colnames(filtered_baselines_decB_table)
  #filtered_baselines_decB_table <- rbind(filtered_baselines_decB_table, new_rows)

  #norm_filtered<-normalize_abundance(filtered_baselines_decB_table)

  for(i in 1:ncol(new_rows)){
    new_rows[nrow(new_rows) - 2,i]<-shannon_index(filtered_baselines_decB_table[,i])
    new_rows[nrow(new_rows) - 1,i]<-simpson_index(filtered_baselines_decB_table[,i])
    new_rows[nrow(new_rows),i]<-EH_index(filtered_baselines_decB_table[,i])
  }

  write.csv(new_rows, 
          file = paste0(output_folder, file_name, "_alpha.csv"), 
          row.names = TRUE)
}

#setPath
base_path <- "Output/MERGED_DAS/"

domain <- c("Bacteria"
#,"Archaea","Eukaryota"
)
gc_treatment <- c("_NO_GC_","_GC_")
tabels<-c("BM_Lesion","Gadolinium","Lesion","Subtentorial")

for (i in 1:length(domain)) {
  for (j in 1:length(gc_treatment)) {
    for (k in 1:length(tabels)) {
      file_name <- paste0(gc_treatment[j], tabels[k])  # Creazione del nome del file

      if (gc_treatment[j] == "_GC_") {
        DAS_alpha(base_path, domain[i], gsub(" ","",paste(domain[i],file_name)), GC_patient, output_folder)
      } else { 
        DAS_alpha(base_path, domain[i], gsub(" ","",paste(domain[i],file_name)), NO_GC_patient, output_folder)
      }
    }
  }
}


##########################from now on only for calculate alpha for das withouth discriminationg gc vs non gc

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

createTab<- function(tab1,tab2,tab3,tab4,dataset,output_mod2,filtered_baselines_decB_table){

  taxa_names1 <- rownames(tab1@tax_table)
  taxa_names2 <- rownames(tab2@tax_table)
  taxa_names3 <- rownames(tab3@tax_table)
  taxa_names4 <- rownames(tab4@tax_table)
    
    
  filtered_baselines_decB_table <- baselines_decB_table[, colnames(baselines_decB_table) %in% NT_HV_patient, drop = FALSE]
  filtered_baselines_decB_table1<- filtered_baselines_decB_table[rownames(baselines_decB_table) %in% taxa_names1, ]
  filtered_baselines_decB_table2<- filtered_baselines_decB_table[rownames(baselines_decB_table) %in% taxa_names2, ]
  filtered_baselines_decB_table3<- filtered_baselines_decB_table[rownames(baselines_decB_table) %in% taxa_names3, ]
  filtered_baselines_decB_table4<- filtered_baselines_decB_table[rownames(baselines_decB_table) %in% taxa_names4, ]

  shannon_rows <- as.data.frame(matrix(0, nrow = 4, ncol = ncol(filtered_baselines_decB_table)))
  colnames(shannon_rows) <- colnames(filtered_baselines_decB_table)
  rownames(shannon_rows)<-c("Lesion","Spinal_Cord","Gadolinium","Subtentorial")

  simpson_rows <- as.data.frame(matrix(0, nrow = 4, ncol = ncol(filtered_baselines_decB_table)))
  colnames(simpson_rows) <- colnames(filtered_baselines_decB_table)
  rownames(simpson_rows)<-c("Lesion","Spinal_Cord","Gadolinium","Subtentorial")

  EH_rows <- as.data.frame(matrix(0, nrow = 4, ncol = ncol(filtered_baselines_decB_table)))
  colnames(EH_rows) <- colnames(filtered_baselines_decB_table)
  rownames(EH_rows)<-c("Lesion","Spinal_Cord","Gadolinium","Subtentorial")


    for(i in 1:ncol(EH_rows)){
        shannon_rows[nrow(shannon_rows) - 3,i]<-shannon_index(filtered_baselines_decB_table1[,i])
        shannon_rows[nrow(shannon_rows) - 2,i]<-shannon_index(filtered_baselines_decB_table2[,i])
        shannon_rows[nrow(shannon_rows)-1,i]<-shannon_index(filtered_baselines_decB_table3[,i])
        shannon_rows[nrow(shannon_rows),i]<-shannon_index(filtered_baselines_decB_table4[,i])

        simpson_rows[nrow(simpson_rows) - 3,i]<-simpson_index(filtered_baselines_decB_table1[,i])
        simpson_rows[nrow(simpson_rows) - 2,i]<-simpson_index(filtered_baselines_decB_table2[,i])
        simpson_rows[nrow(simpson_rows)-1,i]<-simpson_index(filtered_baselines_decB_table3[,i])
        simpson_rows[nrow(simpson_rows),i]<-simpson_index(filtered_baselines_decB_table4[,i])

        EH_rows[nrow(EH_rows) - 3,i]<-EH_index(filtered_baselines_decB_table1[,i])
        EH_rows[nrow(EH_rows) - 2,i]<-EH_index(filtered_baselines_decB_table2[,i])
        EH_rows[nrow(EH_rows)-1,i]<-EH_index(filtered_baselines_decB_table3[,i])
        EH_rows[nrow(EH_rows),i]<-EH_index(filtered_baselines_decB_table4[,i])  
    }
    colnames(shannon_rows) <- paste0(colnames(shannon_rows), "_39")
    colnames(simpson_rows) <- paste0(colnames(simpson_rows), "_39")
    colnames(EH_rows) <- paste0(colnames(EH_rows), "_39")
    write.csv(shannon_rows, file = paste0(output_mod2, "shannon_alpha","_",dataset,".csv"), row.names = TRUE)
    write.csv(simpson_rows, file = paste0(output_mod2, "simpson_alpha","_",dataset,".csv"), row.names = TRUE)
    write.csv(EH_rows, file = paste0(output_mod2, "EH_alpha","_",dataset,".csv"), row.names = TRUE)

}


tab1<- readRDS("Output/DAS_MERGED_HD/Bacteria_hd_39_merged.rds")
tab2<- readRDS("Output/DAS_MERGED_HD/Bacteria_hd_264_merged.rds")
createTabone(tab1,"39",output_folder,filtered_baselines_decB_table)
createTabone(tab2,"264",output_folder,filtered_baselines_decB_table)
createTabone<- function(tab1,dataset,output_mod2,filtered_baselines_decB_table){

  taxa_names1 <- rownames(tab1@tax_table)
    
    
  filtered_baselines_decB_table <- baselines_decB_table
  filtered_baselines_decB_table1<- filtered_baselines_decB_table[rownames(baselines_decB_table) %in% taxa_names1, ]


  shannon_rows <- as.data.frame(matrix(0, nrow = 1, ncol = ncol(filtered_baselines_decB_table)))
  colnames(shannon_rows) <- colnames(filtered_baselines_decB_table)
  rownames(shannon_rows)<-c("HD")

  simpson_rows <- as.data.frame(matrix(0, nrow = 1, ncol = ncol(filtered_baselines_decB_table)))
  colnames(simpson_rows) <- colnames(filtered_baselines_decB_table)
  rownames(simpson_rows)<-c("HD")

  EH_rows <- as.data.frame(matrix(0, nrow = 1, ncol = ncol(filtered_baselines_decB_table)))
  colnames(EH_rows) <- colnames(filtered_baselines_decB_table)
  rownames(EH_rows)<-c("HD")


    for(i in 1:ncol(EH_rows)){
        shannon_rows[nrow(shannon_rows),i]<-shannon_index(filtered_baselines_decB_table1[,i])

        simpson_rows[nrow(simpson_rows) ,i]<-simpson_index(filtered_baselines_decB_table1[,i])

        EH_rows[nrow(EH_rows),i]<-EH_index(filtered_baselines_decB_table1[,i])

    }
  #  colnames(shannon_rows) <- paste0(colnames(shannon_rows), "_39")
  #  colnames(simpson_rows) <- paste0(colnames(simpson_rows), "_39")
  #  colnames(EH_rows) <- paste0(colnames(EH_rows), "_39")
    write.csv(shannon_rows, file = paste0(output_mod2, "shannon_alpha_hd","_",dataset,".csv"), row.names = TRUE)
    write.csv(simpson_rows, file = paste0(output_mod2, "simpson_alpha_hd","_",dataset,".csv"), row.names = TRUE)
    write.csv(EH_rows, file = paste0(output_mod2, "EH_alpha_hd","_",dataset,".csv"), row.names = TRUE)

}
# Example usage
Sh_39 <- read.csv("Output/FuseAlphaTab/shannon_alpha39.csv")
Si_39 <- read.csv("Output/FuseAlphaTab/simpson_alpha39.csv")
EH_39 <- read.csv("Output/FuseAlphaTab/EH_alpha39.csv")

Sh_116 <- read.csv("Output/FuseAlphaTab/shannon_alpha116.csv")
Si_116 <- read.csv("Output/FuseAlphaTab/simpson_alpha116.csv")
EH_116 <- read.csv("Output/FuseAlphaTab/EH_alpha116.csv")

ShAll <- read.csv("Output/FuseAlphaTab/shannon_alphaALL.csv")
SiAll <- read.csv("Output/FuseAlphaTab/simpson_alphaALL.csv")
EHAll <- read.csv("Output/FuseAlphaTab/EH_alphaALL.csv")

# Merge the three tables by row names
# Merge the two data frames by the specified column
merged_shannon <- merge(Sh_39, Sh_116, by = "X", all = TRUE)
merged_shannon <- merge(merged_shannon, ShAll, by="X", all = TRUE)

merged_shannon <- merge(Sh_39, Sh_116, by = "row.names", all = TRUE)
merged_shannon <- merge(merged_shannon, ShAll, by.x = "Row.names", by.y = "row.names", all = TRUE)
rownames(merged_shannon) <- merged_shannon$Row.names
merged_shannon <- merged_shannon[, !colnames(merged_shannon) %in% c("Row.names")]

merged_simpson <- merge(Si_39, Si_116, by = "row.names", all = TRUE)
merged_simpson <- merge(merged_simpson, SiAll, by.x = "Row.names", by.y = "row.names", all = TRUE)
rownames(merged_simpson) <- merged_simpson$Row.names
merged_simpson <- merged_simpson[, !colnames(merged_simpson) %in% c("Row.names")]

merged_EH <- merge(EH_39, EH_116, by = "row.names", all = TRUE)
merged_EH <- merge(merged_EH, EHAll, by.x = "Row.names", by.y = "row.names", all = TRUE)
rownames(merged_EH) <- merged_EH$Row.names
merged_EH <- merged_EH[, !colnames(merged_EH) %in% c("Row.names")]

# Save the merged tables to CSV files
write.csv(merged_shannon, file = paste0(output_mod2, "merged_shannon_alpha.csv"), row.names = TRUE)
write.csv(merged_simpson, file = paste0(output_mod2, "merged_simpson_alpha.csv"), row.names = TRUE)
write.csv(merged_EH, file = paste0(output_mod2, "merged_EH_alpha.csv"), row.names = TRUE)