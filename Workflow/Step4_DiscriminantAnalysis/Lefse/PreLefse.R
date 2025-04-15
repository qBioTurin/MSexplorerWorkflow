source("Settings/utilities.R")
createFolder("Output/LEFSE/")
folder001="Output/LEFSE/001/"
createFolder(folder001)
folder01="Output/LEFSE/01/"
createFolder(folder01)
folder05="Output/LEFSE/05/"
createFolder(folder05)
folderMSHD="Output/LEFSE/MSHD/"
createFolder(folderMSHD)
folderGC="Output/LEFSE/GC/"
createFolder(folderGC)


output_folder_001 = ("Output/LEFSE/001/step0/")
output_folder_01= ("Output/LEFSE/01/step0/")
output_folder_05= ("Output/LEFSE/05/step0/")
output_folderMSHD= ("Output/LEFSE/MSHD/step0/")
output_folderGC=("Output/LEFSE/GC/step0/")

createFolder(output_folder_001)
createFolder(output_folder_01)
createFolder(output_folder_05)
createFolder(output_folderMSHD)
createFolder(output_folderGC)

generate.LEFSE1 <- function(Domain, metadata, column, status, fileName,output_folder ) {
  
  taxatables <- as.data.frame(tax_table(Domain))
  taxatables <- taxatables %>%
    mutate(namesSpecies = paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = "|"))
  taxatables <- rownames_to_column(taxatables, var = "otu_id")
  
  otutables <- data.frame(otu_table(Domain))
  otutables <- as.data.frame(abundances(otutables, transform = "compositional"))
  otutables <- rownames_to_column(otutables, var = "otu_id")
  
  combined_df <- full_join(taxatables, otutables, by = "otu_id")
  combined_df <- combined_df[, 9:80]
  
  comb_t <- t(combined_df)
  colnames(comb_t) <- comb_t[1, ]
  comb_t <- comb_t[-1, ]
  comb_t <- as.data.frame(comb_t)

  metadata$gc_treatment=as.factor(metadata$gc_treatment)
  metadata$category=as.factor(metadata$category)
  metadata$lesion_burden=as.factor(metadata$lesion_burden)
  levels(metadata$lesion_burden) =c("low","high")
  metadata$spinal_cord_lesion=as.factor(metadata$spinal_cord_lesion)
  levels(metadata$spinal_cord_lesion) =c("BM_low","BM_high")
  metadata$gadolinium_contrast=as.factor(metadata$gadolinium_contrast)
  levels(metadata$gadolinium_contrast) =c("NoActive","Active")
  metadata$subtentorial_lesions=as.factor(metadata$subtentorial_lesions)
  levels(metadata$subtentorial_lesions) =c("No","Yes")


  metaSel <- metadata %>%
    select(id, category, gc_treatment, lesion_burden, spinal_cord_lesion, gadolinium_contrast, subtentorial_lesions) %>%
    as.data.frame()
  colnames(metaSel) <- c("rownames", "category", "gc_treatment", "lesion_burden", "spinal_cord_lesion", "gadolinium_contrast", "subtentorial_lesions")

  comb_t <- comb_t %>%
    mutate(rownames = rownames(comb_t))

  joined_df <- comb_t %>%
    inner_join(metaSel, by = "rownames")
  

  total <- c("rownames", "category", "gc_treatment", "lesion_burden", "spinal_cord_lesion", "gadolinium_contrast", "subtentorial_lesions")
  partial <- setdiff(total, column)
  if (column == "category" ) {
    new_joined_df <- joined_df 
  }
  else if(column == "gc_treatment" | status == "both") {
    new_joined_df <- joined_df %>%
      filter(gc_treatment == "positive" | gc_treatment == "negative")
  } 
  else {
     new_joined_df <- joined_df %>%
      filter(gc_treatment == status)
    }  
  
  otuValue <- select(new_joined_df, -all_of(partial))
  otuValue <- t(otuValue)
  otuValue <- as.data.frame(otuValue)
  
  glcLesionFinal <- rbind(new_joined_df$lesion_burden, new_joined_df$rownames, otuValue)
  glcLesionFinal <- rbind(glcLesionFinal[nrow(glcLesionFinal), ], glcLesionFinal[-nrow(glcLesionFinal), ])
  glcLesionFinal <- glcLesionFinal[-2, ]
  

  write.table(glcLesionFinal, file = paste0(output_folder, fileName), sep = "\t", quote = FALSE, col.names = FALSE)
}


ARCH_Supervised_decontam <- readRDS("Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.001.rds")
EUK_Supervised_decontam <- readRDS("Output/SUPERVISED_DEC/Eukaryota_Supervised_decontam0.001.rds")
BACT_Supervised_decontam001 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
BACT_Supervised_decontam01 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds")
BACT_Supervised_decontam05 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds")

metadataB <- read_delim("Output/AlphaMetadata/Bacteria_alpha_metadata.csv", delim = ",", escape_double = FALSE, trim_ws = TRUE)
metadataA <- read_delim("Output/AlphaMetadata/Archaea_alpha_metadata.csv", delim = ",", escape_double = FALSE, trim_ws = TRUE)
metadataE <- read_delim("Output/AlphaMetadata/Eukaryota_alpha_metadata.csv", delim = ",", escape_double = FALSE, trim_ws = TRUE)

generate.LEFSE1(BACT_Supervised_decontam001, metadataB, "category", "both", "BACT_MsHd", output_folderMSHD)
generate.LEFSE1(EUK_Supervised_decontam, metadataE, "category", "both", "EUK_MsHd", output_folderMSHD)
generate.LEFSE1(ARCH_Supervised_decontam, metadataA, "category", "both", "ARCH_MsHd", output_folderMSHD)

generate.LEFSE1(BACT_Supervised_decontam001, metadataB, "gc_treatment", "both", "BACT_GC", output_folderGC)
generate.LEFSE1(EUK_Supervised_decontam, metadataE, "gc_treatment", "both", "EUK_GC", output_folderGC)
generate.LEFSE1(ARCH_Supervised_decontam, metadataA, "gc_treatment", "both", "ARCH_GC", output_folderGC)

analysis <- c("lesion_burden", "spinal_cord_lesion", "gadolinium_contrast", "subtentorial_lesions")

for (i in 1:length(analysis)) {
  name <- paste("BACT_", analysis[i])
  generate.LEFSE1(BACT_Supervised_decontam001, metadataB, analysis[i], "both", paste(name, "_001"), output_folder_001)
  generate.LEFSE1(BACT_Supervised_decontam01, metadataB, analysis[i], "both", paste(name, "_01"), output_folder_01)
  generate.LEFSE1(BACT_Supervised_decontam05, metadataB, analysis[i], "both", paste(name, "_05"), output_folder_05)
}

generate.LEFSE1(BACT_Supervised_decontam001, metadataB, "gc_treatment", "both", paste(name, "_001"), output_folder_001)
generate.LEFSE1(BACT_Supervised_decontam01, metadataB, "gc_treatment", "both", paste(name, "_01"), output_folder_01)
generate.LEFSE1(BACT_Supervised_decontam05, metadataB, "gc_treatment", "both", paste(name, "_05"), output_folder_05)

  system(paste0(
  "docker run -v",normalizePath(folder01),":/input_files/ -it ",
  "fpant/lefse bash Scripts/lefseEx.sh step0/"))

  system(paste0(
  "docker run -v",normalizePath(folder001),":/input_files/ -it ",
  "fpant/lefse bash Scripts/lefseEx.sh step0/"))

    system(paste0(
  "docker run -v",normalizePath(folder05),":/input_files/ -it ",
  "fpant/lefse bash Scripts/lefseEx.sh step0/"))

  system(paste0(
  "docker run -v",normalizePath(folderMSHD),":/input_files/ -it ",
  "fpant/lefse bash Scripts/lefseEx.sh step0/"))

    system(paste0(
  "docker run -v",normalizePath(folderGC),":/input_files/ -it ",
  "fpant/lefse bash Scripts/lefseEx.sh step0/"))

