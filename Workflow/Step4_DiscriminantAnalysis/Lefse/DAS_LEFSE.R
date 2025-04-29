source("Settings/utilities.R")

check_or_pull_image <- function(image_name) {
  # Check if image exists locally
  result <- system(paste("docker image inspect", image_name), intern = TRUE, ignore.stderr = TRUE)
  
  if (length(result) == 0) {
    message("Image not found locally. Pulling from Docker Hub...")
    pull_status <- system(paste("docker pull", image_name))
    
    if (pull_status == 0) {
      message("Image pulled successfully.")
    } else {
      stop("Failed to pull image from Docker Hub.")
    }
  } else {
    message("Image exists locally.")
  }
}

check_or_pull_image("qbioturin/lefse")

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
folderGC_comp="Output/LEFSE/GC_comp/"
createFolder(folderGC_comp)


output_folder_001 = ("Output/LEFSE/001/step0/")
output_folder_01= ("Output/LEFSE/01/step0/")
output_folder_05= ("Output/LEFSE/05/step0/")
output_folderMSHD= ("Output/LEFSE/MSHD/step0/")
output_folderGC=("Output/LEFSE/GC/step0/")
output_folderGC_comp=("Output/LEFSE/GC_comp/step0/")

createFolder(output_folder_001)
createFolder(output_folder_01)
createFolder(output_folder_05)
createFolder(output_folderMSHD)
createFolder(output_folderGC)
createFolder(output_folderGC_comp)

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
  else if(column == "gc_treatment" || status == "both") {
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
  

  write.table(glcLesionFinal, file = gsub(" ","",paste0(output_folder, fileName,"_lefse")), sep = "\t", quote = FALSE, col.names = FALSE)
}


Archaea_Supervised_decontam <- readRDS("Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.001.rds")
Eukaryote_Supervised_decontam <- readRDS("Output/SUPERVISED_DEC/Eukaryote_Supervised_decontam0.001.rds")
Bacteria_Supervised_decontam001 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
Bacteria_Supervised_decontam01 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds")
Bacteria_Supervised_decontam05 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds")

metadataB <- read_delim("Output/AlphaMetadata/Bacteria_alpha_metadata.csv", delim = ",", escape_double = FALSE, trim_ws = TRUE)
metadataA <- read_delim("Output/AlphaMetadata/Archaea_alpha_metadata.csv", delim = ",", escape_double = FALSE, trim_ws = TRUE)
metadataE <- read_delim("Output/AlphaMetadata/Eukaryote_alpha_metadata.csv", delim = ",", escape_double = FALSE, trim_ws = TRUE)

generate.LEFSE1(Bacteria_Supervised_decontam001, metadataB, "category", "both", "Bacteria_MsHd", output_folderMSHD)
generate.LEFSE1(Eukaryote_Supervised_decontam, metadataE, "category", "both", "Eukaryote_MsHd", output_folderMSHD)
generate.LEFSE1(Archaea_Supervised_decontam, metadataA, "category", "both", "Archaea_MsHd", output_folderMSHD)
generate.LEFSE1(Bacteria_Supervised_decontam05, metadataB, "category", "both", "Bacteria_MsHd_05", output_folderMSHD)

generate.LEFSE1(Bacteria_Supervised_decontam001, metadataB, "gc_treatment", "both", "Bacteria_GC", output_folderGC)
generate.LEFSE1(Eukaryote_Supervised_decontam, metadataE, "gc_treatment", "both", "Eukaryote_GC", output_folderGC)
generate.LEFSE1(Archaea_Supervised_decontam, metadataA, "gc_treatment", "both", "Archaea_GC", output_folderGC)

analysis <- c("lesion_burden", "spinal_cord_lesion", "gadolinium_contrast", "subtentorial_lesions","gc_treatment")
status <- c("positive", "negative")
for (i in 1:length(analysis)) {
  name <- gsub(" ","",paste("Bacteria_", analysis[i]))
  generate.LEFSE1(Bacteria_Supervised_decontam001, metadataB, analysis[i], "both", gsub(" ","",paste(name, "_001")), output_folder_001)
  generate.LEFSE1(Bacteria_Supervised_decontam01, metadataB, analysis[i], "both", gsub(" ","",paste(name, "_01")), output_folder_01)
  generate.LEFSE1(Bacteria_Supervised_decontam05, metadataB, analysis[i], "both", gsub(" ","",paste(name, "_05")), output_folder_05)
  for(j in 1:length(status)){
    generate.LEFSE1(Bacteria_Supervised_decontam001, metadataB, analysis[i], status[j],paste0("Bacteria_", analysis[i], "_", status[j]), output_folderGC_comp)
    generate.LEFSE1(Eukaryote_Supervised_decontam, metadataE, analysis[i], status[j], paste0("Eukaryote_", analysis[i], "_", status[j]), output_folderGC_comp)
    generate.LEFSE1(Archaea_Supervised_decontam, metadataA, analysis[i], status[j], paste0("Archaea_", analysis[i], "_", status[j]), output_folderGC_comp)
  }
}


folders <- list(folder01, folder001, folder05, folderMSHD, folderGC, folderGC_comp)
docker_name <- "qbioturin/lefse:latest"
container_ids <- c()

# Avvia i container in background e salva gli ID
for (folder in folders) {
  cmd <- paste0(
    "docker run -v ", normalizePath(folder), ":/input_files/ -d ",
    docker_name, " bash Scripts/lefseEx.sh step0/")
  container_id <- system(cmd, intern = TRUE)
  container_ids <- c(container_ids, container_id)
}

# Aspetta la fine di tutti i container e poi li rimuove
for (id in container_ids) {
  system(paste("docker wait", id))
  system(paste("docker rm", id))
}

print("All containers have finished processing.")
