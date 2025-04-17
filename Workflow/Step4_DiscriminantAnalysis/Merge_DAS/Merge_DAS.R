source("Settings/utilities.R")
output_folderMSHD="Output/merge_DAS/MSHD/"
createFolder(output_folderMSHD)
output_folderGC="Output/merge_DAS/GC/"
createFolder(output_folderGC)
output_folder = "Output/merge_DAS/"
createFolder(output_folder)
output_folder_001 = "Output/merge_DAS/001/"
createFolder(output_folder_001)
output_folder_01 = "Output/merge_DAS/01/"
createFolder(output_folder_01)
output_folder_05 = "Output/merge_DAS/05/"
createFolder(output_folder_05)
output_folderGC_comp="Output/merge_DAS/GC_comp/"
createFolder(output_folderGC_comp)


remove_common_lefse<- function(tab1,tab2){
tab1$Comparison <- paste(tab1$Genus, tab1$Species, sep = " ")
tab2$Comparison <- paste(tab2$Genus, tab2$Species, sep = " ")
# Subtract rows in tab1 that also appear in tab2
tab3 <- tab1[!tab1$Comparison %in% tab2$Comparison, ]
tab3$Comparison <- NULL
  print(paste("Number of lines removed:", nrow(tab1) - nrow(tab3)))
  print(paste("Number of lines remaining:", nrow(tab3)))
return(tab3)
}

remove_common_limma<-function(tab1,tab2){
  tab3 <- tab1[!tab1$X %in% tab2$X, ]
  print(paste("Number of lines removed:", nrow(tab1) - nrow(tab3)))
  print(paste("Number of lines remaining:", nrow(tab3)))
  return(tab3)
}


merge_das <- function(baselines_dec, das_lefse, das_limma, type, name, output_folder) {

  print(paste("Number of rows in das_lefse:", nrow(das_lefse)))
  print(paste("Number of rows in das_limma:", nrow(das_limma)))
  otu <- as.data.frame(abundances(otu_table(baselines_dec), transform = "compositional"))
  baselines_dec <- phyloseq(otu_table(otu, taxa_are_rows = TRUE), tax_table(baselines_dec), sample_data(baselines_dec))
  das_lefse <- das_lefse %>%
    mutate(Genus_Species = paste(Genus, Species, sep = " "))
  taxa <- as.data.frame(baselines_dec@tax_table)
  taxa <- taxa %>%
    mutate(Genus_Species = paste(Genus, Species, sep = " "))
  colnames(das_limma)[colnames(das_limma) == "X"] <- "Genus_Species"
  lefse_list<- das_lefse$Genus_Species
  limma_list<- das_limma$Genus_Species
  if(length(lefse_list)==0 && length(limma_list)==0){
    print(paste("Error: both lefse and limma list are empty. Please check the lefse file", name, "for incorrect species.(could be also both empty)you will not see this rds"))
  }
  else {
    if(length(lefse_list)>0 && length(limma_list)>0){
      DAS <- union(das_lefse$Genus_Species, das_limma$Genus_Species)}
    else if (length(lefse_list)>0) {
      DAS <- das_lefse$Genus_Species
    }
    else{
      DAS <- das_limma$Genus_Species
    }

    DAS_dimension<-length(DAS)
    unique_species<-length(union(lefse_list,limma_list)) 
    if(DAS_dimension!=unique_species){
      print(paste("Error: DAS list is not equal to the union of lefse and limma list. Please check the lefse file", name, "for incorrect species."))
    }
    taxa <- taxa[taxa$Genus_Species %in% DAS, ]
    taxa <- taxa[, -8]
    keep_otu <- rownames(taxa)
    baselines_dec <- prune_taxa(keep_otu, baselines_dec)
    if (type == "GC") {
      samples = as.data.frame(sample_data(baselines_dec))
      samples = samples[samples$gc_treatment %in% c("positive", "negative"),]
      baselines_dec <- prune_samples(samples$id, baselines_dec)
    } else if (type == "negative") {
        samples = data.frame(sample_data(baselines_dec))
        samples = samples %>%
        filter(gc_treatment == "negative") %>%
        droplevels()
        baselines_dec = prune_samples(samples$id, baselines_dec)
    } else if (type == "positive") {
        samples = data.frame(sample_data(baselines_dec))
        samples = samples %>%
        filter(gc_treatment == "positive") %>%
        droplevels()
        baselines_dec = prune_samples(samples$id, baselines_dec)
    }
    final_name <- paste(output_folder, name, "_merged.rds", sep = "")
    print(baselines_dec)
    baselines_dec <- saveRDS(baselines_dec, paste(final_name))
  }
  
}

baselines_dec_001 = readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
baselines_decA = readRDS(file = "Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.001.rds")
baselines_decE = readRDS(file = "Output/SUPERVISED_DEC/Eukaryote_Supervised_decontam0.001.rds")

baselines_dec_01 = readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds")
baselines_dec_05 = readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds")

das_lefseB_MsHd<- read_tsv("Output/LEFSE/MSHD/final_name/Bacteria_MsHd_lefse.res")
das_lefseA_MsHd<- read_tsv("Output/LEFSE/MSHD/final_name/Archaea_MsHd_lefse.res")
das_lefseE_MsHd<- read_tsv("Output/LEFSE/MSHD/final_name/Eukaryote_MsHd_lefse.res")
das_lefseB_GC<- read_tsv("Output/LEFSE/GC/final_name/Bacteria_GC_lefse.res")
das_lefseA_GC<- read_tsv("Output/LEFSE/GC/final_name/Archaea_GC_lefse.res")
das_lefseE_GC<- read_tsv("Output/LEFSE/GC/final_name/Eukaryote_GC_lefse.res")

das_limmaB_MsHd<-read.csv("Output/LIMMA_score/MSHD/Bacteria_category_limma.csv")
das_limmaA_MsHd<-read.csv("Output/LIMMA_score/MSHD/Archaea_category_limma.csv")
das_limmaE_MsHd<-read.csv("Output/LIMMA_score/MSHD/Eukaryote_category_limma.csv")
das_limmaB_GC<-read.csv("Output/LIMMA_score/GC/Bacteria_gc_treatment_limma.csv")
das_limmaA_GC<-read.csv("Output/LIMMA_score/GC/Archaea_gc_treatment_limma.csv")
das_limmaE_GC<-read.csv("Output/LIMMA_score/GC/Eukaryote_gc_treatment_limma.csv")


###### Merge DAS for MsHd
merge_das(baselines_dec_001, das_lefseB_MsHd, das_limmaB_MsHd, "msHd", "Bacteria_MsHd", output_folderMSHD)
merge_das(baselines_decA, das_lefseA_MsHd, das_limmaA_MsHd, "msHd", "Archaea_MsHd", output_folderMSHD)
merge_das(baselines_decE, das_lefseE_MsHd, das_limmaE_MsHd, "msHd", "Eukaryote_MsHd", output_folderMSHD)

##### Merge DAS for GC
merge_das(baselines_dec_001, das_lefseB_GC, das_limmaB_GC, "GC", "Bacteria_GC", output_folderGC)
merge_das(baselines_decA, das_lefseA_GC, das_limmaA_GC, "GC", "Archaea_GC", output_folderGC)
merge_das(baselines_decE, das_lefseE_GC, das_limmaE_GC, "GC", "Eukaryote_GC", output_folderGC)

domain = c("Bacteria", "Archaea", "Eukaryote")
dimension = c("001", "05", "01")
analysis = c("lesion_burden", "spinal_cord_lesion", "gadolinium_contrast", "subtentorial_lesions")
status = c("positive", "negative")

for (i in seq_along(analysis)) {
  for (j in seq_along(dimension)) {
    lefse <- read_tsv(gsub(" ","",paste0("Output/LEFSE/",dimension[j],"/final_name/Bacteria_",analysis[i],"_", dimension[j],"_lefse.res")))
    gc_lefse <-read_tsv(gsub(" ","",paste0("Output/LEFSE/",dimension[j],"/final_name/Bacteria_","gc_treatment_", dimension[j],"_lefse.res")))
    lefse_fin<-remove_common_lefse(lefse,gc_lefse)

    limma <- read.csv(gsub(" ","",paste0("Output/LIMMA_score/",dimension[j],"/Bacteria_",analysis[i],"_",dimension[j],"_limma.csv")))
    gc_limma <- read.csv(gsub(" ","",paste0("Output/LIMMA_score/",dimension[j],"/Bacteria_gc_treatment_",dimension[j],"_limma.csv")))
    limma_fin<-remove_common_limma(limma,gc_limma)

    merge_das(baselines_dec_001, lefse_fin, limma_fin, "msHd", paste0("Bacteria_", analysis[i], "_", dimension[j]),paste0(output_folder,dimension[j],"/"))
  }
}
for (i in seq_along(domain)) {
  for (j in seq_along(analysis)){
    for(k in seq_along(status)){
      lefse <- read_tsv(paste0("Output/LEFSE/GC_comp/final_name/", domain[i], "_",analysis[j], "_", status[k], "_lefse.res"))
      gc_lefse <- read_tsv(paste0("Output/LEFSE/GC_comp/final_name/", domain[i], "_gc_treatment_", status[k], "_lefse.res"))
      lefse_fin <- remove_common_lefse(lefse, gc_lefse)

      limma <- read.csv(paste0("Output/LIMMA_score/GC_comp/", domain[i], "_",analysis[j], "_", status[k], "_limma.csv"))
      gc_limma <- read.csv(paste0("Output/LIMMA_score/GC_comp/", domain[i], "_gc_treatment_", status[k], "_limma.csv"))
      limma_fin <- remove_common_limma(limma, gc_limma)
      if(domain[i] == "Bacteria") {
        merge_das(baselines_dec_001, lefse_fin, limma_fin, domain[i], paste0("Bacteria_", analysis[j], "_", status[k]), output_folderGC_comp)
      } else if (domain[i] == "Archaea") {
        merge_das(baselines_decA, lefse_fin, limma_fin, domain[i], paste0("Archaea_", analysis[j], "_", status[k]),output_folderGC_comp)
      } else if (domain[i] == "Eukaryote") {
        merge_das(baselines_decE, lefse_fin, limma_fin, domain[i], paste0("Eukaryote_", analysis[j], "_", status[k]), output_folderGC_comp)
      }
    }
  }
}

