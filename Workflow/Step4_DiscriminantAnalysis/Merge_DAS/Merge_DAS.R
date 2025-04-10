source("utilities.R")
output_folder <- "Output/MERGED_DAS/"
createFolder(output_folder)
output_folderA <- "Output/MERGED_DAS/Archaea/"
output_folderB <- "Output/MERGED_DAS/Bacteria/"
output_folderE <- "Output/MERGED_DAS/Eukaryota/"
output_folderMod <- "Output/MERGED_DAS/Mod/"
createFolder(output_folderA)
createFolder(output_folderB)
createFolder(output_folderE)
createFolder(output_folderMod)

output_folder_HD<- "Output/DAS_MERGED_HD/"
createFolder(output_folder_HD)


baselines_decB <- readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam.rds")
baselines_decA <- readRDS(file = "Output/SUPERVISED_DEC/Archaea_Supervised_decontam.rds")
baselines_decE <- readRDS(file = "Output/SUPERVISED_DEC/Eukaryota_Supervised_decontam.rds")

das_lefseB_MsHd <- read_tsv("Output/LEFSE/final_output/BACT_MsHd.res")
das_lefseB_GC <- read_tsv("Output/LEFSE/final_output/BACT_GC.res")
das_lefseB_GC_Lesion <- read_tsv("Output/LEFSE/final_output/BACT_positive_lesion_burden.res")
das_lefseB_GC_BM_Lesion <- read_tsv("Output/LEFSE/final_output/BACT_positive_bone_marrow_lesions.res")
das_lefseB_GC_Gadolinium <- read_tsv("Output/LEFSE/final_output/BACT_positive_gadolinium_contrast.res")
das_lefseB_GC_Subtentorial <- read_tsv("Output/LEFSE/final_output/BACT_positive_subtentorial_lesions.res")
das_lefseB_NO_GC_Lesion <- read_tsv("Output/LEFSE/final_output/BACT_negative_lesion_burden.res")
das_lefseB_NO_GC_BM_Lesion <- read_tsv("Output/LEFSE/final_output/BACT_negative_bone_marrow_lesions.res")
das_lefseB_NO_GC_Gadolinium <- read_tsv("Output/LEFSE/final_output/BACT_negative_gadolinium_contrast.res")
das_lefseB_NO_GC_Subtentorial <- read_tsv("Output/LEFSE/final_output/BACT_negative_subtentorial_lesions.res")

das_limmaB_MsHd <- read_csv("Output/LIMMA/Bacteria_limma_category.csv")
das_limmaB_GC <- read_csv("Output/LIMMA/Bacteria_limma_gc_treatment.csv")
das_limmaB_GC_Lesion <- read_csv("Output/LIMMA/Bacteria_limma_lesion_burden_positive.csv")
das_limmaB_GC_BM_Lesion <- read_csv("Output/LIMMA/Bacteria_limma_bone_marrow_lesions_positive.csv")
das_limmaB_GC_Gadolinium <- read_csv("Output/LIMMA/Bacteria_limma_gadolinium_contrast_positive.csv")
das_limmaB_GC_Subtentorial <- read_csv("Output/LIMMA/Bacteria_limma_subtentorial_lesions_positive.csv")
das_limmaB_NO_GC_Lesion <- read_csv("Output/LIMMA/Bacteria_limma_lesion_burden_negative.csv")
das_limmaB_NO_GC_BM_Lesion <- read_csv("Output/LIMMA/Bacteria_limma_bone_marrow_lesions_negative.csv")
das_limmaB_NO_GC_Gadolinium <- read_csv("Output/LIMMA/Bacteria_limma_gadolinium_contrast_negative.csv")
das_limmaB_NO_GC_Subtentorial <- read_csv("Output/LIMMA/Bacteria_limma_subtentorial_lesions_negative.csv")

das_lefseA_MsHd <- read_tsv("Output/LEFSE/final_output/ARCH_MsHd.res")
das_lefseA_GC <- read_tsv("Output/LEFSE/final_output/ARCH_GC.res")
das_lefseA_GC_Lesion <- read_tsv("Output/LEFSE/final_output/ARCH_positive_lesion_burden.res")
das_lefseA_GC_BM_Lesion <- read_tsv("Output/LEFSE/final_output/ARCH_positive_bone_marrow_lesions.res")
das_lefseA_GC_Gadolinium <- read_tsv("Output/LEFSE/final_output/ARCH_positive_gadolinium_contrast.res")
das_lefseA_GC_Subtentorial <- read_tsv("Output/LEFSE/final_output/ARCH_positive_subtentorial_lesions.res")
das_lefseA_NO_GC_Lesion <- read_tsv("Output/LEFSE/final_output/ARCH_negative_lesion_burden.res")
das_lefseA_NO_GC_BM_Lesion <- read_tsv("Output/LEFSE/final_output/ARCH_negative_bone_marrow_lesions.res")
das_lefseA_NO_GC_Gadolinium <- read_tsv("Output/LEFSE/final_output/ARCH_negative_gadolinium_contrast.res")
das_lefseA_NO_GC_Subtentorial <- read_tsv("Output/LEFSE/final_output/ARCH_negative_subtentorial_lesions.res")

das_limmaA_MsHd <- read_csv("Output/LIMMA/Archaea_limma_category.csv")
das_limmaA_GC <- read_csv("Output/LIMMA/Archaea_limma_gc_treatment.csv")
das_limmaA_GC_Lesion <- read_csv("Output/LIMMA/Archaea_limma_lesion_burden_positive.csv")
das_limmaA_GC_BM_Lesion <- read_csv("Output/LIMMA/Archaea_limma_bone_marrow_lesions_positive.csv")
das_limmaA_GC_Gadolinium <- read_csv("Output/LIMMA/Archaea_limma_gadolinium_contrast_positive.csv")
das_limmaA_GC_Subtentorial <- read_csv("Output/LIMMA/Archaea_limma_subtentorial_lesions_positive.csv")
das_limmaA_NO_GC_Lesion <- read_csv("Output/LIMMA/Archaea_limma_lesion_burden_negative.csv")
das_limmaA_NO_GC_BM_Lesion <- read_csv("Output/LIMMA/Archaea_limma_bone_marrow_lesions_negative.csv")
das_limmaA_NO_GC_Gadolinium <- read_csv("Output/LIMMA/Archaea_limma_gadolinium_contrast_negative.csv")
das_limmaA_NO_GC_Subtentorial <- read_csv("Output/LIMMA/Archaea_limma_subtentorial_lesions_negative.csv")

das_lefseE_MsHd <- read_tsv("Output/LEFSE/final_output/EUK_MsHd.res")
das_lefseE_GC <- read_tsv("Output/LEFSE/final_output/EUK_GC.res")
das_lefseE_GC_Lesion <- read_tsv("Output/LEFSE/final_output/EUK_positive_lesion_burden.res")
das_lefseE_GC_BM_Lesion <- read_tsv("Output/LEFSE/final_output/EUK_positive_bone_marrow_lesions.res")
das_lefseE_GC_Gadolinium <- read_tsv("Output/LEFSE/final_output/EUK_positive_gadolinium_contrast.res")
das_lefseE_GC_Subtentorial <- read_tsv("Output/LEFSE/final_output/EUK_positive_subtentorial_lesions.res")
das_lefseE_NO_GC_Lesion <- read_tsv("Output/LEFSE/final_output/EUK_negative_lesion_burden.res")
das_lefseE_NO_GC_BM_Lesion <- read_tsv("Output/LEFSE/final_output/EUK_negative_bone_marrow_lesions.res")
das_lefseE_NO_GC_Gadolinium <- read_tsv("Output/LEFSE/final_output/EUK_negative_gadolinium_contrast.res")
das_lefseE_NO_GC_Subtentorial <- read_tsv("Output/LEFSE/final_output/EUK_negative_subtentorial_lesions.res")


das_limmaE_MsHd <- read_csv("Output/LIMMA/Eukaryota_limma_category.csv")
das_limmaE_GC <- read_csv("Output/LIMMA/Eukaryota_limma_gc_treatment.csv")
das_limmaE_GC_Lesion <- read_csv("Output/LIMMA/Eukaryota_limma_lesion_burden_positive.csv")
das_limmaE_GC_BM_Lesion <- read_csv("Output/LIMMA/Eukaryota_limma_bone_marrow_lesions_positive.csv")
das_limmaE_GC_Gadolinium <- read_csv("Output/LIMMA/Eukaryota_limma_gadolinium_contrast_positive.csv")
das_limmaE_GC_Subtentorial <- read_csv("Output/LIMMA/Eukaryota_limma_subtentorial_lesions_positive.csv")
das_limmaE_NO_GC_Lesion <- read_csv("Output/LIMMA/Eukaryota_limma_lesion_burden_negative.csv")
das_limmaE_NO_GC_BM_Lesion <- read_csv("Output/LIMMA/Eukaryota_limma_bone_marrow_lesions_negative.csv")
das_limmaE_NO_GC_Gadolinium <- read_csv("Output/LIMMA/Eukaryota_limma_gadolinium_contrast_negative.csv")
das_limmaE_NO_GC_Subtentorial <- read_csv("Output/LIMMA/Eukaryota_limma_subtentorial_lesions_negative.csv")




lefse_arrayB_pos <- list(
  das_lefseB_GC_Lesion,
  das_lefseB_GC_BM_Lesion,
  das_lefseB_GC_Gadolinium,
  das_lefseB_GC_Subtentorial
)
lefse_arrayB_neg <- list(
  das_lefseB_NO_GC_Lesion,
  das_lefseB_NO_GC_BM_Lesion,
  das_lefseB_NO_GC_Gadolinium,
  das_lefseB_NO_GC_Subtentorial
)

lefse_arrayA_pos <- list(
  das_lefseA_GC_Lesion,
  das_lefseA_GC_BM_Lesion,
  das_lefseA_GC_Gadolinium,
  das_lefseA_GC_Subtentorial
)
lefse_arrayA_neg <- list(
  das_lefseA_NO_GC_Lesion,
  das_lefseA_NO_GC_BM_Lesion,
  das_lefseA_NO_GC_Gadolinium,
  das_lefseA_NO_GC_Subtentorial
)
lefse_arrayE_pos <- list(
  das_lefseE_GC_Lesion,
  das_lefseE_GC_BM_Lesion,
  das_lefseE_GC_Gadolinium,
  das_lefseE_GC_Subtentorial
)
lefse_arrayE_neg <- list(
  das_lefseE_NO_GC_Lesion,
  das_lefseE_NO_GC_BM_Lesion,
  das_lefseE_NO_GC_Gadolinium,
  das_lefseE_NO_GC_Subtentorial
)

limma_arrayB_pos <- list(
  das_limmaB_GC_Lesion,
  das_limmaB_GC_BM_Lesion,
  das_limmaB_GC_Gadolinium,
  das_limmaB_GC_Subtentorial
)
limma_arrayB_neg <- list(
  das_limmaB_NO_GC_Lesion,
  das_limmaB_NO_GC_BM_Lesion,
  das_limmaB_NO_GC_Gadolinium,
  das_limmaB_NO_GC_Subtentorial
)

limma_arrayA_pos <- list(
  das_limmaA_GC_Lesion,
  das_limmaA_GC_BM_Lesion,
  das_limmaA_GC_Gadolinium,
  das_limmaA_GC_Subtentorial
)
limma_arrayA_neg <- list(
  das_limmaA_NO_GC_Lesion,
  das_limmaA_NO_GC_BM_Lesion,
  das_limmaA_NO_GC_Gadolinium,
  das_limmaA_NO_GC_Subtentorial
)
limma_arrayE_pos <- list(
  das_limmaE_GC_Lesion,
  das_limmaE_GC_BM_Lesion,
  das_limmaE_GC_Gadolinium,
  das_limmaE_GC_Subtentorial
)
limma_arrayE_neg <- list(
  das_limmaE_NO_GC_Lesion,
  das_limmaE_NO_GC_BM_Lesion,
  das_limmaE_NO_GC_Gadolinium,
  das_limmaE_NO_GC_Subtentorial
)

#####################LAST_MOD

das_lefseB_Lesion_mod<-read_tsv("Output/modLefse/BACT_Lesion_mod.rds")
das_lefseB_Gadolinium_mod<-read_tsv("Output/modLefse/BACT_Gadolinium_mod.rds")
das_lefseB_Subtentorial_mod<-read_tsv("Output/modLefse/BACT_Subtentorial_mod.rds")
das_lefseB_BM_Lesion_mod<-read_tsv("Output/modLefse/BACT_BM_Lesion_mod.rds")

dasLimma_Lesion_mod<- read_csv("Output/modLimma/Bacteria_limma_lesion_burden_both_mod.csv")
dasLimma_Gadolinium_mod<- read_csv("Output/modLimma/Bacteria_limma_gadolinium_contrast_both_mod.csv")
dasLimma_Subtentorial_mod<- read_csv("Output/modLimma/Bacteria_limma_subtentorial_lesions_both_mod.csv")
dasLimma_BM_Lesion_mod<- read_csv("Output/modLimma/Bacteria_limma_bone_marrow_lesions_both_mod.csv")

##############DAS_ALE 
empty_limma<- read_csv("Output/DAS_ONLY_LEFSE/Empty_limma.csv")


######DAS_Ale 39

baselines_decB39<- readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_dec_elabundance005.rds")
print(baselines_decB39)

dasLefse39_GC<-read_tsv("Output/DAS_ONLY_LEFSE/final_39/BACT_GC39.res")
dasLefse39_Lesion<-read_tsv("Output/DAS_ONLY_LEFSE/final_39/BACT_Lesion39.res")
dasLefse39_Gadolinium<-read_tsv("Output/DAS_ONLY_LEFSE/final_39/BACT_Gadolinium39.res")
dasLefse39_Subtentorial<-read_tsv("Output/DAS_ONLY_LEFSE/final_39/BACT_Subtentorial39.res")
dasLefse39_BM_Lesion<-read_tsv("Output/DAS_ONLY_LEFSE/final_39/BACT_BM_Lesion39.res")

dasLimma_GC_39<-read.csv("Output/DAS_ONLY_LIMMA/39/Bacteria_39_limma_gc_treatment_both.csv")
dasLimma_Lesion_39<-read.csv("Output/DAS_ONLY_LIMMA/39/Bacteria_39_limma_lesion_burden_both.csv")
dasLimma_Gadolinium_39<-read.csv("Output/DAS_ONLY_LIMMA/39/Bacteria_39_limma_gadolinium_contrast_both.csv")
dasLimma_Subtentorial_39<-read.csv("Output/DAS_ONLY_LIMMA/39/Bacteria_39_limma_subtentorial_lesions_both.csv")
dasLimma_BM_Lesion_39<-read.csv("Output/DAS_ONLY_LIMMA/39/Bacteria_39_limma_bone_marrow_lesions_both.csv")


######DAS_Ale 116
baselines_decB116<- readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam.rds")
print(baselines_decB116)
dasLefse116_GC<-read_tsv("Output/DAS_ONLY_LEFSE/final_116/BACT_GC116.res")

dasLefse116_Lesion<-read_tsv("Output/DAS_ONLY_LEFSE/final_116/BACT_Lesion116.res")
dasLefse116_Gadolinium<-read_tsv("Output/DAS_ONLY_LEFSE/final_116/BACT_Gadolinium116.res")
dasLefse116_Subtentorial<-read_tsv("Output/DAS_ONLY_LEFSE/final_116/BACT_Subtentorial116.res")
dasLefse116_BM_Lesion<-read_tsv("Output/DAS_ONLY_LEFSE/final_116/BACT_BM_Lesion116.res")


dasLimma_GC_116<-read.csv("Output/DAS_ONLY_LIMMA/116/Bacteria_116_limma_gc_treatment_both.csv")
dasLimma_Lesion_116<-read.csv("Output/DAS_ONLY_LIMMA/116/Bacteria_116_limma_lesion_burden_both.csv")
dasLimma_Gadolinium_116<-read.csv("Output/DAS_ONLY_LIMMA/116/Bacteria_116_limma_gadolinium_contrast_both.csv")
dasLimma_Subtentorial_116<-read.csv("Output/DAS_ONLY_LIMMA/116/Bacteria_116_limma_subtentorial_lesions_both.csv")
dasLimma_BM_Lesion_116<-read.csv("Output/DAS_ONLY_LIMMA/116/Bacteria_116_limma_bone_marrow_lesions_both.csv")



######DAS_Ale All
baselines_decBAll<- readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam264.rds")
print(baselines_decBAll)
dasLefseAll_GC<-read_tsv("Output/DAS_ONLY_LEFSE/final_264/BACT_GC264.res")
dasLefseAll_Lesion<-read_tsv("Output/DAS_ONLY_LEFSE/final_264/BACT_Lesion264.res")
dasLefseAll_Gadolinium<-read_tsv("Output/DAS_ONLY_LEFSE/final_264/BACT_Gadolinium264.res")
dasLefseAll_Subtentorial<-read_tsv("Output/DAS_ONLY_LEFSE/final_264/BACT_Subtentorial264.res")
dasLefseAll_BM_Lesion<-read_tsv("Output/DAS_ONLY_LEFSE/final_264/BACT_BM_Lesion264.res")

dasLimma_GC_264<-read.csv("Output/DAS_ONLY_LIMMA/264/Bacteria_264_limma_gc_treatment_both.csv")
dasLimma_Lesion_264<-read.csv("Output/DAS_ONLY_LIMMA/264/Bacteria_264_limma_lesion_burden_both.csv")
dasLimma_Gadolinium_264<-read.csv("Output/DAS_ONLY_LIMMA/264/Bacteria_264_limma_gadolinium_contrast_both.csv")
dasLimma_Subtentorial_264<-read.csv("Output/DAS_ONLY_LIMMA/264/Bacteria_264_limma_subtentorial_lesions_both.csv")
dasLimma_BM_Lesion_264<-read.csv("Output/DAS_ONLY_LIMMA/264/Bacteria_264_limma_bone_marrow_lesions_both.csv")


baselines_dec167<- readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontamMidWay.rds")
print(baselines_dec167)
dasLefse167_GC<-read_tsv("Output/DAS_ONLY_LEFSE/final_167/BACT_GC167.res")
dasLefse167_Lesion<-read_tsv("Output/DAS_ONLY_LEFSE/final_167/BACT_Lesion167.res")
dasLefse167_Gadolinium<-read_tsv("Output/DAS_ONLY_LEFSE/final_167/BACT_Gadolinium167.res")
dasLefse167_Subtentorial<-read_tsv("Output/DAS_ONLY_LEFSE/final_167/BACT_Subtentorial167.res")
dasLefse167_BM_Lesion<-read_tsv("Output/DAS_ONLY_LEFSE/final_167/BACT_BM_Lesion167.res")

dasLimma_GC_167<-read.csv("Output/DAS_ONLY_LIMMA/167/Bacteria_167_limma_gc_treatment_both.csv")
dasLimma_Lesion_167<-read.csv("Output/DAS_ONLY_LIMMA/167/Bacteria_167_limma_lesion_burden_both.csv")
dasLimma_Gadolinium_167<-read.csv("Output/DAS_ONLY_LIMMA/167/Bacteria_167_limma_gadolinium_contrast_both.csv")
dasLimma_Subtentorial_167<-read.csv("Output/DAS_ONLY_LIMMA/167/Bacteria_167_limma_subtentorial_lesions_both.csv")
dasLimma_BM_Lesion_167<-read.csv("Output/DAS_ONLY_LIMMA/167/Bacteria_167_limma_bone_marrow_lesions_both.csv")

output_folder39 <- "Output/DAS_ONLY_LEFSE/39RDS/"
output_folder116 <- "Output/DAS_ONLY_LEFSE/116RDS/"
output_folder167 <- "Output/DAS_ONLY_LEFSE/167RDS/"
output_folder264 <- "Output/DAS_ONLY_LEFSE/264RDS/"

fused_folder39 <- "Output/MERGED_DAS/39/"
fused_folder116 <- "Output/MERGED_DAS/116/"
fused_folder167 <- "Output/MERGED_DAS/167/"
fused_folder264 <- "Output/MERGED_DAS/264/"

createFolder(output_folder39)
createFolder(output_folder116)
createFolder(output_folder264)
createFolder(output_folder167)

createFolder(fused_folder39)
createFolder(fused_folder116)
createFolder(fused_folder264)
createFolder(fused_folder167)
empty_limma<- read_csv("Output/DAS_ONLY_LEFSE/Empty_limma.csv")

merge_das(baselines_decBAll,remove_common(read_tsv("Output/DAS_ONLY_LEFSE_HD/264_Final/BACT_hd_264.res"),dasLefseAll_GC) ,remove_common_limma(read.csv("Output/DAS_ONLY_LIMMA_HD/264/Bacteria_hd_264_limma_category_both.csv"),dasLimma_GC_264) , "All", "Bacteria_hd_264", output_folder_HD)
merge_das(baselines_decB39,remove_common(read_tsv("Output/DAS_ONLY_LEFSE_HD/39_Final/BACT_hd_39.res"),dasLefse39_GC) ,remove_common_limma(read.csv("Output/DAS_ONLY_LIMMA_HD/39/Bacteria_hd_39_limma_category_both.csv"),dasLimma_GC_39), "All", "Bacteria_hd_39", output_folder_HD)
aaaProva<-baselines_decBAll@tax_table
aaaProva<-as.data.frame(aaaProva)
aaaProva<-aaaProva$Species
lenght_array=c("39","116","264")

hello<-remove_common_limma(dasLimma_Lesion_39,dasLimma_GC_39)

remove_common(das)

remove_common<- function(tab1,tab2){
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


bm <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontamOLD.rds")
print(bm)
sub <- readRDS("Output/DAS_ONLY_LEFSE/ALLRDS/Bacteria_SubtentorialAll_merged.rds")

h1 <- bm@tax_table
h2 <- sub@tax_table


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

merge_das(baselines_decB, das_lefseB_MsHd, das_limmaB_MsHd, "msHd", "Bacteria_MsHd", output_folderB)
merge_das(baselines_decA, das_lefseA_MsHd, das_limmaA_MsHd, "msHd", "Archaea_MsHd", output_folderA)
merge_das(baselines_decE, das_lefseE_MsHd, das_limmaE_MsHd, "msHd", "Eukaryota_MsHd", output_folderE)
merge_das(baselines_decB, das_lefseB_GC, das_limmaB_GC, "GC", "Bacteria_GC", output_folderB)
merge_das(baselines_decA, das_lefseA_GC, das_limmaA_GC, "GC", "Archaea_GC", output_folderA)
merge_das(baselines_decE, das_lefseE_GC, das_limmaE_GC, "GC", "Eukaryota_GC", output_folderE)
array_name <- c("Lesion", "BM_Lesion", "Gadolinium", "Subtentorial")
for (i in seq_along(lefse_arrayB_pos)) {
  merge_das(baselines_decB, lefse_arrayB_pos[[i]], limma_arrayB_pos[[i]], "positive", paste0("Bacteria_GC_", array_name[i]), output_folderB)
  merge_das(baselines_decA, lefse_arrayA_pos[[i]], limma_arrayA_pos[[i]], "positive", paste0("Archaea_GC_", array_name[i]), output_folderA)
  merge_das(baselines_decE, lefse_arrayE_pos[[i]], limma_arrayE_pos[[i]], "positive", paste0("Eukaryota_GC_", array_name[i]), output_folderE)
}

for (i in seq_along(lefse_arrayB_neg)) {
  merge_das(baselines_decB, lefse_arrayB_neg[[i]], limma_arrayB_neg[[i]], "negative", paste0("Bacteria_NO_GC_", array_name[i]), output_folderB)
  merge_das(baselines_decA, lefse_arrayA_neg[[i]], limma_arrayA_neg[[i]], "negative", paste0("Archaea_NO_GC_", array_name[i]), output_folderA)
  merge_das(baselines_decE, lefse_arrayE_neg[[i]], limma_arrayE_neg[[i]], "negative", paste0("Eukaryota_NO_GC_", array_name[i]), output_folderE)
}

merge_das(baselines_decBAll, remove_common(dasLefse39_Lesion,dasLefse39_GC), remove_common_limma(dasLimma_Lesion_39,dasLimma_GC_39), "GC", "Bacteria_Lesion39", fused_folder39)
merge_das(baselines_decBAll, remove_common(dasLefse39_Gadolinium,dasLefse39_GC), remove_common_limma(dasLimma_Gadolinium_39,dasLimma_GC_39), "GC", "Bacteria_Gadolinium39", fused_folder39)
merge_das(baselines_decBAll, remove_common(dasLefse39_Subtentorial,dasLefse39_GC), remove_common_limma(dasLimma_Subtentorial_39,dasLimma_GC_39), "GC", "Bacteria_Subtentorial39", fused_folder39)
merge_das(baselines_decBAll, remove_common(dasLefse39_BM_Lesion,dasLefse39_GC), remove_common_limma(dasLimma_BM_Lesion_39,dasLimma_GC_39), "GC", "Bacteria_BM_Lesion39", fused_folder39)

merge_das(baselines_decBAll, remove_common(dasLefse116_Lesion,dasLefse116_GC), remove_common_limma(dasLimma_Lesion_116,dasLimma_GC_116), "GC", "Bacteria_Lesion116", fused_folder116)
merge_das(baselines_decBAll, remove_common(dasLefse116_Gadolinium,dasLefse116_GC), remove_common_limma(dasLimma_Gadolinium_116,dasLimma_GC_116), "GC", "Bacteria_Gadolinium116", fused_folder116)
merge_das(baselines_decBAll, remove_common(dasLefse116_Subtentorial,dasLefse116_GC), remove_common_limma(dasLimma_Subtentorial_116,dasLimma_GC_116), "GC", "Bacteria_Subtentorial116", fused_folder116)
merge_das(baselines_decBAll, remove_common(dasLefse116_BM_Lesion,dasLefse116_GC), remove_common_limma(dasLimma_BM_Lesion_116,dasLimma_GC_116), "GC", "Bacteria_BM_Lesion116", fused_folder116)

merge_das(baselines_decBAll, remove_common(dasLefse167_Lesion,dasLefse167_GC), remove_common_limma(dasLimma_Lesion_167,dasLimma_GC_167), "GC", "Bacteria_Lesion167", fused_folder167)
merge_das(baselines_decBAll, remove_common(dasLefse167_Gadolinium,dasLefse167_GC), remove_common_limma(dasLimma_Gadolinium_167,dasLimma_GC_167), "GC", "Bacteria_Gadolinium167", fused_folder167)
merge_das(baselines_decBAll, remove_common(dasLefse167_Subtentorial,dasLefse167_GC), remove_common_limma(dasLimma_Subtentorial_167,dasLimma_GC_167), "GC", "Bacteria_Subtentorial167", fused_folder167)
merge_das(baselines_decBAll, remove_common(dasLefse167_BM_Lesion,dasLefse167_GC), remove_common_limma(dasLimma_BM_Lesion_167,dasLimma_GC_167), "GC", "Bacteria_BM_Lesion167", fused_folder167)

merge_das(baselines_decBAll, remove_common(dasLefseAll_Lesion,dasLefseAll_GC), remove_common_limma(dasLimma_Lesion_264,dasLimma_GC_264), "GC", "Bacteria_LesionAll", fused_folder264)
merge_das(baselines_decBAll, remove_common(dasLefseAll_Gadolinium,dasLefseAll_GC), remove_common_limma(dasLimma_Gadolinium_264,dasLimma_GC_264), "GC", "Bacteria_GadoliniumAll", fused_folder264)
merge_das(baselines_decBAll, remove_common(dasLefseAll_Subtentorial,dasLefseAll_GC),  remove_common_limma(dasLimma_Subtentorial_264,dasLimma_GC_264), "GC", "Bacteria_SubtentorialAll", fused_folder264)
merge_das(baselines_decBAll, remove_common(dasLefseAll_BM_Lesion,dasLefseAll_GC), remove_common_limma(dasLimma_BM_Lesion_264,dasLimma_GC_264), "GC", "Bacteria_BM_LesionAll", fused_folder264)


merge_das(baselines_decBAll, remove_common(dasLefse39_Lesion,dasLefse39_GC), empty_limma, "GC", "Bacteria_Lesion39", output_folder39)
merge_das(baselines_decBAll, remove_common(dasLefse39_Gadolinium,dasLefse39_GC), empty_limma, "GC", "Bacteria_Gadolinium39", output_folder39)
merge_das(baselines_decBAll, remove_common(dasLefse39_Subtentorial,dasLefse39_GC), empty_limma, "GC", "Bacteria_Subtentorial39", output_folder39)
merge_das(baselines_decBAll, remove_common(dasLefse39_BM_Lesion,dasLefse39_GC), empty_limma, "GC", "Bacteria_BM_Lesion39", output_folder39)

merge_das(baselines_decBAll, remove_common(dasLefse116_Lesion,dasLefse116_GC),empty_limma, "GC", "Bacteria_Lesion116", output_folder116)
merge_das(baselines_decBAll, remove_common(dasLefse116_Gadolinium,dasLefse116_GC),empty_limma, "GC", "Bacteria_Gadolinium116", output_folder116)
merge_das(baselines_decBAll, remove_common(dasLefse116_Subtentorial,dasLefse116_GC), empty_limma, "GC", "Bacteria_Subtentorial116", output_folder116)
merge_das(baselines_decBAll, remove_common(dasLefse116_BM_Lesion,dasLefse116_GC), empty_limma, "GC", "Bacteria_BM_Lesion116", output_folder116)

merge_das(baselines_decBAll, remove_common(dasLefseAll_Lesion,dasLefseAll_GC), empty_limma, "GC", "Bacteria_LesionAll", output_folder264)
merge_das(baselines_decBAll, remove_common(dasLefseAll_Gadolinium,dasLefseAll_GC), empty_limma, "GC", "Bacteria_GadoliniumAll", output_folder264)
merge_das(baselines_decBAll, remove_common(dasLefseAll_Subtentorial,dasLefseAll_GC), empty_limma, "GC", "Bacteria_SubtentorialAll", output_folder264)
merge_das(baselines_decBAll, remove_common(dasLefseAll_BM_Lesion,dasLefseAll_GC), empty_limma, "GC", "Bacteria_BM_LesionAll", output_folder264)


merge_das(baselines_decB, das_lefseB_Lesion_mod, dasLimma_Lesion_mod, "GC", "Bacteria_Lesion_mod", output_folderMod)
merge_das(baselines_decB, das_lefseB_Gadolinium_mod, dasLimma_Gadolinium_mod, "GC", "Bacteria_Gadolinium_mod", output_folderMod)
merge_das(baselines_decB, das_lefseB_Subtentorial_mod, dasLimma_Subtentorial_mod, "GC", "Bacteria_Subtentorial_mod", output_folderMod)
merge_das(baselines_decB, das_lefseB_BM_Lesion_mod, dasLimma_BM_Lesion_mod, "GC", "Bacteria_BM_Lesion_mod", output_folderMod)