source("utilities.R")
output_folder = "Output/LEFSE/step0_score/"
output_folder2= "Output/LEFSE/step0v2_score/"
output_folder3= "Output/LEFSE/step0v116/"
output_folder4= "Output/LEFSE/step0v39/"
output_folder5= "Output/LEFSE/step0v167/"
output_folder6= "Output/LEFSE/step0v264/"
createFolder(output_folder)
createFolder(output_folder2)
createFolder(output_folder3)
createFolder(output_folder4)
createFolder(output_folder5)
createFolder(output_folder6)

output_folder_hd= "Output/DAS_ONLY_LEFSE_HD/"
createFolder(output_folder_hd)
output_folder_HD_256= "Output/DAS_ONLY_LEFSE_HD/256_V0/"
output_folder_HD_39= "Output/DAS_ONLY_LEFSE_HD/39_V0/"

createFolder(output_folder_HD_256)
createFolder(output_folder_HD_39)


ARCH_Supervised_decontam <- readRDS("Output/SUPERVISED_DEC/Archaea_Supervised_decontam.rds")
EUK_Supervised_decontam <- readRDS("Output/SUPERVISED_DEC/Eukaryota_Supervised_decontam.rds")
BACT_Supervised_decontam116<- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam.rds")
BACT_Supervised_decontam39<- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_dec_elabundance005.rds")
print(BACT_Supervised_decontam39)
BACT_Supervised_decontam167<- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontamMidWay.rds")
BACT_Supervised_decontam264<- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam264.rds")
print(BACT_Supervised_decontam264)
metadataB <- read_delim("input/Bacteria_alpha_metadata.csv", delim = ",", escape_double = FALSE,  trim_ws = TRUE)
metadataA <- read_delim("input/Archaea_alpha_metadata.csv", delim = ",", escape_double = FALSE,  trim_ws = TRUE)
metadataE <- read_delim("input/Eukaryota_alpha_metadata.csv", delim = ",", escape_double = FALSE,  trim_ws = TRUE)



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
  metadata$bone_marrow_lesions=as.factor(metadata$bone_marrow_lesions)
  levels(metadata$bone_marrow_lesions) =c("BM_low","BM_high")
  metadata$gadolinium_contrast=as.factor(metadata$gadolinium_contrast)
  levels(metadata$gadolinium_contrast) =c("NoActive","Active")
  metadata$subtentorial_lesions=as.factor(metadata$subtentorial_lesions)
  levels(metadata$subtentorial_lesions) =c("No","Yes")


  metaSel <- metadata %>%
    select(id, category, gc_treatment, lesion_burden, bone_marrow_lesions, gadolinium_contrast, subtentorial_lesions) %>%
    as.data.frame()
  colnames(metaSel) <- c("rownames", "category", "gc_treatment", "lesion_burden", "bone_marrow_lesions", "gadolinium_contrast", "subtentorial_lesions")

  comb_t <- comb_t %>%
    mutate(rownames = rownames(comb_t))

  joined_df <- comb_t %>%
    inner_join(metaSel, by = "rownames")
  

  total <- c("rownames", "category", "gc_treatment", "lesion_burden", "bone_marrow_lesions", "gadolinium_contrast", "subtentorial_lesions")
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

arry1 <- c("positive", "negative")
array2 <- c("lesion_burden", "bone_marrow_lesions", "gadolinium_contrast", "subtentorial_lesions")

for (i in seq_along(arry1)) {
  for (j in seq_along(array2)) {

   # fileName_ARCH <- paste0("ARCH_", arry1[i], "_", array2[j])
   # fileName_EUK <- paste0("EUK_", arry1[i], "_", array2[j])
    fileName_BACT <- paste0("BACT_", arry1[i], "_", array2[j])
    

    #generate.LEFSE1(ARCH_Supervised_decontam, metadataA, array2[j], arry1[i], fileName_ARCH,output_folder )
    #generate.LEFSE1(EUK_Supervised_decontam, metadataE, array2[j], arry1[i], fileName_EUK,output_folder )
    generate.LEFSE1(BACT_Supervised_decontam, metadataB, array2[j], arry1[i], fileName_BACT,output_folder )
  }
}
 
    generate.LEFSE1(BACT_Supervised_decontam, metadataB, "gc_treatment", "both", "BACT_GC",output_folder )
    #generate.LEFSE1(EUK_Supervised_decontam, metadataE, "gc_treatment", both, "EUK_GC",output_folder )
    #generate.LEFSE1(ARCH_Supervised_decontam, metadataA, "gc_treatment", both, "ARCH_GC",output_folder )

    generate.LEFSE1(BACT_Supervised_decontam, metadataB, "category", "both", "BACT_MsHd",output_folder )
    #generate.LEFSE1(EUK_Supervised_decontam, metadataE, "category", both, "EUK_MsHd",output_folder )
    #generate.LEFSE1(ARCH_Supervised_decontam, metadataA, "category", both, "ARCH_MsHd",output_folder )

    generate.LEFSE1(BACT_Supervised_decontam116, metadataB, "lesion_burden", "both", "BACT_Lesion116",output_folder3 )
    generate.LEFSE1(BACT_Supervised_decontam116, metadataB, "bone_marrow_lesions", "both", "BACT_BM_Lesion116",output_folder3 )
    generate.LEFSE1(BACT_Supervised_decontam116, metadataB, "gadolinium_contrast", "both", "BACT_Gadolinium116",output_folder3 )
    generate.LEFSE1(BACT_Supervised_decontam116, metadataB, "subtentorial_lesions", "both", "BACT_Subtentorial116",output_folder3 )
    generate.LEFSE1(BACT_Supervised_decontam116, metadataB, "gc_treatment", "both", "BACT_GC116",output_folder3)

    generate.LEFSE1(BACT_Supervised_decontam39, metadataB, "lesion_burden", "both", "BACT_Lesion39",output_folder4 )
    generate.LEFSE1(BACT_Supervised_decontam39, metadataB, "bone_marrow_lesions", "both", "BACT_BM_Lesion39",output_folder4 )
    generate.LEFSE1(BACT_Supervised_decontam39, metadataB, "gadolinium_contrast", "both", "BACT_Gadolinium39",output_folder4 )
    generate.LEFSE1(BACT_Supervised_decontam39, metadataB, "subtentorial_lesions", "both", "BACT_Subtentorial39",output_folder4)
    generate.LEFSE1(BACT_Supervised_decontam39, metadataB, "gc_treatment", "both", "BACT_GC39",output_folder4)


    generate.LEFSE1(BACT_Supervised_decontam167, metadataB, "lesion_burden", "both", "BACT_Lesion167",output_folder5 )
    generate.LEFSE1(BACT_Supervised_decontam167, metadataB, "bone_marrow_lesions", "both", "BACT_BM_Lesion167",output_folder5 )
    generate.LEFSE1(BACT_Supervised_decontam167, metadataB, "gadolinium_contrast", "both", "BACT_Gadolinium167",output_folder5 )
    generate.LEFSE1(BACT_Supervised_decontam167, metadataB, "subtentorial_lesions", "both", "BACT_Subtentorial167",output_folder5)
    generate.LEFSE1(BACT_Supervised_decontam167, metadataB, "gc_treatment", "both", "BACT_GC167",output_folder5)

  generate.LEFSE1(BACT_Supervised_decontam264, metadataB, "lesion_burden", "both", "BACT_Lesion264",output_folder6)
  generate.LEFSE1(BACT_Supervised_decontam264, metadataB, "bone_marrow_lesions", "both", "BACT_BM_Lesion264",output_folder6)
  generate.LEFSE1(BACT_Supervised_decontam264, metadataB, "gadolinium_contrast", "both", "BACT_Gadolinium264",output_folder6)
  generate.LEFSE1(BACT_Supervised_decontam264, metadataB, "subtentorial_lesions", "both", "BACT_Subtentorial264",output_folder6)
  generate.LEFSE1(BACT_Supervised_decontam264, metadataB, "gc_treatment", "both", "BACT_GC264",output_folder6)

  generate.LEFSE1(BACT_Supervised_decontam264, metadataB, "category", "category", "BACT_hd_264",output_folder_HD_256)
  generate.LEFSE1(BACT_Supervised_decontam39, metadataB, "category", "category", "BACT_hd_39",output_folder_HD_39)
