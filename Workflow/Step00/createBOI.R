
createBOI <- function(rds_path, species_path, folder, name) {
    species_tab <- read.csv(species_path)
    species_list <- species_tab$Species
    rds <- readRDS(rds_path)
    metaData <- as.data.frame(as.matrix(sample_data(rds)))
    keep_samples <- rownames(metaData)[metaData$gc_treatment %in% c("positive", "negative")]
    rds <- prune_samples(keep_samples, rds)
    tax <- as.data.frame(as.matrix(tax_table(rds)))
    tax$GenusSpecies <- paste(tax$Genus, tax$Species, sep = "_")

    keep_taxa <- rownames(tax)[tax$GenusSpecies %in% species_list]
    rds <- prune_taxa(keep_taxa, rds)
    if (!dir.exists(folder)) dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    saveRDS(rds, file = file.path(folder, paste0(name,".rds")))
    print(rds)
}


rds_001 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
species_tab <- read.csv("Output/MAASLIN3/001_disc/001_gadolinium_contrast_onlygc_plus_gc/gadolinium_contrast_MAASLIN3_SR_features_001.csv")

createBOI("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
          "Output/MAASLIN3/001_disc/001_gadolinium_contrast_onlygc_plus_gc/gadolinium_contrast_MAASLIN3_SR_features_001.csv",
          "Output/BOI/001/",
          "Bacteria_gadolinium_001")
createBOI("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
          "Output/MAASLIN3/001_disc/001_lesion_burden_onlygc_plus_gc/lesion_burden_MAASLIN3_SR_features_001.csv",
          "Output/BOI/001/",
          "Bacteria_lesion_001")
createBOI("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
          "Output/MAASLIN3/001_disc/001_spinal_cord_lesion_onlygc_plus_gc/spinal_cord_lesion_MAASLIN3_SR_features_001.csv",
          "Output/BOI/001/",
          "Bacteria_spinal_cord_001")
createBOI("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
          "Output/MAASLIN3/001_disc/001_subtentorial_lesions_onlygc_plus_gc/subtentorial_lesions_MAASLIN3_SR_features_001.csv",
          "Output/BOI/001/",
          "Bacteria_subtentorial_lesions_001")

createBOI("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
          "Output/MAASLIN3/01_disc/01_gadolinium_contrast_onlygc_T/gadolinium_contrast_MAASLIN3_SR_features_01.csv",
          "Output/BOI/01/",
          "Bacteria_gadolinium_01")
createBOI("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
          "Output/MAASLIN3/01_disc/01_lesion_burden_onlygc_T/lesion_burden_MAASLIN3_SR_features_01.csv",
          "Output/BOI/01/",
          "Bacteria_lesion_01")
createBOI("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
          "Output/MAASLIN3/01_disc/01_spinal_cord_lesion_onlygc_T/spinal_cord_lesion_MAASLIN3_SR_features_01.csv",
          "Output/BOI/01/",
          "Bacteria_spinal_cord_01")
createBOI("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
          "Output/MAASLIN3/01_disc/01_subtentorial_lesions_onlygc_T/subtentorial_lesions_MAASLIN3_SR_features_01.csv",
          "Output/BOI/01/",
          "Bacteria_subtentorial_lesions_01")

createBOI("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
          "Output/MAASLIN3/05_disc/05_gadolinium_contrast_onlygc_T/gadolinium_contrast_MAASLIN3_SR_features_05.csv",
          "Output/BOI/05/",
          "Bacteria_gadolinium_05")
createBOI("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
          "Output/MAASLIN3/05_disc/05_lesion_burden_onlygc_T/lesion_burden_MAASLIN3_SR_features_05.csv",
          "Output/BOI/05/",
          "Bacteria_lesion_05")
createBOI("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
          "Output/MAASLIN3/05_disc/05_spinal_cord_lesion_onlygc_T/spinal_cord_lesion_MAASLIN3_SR_features_05.csv",
          "Output/BOI/05/",
          "Bacteria_spinal_cord_05")
createBOI("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
          "Output/MAASLIN3/05_disc/05_subtentorial_lesions_onlygc_T/subtentorial_lesions_MAASLIN3_SR_features_05.csv",
          "Output/BOI/05/",
          "Bacteria_subtentorial_lesions_05")
