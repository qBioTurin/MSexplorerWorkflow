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
    saveRDS(rds, file = file.path(folder, paste0(name, ".rds")))
    print(rds)
}

createBOI2 <- function(rds_path, species_path, folder, name, filter = TRUE) {
    species_tab1 <- read.csv(species_path)
    species_tab <- species_tab1[species_tab1$Presence == "both", ]
    species_list <- species_tab$Species
    rds <- readRDS(rds_path)
    metaData <- as.data.frame(as.matrix(sample_data(rds)))
    if (filter) {
        keep_samples <- rownames(metaData)[metaData$gc_treatment %in% c("positive", "negative")]
    } else {
        keep_samples <- rownames(metaData)[metaData$category %in% c("MS", "HEALTHY")]
    }
    rds <- prune_samples(keep_samples, rds)
    tax <- as.data.frame(as.matrix(tax_table(rds)))
    tax$GenusSpecies <- paste(tax$Genus, tax$Species, sep = "_")

    keep_taxa <- rownames(tax)[tax$GenusSpecies %in% species_list]
    rds <- prune_taxa(keep_taxa, rds)
    if (!dir.exists(folder)) dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    saveRDS(rds, file = file.path(folder, paste0(name, ".rds")))
    print(rds)
}

createBOI_ancomBC2<- function(
    rds_path, species_path, folder, name) {
    species_tab <- read.csv(species_path)
    species_list <- species_tab$Genus_Species
    rds <- readRDS(rds_path)
    metaData <- as.data.frame(as.matrix(sample_data(rds)))
    keep_samples <- rownames(metaData)[metaData$gc_treatment %in% c("positive", "negative")]
    rds <- prune_samples(keep_samples, rds)
    tax <- as.data.frame(as.matrix(tax_table(rds)))
    tax$Genus_Species <- paste(tax$Genus, tax$Species, sep = " ")
    keep_taxa <- rownames(tax)[tax$Genus_Species %in% species_list]
    rds <- prune_taxa(keep_taxa, rds)
    if (!dir.exists(folder)) dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    saveRDS(rds, file = file.path(folder, paste0(name, ".rds")))
    print(rds)
}


rds_001 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
species_tab <- read.csv("Output/MAASLIN3/001_disc/001_gadolinium_contrast_onlygc_plus_gc/gadolinium_contrast_MAASLIN3_SR_features_001.csv")

createBOI(
    "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    "Output/MAASLIN3/001_disc/001_gadolinium_contrast_onlygc_plus_gc/gadolinium_contrast_MAASLIN3_SR_features_001.csv",
    "Output/BOI/001/",
    "Bacteria_gadolinium_001"
)
createBOI(
    "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    "Output/MAASLIN3/001_disc/001_lesion_burden_onlygc_plus_gc/lesion_burden_MAASLIN3_SR_features_001.csv",
    "Output/BOI/001/",
    "Bacteria_lesion_001"
)
createBOI(
    "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    "Output/MAASLIN3/001_disc/001_spinal_cord_lesion_onlygc_plus_gc/spinal_cord_lesion_MAASLIN3_SR_features_001.csv",
    "Output/BOI/001/",
    "Bacteria_spinal_cord_001"
)
createBOI(
    "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    "Output/MAASLIN3/001_disc/001_subtentorial_lesions_onlygc_plus_gc/subtentorial_lesions_MAASLIN3_SR_features_001.csv",
    "Output/BOI/001/",
    "Bacteria_subtentorial_lesions_001"
)

createBOI(
    "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
    "Output/MAASLIN3/01_disc/01_gadolinium_contrast_onlygc_T/gadolinium_contrast_MAASLIN3_SR_features_01.csv",
    "Output/BOI/01/",
    "Bacteria_gadolinium_01"
)
createBOI(
    "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
    "Output/MAASLIN3/01_disc/01_lesion_burden_onlygc_T/lesion_burden_MAASLIN3_SR_features_01.csv",
    "Output/BOI/01/",
    "Bacteria_lesion_01"
)
createBOI(
    "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
    "Output/MAASLIN3/01_disc/01_spinal_cord_lesion_onlygc_T/spinal_cord_lesion_MAASLIN3_SR_features_01.csv",
    "Output/BOI/01/",
    "Bacteria_spinal_cord_01"
)
createBOI(
    "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
    "Output/MAASLIN3/01_disc/01_subtentorial_lesions_onlygc_T/subtentorial_lesions_MAASLIN3_SR_features_01.csv",
    "Output/BOI/01/",
    "Bacteria_subtentorial_lesions_01"
)

createBOI(
    "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
    "Output/MAASLIN3/05_disc/05_gadolinium_contrast_onlygc_T/gadolinium_contrast_MAASLIN3_SR_features_05.csv",
    "Output/BOI/05/",
    "Bacteria_gadolinium_05"
)
createBOI(
    "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
    "Output/MAASLIN3/05_disc/05_lesion_burden_onlygc_T/lesion_burden_MAASLIN3_SR_features_05.csv",
    "Output/BOI/05/",
    "Bacteria_lesion_05"
)
createBOI(
    "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
    "Output/MAASLIN3/05_disc/05_spinal_cord_lesion_onlygc_T/spinal_cord_lesion_MAASLIN3_SR_features_05.csv",
    "Output/BOI/05/",
    "Bacteria_spinal_cord_05"
)
createBOI(
    "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
    "Output/MAASLIN3/05_disc/05_subtentorial_lesions_onlygc_T/subtentorial_lesions_MAASLIN3_SR_features_05.csv",
    "Output/BOI/05/",
    "Bacteria_subtentorial_lesions_05"
)



createBOI2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    species_path = "Output/compare_DAS_MAASLIN3/gadolinium_contrast_001/gadolinium_contrast_001.csv",
    folder = "Output/BOI2/001_unioned/",
    name = "Bacteria_gadolinium_001"
)
createBOI2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
    species_path = "Output/compare_DAS_MAASLIN3/gadolinium_contrast_01/gadolinium_contrast_01.csv",
    folder = "Output/BOI2/01_unioned/",
    name = "Bacteria_gadolinium_01"
)
createBOI2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
    species_path = "Output/compare_DAS_MAASLIN3/gadolinium_contrast_05/gadolinium_contrast_05.csv",
    folder = "Output/BOI2/05_unioned/",
    name = "Bacteria_gadolinium_05"
)

createBOI2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    species_path = "Output/compare_DAS_MAASLIN3/lesion_burden_001/lesion_burden_001.csv",
    folder = "Output/BOI2/001_unioned/",
    name = "Bacteria_lesion_001"
)
createBOI2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
    species_path = "Output/compare_DAS_MAASLIN3/lesion_burden_01/lesion_burden_01.csv",
    folder = "Output/BOI2/01_unioned/",
    name = "Bacteria_lesion_01"
)
createBOI2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
    species_path = "Output/compare_DAS_MAASLIN3/lesion_burden_05/lesion_burden_05.csv",
    folder = "Output/BOI2/05_unioned/",
    name = "Bacteria_lesion_05"
)

createBOI2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    species_path = "Output/compare_DAS_MAASLIN3/spinal_cord_lesion_001/spinal_cord_lesion_001.csv",
    folder = "Output/BOI2/001_unioned/",
    name = "Bacteria_spinal_cord_001"
)
createBOI2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
    species_path = "Output/compare_DAS_MAASLIN3/spinal_cord_lesion_01/spinal_cord_lesion_01.csv",
    folder = "Output/BOI2/01_unioned/",
    name = "Bacteria_spinal_cord_01"
)
createBOI2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
    species_path = "Output/compare_DAS_MAASLIN3/spinal_cord_lesion_05/spinal_cord_lesion_05.csv",
    folder = "Output/BOI2/05_unioned/",
    name = "Bacteria_spinal_cord_05"
)

createBOI2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    species_path = "Output/compare_DAS_MAASLIN3/subtentorial_lesions_001/subtentorial_lesions_001.csv",
    folder = "Output/BOI2/001_unioned/",
    name = "Bacteria_subtentorial_001"
)
createBOI2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
    species_path = "Output/compare_DAS_MAASLIN3/subtentorial_lesions_01/subtentorial_lesions_01.csv",
    folder = "Output/BOI2/01_unioned/",
    name = "Bacteria_subtentorial_01"
)
createBOI2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
    species_path = "Output/compare_DAS_MAASLIN3/subtentorial_lesions_05/subtentorial_lesions_05.csv",
    folder = "Output/BOI2/05_unioned/",
    name = "Bacteria_subtentorial_05"
)



createBOI2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    species_path = "Output/compare_DAS_MAASLIN3/GC_Bacteria/GC_Bacteria.csv",
    folder = "Output/BOI2/GC/",
    name = "Bacteria_GC"
)
createBOI2(
    rds_path = "Output/SUPERVISED_DEC/Eukaryote_Supervised_decontam0.001.rds",
    species_path = "Output/compare_DAS_MAASLIN3/GC_Eukaryote/GC_Eukaryote.csv",
    folder = "Output/BOI2/GC/",
    name = "Eukaryote_GC"
)
createBOI2(
    rds_path = "Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.001.rds",
    species_path = "Output/compare_DAS_MAASLIN3/GC_Archaea/GC_Archaea.csv",
    folder = "Output/BOI2/GC/",
    name = "Archaea_GC"
)

createBOI2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    species_path = "Output/compare_DAS_MAASLIN3/MS_HD_Bacteria/MS_HD_Bacteria.csv",
    folder = "Output/BOI2/MSHD/",
    name = "Bacteria_MSHD",
    filter = FALSE
)
createBOI2(
    rds_path = "Output/SUPERVISED_DEC/Eukaryote_Supervised_decontam0.001.rds",
    species_path = "Output/compare_DAS_MAASLIN3/MS_HD_Eukaryote/MS_HD_Eukaryote.csv",
    folder = "Output/BOI2/MSHD/",
    name = "Eukaryote_MSHD",
    filter = FALSE
)

createBOI2(
    rds_path = "Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.001.rds",
    species_path = "Output/compare_DAS_MAASLIN3/MS_HD_Archaea/MS_HD_Archaea.csv",
    folder = "Output/BOI2/MSHD/",
    name = "Archaea_MSHD",
    filter = FALSE
)


createBOI_ancomBC2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    species_path = "Output/AncomBC2_p/001_disc/001_subtentorial_lesions_onlygc/subtentorial_lesions_ANCOMBC2_significant_results.csv",
    folder = "Output/BOI_AncomBC2/001/",
    name = "Bacteria_subtentorial_lesions_001"
)
createBOI_ancomBC2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    species_path = "Output/AncomBC2_p/001_disc/001_spinal_cord_lesion_onlygc/spinal_cord_lesion_ANCOMBC2_significant_results.csv",
    folder = "Output/BOI_AncomBC2/001/",
    name = "Bacteria_spinal_cord_lesion_001"
)
createBOI_ancomBC2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    species_path = "Output/AncomBC2_p/001_disc/001_lesion_burden_onlygc/lesion_burden_ANCOMBC2_significant_results.csv",
    folder = "Output/BOI_AncomBC2/001/",
    name = "Bacteria_lesion_burden_001"
)
createBOI_ancomBC2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    species_path = "Output/AncomBC2_p/001_disc/001_gadolinium_contrast_onlygc/gadolinium_contrast_ANCOMBC2_significant_results.csv",
    folder = "Output/BOI_AncomBC2/001/",
    name = "Bacteria_gadolinium_contrast_001"
)


createBOI_ancomBC2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
    species_path = "Output/AncomBC2_p/05_disc/05_subtentorial_lesions_onlygc/subtentorial_lesions_ANCOMBC2_significant_results.csv",
    folder = "Output/BOI_AncomBC2/05/",
    name = "Bacteria_subtentorial_lesions_05"
)

createBOI_ancomBC2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
    species_path = "Output/AncomBC2_p/05_disc/05_spinal_cord_lesion_onlygc/spinal_cord_lesion_ANCOMBC2_significant_results.csv",
    folder = "Output/BOI_AncomBC2/05/",
    name = "Bacteria_spinal_cord_lesion_05"
)
createBOI_ancomBC2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
    species_path = "Output/AncomBC2_p/05_disc/05_lesion_burden_onlygc/lesion_burden_ANCOMBC2_significant_results.csv",
    folder = "Output/BOI_AncomBC2/05/",
    name = "Bacteria_lesion_burden_05"
)
createBOI_ancomBC2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
    species_path = "Output/AncomBC2_p/05_disc/05_gadolinium_contrast_onlygc/gadolinium_contrast_ANCOMBC2_significant_results.csv",
    folder = "Output/BOI_AncomBC2/05/",
    name = "Bacteria_gadolinium_contrast_05"
)

createBOI_ancomBC2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
    species_path = "Output/AncomBC2_p/01_disc/01_gadolinium_contrast_onlygc/gadolinium_contrast_ANCOMBC2_significant_results.csv",
    folder = "Output/BOI_AncomBC2/01/",
    name = "Bacteria_gadolinium_contrast_01"
)
createBOI_ancomBC2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
    species_path = "Output/AncomBC2_p/01_disc/01_lesion_burden_onlygc/lesion_burden_ANCOMBC2_significant_results.csv",
    folder = "Output/BOI_AncomBC2/01/",
    name = "Bacteria_lesion_burden_01"
)
createBOI_ancomBC2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
    species_path = "Output/AncomBC2_p/01_disc/01_spinal_cord_lesion_onlygc/spinal_cord_lesion_ANCOMBC2_significant_results.csv",
    folder = "Output/BOI_AncomBC2/01/",
    name = "Bacteria_spinal_cord_lesion_01"
)
createBOI_ancomBC2(
    rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
    species_path = "Output/AncomBC2_p/01_disc/01_subtentorial_lesions_onlygc/subtentorial_lesions_ANCOMBC2_significant_results.csv",
    folder = "Output/BOI_AncomBC2/01/",
    name = "Bacteria_subtentorial_lesions_01"
)   



#####for union

createBOI(
    "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    "Output/compare_DAS_MAASLIN3/gadolinium_contrast_001_union/combined_species_001.csv",
    "Output/BOI/001/",
    "Bacteria_gadolinium_001_union"
)
createBOI(
    "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    "Output/compare_DAS_MAASLIN3/lesion_burden_001_union/combined_species_001.csv",
    "Output/BOI/001/",
    "Bacteria_lesion_001_union"
)
createBOI(
    "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    "Output/compare_DAS_MAASLIN3/spinal_cord_lesion_001_union/combined_species_001.csv",
    "Output/BOI/001/",
    "Bacteria_spinal_cord_001_union"
)
createBOI(
    "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    "Output/compare_DAS_MAASLIN3/subtentorial_lesions_001_union/combined_species_001.csv",
    "Output/BOI/001/",
    "Bacteria_subtentorial_lesions_001_union"
)
