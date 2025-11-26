AncomBC2_run <- function(rds, output_folder, categorySel, elements, discriminant) {
    taxa <- as.data.frame(rds@tax_table)
    taxa$Genus_Species <- paste(taxa$Genus, taxa$Species, sep = " ")
    taxa$otu <- rownames(taxa)
    taxa_tab <- taxa %>% select(otu, Genus_Species)

    # Filtraggio campioni
    otu <- as.data.frame(rds@otu_table)
    metadata_df <- as(rds@sam_data, "data.frame")
    ids_to_retain <- rownames(metadata_df[metadata_df[[categorySel]] %in% elements, ])
    otu <- otu[, colnames(otu) %in% ids_to_retain]
    metadata_df <- metadata_df[rownames(metadata_df) %in% ids_to_retain, ]

    # Variabile discriminante
    metadata_df[[discriminant]] <- as.factor(metadata_df[[discriminant]])

    # Creazione output folder
    if (!dir.exists(output_folder)) {
        dir.create(output_folder, recursive = TRUE)
    }

    # ANCOM-BC2
    res <- ANCOMBC::ancombc2(
        data = otu,
        meta_data = metadata_df,
        fix_formula = discriminant,
        tax_level = NULL,
        p_adj_method = "BH",
        rand_formula = NULL,
        prv_cut = 0.0,
        lib_cut = 0,
        struc_zero = FALSE
    )

    results <- res$res
    discriminant1 <- paste0("q_", gsub(" ", "_", discriminant), "1")
    results <- res$res[, c("taxon", discriminant1)]
    results <- as.data.frame(results)
results <- results %>% dplyr::rename(otu = taxon)

    # Merge con tassonomia
    out2 <- merge(results, taxa_tab, by = "otu")
    out2 <- out2 %>% select(Genus_Species, all_of(discriminant1))
    print(min(out2[[discriminant1]], na.rm = TRUE))
    write.csv(out2, file.path(output_folder, paste0(discriminant, "_ANCOMBC2_results.csv")), row.names = FALSE)
    #saveRDS(out2, file.path(output_folder, "ANCOMBC2_results.rds"))
}
rds_001 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
AncomBC2_run(
    rds = rds_001,
    output_folder = "Output/AncomBC2/001/001_lesion_burden_onlygc",
    categorySel = "gc_treatment", # usato per filtrare
    elements = c("positive", "negative"),
    discriminant = "lesion_burden" # usato per ANCOM-BC2
)
AncomBC2_run(
    rds = rds_001,
    output_folder = "Output/AncomBC2/001/001_spinal_cord lesion_onlygc",
    categorySel = "gc_treatment", # usato per filtrare
    elements = c("positive", "negative"),
    discriminant = "spinal_cord_lesion" # usato per ANCOM-BC2
)
AncomBC2_run(
    rds = rds_001,
    output_folder = "Output/AncomBC2/001/001_gadolinium_contrast_onlygc",
    categorySel = "gc_treatment", # usato per filtrare
    elements = c("positive", "negative"),
    discriminant = "gadolinium_contrast" # usato per ANCOM-BC2
)
AncomBC2_run(
    rds = rds_001,
    output_folder = "Output/AncomBC2/001/001_subtentorial_lesions_onlygc",
    categorySel = "gc_treatment", # usato per filtrare
    elements = c("positive", "negative"),
    discriminant = "subtentorial_lesions" # usato per ANCOM-BC2
)
rds_05 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds")
AncomBC2_run(
    rds = rds_05,
    output_folder = "Output/AncomBC2/05/05_lesion_burden_onlygc",
    categorySel = "gc_treatment", # usato per filtrare
    elements = c("positive", "negative"),
    discriminant = "lesion_burden" # usato per ANCOM-BC2
)
AncomBC2_run(
    rds = rds_05,
    output_folder = "Output/AncomBC2/05/05_spinal_cord_lesion_onlygc",
    categorySel = "gc_treatment", # usato per filtrare
    elements = c("positive", "negative"),
    discriminant = "spinal_cord_lesion" # usato per ANCOM-BC2
)
AncomBC2_run(
    rds = rds_05,
    output_folder = "Output/AncomBC2/05/05_gadolinium_contrast_onlygc",
    categorySel = "gc_treatment", # usato per filtrare
    elements = c("positive", "negative"),
    discriminant = "gadolinium_contrast" # usato per ANCOM-BC2
)
AncomBC2_run(
    rds = rds_05,
    output_folder = "Output/AncomBC2/05/05_subtentorial_lesions_onlygc",
    categorySel = "gc_treatment", # usato per filtrare
    elements = c("positive", "negative"),
    discriminant = "subtentorial_lesions" # usato per ANCOM-BC2
)

rds_01 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds")
AncomBC2_run(
    rds = rds_01,
    output_folder = "Output/AncomBC2/01/01_lesion_burden_onlygc",
    categorySel = "gc_treatment", # usato per filtrare
    elements = c("positive", "negative"),
    discriminant = "lesion_burden" # usato per ANCOM-BC2
)
AncomBC2_run(
    rds = rds_01,
    output_folder = "Output/AncomBC2/01/01_spinal_cord_lesion_onlygc",
    categorySel = "gc_treatment", # usato per filtrare
    elements = c("positive", "negative"),
    discriminant = "spinal_cord_lesion" # usato per ANCOM-BC2
)
AncomBC2_run(
    rds = rds_01,
    output_folder = "Output/AncomBC2/01/01_gadolinium_contrast_onlygc",
    categorySel = "gc_treatment", # usato per filtrare
    elements = c("positive", "negative"),
    discriminant = "gadolinium_contrast" # usato per ANCOM-BC2
)
AncomBC2_run(
    rds = rds_01,
    output_folder = "Output/AncomBC2/01/01_subtentorial_lesions_onlygc",
    categorySel = "gc_treatment", # usato per filtrare
    elements = c("positive", "negative"),
    discriminant = "subtentorial_lesions" # usato per ANCOM-BC2
)
