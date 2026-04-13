Aldex2_run <- function(rds, output_folder, categorySel, elements, discriminant) {

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

    # Variabile discriminante (NON quella del filtraggio)
    conditions <- metadata_df[[discriminant]]

    # Creazione output
    if (!dir.exists(output_folder)) {
        dir.create(output_folder, recursive = TRUE)
    }

    # ALDEx2
    aldex_res <- ALDEx2::aldex(
        otu,
        conditions,
        mc.samples = 128,
        denom = "iqlr",
        test = "t",
        effect = TRUE,
        verbose = FALSE
    )

    # Merge con tassonomia
    aldex_res$otu <- rownames(aldex_res)
    tab <- merge(aldex_res, taxa_tab, by = "otu")
    tab <- tab %>% select(Genus_Species, wi.eBH)

    # Salvataggio
    write.csv(tab, file.path(output_folder, paste0(discriminant,"_aldex2_results.csv")), row.names = FALSE)
    # Salvataggio risultati significativi (p-value < 0.05)
        tab_sig <- tab %>% filter(wi.eBH < 0.05)
        if (nrow(tab_sig) > 0) {
            write.csv(tab_sig, file.path(output_folder, paste0(discriminant,"_aldex2_significant.csv")), row.names = FALSE)
        }
   # saveRDS(aldex_res, file.path(output_folder, "aldex2_results.rds"))
    print(min(aldex_res$wi.eBH, na.rm = TRUE))
   # return(aldex_res)
}

rds_001 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
Aldex2_run(
    rds = rds_001,
    output_folder = "Output/Aldex2/001/001_lesion_burden_onlygc",
    categorySel = "gc_treatment",                  # usato per filtrare
    elements = c("positive", "negative"),
    discriminant = "lesion_burden"                 # usato per ALDEx2
)
Aldex2_run(
    rds = rds_001,
    output_folder = "Output/Aldex2/001/001_spinal_cord_lesion_onlygc",
    categorySel = "gc_treatment",                  # usato per filtrare
    elements = c("positive", "negative"),
    discriminant = "spinal_cord_lesion"           # usato per ALDEx2
)
Aldex2_run(
    rds = rds_001,
    output_folder = "Output/Aldex2/001/001_gadolinium_contrast_onlygc",
    categorySel = "gc_treatment",                  # usato per filtrare
    elements = c("positive", "negative"),
    discriminant = "gadolinium_contrast"          # usato per ALDEx2
)
Aldex2_run(
    rds = rds_001,
    output_folder = "Output/Aldex2/001/001_subtentorial_lesions_onlygc",
    categorySel = "gc_treatment",                  # usato per filtrare
    elements = c("positive", "negative"),
    discriminant = "subtentorial_lesions"         # usato per ALDEx2
)
rds_05 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds")

Aldex2_run(
    rds = rds_05,
    output_folder = "Output/Aldex2/05/05_lesion_burden_onlygc",
    categorySel = "gc_treatment",                  # usato per filtrare
    elements = c("positive", "negative"),
    discriminant = "lesion_burden"                 # usato per ALDEx2
)
Aldex2_run(
    rds = rds_05,
    output_folder = "Output/Aldex2/05/05_spinal_cord_lesion_onlygc",
    categorySel = "gc_treatment",                  # usato per filtrare
    elements = c("positive", "negative"),  
    discriminant = "spinal_cord_lesion"           # usato per ALDEx2
)
Aldex2_run(
    rds = rds_05,
    output_folder = "Output/Aldex2/05/05_gadolinium_contrast_onlygc",
    categorySel = "gc_treatment",                  # usato per filtrare
    elements = c("positive", "negative"),
    discriminant = "gadolinium_contrast"          # usato per ALDEx2
)
Aldex2_run(
    rds = rds_05,
    output_folder = "Output/Aldex2/05/05_subtentorial_lesions_onlygc",
    categorySel = "gc_treatment",                  # usato per filtrare
    elements = c("positive", "negative"),
    discriminant = "subtentorial_lesions"         # usato per ALDEx2
)
rds_01 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds")

Aldex2_run(
    rds = rds_01,
    output_folder = "Output/Aldex2/01/01_lesion_burden_onlygc",
    categorySel = "gc_treatment",                  # usato per filtrare
    elements = c("positive", "negative"),
    discriminant = "lesion_burden"                 # usato per ALDEx2
)
Aldex2_run(
    rds = rds_01,
    output_folder = "Output/Aldex2/01/01_spinal_cord_lesion_onlygc",
    categorySel = "gc_treatment",                  # usato per filtrare
    elements = c("positive", "negative"),  
    discriminant = "spinal_cord_lesion"           # usato per ALDEx2
)
Aldex2_run(
    rds = rds_01,
    output_folder = "Output/Aldex2/01/01_gadolinium_contrast_onlygc",
    categorySel = "gc_treatment",                  # usato per filtrare 
    elements = c("positive", "negative"),
    discriminant = "gadolinium_contrast"          # usato per ALDEx2
)
Aldex2_run(
    rds = rds_01,
    output_folder = "Output/Aldex2/01/01_subtentorial_lesions_onlygc",
    categorySel = "gc_treatment",                  # usato per filtrare
    elements = c("positive", "negative"),     
    discriminant = "subtentorial_lesions"         # usato per ALDEx2
)
