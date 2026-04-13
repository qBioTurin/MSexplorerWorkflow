# # Read two TSV files
# limma <- read.delim("Output/LIMMA_score/05/Bacteria_gadolinium_contrastboth_05_limma_top_table.tsv", header = TRUE, sep = "\t")
# lefse <- read.delim("Output/LEFSE/05/unprocesssed_lefse/Bacteria_gadolinium_contrast_05_lefse.res", header = FALSE, sep = "\t")
# maaslin3 <- read.delim("Output/MAASLIN3/05/05_gadolinium_contrast_onlygc/all_results.tsv", header = TRUE, sep = "\t")
merge_p_value <- function(limma_path, lefse_path, maaslin3_path, analysis_name, aldex_path, ancomb_path, rds_path = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds") {
    rds <- readRDS(rds_path)
    taxa <- as.data.frame(phyloseq::tax_table(rds))
    taxa$genus_species <- paste(taxa$Genus, taxa$Species, sep = " ")
    genus_species_list <- taxa$genus_species
    add_original_names <- function(df, original_names) {
        # normalizzazione robusta
        normalize <- function(x) {
            x |>
                tolower() |>
                gsub("[^a-z0-9]", "", x = _) # elimina tutto tranne lettere e numeri
        }

        # normalizza i nomi originali
        norm_original <- normalize(original_names)

        # crea dizionario: nome_normalizzato → nome_originale
        key <- setNames(original_names, norm_original)

        # normalizza la colonna dei nomi nella tabella
        df$species <- key[normalize(df[[1]])]

        return(df)
    }

    limma <- read.delim(limma_path, header = TRUE, sep = "\t")
    lefse <- read.delim(lefse_path, header = FALSE, sep = "\t")
    maaslin3 <- read.delim(maaslin3_path, header = TRUE, sep = "\t")
    aldex <- read.csv(aldex_path, header = TRUE, sep = ",")
    ancombc2 <- read.csv(ancomb_path, header = TRUE, sep = ",")
    # Rimuovi righe con model = "prevalence" da maaslin3
    # maaslin3 <- maaslin3[maaslin3$model != "prevalence", ]

    # ---- LEFSE: selezione colonne, estrazione specie e normalizzazione + aggregazione ----
    lefse <- lefse[, c("V1", "V5")]
    pval_lefse_col <- paste0("P.val_lefse_", analysis_name)
    colnames(lefse) <- c("species_l", pval_lefse_col)

    # mantieni solo le righe con la struttura attesa
    lefse <- lefse[lengths(gregexpr("\\.", lefse$species_l)) == 6, ]

    # estrai l'ultimo livello con pattern esistente
    lefse$species_l <- sub(".*\\.([^.]+\\.[^.]+)$", "\\1", lefse$species_l)

    # normalizza il nome (fusion)
    lefse$species_l <- gsub("\\.", " ", lefse$species_l)
    lefse <- add_original_names(lefse, genus_species_list)
    lefse$fusion_raw <- lefse$species_l
    # controlla se è presente Ruminococcus sp_FMBCY1 e correggilo
    lefse$species[lefse$species_l == "Ruminococcus sp_FMBCY1"] <- "Ruminococcus sp. FMBCY1"
    lefse <- lefse[, c("species", pval_lefse_col)]
    # ---- LIMMA: selezione colonne e normalizzazione fusion (coerente con LEFSE) ----
    limma <- limma[, c("X", "P.Value")]
    colnames(limma) <- c("species", paste0("P.val_limma_", analysis_name))

    # ---- MAASLIN3: selezione, pulizia e normalizzazione fusion ----
    maaslin3 <- maaslin3[, c("feature", "pval_individual", "model")]
    colnames(maaslin3) <- c("species_m", paste0("P.val_maaslin_", analysis_name), "model_maaslin")

    # assicurati p-value numerico; rimuovi NA convertiti
    #qua bisogna inoltre 
    pval_maas_col <- paste0("P.val_maaslin_", analysis_name)
    maaslin3[[pval_maas_col]] <- suppressWarnings(as.numeric(as.character(maaslin3[[pval_maas_col]])))
    maaslin3 <- maaslin3[!is.na(maaslin3[[pval_maas_col]]), ]

    # rimuovi righe di 'prevalence' se duplicate (correzione: model_maaslin)
    if ("model_maaslin" %in% colnames(maaslin3)) {
        dup_species <- maaslin3$species_m[duplicated(maaslin3$species_m) | duplicated(maaslin3$species_m, fromLast = TRUE)]
        if (length(dup_species) > 0) {
            maaslin3 <- maaslin3[!(maaslin3$species_m %in% dup_species & maaslin3$model_maaslin == "prevalence"), ]
        }
    }
    # imposta p-value a 1 per righe con model_maaslin == "prevalence"
    if ("model_maaslin" %in% colnames(maaslin3)) {
        maaslin3[[pval_maas_col]][maaslin3$model_maaslin == "prevalence"] <- 1
    }

    # normalizza species_m sostituendo underscore con spazi
    maaslin3$species_m <- gsub("_", " ", maaslin3$species_m)
    colnames(maaslin3)[colnames(maaslin3) == "species_m"] <- "species"
    maaslin3 <- maaslin3[, c("species", pval_maas_col)]
    # ---- Merge: limma + lefse + maaslin3 tramite 'fusion' ----
    merged <- merge(limma, lefse, by = "species", all = TRUE)
    merged <- merge(merged, maaslin3, by = "species", all = TRUE)
    merged2 <- as.data.frame(merged)
    rownames(merged2) <- merged2$species

    pval_col_aldex <- paste0("P.val_aldex2_", analysis_name)
    pval_ancombc2 <- paste0("P.val_ancombc2_", analysis_name)
    aldex2 <- aldex %>%
        dplyr::rename(
            species = Genus_Species,
            !!pval_col_aldex := wi.eBH # nota l'uso di !! per il nome dinamico
        )
    ancombc2 <- ancombc2 %>%
        dplyr::rename(
            species = Genus_Species,
            !!pval_ancombc2 := all_of(paste0("q_", gsub(" ", "_", analysis_name), "1"))
        )   


    merged2 <- merge(merged2, aldex2, by = "species", all = TRUE)
    merged2 <- merge(merged2, ancombc2, by = "species", all = TRUE)
    rownames(merged2) <- merged2$species
    merged2 <- merged2[, -1]

    #Converti tutte le colonne in numerico, sostituendo NA e valori non validi con 1
    merged2[] <- lapply(merged2, function(col) {
        num_col <- suppressWarnings(as.numeric(as.character(col)))
        num_col[is.na(num_col)] <- NA
        return(num_col)
     })
    return(merged2)
}

# Gadolinium contrast analysis
merge_p_value(
    limma_path = "Output/LIMMA_score/05/Bacteria_gadolinium_contrastboth_05_limma_top_table.tsv",
    lefse_path = "Output/LEFSE/05/unprocesssed_lefse/Bacteria_gadolinium_contrast_05_lefse.res",
    maaslin3_path = "Output/MAASLIN3/05/05_gadolinium_contrast_onlygc/all_results.tsv",
    aldex_path = "Output/Aldex2/05/05_gadolinium_contrast_onlygc/gadolinium_contrast_aldex2_results.csv",
    ancomb_path = "Output/AncomBC2/05/05_gadolinium_contrast_onlygc/gadolinium_contrast_ANCOMBC2_results.csv",
    analysis_name = "gadolinium_contrast"
) -> merged_df_gadolinium_contrast_05

# Lesion burden analysis
merge_p_value(
    limma_path = "Output/LIMMA_score/05/Bacteria_lesion_burdenboth_05_limma_top_table.tsv",
    lefse_path = "Output/LEFSE/05/unprocesssed_lefse/Bacteria_lesion_burden_05_lefse.res",
    maaslin3_path = "Output/MAASLIN3/05/05_lesion_burden_onlygc/all_results.tsv",
    aldex_path = "Output/Aldex2/05/05_lesion_burden_onlygc/lesion_burden_aldex2_results.csv",
    ancomb_path = "Output/AncomBC2/05/05_lesion_burden_onlygc/lesion_burden_ANCOMBC2_results.csv",
    analysis_name = "lesion_burden"
) -> merged_df_lesion_burden_05

# Spinal cord analysis
merge_p_value(
    limma_path = "Output/LIMMA_score/05/Bacteria_spinal_cord_lesionboth_05_limma_top_table.tsv",
    lefse_path = "Output/LEFSE/05/unprocesssed_lefse/Bacteria_spinal_cord_lesion_05_lefse.res",
    maaslin3_path = "Output/MAASLIN3/05/05_spinal_cord_lesion_onlygc/all_results.tsv",
    aldex_path = "Output/Aldex2/05/05_spinal_cord_lesion_onlygc/spinal_cord_lesion_aldex2_results.csv",
    ancomb_path = "Output/AncomBC2/05/05_spinal_cord_lesion_onlygc/spinal_cord_lesion_ANCOMBC2_results.csv",
    analysis_name = "spinal_cord_lesion"
) -> merged_df_spinal_cord_05

# Subtentorial analysis
merge_p_value(
    limma_path = "Output/LIMMA_score/05/Bacteria_subtentorial_lesionsboth_05_limma_top_table.tsv",
    lefse_path = "Output/LEFSE/05/unprocesssed_lefse/Bacteria_subtentorial_lesions_05_lefse.res",
    maaslin3_path = "Output/MAASLIN3/05/05_subtentorial_lesions_onlygc/all_results.tsv",
    aldex_path = "Output/Aldex2/05/05_subtentorial_lesions_onlygc/subtentorial_lesions_aldex2_results.csv",
    ancomb_path = "Output/AncomBC2/05/05_subtentorial_lesions_onlygc/subtentorial_lesions_ANCOMBC2_results.csv",
    analysis_name = "subtentorial_lesions"
) -> merged_df_subtentorial_05

# Merge all four analysis results using rownames
all_merged <- Reduce(
    function(x, y) {
        merge(x, y, by = "row.names", all = TRUE) |>
            (\(df) {
                rownames(df) <- df$Row.names
                df[, -1, drop = FALSE]
            })()
    },
    list(
        merged_df_gadolinium_contrast_05,
        merged_df_lesion_burden_05,
        merged_df_spinal_cord_05,
        merged_df_subtentorial_05
    )
)

# Create directory if it doesn't exist
if (!dir.exists("Output/p_value")) {
    dir.create("Output/p_value", recursive = TRUE)
}
all_merged$species <- rownames(all_merged)
all_merged <- all_merged[, c("species", setdiff(names(all_merged), "species"))]
write.table(all_merged, "Output/p_value/05_bacteria.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# 01 threshold analysis
merge_p_value(
    limma_path = "Output/LIMMA_score/01/Bacteria_gadolinium_contrastboth_01_limma_top_table.tsv",
    lefse_path = "Output/LEFSE/01/unprocesssed_lefse/Bacteria_gadolinium_contrast_01_lefse.res",
    maaslin3_path = "Output/MAASLIN3/01/01_gadolinium_contrast_onlygc/all_results.tsv",
    aldex_path = "Output/Aldex2/01/01_gadolinium_contrast_onlygc/gadolinium_contrast_aldex2_results.csv",
    ancomb_path = "Output/AncomBC2/01/01_gadolinium_contrast_onlygc/gadolinium_contrast_ANCOMBC2_results.csv",
    analysis_name = "gadolinium_contrast"
) -> merged_df_gadolinium_contrast_01

merge_p_value(
    limma_path = "Output/LIMMA_score/01/Bacteria_lesion_burdenboth_01_limma_top_table.tsv",
    lefse_path = "Output/LEFSE/01/unprocesssed_lefse/Bacteria_lesion_burden_01_lefse.res",
    maaslin3_path = "Output/MAASLIN3/01/01_lesion_burden_onlygc/all_results.tsv",
    aldex_path = "Output/Aldex2/01/01_lesion_burden_onlygc/lesion_burden_aldex2_results.csv",
    ancomb_path = "Output/AncomBC2/01/01_lesion_burden_onlygc/lesion_burden_ANCOMBC2_results.csv",
    analysis_name = "lesion_burden"
) -> merged_df_lesion_burden_01

merge_p_value(
    limma_path = "Output/LIMMA_score/01/Bacteria_spinal_cord_lesionboth_01_limma_top_table.tsv",
    lefse_path = "Output/LEFSE/01/unprocesssed_lefse/Bacteria_spinal_cord_lesion_01_lefse.res",
    maaslin3_path = "Output/MAASLIN3/01/01_spinal_cord_lesion_onlygc/all_results.tsv",
    aldex_path = "Output/Aldex2/01/01_spinal_cord_lesion_onlygc/spinal_cord_lesion_aldex2_results.csv",
    ancomb_path = "Output/AncomBC2/01/01_spinal_cord_lesion_onlygc/spinal_cord_lesion_ANCOMBC2_results.csv",
    analysis_name = "spinal_cord_lesion"
) -> merged_df_spinal_cord_01

merge_p_value(
    limma_path = "Output/LIMMA_score/01/Bacteria_subtentorial_lesionsboth_01_limma_top_table.tsv",
    lefse_path = "Output/LEFSE/01/unprocesssed_lefse/Bacteria_subtentorial_lesions_01_lefse.res",
    maaslin3_path = "Output/MAASLIN3/01/01_subtentorial_lesions_onlygc/all_results.tsv",
    aldex_path = "Output/Aldex2/01/01_subtentorial_lesions_onlygc/subtentorial_lesions_aldex2_results.csv",
    ancomb_path = "Output/AncomBC2/01/01_subtentorial_lesions_onlygc/subtentorial_lesions_ANCOMBC2_results.csv",
    analysis_name = "subtentorial_lesions"
) -> merged_df_subtentorial_01

all_merged <- Reduce(
    function(x, y) {
        merge(x, y, by = "row.names", all = TRUE) |>
            (\(df) {
                rownames(df) <- df$Row.names
                df[, -1, drop = FALSE]
            })()
    },
    list(
        merged_df_gadolinium_contrast_01,
        merged_df_lesion_burden_01,
        merged_df_spinal_cord_01,
        merged_df_subtentorial_01
    )
)
all_merged$species <- rownames(all_merged)
all_merged <- all_merged[, c("species", setdiff(names(all_merged), "species"))]
write.table(all_merged, "Output/p_value/01_bacteria.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
# 001 threshold analysis
merge_p_value(
    limma_path = "Output/LIMMA_score/001/Bacteria_gadolinium_contrastboth_001_limma_top_table.tsv",
    lefse_path = "Output/LEFSE/001/unprocesssed_lefse/Bacteria_gadolinium_contrast_001_lefse.res",
    maaslin3_path = "Output/MAASLIN3/001/001_gadolinium_contrast_onlygc/all_results.tsv",
    aldex_path = "Output/Aldex2/001/001_gadolinium_contrast_onlygc/gadolinium_contrast_aldex2_results.csv",
    ancomb_path = "Output/AncomBC2/001/001_gadolinium_contrast_onlygc/gadolinium_contrast_ANCOMBC2_results.csv",
    analysis_name = "gadolinium_contrast"
) -> merged_df_gadolinium_contrast_001

merge_p_value(
    limma_path = "Output/LIMMA_score/001/Bacteria_lesion_burdenboth_001_limma_top_table.tsv",
    lefse_path = "Output/LEFSE/001/unprocesssed_lefse/Bacteria_lesion_burden_001_lefse.res",
    maaslin3_path = "Output/MAASLIN3/001/001_lesion_burden_onlygc/all_results.tsv",
    aldex_path = "Output/Aldex2/001/001_lesion_burden_onlygc/lesion_burden_aldex2_results.csv",
    ancomb_path = "Output/AncomBC2/001/001_lesion_burden_onlygc/lesion_burden_ANCOMBC2_results.csv",
    analysis_name = "lesion_burden"
) -> merged_df_lesion_burden_001

merge_p_value(
    limma_path = "Output/LIMMA_score/001/Bacteria_spinal_cord_lesionboth_001_limma_top_table.tsv",
    lefse_path = "Output/LEFSE/001/unprocesssed_lefse/Bacteria_spinal_cord_lesion_001_lefse.res",
    maaslin3_path = "Output/MAASLIN3/001/001_spinal_cord_lesion_onlygc/all_results.tsv",
    aldex_path = "Output/Aldex2/001/001_spinal_cord_lesion_onlygc/spinal_cord_lesion_aldex2_results.csv",
    ancomb_path = "Output/AncomBC2/001/001_spinal_cord lesion_onlygc/spinal_cord_lesion_ANCOMBC2_results.csv",
    analysis_name = "spinal_cord_lesion"
) -> merged_df_spinal_cord_001

merge_p_value(
    limma_path = "Output/LIMMA_score/001/Bacteria_subtentorial_lesionsboth_001_limma_top_table.tsv",
    lefse_path = "Output/LEFSE/001/unprocesssed_lefse/Bacteria_subtentorial_lesions_001_lefse.res",
    maaslin3_path = "Output/MAASLIN3/001/001_subtentorial_lesions_onlygc/all_results.tsv",
    aldex_path = "Output/Aldex2/001/001_subtentorial_lesions_onlygc/subtentorial_lesions_aldex2_results.csv",
    ancomb_path = "Output/AncomBC2/001/001_subtentorial_lesions_onlygc/subtentorial_lesions_ANCOMBC2_results.csv",
    analysis_name = "subtentorial_lesions"
) -> merged_df_subtentorial_001

all_merged <- Reduce(
    function(x, y) {
        merge(x, y, by = "row.names", all = TRUE) |>
            (\(df) {
                rownames(df) <- df$Row.names
                df[, -1, drop = FALSE]
            })()
    },
    list(
        merged_df_gadolinium_contrast_001,
        merged_df_lesion_burden_001,
        merged_df_spinal_cord_001,
        merged_df_subtentorial_001
    )
)
all_merged$species <- rownames(all_merged)
all_merged <- all_merged[, c("species", setdiff(names(all_merged), "species"))]
write.table(all_merged, "Output/p_value/001_bacteria.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

