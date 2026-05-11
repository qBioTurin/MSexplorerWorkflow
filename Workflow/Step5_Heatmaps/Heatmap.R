source("Settings/utilities.R")
output_folder <- "Image/Heatmap_Maaslin_intersection_spec/BOI_U_DAS_final/"
createFolder(output_folder)
#### Functions ############
# --- Inizio codice scompattato senza funzioni ---

# Bacteria <- readRDS(file = "Output/merge_DAS/MSHD/Bacteria_MsHd_merged.rds")
# Archaea <- readRDS(file = "Output/merge_DAS/MSHD/Archaea_MsHd_merged.rds")
# Eukaryota <- readRDS(file = "Output/merge_DAS/MSHD/Eukaryote_MsHd_merged.rds")
# filename <- "category"
######### GC_LESION############
create_heatmap <- function(Bacteria, Archaea, Eukaryota, filename, output_folder, order, topBar, status = NULL) {
    # Funzione interna per preparare i dati di un kingdom
    prep_data <- function(kingdom_data) {
        if (is.null(kingdom_data)) {
            return(list(NULL, NULL, NULL))
        }

        taxatables <- as.data.frame(tax_table(kingdom_data)) %>% rownames_to_column("otu_id")
        otutables <- data.frame(otu_table(kingdom_data)) %>% rownames_to_column("otu_id")
        combined_df <- full_join(taxatables, otutables, by = "otu_id")

        combined_df <- combined_df %>% mutate(Genus_species = paste(Genus, Species, sep = " "))
        metadata_kingdom <- combined_df %>% select(Kingdom, Phylum, Genus, Genus_species)
        rownames(metadata_kingdom) <- metadata_kingdom$Genus_species
        metadata_kingdom <- metadata_kingdom %>% select(-Genus_species)

        norm_data <- combined_df[, -c(2:8)]
        rownames(norm_data) <- norm_data$Genus_species
        norm_data <- norm_data[, -c(1, length(norm_data))]

        metadata <- read.csv("InputData/merged_data_cluster_name.csv", header = TRUE, sep = ",", na = c("", " ", "NA"), check.names = TRUE)
        metadata <- metadata %>%
            mutate(
                across(c(id), as.character),
                across(c(
                    sex, category, clinical_presentation, gc_treatment, subtentorial_lesions, spinal_cord_lesion,
                    gadolinium_contrast, sequencing_batch, WORSENING, EDSS_DIAGNOSI, EDSS_PROGRESSIONE, Event,
                    naive, previous_therapy, antibiotic_use, sample_type, lesion_burden
                ), as.factor),
                across(c(age, bmi, EventTime), as.numeric),
                across(c(sample_collection_date), as.Date)
            )
        samples <- colnames(norm_data)
        metadata <- metadata[metadata$id %in% samples, ]
        rownames(metadata) <- metadata$id

        metadata_hm <- metadata %>%
            filter(category %in% c("HEALTHY", "MS")) %>%
            select(
                Sex = sex, Status = category, gc_treatment, Lesion_Burden = lesion_burden,
                Spinal_cord_lesions = spinal_cord_lesion, gadolinium_contrast, subtentorial_lesions
            ) %>%
            mutate(
                Sex = factor(Sex, levels = c("M", "F")),
                Status = factor(Status, levels = c("HEALTHY", "MS"), labels = c("HD", "PwMS")),
                gc_treatment = factor(gc_treatment, levels = c("healthy", "positive", "negative")),
                Lesion_Burden = factor(Lesion_Burden, levels = c(0, 1), labels = c("Low", "High")),
                Spinal_cord_lesions = factor(Spinal_cord_lesions, levels = c(0, 1), labels = c("No", "Yes")),
                subtentorial_lesions = factor(subtentorial_lesions, levels = c(0, 1), labels = c("No", "Yes")),
                gadolinium_contrast = factor(gadolinium_contrast, levels = c(0, 1), labels = c("No", "Yes"))
            )

        norm_data_z <- t(scale(t(norm_data)))
        list(norm_data_z, metadata_hm, metadata_kingdom)
    }

    # Prepara i dati
    data_list <- list(
        Bacteria = prep_data(Bacteria),
        Archaea = prep_data(Archaea),
        Eukaryota = prep_data(Eukaryota)
    )

    # Palette colori
    paletteLength <- 35
    colors_list <- list(
        Bacteria = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(paletteLength + 1)),
        Archaea = colorRampPalette(RColorBrewer::brewer.pal(9, "BrBG"))(paletteLength + 1),
        Eukaryota = colorRampPalette(RColorBrewer::brewer.pal(9, "PuOr"))(paletteLength + 1)
    )

    ann_colors_base <- list(
        Sex = c("M" = "blue", "F" = "pink"),
        Status = c("HD" = "#bad7f2", "PwMS" = "#ef6f6c"),
        gc_treatment = c("positive" = "#4D4D4D", "negative" = "#D7D7D7"),
        Spinal_cord_lesions = c("No" = "#13913B", "Yes" = "brown"),
        Lesion_Burden = c("Low" = "orange", "High" = "purple"),
        subtentorial_lesions = c("No" = "#233E6CFF", "Yes" = "#FFEA46FF"),
        gadolinium_contrast = c("No" = "#0DC9A7", "Yes" = "#E21277")
    )

    heatmap_list <- list()
    colOrder <- NULL
    kingdom <- names(data_list)[1]
    for (kingdom in names(data_list)) {
        dat <- data_list[[kingdom]]
        if (is.null(dat[[1]])) next

        norm_data_z <- dat[[1]]
        metadata_hm <- dat[[2]]
        metadata_kingdom <- dat[[3]]
        norm_data_z[is.na(norm_data_z)] <- 0

        # Ordine colonne
        if (order == "ALL") {
            colOrder <- sort(colnames(norm_data_z))

            # Ordine alfanumerico corretto: prima HD, poi MS, numeri in ordine
            nums <- suppressWarnings(as.numeric(gsub("\\D", "", colOrder)))
            pref <- substr(colOrder, 1, 2)

            colOrder <- colOrder[order(pref, nums)]
        } else {
            if (!order %in% colnames(metadata_hm)) {
                stop(paste("Column", order, "not found in metadata_hm"))
            }

            colOrder <- metadata_hm %>%
                mutate(sample_id = rownames(metadata_hm)) %>%
                arrange(!!sym(order), sample_id) %>%
                pull(sample_id)
            colOrder <- colOrder[colOrder %in% colnames(norm_data_z)]
        }

        # Costruisci break e colori
        # Scala simmetrica centrata sullo zero
        max_abs <- max(abs(range(norm_data_z)))
        myBreaks <- seq(-max_abs, max_abs, length.out = paletteLength + 1)
        col_fun <- circlize::colorRamp2(myBreaks, colors_list[[kingdom]])

        # Annotazioni
        phil_name <- paste0("Phylum ", unique(metadata_kingdom$Kingdom))
        ann_col <- metadata_kingdom %>% select(Phylum)
        colnames(ann_col) <- phil_name

        ann_row <- metadata_hm %>% select(all_of(topBar))
        ann_row <- ann_row[colOrder, , drop = FALSE]


        ann_colors <- ann_colors_base
        ann_colors[[phil_name]] <- setNames(viridisLite::mako(length(unique(metadata_kingdom$Phylum))), levels(as.factor(metadata_kingdom$Phylum)))
        ann_colors <- ann_colors[names(ann_colors) %in% c(names(ann_row), phil_name)]

        # Heatmap
        mat_t <- t(norm_data_z[, colOrder, drop = FALSE])
        cluster_cols_flag <- ncol(mat_t) > 1

        heatmap_list[[kingdom]] <- pheatmap(
            mat_t,
            annotation_col = ann_col,
            annotation_row = ann_row,
            annotation_colors = ann_colors,
            cluster_cols = cluster_cols_flag,
            cluster_rows = FALSE,
            legend = FALSE,
            annotation_legend = FALSE,
            fontsize = 14,
            fontsize_row = 18,
            fontsize_col = 18,
            angle_col = 90,
            col = colors_list[[kingdom]],
            border_color = NA,
            breaks = myBreaks
        )
        ph <- pheatmap(
            mat_t,
            annotation_col = ann_col,
            annotation_row = ann_row,
            annotation_colors = ann_colors,
            cluster_cols = cluster_cols_flag,
            cluster_rows = FALSE,
            legend = TRUE, # necessario
            annotation_legend = TRUE, # necessario
            fontsize = 14,
            fontsize_row = 18,
            fontsize_col = 18,
            angle_col = 45,
            col = colors_list[[kingdom]],
            border_color = NA,
            breaks = myBreaks
        )

        gt <- ph$gtable
        names_gt <- sapply(gt$grobs, function(x) x$name)
        # Controlla quali grobs ci sono
        sapply(gt$grobs, function(x) x$name)
        gtrees <- gt$grobs[sapply(gt$grobs, function(x) inherits(x, "gTree"))]
        lapply(gtrees, function(x) x$children)
        legend_grob_heat <- gtrees[[3]] # legenda heatmap
        legend_grob_legend <- gtrees[[2]] # legenda annotazioni

        # Combina tutto: gtable heatmap + due legende sotto
        legends_row <- arrangeGrob(
            legend_grob_heat,
            legend_grob_legend,
            nrow = 1,
            widths = c(1, 2) # larghezza relativa dei due grob
        )
        # Combina heatmap sopra e legende sotto
        final_merge <- arrangeGrob(
            heatmap_list[[kingdom]]$gtable,
            legends_row,
            ncol = 1,
            heights = c(8, 2) # altezza relativa della heatmap e delle legende
        )

        ggsave(
            plot = final_merge,
            filename = paste0(output_folder, filename, "_", kingdom, "_legends.pdf"),
            width = ceiling(nrow(ann_col) * 0.2 + 5),
            height = 40, limitsize = FALSE
        )
        saveRDS(heatmap_list[[kingdom]], file = paste0(output_folder, filename, "_", kingdom, "_heatmap.rds"))
    }

    # Combina heatmap in grid
    # labels <- names(heatmap_list)
    # text_grobs <- lapply(labels, function(lbl) textGrob(lbl, rot = 90, gp = gpar(fontsize = 14)))

    # g <- grid.arrange(
    #     grobs = c(rbind(text_grobs, lapply(heatmap_list, function(hm) hm$gtable))),
    #     ncol = 2,
    #     heights = rep(1, length(labels)),
    #     widths = c(0.05, 1.5)
    # )

    # # Salva tabelle tassonomiche
    # for (kingdom in names(data_list)) {
    #     if (!is.null(get(kingdom))) {
    #         write.table(as.data.frame(tax_table(get(kingdom))),
    #             file = paste0(output_folder, filename, "_", kingdom, "_taxa_table.tsv"),
    #             sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA
    #         )
    #     }
    # }

    # # Salva PDF
    # ggsave(
    #     plot = g, filename = paste0(output_folder, filename, ".pdf"),
    #     width = 40, height = 80, limitsize = FALSE
    # )
}



# ######################
# ##########ALL###########
# hello <- create_heatmap(
#     Bacteria = readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds"),
#     Archaea = readRDS(file = "Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.001.rds"),
#     Eukaryota = readRDS(file = "Output/SUPERVISED_DEC/Eukaryote_Supervised_decontam0.001.rds"),
#     filename = "ALL", output_folder = output_folder, order = "ALL", topBar = "Status"
# )
# ########################
# ######### MSHD ############
# create_heatmap(
#     Bacteria = readRDS(file = "Output/merge_DAS/MSHD/Bacteria_MsHd_merged.rds"),
#     Archaea = readRDS(file = "Output/merge_DAS/MSHD/Archaea_MsHd_merged.rds"),
#     Eukaryota = readRDS(file = "Output/merge_DAS/MSHD/Eukaryote_MsHd_merged.rds"),
#     filename = "category", output_folder = output_folder, order = "ALL", topBar = "Status"
# )
# # #########GC_LESION############
# create_heatmap(
#     Bacteria = readRDS("Output/merge_DAS/GC/Bacteria_GC_merged.rds"),
#     Archaea = readRDS("Output/merge_DAS/GC/Archaea_GC_merged.rds"),
#     Eukaryota = readRDS("Output/merge_DAS/GC/Eukaryote_GC_merged.rds"),
#     filename = "GC_LESION", output_folder = output_folder, order = "gc_treatment", topBar = "gc_treatment"
# )

# # #########lesion_burden############
# create_heatmap(
#     Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_lesion_burden_positive_merged.rds")),
#     Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_lesion_burden_positive_merged.rds")),
#     Eukaryota = readRDS(paste0("Output/merge_DAS/GC_comp/Eukaryote_lesion_burden_positive_merged.rds")),
#     filename = paste0("lesion_burden_positive"), output_folder = output_folder, order = "Spinal_cord_lesions", topBar = "Spinal_cord_lesions", status = "positive"
# )

# create_heatmap(
#     Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_lesion_burden_negative_merged.rds")),
#     Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_lesion_burden_negative_merged.rds")),
#     Eukaryota = readRDS(paste0("Output/merge_DAS/GC_comp/Eukaryote_lesion_burden_negative_merged.rds")),
#     filename = paste0("lesion_burden_negative"), output_folder = output_folder, order = "Spinal_cord_lesions", topBar = "Spinal_cord_lesions", status = "negative"
# )


# # #########spinal_cord_lesion############
# create_heatmap(
#     Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_spinal_cord_lesion_positive_merged.rds")),
#     Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_spinal_cord_lesion_positive_merged.rds")),
#     Eukaryota = NULL,
#     filename = paste0("spinal_cord_lesion_positive"), output_folder = output_folder, order = "Spinal_cord_lesions", topBar = "Spinal_cord_lesions", status = "positive"
# )


# create_heatmap(
#     Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_spinal_cord_lesion_negative_merged.rds")),
#     Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_spinal_cord_lesion_negative_merged.rds")),
#     Eukaryota = readRDS(paste0("Output/merge_DAS/GC_comp/Eukaryote_spinal_cord_lesion_negative_merged.rds")),
#     filename = paste0("spinal_cord_lesion_negative"), output_folder = output_folder, order = "Spinal_cord_lesions", topBar = "Spinal_cord_lesions", status = "negative"
# )

# ######### subtentorial_lesions############
# create_heatmap(
#     Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_subtentorial_lesions_positive_merged.rds")),
#     Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_subtentorial_lesions_positive_merged.rds")),
#     Eukaryota = readRDS(paste0("Output/merge_DAS/GC_comp/Eukaryote_subtentorial_lesions_positive_merged.rds")),
#     filename = paste0("subtentorial_lesions_positive"), output_folder = output_folder, order = "subtentorial_lesions", topBar = "subtentorial_lesions", status = "positive"
# )

# create_heatmap(
#     Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_subtentorial_lesions_negative_merged.rds")),
#     Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_subtentorial_lesions_negative_merged.rds")),
#     Eukaryota = NULL,
#     filename = paste0("subtentorial_lesions_negative"), output_folder = output_folder, order = "subtentorial_lesions", topBar = "subtentorial_lesions", status = "negative"
# )

# # #########gadolinium_contrast############
# create_heatmap(
#     Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_gadolinium_contrast_positive_merged.rds")),
#     Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_gadolinium_contrast_positive_merged.rds")),
#     Eukaryota = readRDS(paste0("Output/merge_DAS/GC_comp/Eukaryote_gadolinium_contrast_positive_merged.rds")),
#     filename = paste0("gadolinium_contrast_positive"), output_folder = output_folder, order = "gadolinium_contrast", topBar = "gadolinium_contrast", status = "positive"
# )
# create_heatmap(
#     Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_gadolinium_contrast_negative_merged.rds")),
#     Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_gadolinium_contrast_negative_merged.rds")),
#     Eukaryota = readRDS(paste0("Output/merge_DAS/GC_comp/Eukaryote_gadolinium_contrast_negative_merged.rds")),
#     filename = paste0("gadolinium_contrast_negative"), output_folder = output_folder, order = "gadolinium_contrast", topBar = "gadolinium_contrast", status = "negative"
# )
# dim <- c("01", "05", "001")
# for (d in dim) {
#     create_heatmap(
#         Bacteria = readRDS(paste0("Output/merge_DAS/", d, "/Bacteria_gadolinium_contrast_", d, "_merged.rds")),
#         Archaea = NULL,
#         Eukaryota = NULL,
#         filename = paste0("gadolinium_contrast_", d), output_folder = output_folder,
#         order = "gadolinium_contrast", topBar = "gadolinium_contrast"
#     )
#     create_heatmap(
#         Bacteria = readRDS(paste0("Output/merge_DAS/", d, "/Bacteria_spinal_cord_lesion_", d, "_merged.rds")),
#         Archaea = NULL,
#         Eukaryota = NULL,
#         filename = paste0("spinal_cord_lesion_", d), output_folder = output_folder,
#         order = "Spinal_cord_lesions", topBar = "Spinal_cord_lesions"
#     )
#     create_heatmap(
#         Bacteria = readRDS(paste0("Output/merge_DAS/", d, "/Bacteria_subtentorial_lesions_", d, "_merged.rds")),
#         Archaea = NULL,
#         Eukaryota = NULL,
#         filename = paste0("subtentorial_lesions_", d), output_folder = output_folder,
#         order = "subtentorial_lesions", topBar = "subtentorial_lesions"
#     )

#     create_heatmap(
#         Bacteria = readRDS(paste0("Output/merge_DAS/", d, "/Bacteria_lesion_burden_", d, "_merged.rds")),
#         Archaea = NULL,
#         Eukaryota = NULL,
#         filename = paste0("lesion_burden_", d), output_folder = output_folder,
#         order = "Lesion_Burden", topBar = "Lesion_Burden"
#     )
# }

output_folder <- "Image/Rebuttal2/"
create_heatmap(
    Bacteria = NULL,
    Archaea = readRDS("Output/MAASLIN3/Archaea_MSHD_disc/category_MAASLIN3_SR_features_K_T.rds"),
    Eukaryota = NULL,
    filename = paste0("MSHD"), output_folder = output_folder,
    order = "Status", topBar = "Status"
)

create_heatmap(
    Bacteria = NULL,
    Archaea = NULL,
    Eukaryota = readRDS("Output/MAASLIN3/Eukaryote_MSHD_disc/category_MAASLIN3_SR_features_K_T.rds"),
    filename = paste0("MSHD"), output_folder = output_folder,
    order = "Status", topBar = "Status"
)

create_heatmap(
    Bacteria = readRDS("Output/MAASLIN3/Bacteria_MSHD_disc_001/category_MAASLIN3_SR_features_K_T.rds"),
    Archaea = NULL,
    Eukaryota = NULL,
    filename = paste0("MSHD"), output_folder = output_folder,
    order = "Status", topBar = "Status"
)

create_heatmap(
    Bacteria = readRDS("Output/MAASLIN3/Bacteria_gc_treatment_disc/gc_treatment_MAASLIN3_SR_features_K_T.rds"),
    Archaea = NULL,
    Eukaryota = NULL,
    filename = paste0("GC"), output_folder = output_folder,
    order = "gc_treatment", topBar = "gc_treatment"
)


create_heatmap(
    Bacteria = readRDS("Output/MAASLIN3/001_disc/001_gadolinium_contrast_onlygc_plus_gc/gadolinium_contrast_MAASLIN3_SR_features_001_T.rds"),
    Archaea = NULL,
    Eukaryota = NULL,
    filename = paste0("gadolinium_"), output_folder = output_folder,
    order = "gadolinium_contrast", topBar = "gadolinium_contrast"
)


create_heatmap(
    Bacteria = readRDS("Output/MAASLIN3/001_disc/001_spinal_cord_lesion_onlygc_plus_gc/spinal_cord_lesion_MAASLIN3_SR_features_001_T.rds"),
    Archaea = NULL,
    Eukaryota = NULL,
    filename = paste0("spinal_cord_"), output_folder = output_folder,
    order = "Spinal_cord_lesions", topBar = "Spinal_cord_lesions"
)

create_heatmap(
    Bacteria = readRDS("Output/MAASLIN3/001_disc/001_subtentorial_lesions_onlygc_plus_gc/subtentorial_lesions_MAASLIN3_SR_features_001_T.rds"),
    Archaea = NULL,
    Eukaryota = NULL,
    filename = paste0("subtentorial_"), output_folder = output_folder,
    order = "subtentorial_lesions", topBar = "subtentorial_lesions"
)
create_heatmap(
    Bacteria = readRDS("Output/MAASLIN3/001_disc/001_lesion_burden_onlygc_plus_gc/lesion_burden_MAASLIN3_SR_features_001_T.rds"),
    Archaea = NULL,
    Eukaryota = NULL,
    filename = paste0("lesion_burden_"), output_folder = output_folder,
    order = "Lesion_Burden", topBar = "Lesion_Burden"
)
