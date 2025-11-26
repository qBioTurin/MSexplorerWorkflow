source("Settings/utilities.R")
output_folder <- "Image/Heatmap_new_auto/"
createFolder(output_folder)
#### Functions ############
# --- Inizio codice scompattato senza funzioni ---

# Bacteria <- readRDS(file = "Output/merge_DAS/MSHD/Bacteria_MsHd_merged.rds")
# Archaea <- readRDS(file = "Output/merge_DAS/MSHD/Archaea_MsHd_merged.rds")
# Eukaryota <- readRDS(file = "Output/merge_DAS/MSHD/Eukaryote_MsHd_merged.rds")
# filename <- "category"
######### GC_LESION############
create_heatmap <- function(Bacteria, Archaea, Eukaryota, filename, output_folder, order, topBar, status=NULL) {
    if (!is.null(Bacteria)) {
        baselines_dec <- Bacteria
        if (is.null(baselines_dec)) {
            data_bact <- list(NULL, NULL, NULL)
        } else {
            taxatables <- as.data.frame(tax_table(baselines_dec))
            taxatables <- rownames_to_column(taxatables, var = "otu_id")
            otutables <- data.frame(otu_table(baselines_dec))
            otutables <- rownames_to_column(otutables, var = "otu_id")
            combined_df <- full_join(taxatables, otutables, by = "otu_id")
            norm_data <- combined_df %>%
                mutate(Genus_species = paste(Genus, Species, sep = " "))
            metadata_kingdom <- norm_data %>%
                select(Kingdom, Phylum, Genus, Genus_species)
            rownames(metadata_kingdom) <- metadata_kingdom$Genus_species
            metadata_kingdom <- metadata_kingdom %>% select(-Genus_species)
            norm_data <- norm_data[, -c(2:8)]
            rownames(norm_data) <- norm_data$Genus_species
            norm_data <- norm_data[, -c(1, length(norm_data))]
            metadata <- read.csv("InputData/NewMetadataMS_fin.csv",
                header = TRUE,
                sep = ",",
                na = c("", " ", "NA"),
                check.names = TRUE
            )
            metadata <- metadata %>%
                mutate(
                    across(.cols = c(id), .fns = as.character),
                    across(
                        .cols = c(
                            sex, category, clinical_presentation, gc_treatment,
                            subtentorial_lesions, spinal_cord_lesion, gadolinium_contrast,
                            sequencing_batch, WORSENING, EDSS_DIAGNOSI, EDSS_PROGRESSIONE, Event,
                            naive, previous_therapy, antibiotic_use, sample_type, lesion_burden
                        ),
                        .fns = as.factor
                    ),
                    across(.cols = c(age, bmi, EventTime), .fns = as.numeric),
                    across(.cols = c(sample_collection_date), .fns = as.Date)
                )
            samples <- c(colnames(norm_data))
            metadata <- metadata[metadata$id %in% samples, ]
            rownames(metadata) <- metadata$id
            metadata_hm <- metadata %>%
                filter(category %in% c("HEALTHY", "MS")) %>%
                select(
                    Sex = sex,
                    Status = category,
                    gc_treatment = gc_treatment,
                    Lesion_Burden = lesion_burden,
                    Spinal_cord_lesions = spinal_cord_lesion,
                    gadolinium_contrast = gadolinium_contrast,
                    subtentorial_lesions = subtentorial_lesions
                ) %>%
                mutate(
                    Sex = factor(Sex, levels = c("M", "F"), labels = c("M", "F")),
                    Status = factor(Status, levels = c("HEALTHY", "MS"), labels = c("HD", "PwMS")),
                    gc_treatment = factor(gc_treatment, levels = c("healthy", "positive", "negative"), labels = c("healthy", "positive", "negative")),
                    Lesion_Burden = factor(Lesion_Burden, levels = c(0, 1), labels = c("Low", "High")),
                    Spinal_cord_lesions = factor(Spinal_cord_lesions, levels = c(0, 1), labels = c("No", "Yes")),
                    subtentorial_lesions = factor(subtentorial_lesions, levels = c(0, 1), labels = c("No", "Yes")),
                    gadolinium_contrast = factor(gadolinium_contrast, levels = c(0, 1), labels = c("No", "Yes"))
                )
            norm_data_z <- t(scale(t(norm_data)))
            data_bact <- list(norm_data_z, metadata_hm, metadata_kingdom)
        }
    }

    if (!is.null(Archaea)) {
        baselines_dec <- Archaea
        if (is.null(baselines_dec)) {
            data_arch <- list(NULL, NULL, NULL)
        } else {
            taxatables <- as.data.frame(tax_table(baselines_dec))
            taxatables <- rownames_to_column(taxatables, var = "otu_id")
            otutables <- data.frame(otu_table(baselines_dec))
            otutables <- rownames_to_column(otutables, var = "otu_id")
            combined_df <- full_join(taxatables, otutables, by = "otu_id")
            norm_data <- combined_df %>%
                mutate(Genus_species = paste(Genus, Species, sep = " "))
            metadata_kingdom <- norm_data %>%
                select(Kingdom, Phylum, Genus, Genus_species)
            rownames(metadata_kingdom) <- metadata_kingdom$Genus_species
            metadata_kingdom <- metadata_kingdom %>% select(-Genus_species)
            norm_data <- norm_data[, -c(2:8)]
            rownames(norm_data) <- norm_data$Genus_species
            norm_data <- norm_data[, -c(1, length(norm_data))]
            metadata <- read.csv("InputData/NewMetadataMS_fin.csv",
                header = TRUE,
                sep = ",",
                na = c("", " ", "NA"),
                check.names = TRUE
            )
            metadata <- metadata %>%
                mutate(
                    across(.cols = c(id), .fns = as.character),
                    across(
                        .cols = c(
                            sex, category, clinical_presentation, gc_treatment,
                            subtentorial_lesions, spinal_cord_lesion, gadolinium_contrast,
                            sequencing_batch, WORSENING, EDSS_DIAGNOSI, EDSS_PROGRESSIONE, Event,
                            naive, previous_therapy, antibiotic_use, sample_type, lesion_burden
                        ),
                        .fns = as.factor
                    ),
                    across(.cols = c(age, bmi, EventTime), .fns = as.numeric),
                    across(.cols = c(sample_collection_date), .fns = as.Date)
                )
            samples <- c(colnames(norm_data))
            metadata <- metadata[metadata$id %in% samples, ]
            rownames(metadata) <- metadata$id
            metadata_hm <- metadata %>%
                filter(category %in% c("HEALTHY", "MS")) %>%
                select(
                    Sex = sex,
                    Status = category,
                    gc_treatment = gc_treatment,
                    Lesion_Burden = lesion_burden,
                    Spinal_cord_lesions = spinal_cord_lesion,
                    gadolinium_contrast = gadolinium_contrast,
                    subtentorial_lesions = subtentorial_lesions
                ) %>%
                mutate(
                    Sex = factor(Sex, levels = c("M", "F"), labels = c("M", "F")),
                    Status = factor(Status, levels = c("HEALTHY", "MS"), labels = c("HD", "PwMS")),
                    gc_treatment = factor(gc_treatment, levels = c("healthy", "positive", "negative"), labels = c("healthy", "positive", "negative")),
                    Lesion_Burden = factor(Lesion_Burden, levels = c(0, 1), labels = c("Low", "High")),
                    Spinal_cord_lesions = factor(Spinal_cord_lesions, levels = c(0, 1), labels = c("No", "Yes")),
                    subtentorial_lesions = factor(subtentorial_lesions, levels = c(0, 1), labels = c("No", "Yes")),
                    gadolinium_contrast = factor(gadolinium_contrast, levels = c(0, 1), labels = c("No", "Yes"))
                )
            norm_data_z <- t(scale(t(norm_data)))
            data_arch <- list(norm_data_z, metadata_hm, metadata_kingdom)
        }
    }

    if (!is.null(Eukaryota)) {
        baselines_dec <- Eukaryota
        if (is.null(baselines_dec)) {
            data_euk <- list(NULL, NULL, NULL)
        } else {
            taxatables <- as.data.frame(tax_table(baselines_dec))
            taxatables <- rownames_to_column(taxatables, var = "otu_id")
            otutables <- data.frame(otu_table(baselines_dec))
            otutables <- rownames_to_column(otutables, var = "otu_id")
            combined_df <- full_join(taxatables, otutables, by = "otu_id")
            norm_data <- combined_df %>%
                mutate(Genus_species = paste(Genus, Species, sep = " "))
            metadata_kingdom <- norm_data %>%
                select(Kingdom, Phylum, Genus, Genus_species)
            rownames(metadata_kingdom) <- metadata_kingdom$Genus_species
            metadata_kingdom <- metadata_kingdom %>% select(-Genus_species)
            norm_data <- norm_data[, -c(2:8)]
            rownames(norm_data) <- norm_data$Genus_species
            norm_data <- norm_data[, -c(1, length(norm_data))]
            metadata <- read.csv("InputData/NewMetadataMS_fin.csv",
                header = TRUE,
                sep = ",",
                na = c("", " ", "NA"),
                check.names = TRUE
            )
            metadata <- metadata %>%
                mutate(
                    across(.cols = c(id), .fns = as.character),
                    across(
                        .cols = c(
                            sex, category, clinical_presentation, gc_treatment,
                            subtentorial_lesions, spinal_cord_lesion, gadolinium_contrast,
                            sequencing_batch, WORSENING, EDSS_DIAGNOSI, EDSS_PROGRESSIONE, Event,
                            naive, previous_therapy, antibiotic_use, sample_type, lesion_burden
                        ),
                        .fns = as.factor
                    ),
                    across(.cols = c(age, bmi, EventTime), .fns = as.numeric),
                    across(.cols = c(sample_collection_date), .fns = as.Date)
                )
            samples <- c(colnames(norm_data))
            metadata <- metadata[metadata$id %in% samples, ]
            rownames(metadata) <- metadata$id
            metadata_hm <- metadata %>%
                filter(category %in% c("HEALTHY", "MS")) %>%
                select(
                    Sex = sex,
                    Status = category,
                    gc_treatment = gc_treatment,
                    Lesion_Burden = lesion_burden,
                    Spinal_cord_lesions = spinal_cord_lesion,
                    gadolinium_contrast = gadolinium_contrast,
                    subtentorial_lesions = subtentorial_lesions
                ) %>%
                mutate(
                    Sex = factor(Sex, levels = c("M", "F"), labels = c("M", "F")),
                    Status = factor(Status, levels = c("HEALTHY", "MS"), labels = c("HD", "PwMS")),
                    gc_treatment = factor(gc_treatment, levels = c("healthy", "positive", "negative"), labels = c("healthy", "positive", "negative")),
                    Lesion_Burden = factor(Lesion_Burden, levels = c(0, 1), labels = c("Low", "High")),
                    Spinal_cord_lesions = factor(Spinal_cord_lesions, levels = c(0, 1), labels = c("No", "Yes")),
                    subtentorial_lesions = factor(subtentorial_lesions, levels = c(0, 1), labels = c("No", "Yes")),
                    gadolinium_contrast = factor(gadolinium_contrast, levels = c(0, 1), labels = c("No", "Yes"))
                )
            norm_data_z <- t(scale(t(norm_data)))
            data_euk <- list(norm_data_z, metadata_hm, metadata_kingdom)
        }
    }

    # --- Costruzione heatmap per ogni kingdom ---

    paletteLength <- 35
    myColors1 <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(paletteLength + 1))
    myColors2 <- colorRampPalette(RColorBrewer::brewer.pal(9, "PuOr"))(paletteLength + 1)
    myColors3 <- colorRampPalette(RColorBrewer::brewer.pal(9, "BrBG"))(paletteLength + 1)

    ann_colors <- list(
        Sex = setNames(c("blue", "pink"), levels(as.factor(data_bact[[2]]$Sex))),
        Status = setNames(c("#bad7f2", "#ef6f6c"), levels(as.factor(data_bact[[2]]$Status))),
        gc_treatment = setNames(c("#4D4D4D", "#D7D7D7"), c("positive", "negative")),
        Spinal_cord_lesions = setNames(c("brown", "#13913B"), levels(as.factor(data_bact[[2]]$Spinal_cord_lesions))),
        Lesion_Burden = setNames(c("orange", "purple"), levels(as.factor(data_bact[[2]]$Lesion_Burden))),
        subtentorial_lesions = setNames(c("#233E6CFF", "#FFEA46FF"), levels(as.factor(data_bact[[2]]$subtentorial_lesions))),
        gadolinium_contrast = setNames(c("#E21277", "#0DC9A7"), levels(as.factor(data_bact[[2]]$gadolinium_contrast)))
    )

    heatmap_list <- list()
    colOrder_bact <- NULL

    if (!is.null(data_bact[[1]])) {
        if (order == "ALL") {
            colOrder_bact <- sort(colnames(data_bact[[1]])) ############ qua devo modificare per il tipo di ordinamento
        } else {
            meta_tmp <- data_bact[[2]]
            meta_tmp$sample_id <- rownames(meta_tmp)
            colOrder_bact <- meta_tmp %>%
                arrange(!!sym(order)) %>%
                pull(sample_id)
            # Tieni solo quelli presenti nella matrice
            colOrder_bact <- colOrder_bact[colOrder_bact %in% colnames(data_bact[[1]])]

        }
        if (!missing(status) && !is.null(status)) {
           metadata_hm <- metadata_hm %>% filter(gc_treatment == status)
            norm_data_z <- norm_data_z[, rownames(metadata_hm), drop = FALSE]
        }
        nrow(metadata_hm)
        norm_data_z <- data_bact[[1]]
        metadata_hm <- data_bact[[2]]
        metadata_kingdom <- data_bact[[3]]
        norm_data_z[is.na(norm_data_z) | is.nan(norm_data_z)] <- 0
        myBreaks <- c(
            seq(min(norm_data_z), 0, length.out = ceiling(paletteLength / 2) + 1),
            seq(-1 / paletteLength, max(norm_data_z), length.out = floor(paletteLength / 2))
        )
        col_fun <- circlize::colorRamp2(myBreaks, rev(myColors1))
        ann_colors$Phylum <- setNames(
            viridisLite::mako(n = length(unique(metadata_kingdom$Phylum))),
            levels(as.factor(metadata_kingdom$Phylum))
        )
        metadata_hm <- metadata_hm %>% select(
            subtentorial_lesions, gadolinium_contrast,
            Spinal_cord_lesions, Lesion_Burden,
            gc_treatment, Sex, Status
        )
        metadata_hm$gc_treatment[metadata_hm$gc_treatment == "healthy"] <- NA
        if (length(unique(metadata_hm$Status)) == 1) {
            metadata_hm <- metadata_hm %>% select(-Status)
        }
        if (length(unique(metadata_hm$Lesion_Burden)) == 1) {
            metadata_hm <- metadata_hm %>% select(-Lesion_Burden)
        }
        if (length(unique(metadata_hm$gc_treatment)) == 1) {
            metadata_hm <- metadata_hm %>% select(-gc_treatment)
        }
        phil_name <- paste0("Phylum ", unique(metadata_kingdom$Kingdom))
        names(ann_colors)[names(ann_colors) == "Phylum"] <- phil_name
        metadata_kingdom <- metadata_kingdom %>% select(Phylum)
        colnames(metadata_kingdom) <- phil_name
        metadata_hm_nw <- metadata_hm %>%
            select(all_of(topBar)) ########### per poterne mettere solo una
        heatmap_bact <- pheatmap::pheatmap(norm_data_z[, colOrder_bact],
            annotation_col = metadata_hm_nw,
            annotation_colors = ann_colors,
            annotation_row = metadata_kingdom,
            cluster_cols = FALSE,
            cluster_rows = TRUE,
            fontsize = 14,          # testo generale
            fontsize.legend = 18,
            fontsize_row = 18,
            fontsize_col = 12,
            col = myColors1,
            border_color = NA,
            breaks = myBreaks,
            # silent = T
        )
        heatmap_list$Bacteria <- heatmap_bact
    }

    if (exists("data_arch") && !is.null(data_arch[[1]])) {
        norm_data_z <- data_arch[[1]]
        metadata_hm <- data_arch[[2]]
        metadata_kingdom <- data_arch[[3]]
        norm_data_z[is.na(norm_data_z) | is.nan(norm_data_z)] <- 0
        myBreaks <- c(
            seq(min(norm_data_z), 0, length.out = ceiling(paletteLength / 2) + 1),
            seq(-1 / paletteLength, max(norm_data_z), length.out = floor(paletteLength / 2))
        )
        ann_colors_arch <- list()
        ann_colors_arch$Phylum <- setNames(
            viridisLite::cividis(n = length(unique(metadata_kingdom$Phylum))),
            levels(as.factor(metadata_kingdom$Phylum))
        )
        phil_name <- paste0("Phylum ", unique(metadata_kingdom$Kingdom))
        names(ann_colors_arch)[names(ann_colors_arch) == "Phylum"] <- phil_name
        metadata_kingdom <- metadata_kingdom %>% select(Phylum)
        colnames(metadata_kingdom) <- phil_name
        heatmap_arch <- pheatmap::pheatmap(norm_data_z[, colOrder_bact],
            annotation_col = NA,
            annotation_colors = ann_colors_arch,
            annotation_row = metadata_kingdom,
            cluster_cols = FALSE,
            cluster_rows = TRUE,
            fontsize = 14,          # testo generale
            fontsize.legend = 18,
            fontsize_row = 18,
            fontsize_col = 12,
            col = myColors3,
            border_color = NA,
            breaks = myBreaks,
            silent = T
        )
        if (!is.null(data_bact[[1]])) heatmap_arch$gtable$widths <- heatmap_list$Bacteria$gtable$widths
        heatmap_list$Archaea <- heatmap_arch
    }

    if (exists("data_euk") && !is.null(data_euk[[1]])) {
        norm_data_z <- data_euk[[1]]
        metadata_hm <- data_euk[[2]]
        metadata_kingdom <- data_euk[[3]]
        norm_data_z[is.na(norm_data_z) | is.nan(norm_data_z)] <- 0
        myBreaks <- c(
            seq(min(norm_data_z), 0, length.out = ceiling(paletteLength / 2) + 1),
            seq(-1 / paletteLength, max(norm_data_z), length.out = floor(paletteLength / 2))
        )
        ann_colors_euk <- list()
        ann_colors_euk$Phylum <- setNames(
            viridisLite::magma(n = length(unique(metadata_kingdom$Phylum))),
            levels(as.factor(metadata_kingdom$Phylum))
        )
        phil_name <- paste0("Phylum ", unique(metadata_kingdom$Kingdom))
        names(ann_colors_euk)[names(ann_colors_euk) == "Phylum"] <- phil_name
        metadata_kingdom <- metadata_kingdom %>% select(Phylum)
        colnames(metadata_kingdom) <- phil_name
        heatmap_euk <- pheatmap::pheatmap(norm_data_z[, colOrder_bact],
            annotation_col = NA,
            annotation_colors = ann_colors_euk,
            annotation_row = metadata_kingdom,
            cluster_cols = FALSE,
            cluster_rows = TRUE,
            col = myColors2,
            fontsize = 14,          # testo generale
            fontsize.legend = 18,
            fontsize_row = 18,
            fontsize_col = 12,
            border_color = NA,
            breaks = myBreaks,
            silent = T
        )
        if (!is.null(data_bact[[1]])) {
            heatmap_euk$gtable$widths <- heatmap_list$Bacteria$gtable$widths
        } else if (!is.null(data_arch[[1]])) heatmap_euk$gtable$widths <- heatmap_list$Archaea$gtable$widths
        heatmap_list$Eukaryota <- heatmap_euk
    }
    ##################### àprova
    # Numero di righe effettive della matrice (non solo del gtable)
    bactrow <- if (!is.null(data_bact[[1]])) nrow(data_bact[[1]]) else 0
    archrow <- if (exists("data_arch") && !is.null(data_arch[[1]])) nrow(data_arch[[1]]) else 0
    eukrow <- if (exists("data_euk") && !is.null(data_euk[[1]])) nrow(data_euk[[1]]) else 0
    sum <- bactrow + archrow + eukrow
    bactrow1 <- 7 + ((bactrow / sum) * 85)
    archrow1 <- if (exists("data_arch") && !is.null(data_arch[[1]])) 4 + ((archrow / sum) * 85) else 0
    eukrow1 <- if (exists("data_euk") && !is.null(data_euk[[1]])) 4 + ((eukrow / sum) * 85) else 0

    labels <- names(heatmap_list)
    text_grobs <- lapply(labels, function(lbl) textGrob(lbl, rot = 90, gp = gpar(fontsize = 14)))

    g <- grid.arrange(
        grobs = c(rbind(text_grobs, lapply(heatmap_list, function(hm) hm$gtable))),
        ncol = 2,
        heights = c(bactrow1, archrow1, eukrow1),
        widths = c(0.05, 1.5)
    )
    # Salva le taxa tables in un unico PDF
    taxa_tables <- list(
        Bacteria = if (!is.null(Bacteria)) as.data.frame(Bacteria@tax_table) else NULL,
        Archaea = if (!is.null(Archaea)) as.data.frame(Archaea@tax_table) else NULL,
        Eukaryota = if (exists("Eukaryota") && !is.null(Eukaryota)) as.data.frame(Eukaryota@tax_table) else NULL
    )
    taxa_tables <- taxa_tables[!sapply(taxa_tables, is.null)]
    if (length(taxa_tables) > 0) {
        for (kingdom in names(taxa_tables)) {
            write.table(
                taxa_tables[[kingdom]],
                file = gsub(" ", "", paste0(output_folder, filename, "_", kingdom, "_taxa_table.tsv")),
                sep = "\t",
                quote = FALSE,
                row.names = TRUE,
                col.names = NA
            )
        }
    }
    # ggsave(
    #     plot = g,
    #     filename = gsub(" ", "", paste0(output_folder, filename, ".pdf")),
    #     height = 15+((bactrow + archrow + eukrow)/5), width = 25, limitsize = FALSE
    # )
}
# --- Fine codice scompattato ---


# ######################
# ##########ALL###########
create_heatmap(
    Bacteria = readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds"),
    Archaea = readRDS(file = "Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.001.rds"),
    Eukaryota = readRDS(file = "Output/SUPERVISED_DEC/Eukaryote_Supervised_decontam0.001.rds"),
    filename = "ALL", output_folder = output_folder, order = "ALL", topBar = "Status"
)

# ########################
# ######### MSHD ############
create_heatmap(
   Bacteria = readRDS(file = "Output/merge_DAS/MSHD/Bacteria_MsHd_merged.rds"),
   Archaea = readRDS(file = "Output/merge_DAS/MSHD/Archaea_MsHd_merged.rds"),
   Eukaryota = readRDS(file = "Output/merge_DAS/MSHD/Eukaryote_MsHd_merged.rds"),
   filename = "category", output_folder = output_folder, order = "ALL", topBar = "Status"
)
# #########GC_LESION############
 create_heatmap(
   Bacteria = readRDS("Output/merge_DAS/GC/Bacteria_GC_merged.rds"),
   Archaea = readRDS("Output/merge_DAS/GC/Archaea_GC_merged.rds"),
   Eukaryota = readRDS("Output/merge_DAS/GC/Eukaryote_GC_merged.rds"),
   filename = "GC_LESION", output_folder = output_folder, order = "gc_treatment", topBar = "gc_treatment"
)

# #########lesion_burden############
create_heatmap(
         Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_lesion_burden_positive_merged.rds")) ,
         Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_lesion_burden_positive_merged.rds")),
         Eukaryota = readRDS(paste0("Output/merge_DAS/GC_comp/Eukaryote_lesion_burden_positive_merged.rds")),
         filename = paste0("lesion_burden_positive"), output_folder = output_folder, order = "Spinal_cord_lesions", topBar = "Spinal_cord_lesions", status = "positive"
)

  create_heatmap(
         Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_lesion_burden_negative_merged.rds")) ,
         Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_lesion_burden_negative_merged.rds")),
         Eukaryota = readRDS(paste0("Output/merge_DAS/GC_comp/Eukaryote_lesion_burden_negative_merged.rds")),
         filename = paste0("lesion_burden_negative"), output_folder = output_folder, order = "Spinal_cord_lesions", topBar = "Spinal_cord_lesions", status = "negative"
)


# #########spinal_cord_lesion############
 create_heatmap(
         Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_spinal_cord_lesion_positive_merged.rds")) ,
         Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_spinal_cord_lesion_positive_merged.rds")),
        Eukaryota =NULL,
       filename = paste0("spinal_cord_lesion_positive"), output_folder = output_folder, order = "Spinal_cord_lesions", topBar = "Spinal_cord_lesions", status = "positive"
)


 create_heatmap(
        Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_spinal_cord_lesion_negative_merged.rds")) ,
         Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_spinal_cord_lesion_negative_merged.rds")),
         Eukaryota = readRDS(paste0("Output/merge_DAS/GC_comp/Eukaryote_spinal_cord_lesion_negative_merged.rds")),
         filename = paste0("spinal_cord_lesion_negative"), output_folder = output_folder, order = "Spinal_cord_lesions", topBar = "Spinal_cord_lesions", status = "negative"
)

 #########subtentorial_lesions############
create_heatmap(
         Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_subtentorial_lesions_positive_merged.rds")) ,
         Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_subtentorial_lesions_positive_merged.rds")),
         Eukaryota = readRDS(paste0("Output/merge_DAS/GC_comp/Eukaryote_subtentorial_lesions_positive_merged.rds")),
         filename = paste0("subtentorial_lesions_positive"), output_folder = output_folder, order = "subtentorial_lesions", topBar = "subtentorial_lesions", status = "positive"
)

 create_heatmap(
         Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_subtentorial_lesions_negative_merged.rds")) ,
         Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_subtentorial_lesions_negative_merged.rds")),
         Eukaryota = NULL,
         filename = paste0("subtentorial_lesions_negative"), output_folder = output_folder, order = "subtentorial_lesions", topBar = "subtentorial_lesions", status = "negative"
)

# #########gadolinium_contrast############
 create_heatmap(
         Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_gadolinium_contrast_positive_merged.rds")) ,
         Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_gadolinium_contrast_positive_merged.rds")),
         Eukaryota = readRDS(paste0("Output/merge_DAS/GC_comp/Eukaryote_gadolinium_contrast_positive_merged.rds")),
         filename = paste0("gadolinium_contrast_positive"), output_folder = output_folder, order = "gadolinium_contrast", topBar = "gadolinium_contrast", status = "positive"
 )
 create_heatmap(
         Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_gadolinium_contrast_negative_merged.rds")) ,
         Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_gadolinium_contrast_negative_merged.rds")),
         Eukaryota = readRDS(paste0("Output/merge_DAS/GC_comp/Eukaryote_gadolinium_contrast_negative_merged.rds")),
         filename = paste0("gadolinium_contrast_negative"), output_folder = output_folder, order = "gadolinium_contrast", topBar = "gadolinium_contrast", status = "negative"
 )
dim=c("01","05","001")
d="05"
for(d in dim){
  create_heatmap(
         Bacteria = readRDS(paste0("Output/merge_DAS/",d,"/Bacteria_gadolinium_contrast_",d,"_merged.rds")) ,
         Archaea = NULL,
         Eukaryota = NULL,
         filename = paste0("gadolinium_contrast_",d), output_folder = output_folder,
          order = "gadolinium_contrast", topBar = "gadolinium_contrast"
 )
 create_heatmap(
         Bacteria = readRDS(paste0("Output/merge_DAS/",d,"/Bacteria_spinal_cord_lesion_",d,"_merged.rds")) ,
         Archaea = NULL,
         Eukaryota = NULL,
         filename = paste0("spinal_cord_lesion_",d), output_folder = output_folder, 
         order = "Spinal_cord_lesions", topBar = "Spinal_cord_lesions"
 )
    create_heatmap(
         Bacteria = readRDS(paste0("Output/merge_DAS/",d,"/Bacteria_subtentorial_lesions_",d,"_merged.rds")) ,
         Archaea = NULL,
         Eukaryota = NULL,
         filename = paste0("subtentorial_lesions_",d), output_folder = output_folder,
          order = "subtentorial_lesions", topBar = "subtentorial_lesions")
 
    create_heatmap(
            Bacteria = readRDS(paste0("Output/merge_DAS/",d,"/Bacteria_lesion_burden_",d,"_merged.rds")) ,
            Archaea = NULL,
            Eukaryota = NULL,
            filename = paste0("lesion_burden_",d), output_folder = output_folder, 
            order = "Lesion_Burden", topBar = "Lesion_Burden")
    
}

