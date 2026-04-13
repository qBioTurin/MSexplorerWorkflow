compare_species_lists <- function(rds_path, csv_path) {
    ll_sp <- readRDS(rds_path)
    species <- as.data.frame(tax_table(ll_sp))
    species$Species <- paste0(species$Genus, "_", species$Species)
    species_list_ll <- species$Species

    maaslin <- read_csv(csv_path)
    species_list_maaslin <- maaslin$Species

    return(list(species_list_ll = species_list_ll, species_list_maaslin = species_list_maaslin))
}
compare_and_save_lists <- function(list1, list2, outdir, nameFile,
                                   name1 = "limma_lefse", name2 = "maaslin3") {
    library(VennDiagram)
    library(gridExtra)
    library(grid)

    vec1 <- unique(as.character(na.omit(unlist(list1))))
    vec2 <- unique(as.character(na.omit(unlist(list2))))

    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    only1 <- setdiff(vec1, vec2)
    only2 <- setdiff(vec2, vec1)
    common <- intersect(vec1, vec2)

    # salva tabella dei risultati
    all_species <- sort(unique(c(vec1, vec2)))
    presence <- vapply(all_species, function(sp) {
        in1 <- sp %in% vec1
        in2 <- sp %in% vec2
        if (in1 && in2) {
            "both"
        } else if (in1) {
            paste0("only_", name1)
        } else {
            paste0("only_", name2)
        }
    }, character(1))

    mat <- data.frame(
        Species = all_species,
        Presence = presence,
        stringsAsFactors = FALSE
    )

    write.table(mat,
        file = file.path(outdir, paste0(nameFile, ".csv")),
        sep = ",", col.names = TRUE, row.names = FALSE, quote = TRUE
    )

    # --- genera il Venn diagram con fill uguale al contorno ---
    venn_list <- setNames(list(vec1, vec2), c(name1, name2))

    vd <- venn.diagram(
        x = venn_list,
        filename = NULL,
        lwd = 2,
        lty = "solid",
        fill = c("gold", "seagreen"), # fill uguale al contorno
        col = c("gold", "seagreen"), # contorni
        cex = 1.2, # dimensione numeri dentro gli insiemi
        fontface = "bold",
        fontfamily = "sans",
        cat.cex = 0.8, # dimensione dei nomi dei cerchi
        cat.fontface = "bold",
        cat.default.pos = "outer",
        cat.pos = c(-27, 27),
        cat.dist = c(0.055, 0.055),
        cat.fontfamily = "sans"
    )
    venn_grob <- grid.grabExpr(grid.draw(vd))

    # --- prepara le tabelle di specie sotto con nomi colonne personalizzate ---
    df_ll <- if (length(only1) == 0) data.frame(tmp = "(nessuna)") else data.frame(only1)
    colnames(df_ll) <- name1
    df_common <- if (length(common) == 0) data.frame(tmp = "(nessuna)") else data.frame(common)
    colnames(df_common) <- "both"
    df_ma <- if (length(only2) == 0) data.frame(tmp = "(nessuna)") else data.frame(only2)
    colnames(df_ma) <- name2

    p_ll <- tableGrob(df_ll, rows = NULL, theme = ttheme_default(core = list(fg_params = list(cex = 0.8))))
    p_common <- tableGrob(df_common, rows = NULL, theme = ttheme_default(core = list(fg_params = list(cex = 0.8))))
    p_ma <- tableGrob(df_ma, rows = NULL, theme = ttheme_default(core = list(fg_params = list(cex = 0.8))))

    bottom_table <- arrangeGrob(p_ll, p_common, p_ma, ncol = 3)

    # --- combina diagramma + specie sotto ---
    final_plot <- grid.arrange(
        venn_grob,
        bottom_table,
        nrow = 2,
        heights = c(1, 2),
        top = textGrob(nameFile, gp = gpar(fontsize = 16, fontface = "bold"))
    )

    # salva PDF
    pdf(file.path(outdir, paste0(nameFile, "_venn.pdf")), width = 12, height = 10 + (max(c(length(only1), length(only2), length(common))) * 0.25))
    grid.draw(final_plot)
    dev.off()
    invisible(list(only1 = only1, only2 = only2, common = common))
}

gado_01 <- compare_species_lists(
    rds_path = "Output/merge_DAS/01/Bacteria_gadolinium_contrast_01_merged.rds",
    csv_path = "Output/MAASLIN3/01_disc/01_gadolinium_contrast_onlygc_T/gadolinium_contrast_MAASLIN3_SR_features_01.csv"
)
compare_and_save_lists(
    list1 = gado_01$species_list_ll,
    list2 = gado_01$species_list_maaslin,
    outdir = "Output/compare_DAS_MAASLIN3/gadolinium_contrast_01",
    nameFile = "gadolinium_contrast_01"
)


lesion_01 <- compare_species_lists(
    rds_path = "Output/merge_DAS/01/Bacteria_lesion_burden_01_merged.rds",
    csv_path = "Output/MAASLIN3/01_disc/01_lesion_burden_onlygc_T/lesion_burden_MAASLIN3_SR_features_01.csv"
)
compare_and_save_lists(
    list1 = lesion_01$species_list_ll,
    list2 = lesion_01$species_list_maaslin,
    outdir = "Output/compare_DAS_MAASLIN3/lesion_burden_01",
    nameFile = "lesion_burden_01"
)

spinal_01 <- compare_species_lists(
    rds_path = "Output/merge_DAS/01/Bacteria_spinal_cord_lesion_01_merged.rds",
    csv_path = "Output/MAASLIN3/01_disc/01_spinal_cord_lesion_onlygc_T/spinal_cord_lesion_MAASLIN3_SR_features_01.csv"
)
compare_and_save_lists(
    list1 = spinal_01$species_list_ll,
    list2 = spinal_01$species_list_maaslin,
    outdir = "Output/compare_DAS_MAASLIN3/spinal_cord_lesion_01",
    nameFile = "spinal_cord_lesion_01"
)
subtentorial_01 <- compare_species_lists(
    rds_path = "Output/merge_DAS/01/Bacteria_subtentorial_lesions_01_merged.rds",
    csv_path = "Output/MAASLIN3/01_disc/01_subtentorial_lesions_onlygc_T/subtentorial_lesions_MAASLIN3_SR_features_01.csv"
)
compare_and_save_lists(
    list1 = subtentorial_01$species_list_ll,
    list2 = subtentorial_01$species_list_maaslin,
    outdir = "Output/compare_DAS_MAASLIN3/subtentorial_lesions_01",
    nameFile = "subtentorial_lesions_01"
)

# combina le liste di specie dei 4 confronti e genera un confronto finale
combined_ll <- unique(as.character(na.omit(unlist(c(
    gado_01$species_list_ll,
    lesion_01$species_list_ll,
    spinal_01$species_list_ll,
    subtentorial_01$species_list_ll
)))))

combined_ma <- unique(as.character(na.omit(unlist(c(
    gado_01$species_list_maaslin,
    lesion_01$species_list_maaslin,
    spinal_01$species_list_maaslin,
    subtentorial_01$species_list_maaslin
)))))

dir_final <- "Output/compare_DAS_MAASLIN3/final_01"
dir.create(dir_final, recursive = TRUE, showWarnings = FALSE)

# salva le liste combinate
write.table(data.frame(Species = sort(combined_ll)),
    file = file.path(dir_final, "combined_limma_lefse_01.csv"),
    sep = ",", row.names = FALSE, quote = TRUE
)
write.table(data.frame(Species = sort(combined_ma)),
    file = file.path(dir_final, "combined_maaslin3_01.csv"),
    sep = ",", row.names = FALSE, quote = TRUE
)

# esegue il confronto finale
compare_and_save_lists(
    list1 = combined_ll,
    list2 = combined_ma,
    outdir = dir_final,
    nameFile = "comparison_01",
    name1 = "limma_lefse_all",
    name2 = "maaslin3_all"
)


gado_001 <- compare_species_lists(
    rds_path = "Output/merge_DAS/001/Bacteria_gadolinium_contrast_001_merged.rds",
    csv_path = "Output/MAASLIN3/001_disc/001_gadolinium_contrast_onlygc_plus_gc/gadolinium_contrast_MAASLIN3_SR_features_001.csv"
)

compare_and_save_lists(
    list1 = gado_001$species_list_ll,
    list2 = gado_001$species_list_maaslin,
    outdir = "Output/compare_DAS_MAASLIN3/gadolinium_contrast_001",
    nameFile = "gadolinium_contrast_001"
)
lesion_001 <- compare_species_lists(
    rds_path = "Output/merge_DAS/001/Bacteria_lesion_burden_001_merged.rds",
    csv_path = "Output/MAASLIN3/001_disc/001_lesion_burden_onlygc_plus_gc/lesion_burden_MAASLIN3_SR_features_001.csv"
)
compare_and_save_lists(
    list1 = lesion_001$species_list_ll,
    list2 = lesion_001$species_list_maaslin,
    outdir = "Output/compare_DAS_MAASLIN3/lesion_burden_001",
    nameFile = "lesion_burden_001"
)
spinal_001 <- compare_species_lists(
    rds_path = "Output/merge_DAS/001/Bacteria_spinal_cord_lesion_001_merged.rds",
    csv_path = "Output/MAASLIN3/001_disc/001_spinal_cord_lesion_onlygc_plus_gc/spinal_cord_lesion_MAASLIN3_SR_features_001.csv"
)
compare_and_save_lists(
    list1 = spinal_001$species_list_ll,
    list2 = spinal_001$species_list_maaslin,
    outdir = "Output/compare_DAS_MAASLIN3/spinal_cord_lesion_001",
    nameFile = "spinal_cord_lesion_001"
)
subtentorial_001 <- compare_species_lists(
    rds_path = "Output/merge_DAS/001/Bacteria_subtentorial_lesions_001_merged.rds",
    csv_path = "Output/MAASLIN3/001_disc/001_subtentorial_lesions_onlygc_plus_gc/subtentorial_lesions_MAASLIN3_SR_features_001.csv"
)
compare_and_save_lists(
    list1 = subtentorial_001$species_list_ll,
    list2 = subtentorial_001$species_list_maaslin,
    outdir = "Output/compare_DAS_MAASLIN3/subtentorial_lesions_001",
    nameFile = "subtentorial_lesions_001"
)

# combina le liste di specie dei 4 confronti e genera un confronto finale
combined_ll_001 <- unique(as.character(na.omit(unlist(c(
    gado_001$species_list_ll,
    lesion_001$species_list_ll,
    spinal_001$species_list_ll,
    subtentorial_001$species_list_ll
)))))
combined_ma_001 <- unique(as.character(na.omit(unlist(c(
    gado_001$species_list_maaslin,
    lesion_001$species_list_maaslin,
    spinal_001$species_list_maaslin,
    subtentorial_001$species_list_maaslin
)))))
dir_final_001 <- "Output/compare_DAS_MAASLIN3/final_001"
dir.create(dir_final_001, recursive = TRUE, showWarnings = FALSE)
# salva le liste combinate
write.table(data.frame(Species = sort(combined_ll_001)),
    file = file.path(dir_final_001, "combined_limma_lefse_001.csv"),
    sep = ",", row.names = FALSE, quote = TRUE
)
write.table(data.frame(Species = sort(combined_ma_001)),
    file = file.path(dir_final_001, "combined_maaslin3_001.csv"),
    sep = ",", row.names = FALSE, quote = TRUE
)
# esegue il confronto finale
compare_and_save_lists(
    list1 = combined_ll_001,
    list2 = combined_ma_001,
    outdir = dir_final_001,
    nameFile = "comparison_001",
    name1 = "limma_lefse_all",
    name2 = "maaslin3_all"
)

gado_05 <- compare_species_lists(
    rds_path = "Output/merge_DAS/05/Bacteria_gadolinium_contrast_05_merged.rds",
    csv_path = "Output/MAASLIN3/05_disc/05_gadolinium_contrast_onlygc_T/gadolinium_contrast_MAASLIN3_SR_features_05.csv"
)
compare_and_save_lists(
    list1 = gado_05$species_list_ll,
    list2 = gado_05$species_list_maaslin,
    outdir = "Output/compare_DAS_MAASLIN3/gadolinium_contrast_05",
    nameFile = "gadolinium_contrast_05"
)
lesion_05 <- compare_species_lists(
    rds_path = "Output/merge_DAS/05/Bacteria_lesion_burden_05_merged.rds",
    csv_path = "Output/MAASLIN3/05_disc/05_lesion_burden_onlygc_T/lesion_burden_MAASLIN3_SR_features_05.csv"
)
compare_and_save_lists(
    list1 = lesion_05$species_list_ll,
    list2 = lesion_05$species_list_maaslin,
    outdir = "Output/compare_DAS_MAASLIN3/lesion_burden_05",
    nameFile = "lesion_burden_05"
)
spinal_05 <- compare_species_lists(
    rds_path = "Output/merge_DAS/05/Bacteria_spinal_cord_lesion_05_merged.rds",
    csv_path = "Output/MAASLIN3/05_disc/05_spinal_cord_lesion_onlygc_T/spinal_cord_lesion_MAASLIN3_SR_features_05.csv"
)
compare_and_save_lists(
    list1 = spinal_05$species_list_ll,
    list2 = spinal_05$species_list_maaslin,
    outdir = "Output/compare_DAS_MAASLIN3/spinal_cord_lesion_05",
    nameFile = "spinal_cord_lesion_05"
)
subtentorial_05 <- compare_species_lists(
    rds_path = "Output/merge_DAS/05/Bacteria_subtentorial_lesions_05_merged.rds",
    csv_path = "Output/MAASLIN3/05_disc/05_subtentorial_lesions_onlygc_T/subtentorial_lesions_MAASLIN3_SR_features_05.csv"
)
compare_and_save_lists(
    list1 = subtentorial_05$species_list_ll,
    list2 = subtentorial_05$species_list_maaslin,
    outdir = "Output/compare_DAS_MAASLIN3/subtentorial_lesions_05",
    nameFile = "subtentorial_lesions_05"
)
# combina le liste di specie dei 4 confronti e genera un confronto finale
combined_ll_05 <- unique(as.character(na.omit(unlist(c(
    gado_05$species_list_ll,
    lesion_05$species_list_ll,
    spinal_05$species_list_ll,
    subtentorial_05$species_list_ll
)))))
combined_ma_05 <- unique(as.character(na.omit(unlist(c(
    gado_05$species_list_maaslin,
    lesion_05$species_list_maaslin,
    spinal_05$species_list_maaslin,
    subtentorial_05$species_list_maaslin
)))))
dir_final_05 <- "Output/compare_DAS_MAASLIN3/final_05"
dir.create(dir_final_05, recursive = TRUE, showWarnings = FALSE)
# salva le liste combinate
write.table(data.frame(Species = sort(combined_ll_05)),
    file = file.path(dir_final_05, "combined_limma_lefse_05.csv"),
    sep = ",", row.names = FALSE, quote = TRUE
)
write.table(data.frame(Species = sort(combined_ma_05)),
    file = file.path(dir_final_05, "combined_maaslin3_05.csv"),
    sep = ",", row.names = FALSE, quote = TRUE
)
# esegue il confronto finale
compare_and_save_lists(
    list1 = combined_ll_05,
    list2 = combined_ma_05,
    outdir = dir_final_05,
    nameFile = "comparison_05",
    name1 = "limma_lefse_all",
    name2 = "maaslin3_all"
)


MS_HDBact <- compare_species_lists(
    rds_path = "Output/merge_DAS/MSHD/Bacteria_MsHd_merged.rds",
    csv_path = "Output/MAASLIN3/Bacteria_MSHD_disc_T/category_MAASLIN3_SR_features_K.csv"
)
compare_and_save_lists(
    list1 = MS_HDBact$species_list_ll,
    list2 = MS_HDBact$species_list_maaslin,
    outdir = "Output/compare_DAS_MAASLIN3/MS_HD_Bacteria",
    nameFile = "MS_HD_Bacteria"
)

MS_HDArch <- compare_species_lists(
    rds_path = "Output/merge_DAS/MSHD/Archaea_MsHd_merged.rds",
    csv_path = "Output/MAASLIN3/Archaea_MSHD_disc_T/category_MAASLIN3_SR_features_K.csv"
)
compare_and_save_lists(
    list1 = MS_HDArch$species_list_ll,
    list2 = MS_HDArch$species_list_maaslin,
    outdir = "Output/compare_DAS_MAASLIN3/MS_HD_Archaea",
    nameFile = "MS_HD_Archaea"
)
MS_HDEuk <- compare_species_lists(
    rds_path = "Output/merge_DAS/MSHD/Eukaryote_MsHd_merged.rds",
    csv_path = "Output/MAASLIN3/Eukaryote_MSHD_disc_T/category_MAASLIN3_SR_features_K.csv"
)
compare_and_save_lists(
    list1 = MS_HDEuk$species_list_ll,
    list2 = MS_HDEuk$species_list_maaslin,
    outdir = "Output/compare_DAS_MAASLIN3/MS_HD_Eukaryote",
    nameFile = "MS_HD_Eukaryote"
)

GC_Bact <- compare_species_lists(
    rds_path = "Output/merge_DAS/GC/Bacteria_GC_merged.rds",
    csv_path = "Output/MAASLIN3/Bacteria_GC_TREATMENT_disc_T/gc_treatment_MAASLIN3_SR_features_K.csv"
)
compare_and_save_lists(
    list1 = GC_Bact$species_list_ll,
    list2 = GC_Bact$species_list_maaslin,
    outdir = "Output/compare_DAS_MAASLIN3/GC_Bacteria",
    nameFile = "GC_Bacteria"
)
GC_euk <- compare_species_lists(
    rds_path = "Output/merge_DAS/GC/Eukaryote_GC_merged.rds",
    csv_path = "Output/MAASLIN3/Eukaryote_GC_TREATMENT_disc_T/gc_treatment_MAASLIN3_SR_features_K.csv"
)
compare_and_save_lists(
    list1 = GC_euk$species_list_ll,
    list2 = GC_euk$species_list_maaslin,
    outdir = "Output/compare_DAS_MAASLIN3/GC_Eukaryote",
    nameFile = "GC_Eukaryote"
)
GC_Arch <- compare_species_lists(
    rds_path = "Output/merge_DAS/GC/Archaea_GC_merged.rds",
    csv_path = "Output/MAASLIN3/Archaea_GC_TREATMENT_disc_T/gc_treatment_MAASLIN3_SR_features_K.csv"
)
compare_and_save_lists(
    list1 = GC_Arch$species_list_ll,
    list2 = GC_Arch$species_list_maaslin,
    outdir = "Output/compare_DAS_MAASLIN3/GC_Archaea",
    nameFile = "GC_Archaea"
)

combined_gado_001 <- unique(c(gado_001$species_list_ll, gado_001$species_list_maaslin))
dir_gado_001 <- "Output/compare_DAS_MAASLIN3/gadolinium_contrast_001_union"
dir.create(dir_gado_001, recursive = TRUE, showWarnings = FALSE)
write.table(data.frame(Species = sort(combined_gado_001)),
    file = file.path(dir_gado_001, "combined_species_001.csv"),
    sep = ",", row.names = FALSE, quote = TRUE
)
combined_lesion_001 <- unique(c(lesion_001$species_list_ll, lesion_001$species_list_maaslin))
dir_lesion_001 <- "Output/compare_DAS_MAASLIN3/lesion_burden_001_union"
dir.create(dir_lesion_001, recursive = TRUE, showWarnings = FALSE)
write.table(data.frame(Species = sort(combined_lesion_001)),
    file = file.path(dir_lesion_001, "combined_species_001.csv"),
    sep = ",", row.names = FALSE, quote = TRUE
)

combined_spinal_001 <- unique(c(spinal_001$species_list_ll, spinal_001$species_list_maaslin))
dir_spinal_001 <- "Output/compare_DAS_MAASLIN3/spinal_cord_lesion_001_union"
dir.create(dir_spinal_001, recursive = TRUE, showWarnings = FALSE)
write.table(data.frame(Species = sort(combined_spinal_001)),
    file = file.path(dir_spinal_001, "combined_species_001.csv"),
    sep = ",", row.names = FALSE, quote = TRUE
)

combined_subtentorial_001 <- unique(c(subtentorial_001$species_list_ll, subtentorial_001$species_list_maaslin))
dir_subtentorial_001 <- "Output/compare_DAS_MAASLIN3/subtentorial_lesions_001_union"
dir.create(dir_subtentorial_001, recursive = TRUE, showWarnings = FALSE)
write.table(data.frame(Species = sort(combined_subtentorial_001)),
    file = file.path(dir_subtentorial_001, "combined_species_001.csv"),
    sep = ",", row.names = FALSE, quote = TRUE
)
