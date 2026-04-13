merge_p_value_maaslin <- function(maaslin3_path, analysis_name) {
    maaslin3 <- read.delim(maaslin3_path, header = TRUE, sep = "\t")

    if (!"metadata" %in% colnames(maaslin3)) {
        stop(paste0("Colonna 'metadata' non trovata in: ", maaslin3_path))
    }

    maaslin3$metadata <- trimws(as.character(maaslin3$metadata))
    maaslin3 <- maaslin3[tolower(maaslin3$metadata) %in% tolower(analysis_name), ]

    maaslin3 <- maaslin3[, c("feature", "pval_individual","coef", "model")]
    rownames(maaslin3)<-paste(maaslin3$feature, maaslin3$model, sep = "_")
    colnames(maaslin3)<-paste0(colnames(maaslin3), "_", analysis_name)
    colnames(maaslin3)[2] <- gsub("pval", "P.val", colnames(maaslin3)[2])

    return(as.data.frame(maaslin3))
}

merge_analysis_group <- function(threshold_label) {
    merged_df_gadolinium_contrast <- merge_p_value_maaslin(
        maaslin3_path = paste0("Output/MAASLIN3/", threshold_label, "_disc_cluster/", threshold_label, "_gadolinium_contrast_onlygc_plus_gc/all_results.tsv"),
        analysis_name = "gadolinium_contrast"
    )

    merged_df_lesion_burden <- merge_p_value_maaslin(
        maaslin3_path = paste0("Output/MAASLIN3/", threshold_label, "_disc_cluster/", threshold_label, "_lesion_burden_onlygc_plus_gc/all_results.tsv"),
        analysis_name = "lesion_burden"
    )

    merged_df_spinal_cord <- merge_p_value_maaslin(
        maaslin3_path = paste0("Output/MAASLIN3/", threshold_label, "_disc_cluster/", threshold_label, "_spinal_cord_lesion_onlygc_plus_gc/all_results.tsv"),
        analysis_name = "spinal_cord_lesion"
    )

    merged_df_subtentorial <- merge_p_value_maaslin(
        maaslin3_path = paste0("Output/MAASLIN3/", threshold_label, "_disc_cluster/", threshold_label, "_subtentorial_lesions_onlygc_plus_gc/all_results.tsv"),
        analysis_name = "subtentorial_lesions"
    )

    all_merged <- Reduce(
        function(x, y) {
            merge(x, y, by = "row.names", all = TRUE) |>
                (\(df) {
                    rownames(df) <- df$Row.names
                    df[, -1, drop = FALSE]
                })()
        },
        list(
            merged_df_gadolinium_contrast,
            merged_df_lesion_burden,
            merged_df_spinal_cord,
            merged_df_subtentorial
        )
    )

    all_merged$species <- rownames(all_merged)
    all_merged <- all_merged[, c("species", setdiff(names(all_merged), "species"))]

    return(all_merged)
}

if (!dir.exists("Output/p_value_Cluster")) {
    dir.create("Output/p_value_Cluster", recursive = TRUE)
}

all_merged_0 <- merge_analysis_group(threshold_label = "0")
write.table(all_merged_0, "Output/p_value_Cluster/0_bacteria_maaslin.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

all_merged_05 <- merge_analysis_group("05")
write.table(all_merged_05, "Output/p_value_Cluster/05_bacteria_maaslin.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

all_merged_01 <- merge_analysis_group("01")
write.table(all_merged_01, "Output/p_value_Cluster/01_bacteria_maaslin.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

all_merged_001 <- merge_analysis_group("001")
write.table(all_merged_001, "Output/p_value_Cluster/001_bacteria_maaslin.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


category_0=merge_p_value_maaslin(
    maaslin3_path = "Output/MAASLIN3/Bacteria_MSHD_disc/all_results.tsv",
    analysis_name = "category"
)

write.table(
    category_0,
    file = "Output/p_value_Cluster/0_bacteria_maaslin_category.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
category_001=merge_p_value_maaslin(
    maaslin3_path = "Output/MAASLIN3/Bacteria_MSHD_disc_001/all_results.tsv",
    analysis_name = "category"
)
write.table(
    category_001,
    file = "Output/p_value_Cluster/001_bacteria_maaslin_category.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
