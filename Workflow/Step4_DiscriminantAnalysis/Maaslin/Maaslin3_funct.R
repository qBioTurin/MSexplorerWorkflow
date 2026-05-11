Maaslin3_funct <- function(
    rds, output_folder, main_el, categorySel, elements, discr_els, perc = "K",model_to_use="both", median_comparison_abundance_in = TRUE,
    fx_effects = NULL, random_effects = NULL, group_effects = NULL,
    ordered_effects = NULL, strata_effects = NULL, heatmap_vars = NULL) {
    # Create data for maaslin3
    otu <- rds@otu_table
    taxa <- rds@tax_table
    metadata_df <- as(rds@sam_data, "data.frame")
    ids_to_retain <- rownames(metadata_df[metadata_df[[categorySel]] %in% elements, ])
    otu <- otu[, colnames(otu) %in% ids_to_retain]
    metadata_df <- metadata_df[rownames(metadata_df) %in% ids_to_retain, ]
    taxa <- as.data.frame(taxa)
    taxa$Genus_Species <- paste(taxa$Genus, taxa$Species, sep = "_")
    Tnew <- taxa %>% select(Genus_Species)
    Tnew$taxid <- rownames(Tnew)
    if (discr_els %>% length() == 0) {
        formulaInUse <- main_el
    } else {
        formulaInUse <- paste(main_el, paste(discr_els, collapse = "+"), sep = "+")
    }
    otu <- as.data.frame(otu)
    otu$taxid <- rownames(otu)
    # Merge the two data frames by taxid
    merged_df <- merge(otu, Tnew, by = "taxid")
    rownames(merged_df) <- merged_df$Genus_Species
    merged_df <- merged_df %>% select(-c(taxid, Genus_Species))
    merged_df <- t(merged_df)

    # Create metadata for masslin3

    class(metadata_df)
    metadata_df <- metadata_df %>%
        select(category, gc_treatment, lesion_burden, spinal_cord_lesion, gadolinium_contrast, subtentorial_lesions, sex, age, bmi,Cluster)
    if (main_el != "category") {
        metadata_df <- metadata_df %>%
            filter(category == "MS")
        metadata_df$gc_treatment <- factor(metadata_df$gc_treatment,
            levels = c("positive", "negative")
        )
    }
    metadata_df$category <- as.factor(metadata_df$category)
    levels(metadata_df$category) <- c("HEALTHY", "MS")
    metadata_df$sex <- as.factor(metadata_df$sex)
    levels(metadata_df$sex) <- c("F", "M")
    metadata_df$age <- as.numeric(metadata_df$age)
    metadata_df$bmi <- as.numeric(metadata_df$bmi)
    metadata_df$gc_treatment <- as.factor(metadata_df$gc_treatment)
    metadata_df$category <- as.factor(metadata_df$category)
    metadata_df$lesion_burden <- as.factor(metadata_df$lesion_burden)
    levels(metadata_df$lesion_burden) <- c("low", "high")
    metadata_df$spinal_cord_lesion <- as.factor(metadata_df$spinal_cord_lesion)
    levels(metadata_df$spinal_cord_lesion) <- c("BM_low", "BM_high")
    metadata_df$gadolinium_contrast <- as.factor(metadata_df$gadolinium_contrast)
    levels(metadata_df$gadolinium_contrast) <- c("NoActive", "Active")
    metadata_df$subtentorial_lesions <- as.factor(metadata_df$subtentorial_lesions)
    levels(metadata_df$subtentorial_lesions) <- c("No", "Yes")
    metadata_df$Cluster <- as.factor(metadata_df$Cluster)
    levels(metadata_df$Cluster) <- c("1","2")
    metadata_df$sex <- as.factor(metadata_df$sex)
    metadata_df$age <- as.numeric(metadata_df$age)
    metadata_df$bmi <- as.numeric(metadata_df$bmi)


    if (!dir.exists(output_folder)) {
        dir.create(output_folder, recursive = TRUE)
    }
    # el_plot <- levels(metadata_df[[main_el]])[1]
    el_plot2 <- levels(metadata_df[[main_el]])[2]
    # main_el_plot <- paste(main_el, el_plot)
    second_el_plot <- paste(main_el, el_plot2)
    coef_plot_vars <- c(second_el_plot)
    formulaInUse <- paste("~", formulaInUse, sep = "")
    discr_els_formatted <- unlist(lapply(discr_els, function(d) {
        if (is.factor(metadata_df[[d]])) {
            paste(d, levels(metadata_df[[d]])) # tutti i livelli
        } else {
            d # continuous
        }
    }))
    ##for stability of results across runs
    set.seed(2204)
    Sys.setenv(OMP_NUM_THREADS=1, MKL_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1)
    fit_out <- maaslin3(
        input_data = merged_df,
        input_metadata = metadata_df,
        output = output_folder,
        formula = formulaInUse,
        fixed_effects = fx_effects,
        plot_summary_plot = TRUE,
        coef_plot_vars = coef_plot_vars,
        heatmap_vars = heatmap_vars,
        random_effects = random_effects,
        group_effects = group_effects,
        ordered_effects = ordered_effects,
        strata_effects = strata_effects,
        normalization = "TSS",
        transform = "LOG",
        augment = TRUE,
        standardize = TRUE,
        max_significance = 0.1,
        median_comparison_abundance = median_comparison_abundance_in,
        median_comparison_prevalence = FALSE,
        max_pngs = 30,
        summary_plot_first_n = 200,
       # pval_filter = 0.05,
        cores = 1
    )
    summary_plot_file <- file.path(output_folder, "figures", "summary_plot.pdf")
    # if (file.exists(summary_plot_file)) {
    #     try(system2("xdg-open", summary_plot_file), silent = TRUE)
    # } else {
    #     message("summary_plot.pdf not generated for this run (no plottable/significant associations).")
    # }
    if (median_comparison_abundance_in) {
        last <- "T"
    } else {
        last <- "F"
    }
    all_results <- read_tsv(file.path(output_folder, "all_results.tsv"), show_col_types = FALSE)
    all_results_sig <- all_results %>% dplyr::filter(metadata == main_el & pval_individual < 0.05 & (model_to_use=="both" | model == model_to_use))
    all_results_sig <- all_results_sig %>%
        select(feature, metadata, value, coef, null_hypothesis, pval_individual, qval_individual, pval_joint, pval_individual, model, stderr)

    if (model_to_use == "both") {
        combined_pvals <- all_results %>%
            dplyr::filter(metadata == main_el, model %in% c("abundance", "prevalence")) %>%
            dplyr::select(feature, model, pval_individual) %>%
            dplyr::group_by(feature, model) %>%
            dplyr::summarise(
                pval_individual = if (all(is.na(pval_individual))) NA_real_ else min(pval_individual, na.rm = TRUE),
                .groups = "drop"
            ) %>%
            tidyr::pivot_wider(
                names_from = model,
                values_from = pval_individual,
                names_glue = "pval_individual_{model}"
            ) %>%
            dplyr::filter(feature %in% unique(all_results_sig$feature)) %>%
            dplyr::arrange(feature)

        combined_pvals_csv <- file.path(
            output_folder,
            paste0(main_el, "_MAASLIN3_SR_species_pvals_abundance_prevalence_", perc, "_", last, ".csv")
        )
        write.csv(combined_pvals, combined_pvals_csv, row.names = FALSE)
    }

    results_csv <- file.path(output_folder, paste0(main_el, "_MAASLIN3_SR_", perc, "_", last, ".csv"))
    write.csv(all_results_sig, results_csv, row.names = FALSE)
    features_unique <- unique(all_results_sig$feature)
    features_unique <- data.frame(Species = features_unique)
    features_csv <- file.path(output_folder, paste0(main_el, "_MAASLIN3_SR_features_", perc, "_", last, ".csv"))
    write.csv(features_unique, features_csv, row.names = FALSE, quote = FALSE)

    features_from_csv <- read.csv(features_csv, stringsAsFactors = FALSE)
    selected_species <- unique(stats::na.omit(features_from_csv$Species))

    rds_filtered <- rds
    rds_filtered <- prune_samples(sample_names(rds_filtered) %in% ids_to_retain, rds_filtered)

    if (main_el != "category") {
        rds_filtered <- subset_samples(rds_filtered, category == "MS")
    }

    taxa_df <- as.data.frame(tax_table(rds_filtered))
    taxa_df$Genus_Species <- paste(taxa_df$Genus, taxa_df$Species, sep = "_")
    taxa_to_keep <- rownames(taxa_df)[taxa_df$Genus_Species %in% selected_species]
    rds_filtered <- prune_taxa(taxa_to_keep, rds_filtered)

    rds_file <- file.path(output_folder, paste0(main_el, "_MAASLIN3_SR_features_", perc, "_", last, ".rds"))
    saveRDS(rds_filtered, rds_file)
    print(rds_filtered)
}
