source("Settings/utilities.R")
output_folder = "Output/NEW_PLOT_RDS/"
seed_beta <- 123
createFolder(output_folder)

format_permanova_pvalue <- function(p_value) {
    if (is.na(p_value)) {
        return("NA")
    }

    if (p_value < 0.001) {
        return("< 0.001")
    }

    formatC(signif(p_value, 3), format = "fg", flag = "#")
}

Beta <- function(baselines_dec, kingdom, output_folder) {
    baselines_dec <- prune_taxa(taxa_sums(baselines_dec) > 0, baselines_dec)
    baselines_dec <- prune_samples(sample_sums(baselines_dec) > 0, baselines_dec)

    ## CATEGORY
    ############
    # New labels for category
    sample_data(baselines_dec)$category <- factor(
        sample_data(baselines_dec)$category, 
        levels = c("HEALTHY", "MS"), 
        labels = c("HD", "MS")
    )

    # Calculate Bray
    wUF.ordu = ordinate(baselines_dec, method = "PCoA", distance = "bray")
    baselines_dec = phyloseq(
        otu_table(baselines_dec), 
        tax_table(baselines_dec), 
        sample_data(baselines_dec)
    )

    # Plot PCoA
    # Category
    category = plot_ordination(
        baselines_dec, 
        wUF.ordu, 
        type = "sample_type", 
        color = "category"
    ) + 
        theme_classic() + 
        stat_ellipse() + 
        ggtitle("PCoA of bray Curtis distance") + 
        guides(color = guide_legend(title = "Disease Status")) +
        geom_point(size = 4) +
        scale_color_manual(values = rev(colors_venn)) +
        theme(
            axis.title.x = element_text(size = 15),  
            axis.title.y = element_text(size = 15),  
            axis.text.x = element_text(size = 5),  
            axis.text.y = element_text(size = 15),   
            plot.title = element_text(size = 20),    
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15)
        )
    bray_dist <- phyloseq::distance(baselines_dec, method = "bray")
    sample_df <- data.frame(sample_data(baselines_dec))
    set.seed(seed_beta)
    permanova_cat <- vegan::adonis2(bray_dist ~ category,
                                data = sample_df,
                                permutations = 9999)
    pval_cat <- permanova_cat$`Pr(>F)`[1]
    pval_cat_label <- paste0("PERMANOVA p = ", format_permanova_pvalue(pval_cat))
    category <- category +
        labs(subtitle = pval_cat_label) +
        theme(plot.subtitle = element_text(size = 14))
    print(paste("PERMANOVA p-value for category in", kingdom, ":", pval_cat))
  
    saveRDS(category, gsub(" ", "", paste(output_folder, kingdom, "_beta_cat.rds")))

    ## GC_TREATMENT
    ###############
    # Remove NAs from gc_treatment and subset only MS samples
    MS = subset_samples(baselines_dec, category == "MS")

    sample_data_df <- data.frame(sample_data(baselines_dec))
    filtered_sample_data_df <- sample_data_df %>%
        filter(!is.na(gc_treatment)) %>%
        filter(category == "MS")

    # Update the phyloseq object with the filtered sample data
    sample_data(MS) <- sample_data(filtered_sample_data_df)
    MS <- prune_taxa(taxa_sums(MS) > 0, MS)
    MS <- prune_samples(sample_sums(MS) > 0, MS)

    # Change legend label
    sample_data(MS)$gc_treatment <- factor(
        sample_data(MS)$gc_treatment, 
        levels = c("negative", "positive")
    )

    # Calculate Bray
    ordination <- ordinate(MS, method = "PCoA", distance = "bray")
    gc_treatment = plot_ordination(
        MS, 
        ordination, 
        type = "sample_type", 
        color = "gc_treatment"
    ) + 
        theme_classic() + 
        stat_ellipse() + 
        ggtitle("PCoA of bray Curtis distance") + 
        guides(color = guide_legend(title = "Glucocorticoid Treatment")) +
        geom_point(size = 4) +
        scale_color_manual(values=c("negative" = "#D7D7D7", "positive" = "#4D4D4D")) +
        theme(
            axis.title.x = element_text(size = 15),  
            axis.title.y = element_text(size = 15),  
            axis.text.x = element_text(size = 5),   
            axis.text.y = element_text(size = 15),   
            plot.title = element_text(size = 20),    
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15)
        )
    bray_dist <- phyloseq::distance(MS, method = "bray")
    sample_df <- data.frame(sample_data(MS))
    set.seed(seed_beta)
    permanova_cat <- vegan::adonis2(bray_dist ~ gc_treatment, data = sample_df)
    pval_cat <- permanova_cat$`Pr(>F)`[1]
    pval_cat_label <- paste0("PERMANOVA p = ", format_permanova_pvalue(pval_cat))
    gc_treatment <- gc_treatment +
        labs(subtitle = pval_cat_label) +
        theme(plot.subtitle = element_text(size = 14))
    print(paste("PERMANOVA p-value for gc_treatment in", kingdom, ":", pval_cat))  # Print the p-value    
    saveRDS(gc_treatment, gsub(" ", "", paste(output_folder, kingdom, "_beta_gc.rds")))
    saveRDS(gc_treatment, gsub(" ", "", paste(output_folder, kingdom, "_beta_gc.rds")))
}

execute_beta <- function() {
    baselines_decB = readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
    baselines_decA = readRDS(file = "Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.001.rds")
    baselines_decE = readRDS(file = "Output/SUPERVISED_DEC/Eukaryote_Supervised_decontam0.001.rds")

    Beta(baselines_decB, "Bacteria", output_folder)
    Beta(baselines_decA, "Archaea", output_folder)
    Beta(baselines_decE, "Eukaryote", output_folder)    
}

execute_beta()
