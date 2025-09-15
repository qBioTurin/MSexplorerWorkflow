source("Settings/utilities.R")
output_folder = "Output/NEW_PLOT_RDS/"
createFolder(output_folder)

Beta <- function(baselines_dec, kingdom, output_folder) {
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
            axis.text.x = element_text(size = 15),  
            axis.text.y = element_text(size = 15),   
            plot.title = element_text(size = 20),    
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15)
        )
    bray_dist <- phyloseq::distance(baselines_dec, method = "bray")
    sample_df <- data.frame(sample_data(baselines_dec))
    permanova_cat <- vegan::adonis2(bray_dist ~ category, data = sample_df)
    pval_cat <- signif(permanova_cat$`Pr(>F)`[1], 3) 
    print(paste("PERMANOVA p-value for category in", kingdom, ":", pval_cat))  # Print the p-value
    ####pval_cat 0.001
    saveRDS(category, gsub(" ", "", paste(output_folder, kingdom, "_beta_cat.rds")))

    ## GC_TREATMENT
    ###############
    # Remove NAs from gc_treatment and subset only MS samples
    #MS = subset_samples(baselines_dec, category == c("MS","HEALTHY"))
    MS = baselines_dec
    sample_data_df <- data.frame(sample_data(baselines_dec))
    filtered_sample_data_df <- sample_data_df %>%
        filter(!is.na(gc_treatment)) 

    # Update the phyloseq object with the filtered sample data
    sample_data(MS) <- sample_data(filtered_sample_data_df)

    # Change legend label
    sample_data(MS)$gc_treatment <- factor(
        sample_data(MS)$gc_treatment, 
        levels = c("negative", "positive", "healthy")
    )
    unique(sample_data(MS)$gc_treatment)

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
        scale_color_manual(values = c("healthy" = "#6EE2FF99", "positive" = "#4D4D4D", "negative" = "#D7D7D7")) +
        theme(
            axis.title.x = element_text(size = 15),  
            axis.title.y = element_text(size = 15),  
            axis.text.x = element_text(size = 15),   
            axis.text.y = element_text(size = 15),   
            plot.title = element_text(size = 20),    
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15)
        )
    bray_dist <- phyloseq::distance(MS, method = "bray")
    sample_df <- data.frame(sample_data(MS))
    permanova_cat <- vegan::adonis2(bray_dist ~ gc_treatment, data = sample_df)
    pval_cat <- signif(permanova_cat$`Pr(>F)`[1], 3) 
    print(paste("PERMANOVA p-value for gc_treatment in", kingdom, ":", pval_cat))  # Print the p-value    
    saveRDS(gc_treatment, gsub(" ", "", paste(output_folder, kingdom, "_beta_gc.rds")))
}

execute_beta <- function() {
    baselines_decB = readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
    baselines_decA = readRDS(file = "Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.001.rds")
    baselines_decE = readRDS(file = "Output/SUPERVISED_DEC/Eukaryote_Supervised_decontam0.001.rds")

    Beta(baselines_dec = baselines_decB, kingdom = "Bacteria", output_folder)
    Beta(baselines_dec = baselines_decA, kingdom = "Archaea", output_folder)
    Beta(baselines_dec = baselines_decE, kingdom = "Eukaryote", output_folder)
}

execute_beta()
