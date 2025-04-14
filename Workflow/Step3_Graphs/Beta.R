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
    beta1 = plot_ordination(
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
    beta1
    saveRDS(beta1, gsub(" ", "", paste(output_folder, kingdom, "_beta_cat.rds")))

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

    # Change legend label
    sample_data(MS)$gc_treatment <- factor(
        sample_data(MS)$gc_treatment, 
        levels = c("negative", "positive")
    )

    # Calculate Bray
    ordination <- ordinate(MS, method = "PCoA", distance = "bray")
    beta2 = plot_ordination(
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
        scale_color_manual(values = rev(colors_venn4)) +
        theme(
            axis.title.x = element_text(size = 15),  
            axis.title.y = element_text(size = 15),  
            axis.text.x = element_text(size = 15),   
            axis.text.y = element_text(size = 15),   
            plot.title = element_text(size = 20),    
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15)
        )
    saveRDS(beta1, gsub(" ", "", paste(output_folder, kingdom, "_beta_gc.rds")))
}

execute_beta <- function() {
    baselines_decB = readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
    baselines_decA = readRDS(file = "Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.001.rds")
    baselines_decE = readRDS(file = "Output/SUPERVISED_DEC/Eukaryota_Supervised_decontam0.001.rds")

    Beta(baselines_decB, "Bacteria", output_folder)
    Beta(baselines_decA, "Archaea", output_folder)
    Beta(baselines_decE, "Eukaryota", output_folder)    
}

execute_beta()