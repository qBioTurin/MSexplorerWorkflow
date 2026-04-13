library(phyloseq)
# install.packages("randomForest")
library(randomForest)
# install.packages("caret")
library(caret)
library(pROC)

# phyloseq is the phyloseq object containing the data
# output_folder is the folder where the results will be saved
# discriminant is the column used to discriminate the samples, e.g. "category"
# categorySel is the column use to select the patient to use and elements
# are the one witch to run the analysis es positive/negative for gc_treatment
#
RandomForest <- function(physeq, output_folder, discriminant, categorySel, elements) {
    if (!dir.exists(output_folder)) {
        dir.create(output_folder, recursive = TRUE)
    }

    # Sample metadata
    metadata <- as(sample_data(physeq), "data.frame")

    # Taxonomy
    taxa <- as.data.frame(tax_table(physeq))
    taxa$Genus_Species <- paste(taxa$Genus, taxa$Species, sep = "_")
    taxa$id <- rownames(taxa)
    taxa <- taxa %>% select(id, Genus_Species)

    # Filtra i campioni desiderati
    ids_to_retain <- rownames(
        metadata[
            metadata[[categorySel]] %in% elements &
                !is.na(metadata[[discriminant]]),
        ]
    )

    physeq <- prune_samples(ids_to_retain, physeq)

    # Trasformazione in abbondanze relative
    physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

    # OTU table normalizzata
    otu_df <- as.data.frame(otu_table(physeq_rel))
    otu_df$id <- rownames(otu_df)
    otu_df <- left_join(otu_df, taxa, by = "id")
    rownames(otu_df) <- otu_df$Genus_Species
    otu_df <- otu_df %>% select(-id, -Genus_Species)
    otu_df <- t(otu_df)

    # Metadata aggiornato
    metadata <- as(sample_data(physeq_rel), "data.frame")
    target_var <- metadata[[discriminant]]
    target_var <- as.factor(target_var)

    # Random Forest
    set.seed(123)
    rf_model <- randomForest(x = otu_df, y = target_var, importance = TRUE, ntree = 500)

    # Salva grafico importanza
    png(file.path(output_folder, "Top_20_Taxa_Predictive_category.png"), width = 1000, height = 800)
    varImpPlot(rf_model, n.var = 20, main = "Top 20 taxa (species) predictive")
    dev.off()

    # Tabella importanza
    importance_df <- as.data.frame(importance(rf_model))
    importance_df$TaxaID <- rownames(importance_df)
    importance_df <- importance_df[order(importance_df$MeanDecreaseGini, decreasing = TRUE), ]
    write.csv(importance_df, file = file.path(output_folder, "importance_df_category.csv"), row.names = FALSE)

    # ROC curve se binaria
    if (length(levels(target_var)) == 2) {
        prob <- predict(rf_model, type = "prob")[, 2]
        roc_obj <- roc(target_var, prob)
        png(file.path(output_folder, "ROC_curve_category.png"), width = 1000, height = 800)
        plot(roc_obj, main = paste("AUC =", round(auc(roc_obj), 2)))
        dev.off()
    }
}

# 📦 Carica i phyloseq decontaminati
physeq001 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
physeq01 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds")
physeq05 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds")

# 📋 Lista dei phyloseq con etichette
physeq_list <- list(
    "001" = physeq001,
    "01"  = physeq01,
    "05"  = physeq05
)

analysis_list <- list(
    list(name = "category", discrim = "category", filter_col = "category", filter_vals = c("MS", "HEALTHY")),
    list(name = "gc_treatment", discrim = "gc_treatment", filter_col = "category", filter_vals = c("MS")),
    list(name = "lesion_burden", discrim = "lesion_burden", filter_col = "gc_treatment", filter_vals = c("positive", "negative")),
    list(name = "spinal_cord_lesion", discrim = "spinal_cord_lesion", filter_col = "gc_treatment", filter_vals = c("positive", "negative")),
    list(name = "gadolinium_contrast", discrim = "gadolinium_contrast", filter_col = "gc_treatment", filter_vals = c("positive", "negative")),
    list(name = "subtentorial_lesions", discrim = "subtentorial_lesions", filter_col = "gc_treatment", filter_vals = c("positive", "negative"))
)

for (dec_level in names(physeq_list)) {
    phy <- physeq_list[[dec_level]]

    for (analysis in analysis_list) {
        output_dir <- file.path("Output", "RandomForest", paste0(analysis$name, "_", dec_level))

        cat("Running:", analysis$name, "on dec =", dec_level, "\n")

        tryCatch(
            {
                RandomForest(
                    physeq = phy,
                    output_folder = output_dir,
                    discriminant = analysis$discrim,
                    categorySel = analysis$filter_col,
                    elements = analysis$filter_vals
                )
            },
            error = function(e) {
                cat("Error in", analysis$name, "on", dec_level, ":", e$message, "\n")
            }
        )
    }
}
