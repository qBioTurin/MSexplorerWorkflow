
maaslin_plot_results_from_output22<-function (output, metadata, normalization, transform, feature_specific_covariate = NULL, 
    feature_specific_covariate_name = NULL, feature_specific_covariate_record = NULL, 
    median_comparison_abundance = TRUE, median_comparison_prevalence = FALSE, 
    max_significance = 0.1, plot_summary_plot = TRUE, summary_plot_first_n = 25, 
    coef_plot_vars = NULL, heatmap_vars = NULL, plot_associations = TRUE, 
    max_pngs = 30, balanced = FALSE, save_plots_rds = FALSE) 
{
    ret_plots <- list()
    if (!file.exists(output)) {
        logging::loginfo("Creating output folder")
        dir.create(output)
    }
    if (plot_summary_plot || plot_associations) {
        figures_folder <- file.path(output, "figures")
        if (!file.exists(figures_folder)) {
            logging::loginfo("Creating output figures folder")
            dir.create(figures_folder)
        }
    }
    all_results_file <- paste0(gsub("/$", "", output), "/", "all_results.tsv")
    if (!file.exists(all_results_file)) {
        stop(sprintf("Please generate the results file first: %s", 
            all_results_file))
    }
    merged_results <- utils::read.csv(all_results_file, sep = "\t")
    merged_results$model[merged_results$model == "abundance"] <- "linear"
    merged_results$model[merged_results$model == "prevalence"] <- "logistic"
    if (plot_summary_plot) {
        summary_plot_file <- file.path(figures_folder, "summary_plot.pdf")
        logging::loginfo("Writing summary plot of\n                        significant results to file: %s", 
            summary_plot_file)
        if (!is.null(coef_plot_vars) & length(coef_plot_vars) == 
            1) {
            coef_plot_vars <- trimws(unlist(strsplit(coef_plot_vars, 
                ",")))
        }
        if (!is.null(heatmap_vars) & length(heatmap_vars) == 
            1) {
            heatmap_vars <- trimws(unlist(strsplit(heatmap_vars, 
                ",")))
        }
        summary_plot <- maaslin3_summary_plot(merged_results, 
            summary_plot_file, figures_folder, first_n = summary_plot_first_n, 
            max_significance = max_significance, coef_plot_vars = coef_plot_vars, 
            heatmap_vars = heatmap_vars, median_comparison_abundance = median_comparison_abundance, 
            median_comparison_prevalence = median_comparison_prevalence, 
            balanced = balanced, save_plots_rds = save_plots_rds)
        ret_plots[["summary_plot"]] <- summary_plot
    }
    if (plot_associations) {
        features_file <- paste0(gsub("/$", "", output), "/", 
            "features/data_transformed.tsv")
        if (!file.exists(features_file)) {
            stop(sprintf("Please generate the results file first: %s", 
                features_file))
        }
        transformed_data <- utils::read.csv(features_file, sep = "\t", 
            row.names = 1, check.names = FALSE)
        logging::loginfo(paste("Writing association plots", "(one for each significant association)", 
            "to output folder: %s"), figures_folder)
        if (!is.null(feature_specific_covariate)) {
            tryCatch({
                feature_specific_covariate <- feature_specific_covariate[rownames(transformed_data), 
                  colnames(transformed_data)]
            }, error = function(e) {
                stop("feature_specific_covariate does not contain the features\n                    and samples of the filtered data.")
            })
        }
        if (missing("normalization") | missing("transform") | 
            missing("metadata")) {
            stop("Missing normalization, transform, or metadata argument to\n                maaslin_plot_results_from_output")
        }
        plots_out <- tryCatch({
            withCallingHandlers({
                maaslin3_association_plots(merged_results = merged_results, 
                  metadata = metadata, features = transformed_data, 
                  max_significance = max_significance, figures_folder = figures_folder, 
                  max_pngs = max_pngs, normalization = normalization, 
                  transform = transform, feature_specific_covariate = feature_specific_covariate, 
                  feature_specific_covariate_name = feature_specific_covariate_name, 
                  feature_specific_covariate_record = feature_specific_covariate_record, 
                  save_plots_rds = save_plots_rds)
            }, warning = function(w) {
                invokeRestart("muffleWarning")
            })
        })
        ret_plots[["assocation_plots"]] <- plots_out
    }
    if ("logging::writeToFile" %in% names(logging::getLogger()[["handlers"]])) {
        logging::removeHandler("logging::writeToFile")
    }
    return(ret_plots)
}
all_results <- read.csv('hmp2_output/all_results.tsv', sep='\t')
all_results <- all_results %>%
    mutate(metadata = case_when(metadata == 'age' ~ 'Age',
                                metadata == 'antibiotics' ~ 'Abx',
                                metadata == 'diagnosis' ~ 'Diagnosis',
                                metadata == 'dysbiosis_state' ~ 'Dysbiosis',
                                metadata == 'reads' ~ 'Read depth'),
        value = case_when(value == 'dysbiosis_CD' ~ 'CD',
                            value == 'dysbiosis_UC' ~ 'UC',
                            value == 'Yes' ~ 'Used', # Antibiotics
                            value == 'age' ~ 'Age',
                            value == 'reads' ~ 'Read depth',
                            TRUE ~ value),
        feature = gsub('_', ' ', feature) %>%
            gsub(pattern = 'sp ', replacement = 'sp. '))

# Write results
write.table(all_results, 'hmp2_output/all_results.tsv', sep='\t')

# Set the new heatmap and coefficient plot variables and order them
heatmap_vars = c('Dysbiosis UC', 'Diagnosis UC',
                'Abx Used', 'Age', 'Read depth')
coef_plot_vars = c('Dysbiosis CD', 'Diagnosis CD')

# This section is necessary for updating the association plots
taxa_table_copy <- taxa_table
colnames(taxa_table_copy) <- gsub('_', ' ', colnames(taxa_table_copy)) %>%
    gsub(pattern = 'sp ', replacement = 'sp. ')

# Rename the features in the norm transformed data file
data_transformed <-
    read.csv('hmp2_output/features/data_transformed.tsv', sep='\t')
colnames(data_transformed) <-
    gsub('_', ' ', colnames(data_transformed)) %>%
    gsub(pattern = 'sp ', replacement = 'sp. ')
write.table(data_transformed,
            'hmp2_output/features/data_transformed.tsv',
            sep='\t', row.names = FALSE)

# Rename the metadata like in the outputs table
metadata_copy <- metadata
colnames(metadata_copy) <-
    case_when(colnames(metadata_copy) == 'age' ~ 'Age',
            colnames(metadata_copy) == 'antibiotics' ~ 'Abx',
            colnames(metadata_copy) == 'diagnosis' ~ 'Diagnosis',
            colnames(metadata_copy) == 'dysbiosis_state' ~ 'Dysbiosis',
            colnames(metadata_copy) == 'reads' ~ 'Read depth',
            TRUE ~ colnames(metadata_copy))
metadata_copy <- metadata_copy %>%
    mutate(Dysbiosis = case_when(Dysbiosis == 'dysbiosis_UC' ~ 'UC',
                                Dysbiosis == 'dysbiosis_CD' ~ 'CD',
                                Dysbiosis == 'none' ~ 'None') %>%
            factor(levels = c('None', 'UC', 'CD')),
        Abx = case_when(Abx == 'Yes' ~ 'Used',
                        Abx == 'No' ~ 'Not used') %>%
            factor(levels = c('Not used', 'Used')),
        Diagnosis = case_when(Diagnosis == 'nonIBD' ~ 'non-IBD',
                                TRUE ~ Diagnosis) %>%
            factor(levels = c('non-IBD', 'UC', 'CD')))

# Recreate the plots
scatter_plots <- maaslin_plot_results_from_output22(
    output = 'hmp2_output',
    metadata = metadata_copy,
    normalization = "TSS",
    transform = "LOG",
    median_comparison_abundance = TRUE,
    median_comparison_prevalence = FALSE,
    max_significance = 0.1,
    max_pngs = 20)