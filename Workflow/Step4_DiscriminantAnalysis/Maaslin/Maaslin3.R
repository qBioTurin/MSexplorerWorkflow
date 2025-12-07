input_data <- merged_df
input_metadata <- metadata_df
output <- output_folder
formula <- formulaInUse
fixed_effects <- fx_effects
reference <- NULL
random_effects <- random_effects
group_effects <- group_effects
ordered_effects <- ordered_effects
strata_effects <- strata_effects
feature_specific_covariate <- NULL
feature_specific_covariate_name <- NULL
feature_specific_covariate_record <- NULL
min_abundance <- 0
min_prevalence <- 0
zero_threshold <- 0
min_variance <- 0
max_significance <- 0.1
normalization <- "TSS"
transform <- "LOG"
correction <- "BH"
standardize <- TRUE
unscaled_abundance <- NULL
median_comparison_abundance <- median_comparison_abundance_in
median_comparison_prevalence <- FALSE
median_comparison_abundance_threshold <- 0
median_comparison_prevalence_threshold <- 0
subtract_median <- FALSE
warn_prevalence <- TRUE
small_random_effects <- FALSE
augment <- TRUE
evaluate_only <- NULL
plot_summary_plot <- TRUE
summary_plot_first_n <- 50
coef_plot_vars <- coef_plot_vars
heatmap_vars <- NULL
plot_associations <- TRUE
max_pngs <- 30
cores <- 1
save_models <- FALSE
save_plots_rds <- FALSE
verbosity <- "FINEST"
summary_plot_balanced <- FALSE
assay.type <- 1
# --- FINE VARIABILI MAASLIN22 ---
library(maaslin3)
## cambio nomi per creazione plot
output <- output
transformed_data <- transformed_data
unstandardized_metadata <- metadata
fit_data_abundance <- maaslin_results$fit_data_abundance
fit_data_prevalence <- maaslin_results$fit_data_prevalence
normalization <- normalization
transform <- transform
feature_specific_covariate <- feature_specific_covariate
feature_specific_covariate_name <- feature_specific_covariate_name
feature_specific_covariate_record <- feature_specific_covariate_record
median_comparison_abundance <- median_comparison_abundance
median_comparison_prevalence <- median_comparison_prevalence
max_significance <- max_significance
plot_summary_plot <- plot_summary_plot
summary_plot_first_n <- summary_plot_first_n
coef_plot_vars <- coef_plot_vars
heatmap_vars <- heatmap_vars
plot_associations <- plot_associations
max_pngs <- max_pngs
balanced <- summary_plot_balanced
save_plots_rds <- save_plots_rds



maaslin22 <- function(
    input_data, input_metadata = NULL, output, formula = NULL,
    fixed_effects = NULL, reference = NULL, random_effects = NULL,
    group_effects = NULL, ordered_effects = NULL, strata_effects = NULL,
    feature_specific_covariate = NULL, feature_specific_covariate_name = NULL,
    feature_specific_covariate_record = NULL, min_abundance = 0,
    min_prevalence = 0, zero_threshold = 0, min_variance = 0,
    max_significance = 0.1, normalization = "TSS", transform = "LOG",
    correction = "BH", standardize = TRUE, unscaled_abundance = NULL,
    median_comparison_abundance = TRUE, median_comparison_prevalence = FALSE,
    median_comparison_abundance_threshold = 0, median_comparison_prevalence_threshold = 0,
    subtract_median = FALSE, warn_prevalence = TRUE, small_random_effects = FALSE,
    augment = TRUE, evaluate_only = NULL, plot_summary_plot = TRUE,
    summary_plot_first_n = 25, coef_plot_vars = NULL, heatmap_vars = NULL,
    plot_associations = TRUE, max_pngs = 30, cores = 1, save_models = FALSE,
    save_plots_rds = FALSE, verbosity = "FINEST", summary_plot_balanced = FALSE,
    assay.type = 1) {
    match.arg(verbosity, c(
        "FINEST", "FINER", "FINE", "DEBUG",
        "INFO", "WARN", "ERROR"
    ))
    print("Starting Maaslin3 analysis...")
    logging::logReset()
    normalization <- toupper(normalization)
    transform <- toupper(transform)
    correction_choices <- c("BH", "bonferroni", "holm", "hochberg", "hommel", "BY", "fdr", "none")
    correction <- correction_choices[match(
        toupper(correction),
        toupper(correction_choices)
    )]
    if (methods::is(formula, "formula")) {
        formula <- paste0(trimws(deparse(formula)), collapse = " ")
    }
    if (inherits(input_data, "SummarizedExperiment")) {
        summarized_experiment_out <- maaslin_read_summarized_experiment_data(
            input_data,
            assay.type
        )
        input_data <- summarized_experiment_out[["data"]]
        input_metadata <- summarized_experiment_out[["metadata"]]
    }
    print("Starting Maaslin3 analysis...")
    maaslin_log_arguments(
        input_data, input_metadata, output,
        formula, fixed_effects, reference, random_effects, group_effects,
        ordered_effects, strata_effects, feature_specific_covariate,
        feature_specific_covariate_name, feature_specific_covariate_record,
        min_abundance, min_prevalence, zero_threshold, min_variance,
        max_significance, normalization, transform, correction,
        standardize, unscaled_abundance, median_comparison_abundance,
        median_comparison_prevalence, median_comparison_abundance_threshold,
        median_comparison_prevalence_threshold, subtract_median,
        warn_prevalence, small_random_effects, augment, evaluate_only,
        plot_summary_plot, summary_plot_first_n, coef_plot_vars,
        heatmap_vars, plot_associations, max_pngs, cores, save_models,
        save_plots_rds, verbosity, summary_plot_balanced
    )
    read_data_list <- maaslin_read_data(
        input_data, input_metadata,
        feature_specific_covariate, unscaled_abundance
    )
    read_data_list <- maaslin_reorder_data(
        read_data_list$data,
        read_data_list$metadata, read_data_list$feature_specific_covariate,
        read_data_list$unscaled_abundance
    )
    data <- read_data_list$data # tabella di otu con pazienti per righe e specie per colonne
    metadata <- read_data_list$metadata
    unscaled_abundance <- read_data_list$unscaled_abundance
    feature_specific_covariate <- read_data_list$feature_specific_covariate
    if (is.null(formula)) {
        formulas <- maaslin_compute_formula(
            data, metadata, fixed_effects,
            random_effects, group_effects, ordered_effects, strata_effects,
            feature_specific_covariate_name
        )
    } else {
        if (!is.null(fixed_effects) | !is.null(random_effects) |
            !is.null(group_effects) | !is.null(ordered_effects) |
            !is.null(strata_effects)) {
            check_null_error <- paste0(
                "random_effects, group_effects, ",
                "ordered_effects, and strata_effects must be NULL ",
                "when formula is not NULL"
            )
            stop(check_null_error)
        }
        formulas <- maaslin_check_formula(
            data, metadata, formula,
            feature_specific_covariate_name
        )
    }
    formula <- formulas$formula
    random_effects_formula <- formulas$random_effects_formula
    normalized_data <- maaslin_normalize(
        data, output, zero_threshold,
        normalization, unscaled_abundance
    ) ## otu table normalizata
    print("Data normalization completed.")
    filtered_data <- maaslin_filter(
        normalized_data, output,
        min_abundance, min_prevalence, zero_threshold, min_variance
    ) ## otu filtrata
    transformed_data <- maaslin_transform(
        filtered_data, output,
        transform
    )
    standardized_metadata <- maaslin_process_metadata(
        metadata,
        formula, fixed_effects, reference, feature_specific_covariate_name,
        standardize
    )
    print("Metadata standardization completed.")
    maaslin_results <- maaslin3::maaslin_fit(filtered_data, transformed_data,
        standardized_metadata, formula, random_effects_formula,
        feature_specific_covariate, feature_specific_covariate_name,
        feature_specific_covariate_record, zero_threshold, max_significance,
        correction, median_comparison_abundance, median_comparison_prevalence,
        median_comparison_abundance_threshold, median_comparison_prevalence_threshold,
        subtract_median, warn_prevalence, small_random_effects,
        augment, evaluate_only, cores,
        save_models = TRUE, data,
        min_abundance, min_prevalence, min_variance
    )
    maaslin_write_results(
        output, maaslin_results$fit_data_abundance,
        maaslin_results$fit_data_prevalence, random_effects_formula,
        max_significance, save_models
    )
    print("Maaslin3 plotting started.")
    if (plot_summary_plot | plot_associations) {
        tryCatch({
            withCallingHandlers(
                {
                    maaslin_plot_results22(
                        output = output,
                        transformed_data = transformed_data,
                        unstandardized_metadata = unstandardized_metadata,
                        fit_data_abundance = fit_data_abundance,
                        fit_data_prevalence = fit_data_prevalence,
                        normalization = normalization,
                        transform = transform,
                        feature_specific_covariate = feature_specific_covariate,
                        feature_specific_covariate_name = feature_specific_covariate_name,
                        feature_specific_covariate_record = feature_specific_covariate_record,
                        median_comparison_abundance = median_comparison_abundance,
                        median_comparison_prevalence = median_comparison_prevalence,
                        max_significance = max_significance,
                        plot_summary_plot = plot_summary_plot,
                        summary_plot_first_n = summary_plot_first_n,
                        coef_plot_vars = coef_plot_vars,
                        heatmap_vars = heatmap_vars,
                        plot_associations = plot_associations,
                        max_pngs = max_pngs,
                        balanced = balanced,
                        save_plots_rds = save_plots_rds
                    )
                },
                warning = function(w) {
                    invokeRestart("muffleWarning")
                }
            )
        })
    }
    if ("logging::writeToFile" %in% names(logging::getLogger()[["handlers"]])) {
        logging::removeHandler("logging::writeToFile")
    }
    return(list(
        data = data, normalized_data = normalized_data,
        filtered_data = filtered_data, transformed_data = transformed_data,
        metadata = metadata, standardized_metadata = standardized_metadata,
        formula = formulas, fit_data_abundance = maaslin_results$fit_data_abundance,
        fit_data_prevalence = maaslin_results$fit_data_prevalence
    ))
}

maaslin_plot_results22 <- function(
    output, transformed_data, unstandardized_metadata,
    fit_data_abundance, fit_data_prevalence, normalization, transform,
    feature_specific_covariate = NULL, feature_specific_covariate_name = NULL,
    feature_specific_covariate_record = NULL, median_comparison_abundance = TRUE,
    median_comparison_prevalence = FALSE, max_significance = 0.1,
    plot_summary_plot = TRUE, summary_plot_first_n = 25, coef_plot_vars = NULL,
    heatmap_vars = NULL, plot_associations = TRUE, max_pngs = 30,
    balanced = FALSE, save_plots_rds = FALSE) {
    print("Starting Maaslin3 plotting...")
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
    print("Merging results for plotting...")
    if (is.null(fit_data_abundance$results)) {
        merged_results <- fit_data_prevalence$results
        print(merged_results)
    } else if (is.null(fit_data_prevalence$results)) {
        merged_results <- fit_data_abundance$results
    } else {
        merged_results <- rbind(fit_data_abundance$results, fit_data_prevalence$results)
    }
    if (plot_summary_plot) {
        summary_plot_file <- file.path(figures_folder, "summary_plot.pdf")
        logging::loginfo(
            "Writing summary plot of significant\n                        results to file: %s",
            summary_plot_file
        )
        if (!is.null(coef_plot_vars) & length(coef_plot_vars) ==
            1) {
            coef_plot_vars <- trimws(unlist(strsplit(
                coef_plot_vars,
                ","
            )))
        }
        if (!is.null(heatmap_vars) & length(heatmap_vars) ==
            1) {
            heatmap_vars <- trimws(unlist(strsplit(
                heatmap_vars,
                ","
            )))
        }
        print(paste("Generating summary plot...", merged_results))
        summary_plot <- maaslin3_summary_plot22(merged_results,
            summary_plot_file, figures_folder,
            first_n = summary_plot_first_n,
            max_significance = max_significance, coef_plot_vars = coef_plot_vars,
            heatmap_vars = heatmap_vars, median_comparison_abundance = median_comparison_abundance,
            median_comparison_prevalence = median_comparison_prevalence,
            balanced = balanced, save_plots_rds = save_plots_rds
        )
        ret_plots[["summary_plot"]] <- summary_plot
    }
    if (plot_associations) {
        logging::loginfo(paste(
            "Writing association plots", "(one for each significant association)",
            "to output folder: %s"
        ), figures_folder)
        plots_out <- tryCatch({
            withCallingHandlers(
                {
                    maaslin3_association_plots(
                        merged_results = merged_results,
                        metadata = unstandardized_metadata, features = transformed_data,
                        max_significance = max_significance, figures_folder = figures_folder,
                        max_pngs = max_pngs, normalization = normalization,
                        transform = transform, feature_specific_covariate = feature_specific_covariate,
                        feature_specific_covariate_name = feature_specific_covariate_name,
                        feature_specific_covariate_record = feature_specific_covariate_record,
                        save_plots_rds = save_plots_rds
                    )
                },
                warning = function(w) {
                    invokeRestart("muffleWarning")
                }
            )
        })
        print("Generating association plots...")
        ret_plots[["assocation_plots"]] <- plots_out
    }
    return(ret_plots)
}

maaslin3_summary_plot22 <- function(
    merged_results, summary_plot_file, figures_folder,
    first_n = 30, max_significance = 0.1, coef_plot_vars = NULL,
    heatmap_vars = NULL, median_comparison_abundance = FALSE,
    median_comparison_prevalence = FALSE, balanced = FALSE, save_plots_rds = FALSE) {
    ret_plots <- list()
    if (first_n > 200) {
        logging::logerror(paste("At most 200 features can be plotted in the heatmap. \n                    Please choose a smaller first_n."))
        return()
    }
    print("Preprocessing merged results...")
    merged_results <- preprocess_merged_results22(merged_results)
    if (is.null(merged_results) || nrow(merged_results) == 0) {
        return()
    }
    median_df <- merged_results %>%
        dplyr::group_by(
            .data$full_metadata_name,
            .data$model
        ) %>%
        dplyr::summarize(
            median_val = unique(.data$null_hypothesis),
            .groups = "drop"
        )
    if (!is.null(coef_plot_vars) | !is.null(heatmap_vars)) {
        if (any(!c(coef_plot_vars, heatmap_vars) %in% unique(merged_results$full_metadata_name))) {
            logging::loginfo(paste0("The following specified variables were not \n                        found in the associations: ",
                paste0(
                    setdiff(
                        c(coef_plot_vars, heatmap_vars),
                        unique(merged_results$full_metadata_name)
                    ),
                    collapse = ", "
                ),
                collapse = ""
            ))
            logging::loginfo(paste0("Available associations: ",
                paste0(unique(merged_results$full_metadata_name),
                    collapse = ", "
                ),
                collapse = ""
            ))
            print("Filtering merged results based on specified variables...")
            return(NULL)
        }
        merged_results <- merged_results[merged_results$full_metadata_name %in%
            c(coef_plot_vars, heatmap_vars), ]
    }
    merged_results_joint_only <- unique(merged_results[, c(
        "feature",
        "pval_individual", "full_metadata_name"
    )])
    merged_results_joint_only <- merged_results_joint_only[order(merged_results_joint_only$pval_individual), ]
    if (length(unique(merged_results_joint_only$feature)) < first_n) {
        first_n <- length(unique(merged_results_joint_only$feature))
        signif_taxa <- unique(merged_results_joint_only$feature)[seq(first_n)]
    } else {
        if (balanced) {
            if (is.null(coef_plot_vars)) {
                logging::logerror(paste("Balanced plotting requires \n                            you set the variables you \n                            want to plot using\n                            the parameter coef_plot_vars"))
                return()
            } else {
                first_n_per <- first_n / length(coef_plot_vars)
                signif_taxa <- merged_results_joint_only %>%
                    dplyr::group_by(.data$full_metadata_name) %>%
                    dplyr::arrange(dplyr::desc(-.data$pval_individual),
                        .by_group = TRUE
                    ) %>%
                    dplyr::slice_head(n = ceiling(first_n_per)) %>%
                    dplyr::pull(.data$feature) %>%
                    unique()
            }
        } else {
            signif_taxa <- unique(merged_results_joint_only$feature)[seq(first_n)]
        }
    }
    merged_results_sig <- merged_results %>% dplyr::filter(.data$feature %in%
        signif_taxa)
    ord_feature <- with(merged_results_sig, reorder(
        feature,
        pval_individual
    ))
    ord_feature <- levels(ord_feature)
    merged_results_sig$feature <- factor(merged_results_sig$feature,
        levels = ord_feature
    )
    if (is.null(coef_plot_vars)) {
        mean_log_qval <- merged_results_sig %>%
            dplyr::group_by(.data$full_metadata_name) %>%
            dplyr::summarise(mean_value = mean(log(.data$pval_individual),
                na.rm = TRUE
            ))
        coef_plot_vars <- mean_log_qval$full_metadata_name[order(mean_log_qval$mean_value)]
        coef_plot_vars <- setdiff(coef_plot_vars, heatmap_vars)
        if (length(coef_plot_vars) > 0) {
            coef_plot_vars <- coef_plot_vars[seq(min(2, length(coef_plot_vars)))]
        }
    }
    if (is.null(heatmap_vars)) {
        mean_log_qval <- merged_results_sig %>%
            dplyr::group_by(.data$full_metadata_name) %>%
            dplyr::summarise(mean_value = mean(log(.data$pval_individual),
                na.rm = TRUE
            ))
        heatmap_vars <- mean_log_qval$full_metadata_name[order(mean_log_qval$mean_value)]
        heatmap_vars <- setdiff(heatmap_vars, coef_plot_vars)
    }
    if (length(coef_plot_vars) > 0 & sum(merged_results_sig$full_metadata_name %in%
        coef_plot_vars) >= 1) {
        if (balanced) {
            plot_thres <- 5
        } else {
            plot_thres <- 10
        }
        p1 <- make_coef_plot(merged_results_sig, coef_plot_vars,
            max_significance, median_comparison_prevalence, median_comparison_abundance,
            median_df,
            plot_threshold = plot_thres
        )
    } else {
        p1 <- NULL
    }
    if (length(heatmap_vars) > 0 & sum(merged_results_sig$full_metadata_name %in%
        heatmap_vars) >= 1) {
        p2 <- make_heatmap_plot(
            merged_results_sig, heatmap_vars,
            max_significance, median_comparison_prevalence, median_comparison_abundance
        )
        if (!is.null(p1)) {
            p2 <- p2 + ggplot2::theme(
                axis.text.y = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
            )
        }
    } else {
        p2 <- NULL
    }
    if (!is.null(p1) & !is.null(p2)) {
        final_plot <- patchwork::wrap_plots(p1, p2,
            ncol = 3,
            widths = c(max(0, length(coef_plot_vars) * (max(
                15,
                max(nchar(as.character(coef_plot_vars)))
            )) / 15 -
                2) + 2, max(0, length(heatmap_vars) / 4 - 2) +
                2, 0.5), guides = "collect"
        )
    } else if (is.null(p1) & !is.null(p2)) {
        final_plot <- p2
    } else if (!is.null(p1) & is.null(p2)) {
        final_plot <- p1
    } else {
        final_plot <- NULL
    }
    ret_plots[["coefficient"]] <- p1
    ret_plots[["heat"]] <- p2
    ret_plots[["final"]] <- final_plot
    if (!is.null(final_plot)) {
        height_out <- 9.5 + max(first_n / 5 - 5, 0) + max(nchar(c(
            as.character(coef_plot_vars),
            as.character(heatmap_vars)
        ))) / 10
        width_out <- 5 + max(nchar(merged_results$feature)) / 12 +
            (length(coef_plot_vars) * (max(20, max(nchar(as.character(coef_plot_vars))))) / 20) *
                2.5 + length(heatmap_vars) * 0.25
        tryCatch({
            withCallingHandlers(
                {
                    ggplot2::ggsave(summary_plot_file,
                        plot = final_plot,
                        height = height_out, width = width_out
                    )
                },
                warning = function(w) {
                    invokeRestart("muffleWarning")
                }
            )
        })
        png_file <- file.path(figures_folder, "summary_plot.png")
        tryCatch({
            withCallingHandlers(
                {
                    ggplot2::ggsave(png_file,
                        plot = final_plot,
                        height = height_out, width = width_out
                    )
                },
                warning = function(w) {
                    invokeRestart("muffleWarning")
                }
            )
        })
        if (save_plots_rds) {
            saveRDS(final_plot, file = paste(figures_folder,
                "/", "summary_plot_gg.RDS",
                sep = ""
            ))
        }
    }
    return(ret_plots)
}

preprocess_merged_results22 <- function(merged_results) {
    merged_results$error[grepl(
        "Prevalence association possibly induced",
        merged_results$error
    )] <- NA
    merged_results$error[grepl(
        "<4 average observations per random effect",
        merged_results$error
    )] <- NA
    merged_results <- merged_results[is.na(merged_results$error) &
        !is.na(merged_results$qval_individual) & !is.na(merged_results$coef), ]
    if (nrow(merged_results) == 0) {
        logging::loginfo(paste("No associations were without errors. \n                No summary plot generated."))
        return(NULL)
    }
    merged_results$model <- ifelse(merged_results$model == "linear",
        "Abundance", "Prevalence"
    )
    merged_results$full_metadata_name <- ifelse(merged_results$metadata ==
        merged_results$value, merged_results$metadata, paste0(
        merged_results$metadata,
        " ", merged_results$value
    ))
    return(merged_results)
}
make_coef_plot <- function(
    merged_results_sig, coef_plot_vars, max_significance,
    median_comparison_prevalence, median_comparison_abundance,
    median_df, plot_threshold = 10) {
    coef_plot_data <- merged_results_sig[merged_results_sig$full_metadata_name %in%
        coef_plot_vars, ]
    quantile_df <- coef_plot_data %>%
        dplyr::group_by(.data$full_metadata_name) %>%
        dplyr::summarise(
            lower_q = median(.data$coef) - plot_threshold *
                (median(.data$coef) - quantile(.data$coef, 0.25)),
            upper_q = median(.data$coef) + plot_threshold * (quantile(
                .data$coef,
                0.75
            ) - median(.data$coef))
        ) %>%
        data.frame(check.names = FALSE)
    rownames(quantile_df) <- quantile_df$full_metadata_name
    coef_plot_data <- coef_plot_data[coef_plot_data$qval_individual <
        max_significance | (coef_plot_data$coef > quantile_df[
        coef_plot_data$full_metadata_name,
        "lower_q"
    ] & coef_plot_data$coef < quantile_df[
        coef_plot_data$full_metadata_name,
        "upper_q"
    ]), ]
    custom_break_fun <- function(n) {
        return(function(x) {
            extended_breaks <- (scales::breaks_extended(n))(x)
            if (max(x) > 0) {
                extended_breaks <- extended_breaks[extended_breaks <=
                    max(x) * 0.9]
            } else {
                extended_breaks <- extended_breaks[extended_breaks <=
                    max(x) * 1.1]
            }
            if (min(x) > 0) {
                extended_breaks <- extended_breaks[extended_breaks >=
                    min(x) * 1.1]
            } else {
                extended_breaks <- extended_breaks[extended_breaks >=
                    min(x) * 0.9]
            }
            extended_breaks
        })
    }
    p1 <- ggplot2::ggplot(coef_plot_data, ggplot2::aes(
        x = .data$coef,
        y = .data$feature
    ))
    if (median_comparison_prevalence | median_comparison_abundance) {
        p1 <- p1 + ggplot2::guides(linetype = ggplot2::guide_legend(
            title = "Null hypothesis",
            order = 1
        ), ) + ggplot2::geom_vline(data = median_df[median_df$full_metadata_name %in%
            coef_plot_vars, ], ggplot2::aes(
            xintercept = .data$median_val,
            linetype = .data$model
        ), color = "darkgray") + ggplot2::scale_linetype_manual(values = c(
            Prevalence = "dashed",
            Abundance = "solid"
        ))
    } else {
        p1 <- p1 + ggplot2::geom_vline(ggplot2::aes(xintercept = 0),
            color = "darkgray", linetype = "dashed"
        )
    }
    scale_fill_gradient_limits <- c(
        min(max_significance, 10^floor(log10(min(coef_plot_data$qval_individual)))),
        1
    )
    if (min(coef_plot_data$qval_individual) < max_significance) {
        scale_fill_gradient_breaks <- c(
            10^floor(log10(min(coef_plot_data$qval_individual))),
            max_significance, 1
        )
    } else {
        scale_fill_gradient_breaks <- c(max_significance, 1)
    }
    if (min(coef_plot_data$qval_individual) < max_significance) {
        scale_fill_gradient_labels <- c(
            paste0("1e", floor(log10(min(coef_plot_data$qval_individual)))),
            paste0(max_significance), "1"
        )
    } else {
        scale_fill_gradient_labels <- c(
            paste0(max_significance),
            "1"
        )
    }
    p1 <- p1 + ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$coef -
        .data$stderr, xmax = .data$coef + .data$stderr), width = 0.2) +
        ggplot2::geom_point(data = coef_plot_data[coef_plot_data$model ==
            "Prevalence", ], ggplot2::aes(
            shape = .data$model,
            fill = .data$pval_individual
        ), size = 4.5, color = "black") +
        ggplot2::scale_fill_gradient(
            low = "#008B8B", high = "white",
            limits = scale_fill_gradient_limits, breaks = scale_fill_gradient_breaks,
            labels = scale_fill_gradient_labels, transform = scales::pseudo_log_trans(sigma = 0.001),
            name = bquote("Prevalence" ~ P["value"])
        ) + ggnewscale::new_scale_fill() +
        ggplot2::geom_point(data = coef_plot_data[coef_plot_data$model ==
            "Abundance", ], ggplot2::aes(
            shape = .data$model,
            fill = .data$pval_individual
        ), size = 4.5, color = "black") +
        ggplot2::scale_fill_gradient(
            low = "#8B008B", high = "white",
            limits = scale_fill_gradient_limits, breaks = scale_fill_gradient_breaks,
            labels = scale_fill_gradient_labels, transform = scales::pseudo_log_trans(sigma = 0.001),
            name = bquote("Abundance" ~ P["value"])
        ) + ggplot2::scale_x_continuous(
            breaks = custom_break_fun(n = 6),
            limits = c(min(coef_plot_data$coef) - quantile(
                coef_plot_data$stderr,
                0.8
            ), max(coef_plot_data$coef) + quantile(
                coef_plot_data$stderr,
                0.8
            ))
        ) + ggplot2::scale_shape_manual(
            name = "Association",
            values = c(21, 24)
        ) + ggplot2::guides(shape = ggplot2::guide_legend(order = 2), ) + ggplot2::labs(
            x = expression(paste(beta, " coefficient")),
            y = "Feature"
        ) + ggplot2::theme_bw() + ggplot2::theme(
            axis.title = ggplot2::element_text(size = 16),
            axis.text.y = ggplot2::element_text(size = 14), axis.text.x = ggplot2::element_text(size = 14),
            legend.title = ggplot2::element_text(size = 16), legend.text = ggplot2::element_text(
                size = 14,
                face = "plain"
            ), legend.position = "right", legend.background = ggplot2::element_rect(fill = "transparent"),
            panel.spacing = ggplot2::unit(0, "lines"), panel.grid.minor = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(size = 14), strip.background = ggplot2::element_rect(fill = "transparent")
        ) +
        ggplot2::facet_wrap(~ factor(full_metadata_name, levels = unique(coef_plot_vars)),
            scales = "free_x", ncol = length(coef_plot_vars)
        )
    return(p1)
}


maaslin3_association_plots <- function(
    merged_results, metadata, features, max_significance = 0.1,
    figures_folder, max_pngs = 10, normalization, transform,
    feature_specific_covariate = NULL, feature_specific_covariate_name = NULL,
    feature_specific_covariate_record = NULL, save_plots_rds = FALSE) {
    match.arg(normalization, c("TSS", "CLR", "NONE"))
    match.arg(transform, c("LOG", "PLOG", "NONE"))
    merged_results$error[grepl(
        "Prevalence association possibly induced",
        merged_results$error
    )] <- NA
    merged_results$error[grepl(
        "<4 average observations per random effect",
        merged_results$error
    )] <- NA
    new_name_normalization <- c(
        "Total sum scaling", "Center log ratio",
        "None"
    )
    names(new_name_normalization) <- c("TSS", "CLR", "NONE")
    normalization <- new_name_normalization[normalization]
    new_name_transformation <- c(
        "Log base 2", "Pseudo-log base 2",
        "None"
    )
    names(new_name_transformation) <- c("LOG", "PLOG", "NONE")
    transformation <- new_name_transformation[transform]
    merged_results <- merged_results[is.na(merged_results$error) &
        !is.na(merged_results$qval_individual) & merged_results$qval_individual <
        max_significance, ]
    if (nrow(merged_results) == 0) {
        logging::loginfo(paste("All associations had errors \n                                or were insignificant."))
        return(NULL)
    }
    merged_results <- merged_results[order(merged_results$qval_individual), ]
    logging::loginfo(paste(
        "Plotting associations from most",
        "to least significant,", "grouped by metadata"
    ))
    saved_plots <- list()
    features_by_metadata <- unique(merged_results[, c(
        "feature",
        "metadata", "model"
    )])
    for (row_num in seq(min(nrow(features_by_metadata), max_pngs))) {
        feature_name <- features_by_metadata[row_num, "feature"]
        feature_abun <- data.frame(
            sample = rownames(features),
            feature_abun = features[, feature_name], check.names = FALSE
        )
        metadata_name <- features_by_metadata[row_num, "metadata"]
        if (!is.null(feature_specific_covariate_name)) {
            if (metadata_name == feature_specific_covariate_name) {
                metadata_sub <- data.frame(
                    sample = rownames(feature_specific_covariate),
                    metadata = feature_specific_covariate[, feature_name],
                    check.names = FALSE
                )
            } else {
                metadata_sub <- data.frame(
                    sample = rownames(metadata),
                    metadata = metadata[, metadata_name], check.names = FALSE
                )
            }
        } else {
            metadata_sub <- data.frame(
                sample = rownames(metadata),
                metadata = metadata[, metadata_name], check.names = FALSE
            )
        }
        joined_features_metadata <- dplyr::inner_join(feature_abun,
            metadata_sub,
            by = c("sample")
        )
        model_name <- features_by_metadata[row_num, "model"]
        this_signif_association <- merged_results[merged_results$feature ==
            feature_name & merged_results$metadata == metadata_name &
            merged_results$model == model_name, ]
        if ("linear" == model_name) {
            temp_plot <- make_lm_plot(
                this_signif_association,
                joined_features_metadata, metadata, metadata_name,
                feature_name, normalization, transformation,
                feature_specific_covariate_name, feature_specific_covariate
            )
        }
        if ("logistic" == model_name) {
            temp_plot <- make_logistic_plot(
                this_signif_association,
                joined_features_metadata, metadata, metadata_name,
                feature_name, normalization, transformation,
                feature_specific_covariate_name, feature_specific_covariate
            )
        }
        saved_plots[[metadata_name]][[feature_name]][[model_name]] <- temp_plot
    }
    association_plots_folder <- file.path(figures_folder, "association_plots")
    if (!file.exists(association_plots_folder)) {
        dir.create(association_plots_folder)
    }
    vapply(names(saved_plots), function(metadata_variable) {
        if (save_plots_rds) {
            saveRDS(saved_plots[[metadata_variable]], file = file.path(
                association_plots_folder,
                paste0(make.names(metadata_variable), "_gg_associations.RDS")
            ))
        }
        vapply(names(saved_plots[[metadata_variable]]), function(feature) {
            vapply(
                names(saved_plots[[metadata_variable]][[feature]]),
                function(model_name) {
                    this_plot <- saved_plots[[metadata_variable]][[feature]][[model_name]]
                    association_plots_sub_folder <- file.path(
                        association_plots_folder,
                        make.names(metadata_variable), model_name
                    )
                    if (!file.exists(association_plots_sub_folder)) {
                        dir.create(association_plots_sub_folder,
                            recursive = TRUE
                        )
                    }
                    png_file <- file.path(
                        association_plots_sub_folder,
                        paste0(
                            make.names(metadata_variable), "_",
                            make.names(feature), "_", model_name, ".png"
                        )
                    )
                    height <- max(960, 18 * max(nchar(unlist(strsplit(
                        this_plot$labels$y,
                        "\n"
                    )))))
                    tryCatch({
                        withCallingHandlers(
                            {
                                ggplot2::ggsave(
                                    filename = png_file, plot = this_plot,
                                    dpi = 600, width = 960 / 300, height = height / 300
                                )
                            },
                            warning = function(w) {
                                invokeRestart("muffleWarning")
                            }
                        )
                    })
                    return(0)
                }, numeric(1)
            )
            return(0)
        }, numeric(1))
        return(0)
    }, numeric(1))
    return(saved_plots)
}

# BiocManager::install("biobakery/maaslin3")
library(maaslin3)
# Load an RDS file in R



Maaslin3 <- function(
    rds, output_folder, main_el, categorySel, elements, discr_els, perc = "K", median_comparison_abundance_in = TRUE,
    fx_effects = NULL, random_effects = NULL, group_effects = NULL,
    ordered_effects = NULL, strata_effects = NULL) {
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
        select(category, gc_treatment, lesion_burden, spinal_cord_lesion, gadolinium_contrast, subtentorial_lesions, sex, age, bmi)
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
    fit_out <- maaslin3(
        input_data = merged_df,
        input_metadata = metadata_df,
        output = output_folder,
        formula = formulaInUse,
        fixed_effects = fx_effects,
        plot_summary_plot = TRUE,
        coef_plot_vars = coef_plot_vars,
        # heatmap_vars =  discr_els_formatted,
        random_effects = random_effects,
        group_effects = group_effects,
        ordered_effects = ordered_effects,
        strata_effects = strata_effects,
        normalization = "TSS",
        transform = "LOG",
        augment = TRUE,
        standardize = TRUE,
        max_significance = 0.1,
        median_comparison_abundance = TRUE,
        median_comparison_prevalence = FALSE,
        max_pngs = 1,
        summary_plot_first_n = 50,
        cores = 1
    )
    all_results <- read_tsv(file.path(output_folder, "all_results.tsv"))
    all_results_sig <- all_results %>% dplyr::filter(metadata == main_el & pval_individual < 0.05)
    all_results_sig <- all_results_sig %>%
        select(feature, metadata, value, coef, null_hypothesis, pval_individual, qval_individual, pval_joint, pval_individual, model, stderr)
    write.csv(all_results_sig, file.path(output_folder, paste0(main_el, "_MAASLIN3_SR_", perc, ".csv")), row.names = FALSE)
    features_unique <- unique(all_results_sig$feature)
    features_unique <- data.frame(Species = features_unique)
    write.csv(features_unique, file.path(output_folder, paste0(main_el, "_MAASLIN3_SR_features_", perc, ".csv")), row.names = FALSE, quote=FALSE)
}
rds_001 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
library(dplyr)
Maaslin3(
    rds = rds_001, output_folder = "Output/MAASLIN3/001_disc/001_lesion_burden_onlygc_plus_gc",
    main_el = "lesion_burden", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "001"
)

Maaslin3(
    rds = rds_001, output_folder = "Output/MAASLIN3/001_disc/001_spinal_cord_lesion_onlygc_plus_gc",
    main_el = "spinal_cord_lesion", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "001"
)

Maaslin3(
    rds = rds_001, output_folder = "Output/MAASLIN3/001_disc/001_gadolinium_contrast_onlygc_plus_gc",
    main_el = "gadolinium_contrast", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "001"
)

Maaslin3(
    rds = rds_001, output_folder = "Output/MAASLIN3/001_disc/001_subtentorial_lesions_onlygc_plus_gc",
    main_el = "subtentorial_lesions", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "001"
)







rds_01 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds")

Maaslin3(
    rds = rds_01, output_folder = "Output/MAASLIN3/01_disc/01_lesion_burden_onlygc_T",
    main_el = "lesion_burden", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "01"
)

Maaslin3(
    rds = rds_01, output_folder = "Output/MAASLIN3/01_disc/01_spinal_cord_lesion_onlygc_T",
    main_el = "spinal_cord_lesion", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "01"
)

Maaslin3(
    rds = rds_01, output_folder = "Output/MAASLIN3/01_disc/01_gadolinium_contrast_onlygc_T",
    main_el = "gadolinium_contrast", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "01"
)

Maaslin3(
    rds = rds_01, output_folder = "Output/MAASLIN3/01_disc/01_subtentorial_lesions_onlygc_T",
    main_el = "subtentorial_lesions", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "01"
)

rds_05 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds")

Maaslin3(
    rds = rds_05, output_folder = "Output/MAASLIN3/05_disc/05_lesion_burden_onlygc_T",
    main_el = "lesion_burden", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "05"
)

Maaslin3(
    rds = rds_05, output_folder = "Output/MAASLIN3/05_disc/05_spinal_cord_lesion_onlygc_T",
    main_el = "spinal_cord_lesion", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "05"
)

Maaslin3(
    rds = rds_05, output_folder = "Output/MAASLIN3/05_disc/05_gadolinium_contrast_onlygc_T",
    main_el = "gadolinium_contrast", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "05"
)

Maaslin3(
    rds = rds_05, output_folder = "Output/MAASLIN3/05_disc/05_subtentorial_lesions_onlygc_T",
    main_el = "subtentorial_lesions", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "05"
)


rds_Euk <- readRDS("Output/SUPERVISED_DEC/Eukaryote_Supervised_decontam0.001.rds")
Maaslin3(
    rds = rds_Euk, output_folder = "Output/MAASLIN3/Eukaryote_MSHD_disc_T",
    main_el = "category", categorySel = "category", elements = c("HEALTHY", "MS"), discr_els = c("sex", "age", "bmi")
)
rds_Arch <- readRDS("Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.001.rds")
Maaslin3(
    rds = rds_Arch, output_folder = "Output/MAASLIN3/Archaea_MSHD_disc_T",
    main_el = "category", categorySel = "category", elements = c("HEALTHY", "MS"), discr_els = c("sex", "age", "bmi")
)
Maaslin3(
    rds = rds_001, output_folder = "Output/MAASLIN3/Bacteria_MSHD_disc_T",
    main_el = "category", categorySel = "category", elements = c("HEALTHY", "MS"), discr_els = c("sex", "age", "bmi")
)


Maaslin3(
    rds = rds_Euk, output_folder = "Output/MAASLIN3/Eukaryote_GC_TREATMENT_disc_T",
    main_el = "gc_treatment", categorySel = "category", elements = c("MS"), discr_els = c("sex", "age", "bmi")
)
Maaslin3(
    rds = rds_Arch, output_folder = "Output/MAASLIN3/Archaea_GC_TREATMENT_disc_T",
    main_el = "gc_treatment", categorySel = "category", elements = c("MS"), discr_els = c("sex", "age", "bmi")
)
Maaslin3(
    rds = rds_001, output_folder = "Output/MAASLIN3/Bacteria_GC_TREATMENT_disc_T",
    main_el = "gc_treatment", categorySel = "category", elements = c("MS"), discr_els = c("sex", "age", "bmi")
)





########################### false


# rds_001 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
# library(dplyr)
# Maaslin3(
#     rds = rds_001, output_folder = "Output/MAASLIN3/001_disc/001_lesion_burden_onlygc_F",
#     main_el = "lesion_burden", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
#     perc = "001", median_comparison_abundance_in = FALSE
# )

# Maaslin3(
#     rds = rds_001, output_folder = "Output/MAASLIN3/001_disc/001_spinal_cord_lesion_onlygc_F",
#     main_el = "spinal_cord_lesion", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
#     perc = "001", median_comparison_abundance_in = FALSE
# )

# Maaslin3(
#     rds = rds_001, output_folder = "Output/MAASLIN3/001_disc/001_gadolinium_contrast_onlygc_F",
#     main_el = "gadolinium_contrast", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
#     perc = "001", median_comparison_abundance_in = FALSE
# )
# Maaslin3(
#     rds = rds_001, output_folder = "Output/MAASLIN3/001_disc/001_subtentorial_lesions_onlygc_F",
#     main_el = "subtentorial_lesions", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
#     perc = "001", median_comparison_abundance_in = FALSE
# )
# rds_01 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds")

# Maaslin3(
#     rds = rds_01, output_folder = "Output/MAASLIN3/01_disc/01_lesion_burden_onlygc_F",
#     main_el = "lesion_burden", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
#     perc = "01", median_comparison_abundance_in = FALSE
# )

# Maaslin3(
#     rds = rds_01, output_folder = "Output/MAASLIN3/01_disc/01_spinal_cord_lesion_onlygc_F",
#     main_el = "spinal_cord_lesion", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
#     perc = "01", median_comparison_abundance_in = FALSE
# )
# Maaslin3(
#     rds = rds_01, output_folder = "Output/MAASLIN3/01_disc/01_gadolinium_contrast_onlygc_F",
#     main_el = "gadolinium_contrast", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
#     perc = "01", median_comparison_abundance_in = FALSE
# )
# Maaslin3(
#     rds = rds_01, output_folder = "Output/MAASLIN3/01_disc/01_subtentorial_lesions_onlygc_F",
#     main_el = "subtentorial_lesions", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
#     perc = "01", median_comparison_abundance_in = FALSE
# )
# rds_05 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds")

# Maaslin3(
#     rds = rds_05, output_folder = "Output/MAASLIN3/05_disc/05_lesion_burden_onlygc_F",
#     main_el = "lesion_burden", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
#     perc = "05", median_comparison_abundance_in = FALSE
# )

# Maaslin3(
#     rds = rds_05, output_folder = "Output/MAASLIN3/05_disc/05_spinal_cord_lesion_onlygc_F",
#     main_el = "spinal_cord_lesion", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
#     perc = "05", median_comparison_abundance_in = FALSE
# )
# Maaslin3(
#     rds = rds_05, output_folder = "Output/MAASLIN3/05_disc/05_gadolinium_contrast_onlygc_F",
#     main_el = "gadolinium_contrast", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
#     perc = "05", median_comparison_abundance_in = FALSE
# )
# Maaslin3(
#     rds = rds_05, output_folder = "Output/MAASLIN3/05_disc/05_subtentorial_lesions_onlygc_F",
#     main_el = "subtentorial_lesions", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
#     perc = "05", median_comparison_abundance_in = FALSE
# )

# rds_Euk <- readRDS("Output/SUPERVISED_DEC/Eukaryote_Supervised_decontam0.001.rds")
# Maaslin3(
#     rds = rds_Euk, output_folder = "Output/MAASLIN3/Eukaryote_MSHD_disc_F",
#     main_el = "category", categorySel = "category", elements = c("HEALTHY", "MS"), discr_els = c("sex", "age", "bmi"),
#     median_comparison_abundance_in = FALSE
# )

# rds_Arch <- readRDS("Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.001.rds")

# Maaslin3(
#     rds = rds_Arch, output_folder = "Output/MAASLIN3/Archaea_MSHD_disc_F",
#     main_el = "category", categorySel = "category", elements = c("HEALTHY", "MS"), discr_els = c("sex", "age", "bmi"),
#     median_comparison_abundance_in = FALSE
# )
# Maaslin3(
#     rds = rds_001, output_folder = "Output/MAASLIN3/Bacteria_MSHD_disc_F",
#     main_el = "category", categorySel = "category", elements = c("HEALTHY", "MS"), discr_els = c("sex", "age", "bmi"),
#     median_comparison_abundance_in = FALSE
# )
