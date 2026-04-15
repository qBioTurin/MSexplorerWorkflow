source("Settings/utilities.R")
output_folder <- "Output/NEW_PLOT_RDS/"
meta_alpha <- "Output/AlphaMetadata/"

createFolder(output_folder)
createFolder(meta_alpha)
# Calculate alfa diversity
######################

Alpha <- function(baselines_dec, Domain, output_folder) {
  baselines_dec_richness <- estimate_richness(baselines_dec, split = TRUE, measures = c("Shannon", "Simpson"))
  baselines_dec_richness <- data.frame(id = row.names(baselines_dec_richness), baselines_dec_richness)

  baselines_dec_metadata <- data.frame(baselines_dec@sam_data)
  baselines_dec_metadata <- merge(baselines_dec_richness, baselines_dec_metadata, all.y = TRUE, by = "id")



  write.csv(baselines_dec_metadata, file = gsub(" ", "", paste(meta_alpha, Domain, "_alpha_metadata.csv")), row.names = FALSE)

  baselines_dec_metadata1 <- baselines_dec_metadata %>%
    tidyr::gather(Shannon, Simpson, key = "index", value = "value")
  ###### healty vs ms

  custom_labels <- c("negative" = "Untreated", "positive" = "Treated")

  category <- ggplot(baselines_dec_metadata1, aes(x = category, y = value, na.rm = TRUE)) +
    geom_boxplot(aes(fill = category)) +
    facet_wrap(~index, scales = "free_y", nrow = 1) +
    labs(x = "", y = "Diversity Indexes", fill = "") +
    geom_signif(
      comparisons = list(c("HEALTHY", "MS")),
      map_signif_level = function(p) {
        sapply(p, function(x) {
          if (x < 0.001) {
            paste0("*** (", signif(x, 2), ")")
          } else if (x < 0.01) {
            paste0("** (", signif(x, 2), ")")
          } else if (x < 0.05) {
            paste0("* (", signif(x, 2), ")")
          } else if (x < 0.1) {
            paste0("NS. (", signif(x, 2), ")")
          } else {
            "NS."
          }
        })
      },
      textsize = 5
    ) +
    theme_classic() +
    scale_x_discrete(labels = custom_labels) +
    scale_fill_manual(values = c("#6EE2FF99", "#ff410D99"), labels = c("HD", "MS")) +
    theme(
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      plot.title = element_text(size = 20),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15),
      strip.text = element_text(size = 20)
    )
  category
  saveRDS(category, gsub(" ", "", paste(output_folder, Domain, "_alpha_category.rds")))


  ###### Glucorticoid_treatment
  sum(is.na(baselines_dec_metadata1$gc_treatment))
  clean_baselines_dec_metadata1 <- baselines_dec_metadata1[!is.na(baselines_dec_metadata1$gc_treatment), ]

  # MS1 <- clean_baselines_dec_metadata1 %>%
  #  filter(category == "MS") %>%
  #  droplevels()
  MS1 <- clean_baselines_dec_metadata1 %>%
    filter(gc_treatment %in% c("positive", "negative", "healthy")) %>%
    droplevels()
  custom_labels_gc <- c("healthy" = "Healthy", "positive" = "Treated", "negative" = "Untreated")


  create_plot_gc_treatment <- function(plot_MS, y_name) {
    plot_MS$gc_treatment <- factor(plot_MS$gc_treatment,
      levels = c("healthy", "positive", "negative")
    )
    y_max <- max(plot_MS$value, na.rm = TRUE)
    offset <- 0.1 * diff(range(plot_MS$value, na.rm = TRUE))

    comparisons <- list(
      c("healthy", "positive"),
      c("positive", "negative"),
      c("healthy", "negative")
    )


    y_positions <- y_max + (1:length(comparisons)) * offset

    gc_treatment <- ggplot(plot_MS, aes(x = gc_treatment, y = value)) +
      geom_boxplot(aes(fill = gc_treatment)) +
      facet_wrap(~index, scales = "free_y", nrow = 1) +
      labs(x = "", y = y_name, fill = "") +
      geom_signif(
        comparisons = comparisons,
        y_position = y_positions,
        map_signif_level = function(p) {
          sapply(p, function(x) {
            if (x < 0.001) {
              paste0("*** (", signif(x, 2), ")")
            } else if (x < 0.01) {
              paste0("** (", signif(x, 2), ")")
            } else if (x < 0.05) {
              paste0("* (", signif(x, 2), ")")
            } else if (x < 0.1) {
              paste0("NS. (", signif(x, 2), ")")
            } else {
              "NS."
            }
          })
        },
        textsize = 6
      ) +
      theme_classic() +
      scale_x_discrete(labels = custom_labels_gc) +
      scale_fill_manual(values = c(
        "healthy" = "#6EE2FF99",
        "positive" = "#4D4D4D",
        "negative" = "#D7D7D7"
      )) +
      theme(
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 20),
        legend.position = "none"
      )
  }
  shannon <- create_plot_gc_treatment(MS1 %>% filter(index == "Shannon"), "Diversity Indexes") + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  simpson <- create_plot_gc_treatment(MS1 %>% filter(index == "Simpson"), "") + theme(axis.title.x = element_blank(), axis.text.x = element_blank())

  gc_treatment <- shannon + simpson

  saveRDS(gc_treatment, gsub(" ", "", paste(output_folder, Domain, "_alpha_gc_treatment.rds")))
}


execute_alpha <- function() {
  baselines_decB <- readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.rds")
  baselines_decA <- readRDS(file = "Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.rds")
  baselines_decE <- readRDS(file = "Output/SUPERVISED_DEC/Eukaryote_Supervised_decontam0.rds")

  Alpha(baselines_dec = baselines_decB, Domain = "Bacteria", output_folder)
  Alpha(baselines_dec = baselines_decA, Domain = "Archaea", output_folder)
  Alpha(baselines_dec = baselines_decE, Domain = "Eukaryote", output_folder)
}

execute_alpha()
