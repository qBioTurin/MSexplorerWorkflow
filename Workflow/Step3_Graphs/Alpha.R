source("Settings/utilities.R")
output_folder <- "Output/NEW_PLOT_RDS/"
meta_alpha <- "Output/AlphaMetadata/"

createFolder(output_folder)
createFolder(meta_alpha)
# Calculate alfa diversity
######################

Alpha <- function(baselines_dec, Domain, output_folder) {
  baselines_dec_richness <- estimate_richness(baselines_dec, split = TRUE, measures = c("Observed", "Shannon", "Simpson", "Chao1"))
  baselines_dec_richness <- data.frame(id = row.names(baselines_dec_richness), baselines_dec_richness)

  baselines_dec_metadata <- data.frame(baselines_dec@sam_data)
  baselines_dec_metadata <- merge(baselines_dec_richness, baselines_dec_metadata, all.y = TRUE, by = "id")



  write.csv(baselines_dec_metadata, file = gsub(" ", "", paste(meta_alpha, Domain, "_alpha_metadata.csv")), row.names = FALSE)

  baselines_dec_metadata1 <- baselines_dec_metadata %>%
    tidyr::gather(Observed, Simpson, Shannon, key = "index", value = "value")
  ###### healty vs ms

  custom_labels <- c("negative" = "Untreated", "positive" = "Treated")

  category <- ggplot(baselines_dec_metadata1, aes(x = category, y = value, na.rm = TRUE)) +
    geom_boxplot(aes(fill = category)) +
    facet_wrap(~index, scales = "free_y", nrow = 1) +
    labs(x = "", y = "Diversity Indexes", fill = "") +
    geom_signif(
      comparisons = list(c("HEALTHY", "MS")),
      map_signif_level = function(p) {
        if (p < 0.001) {
          return(paste0("*** (", signif(p, 2), ")"))
        } else if (p < 0.01) {
          return(paste0("** (", signif(p, 2), ")"))
        } else if (p < 0.049) {
          return(paste0("* (", signif(p, 2), ")"))
        } else if (p < 0.055) {
          return(paste0("NS. (", signif(p, 2), ")"))
        } else {
          return(paste0("NS."))
        }
      }
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

  saveRDS(category, gsub(" ", "", paste(output_folder, Domain, "_alpha_category.rds")))


  ###### Glucorticoid_treatment
  sum(is.na(baselines_dec_metadata1$gc_treatment))
  clean_baselines_dec_metadata1 <- baselines_dec_metadata1[!is.na(baselines_dec_metadata1$gc_treatment), ]

  #MS1 <- clean_baselines_dec_metadata1 %>%
  #  filter(category == "MS") %>%
  #  droplevels()
  MS1 <- clean_baselines_dec_metadata1 %>%
    filter(gc_treatment %in% c("positive", "negative","healthy")) %>%
    droplevels()
  custom_labels_gc <- c("healthy" = "Healthy", "positive" = "Treated", "negative" = "Untreated")

  # Calcola y_position dinamicamente per ogni indice
  library(dplyr)
  y_positions <- MS1 %>%
    group_by(index) %>%
    summarise(max_value = max(value, na.rm = TRUE)) %>%
    mutate(
      y_pos_1 = max_value * 1.05,
      y_pos_2 = max_value * 1.15,
      y_pos_3 = max_value * 1.25
    )

  # Crea una lista di y_position per ogni indice
  y_position_list <- split(y_positions[, c("y_pos_1", "y_pos_2", "y_pos_3")], y_positions$index)

  gc_treatment <- ggplot(MS1, aes(x = factor(gc_treatment, levels = c("healthy", "positive", "negative")), y = value, na.rm = TRUE)) +
    geom_boxplot(aes(fill = gc_treatment)) +
    facet_wrap(~index, scales = "free_y", nrow = 1) +
    labs(x = "", y = "Diversity Indexes", fill = "") +
    geom_signif(
      comparisons = list(c("healthy", "positive"), c("positive", "negative"), c("healthy", "negative")),
      map_signif_level = function(p) {
        if (p < 0.001) {
          return(paste0("*** (", signif(p, 2), ")"))
        } else if (p < 0.01) {
          return(paste0("** (", signif(p, 2), ")"))
        } else if (p < 0.049) {
          return(paste0("* (", signif(p, 2), ")"))
        } else if (p < 0.055) {
          return(paste0("NS. (", signif(p, 2), ")"))
        } else {
          return(paste0("NS."))
        }
      },
      y_position = unlist(y_position_list) # valori dinamici per ogni indice
    ) +
    theme_classic() +
    scale_x_discrete(labels = custom_labels_gc) +
    scale_fill_manual(values = c("healthy" = "#6EE2FF99", "positive" = "#4D4D4D", "negative" = "#D7D7D7")) +
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
  saveRDS(gc_treatment, gsub(" ", "", paste(output_folder, Domain, "_alpha_gc_treatment.rds")))
}


execute_alpha <- function() {
  baselines_decB <- readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
  baselines_decA <- readRDS(file = "Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.001.rds")
  baselines_decE <- readRDS(file = "Output/SUPERVISED_DEC/Eukaryote_Supervised_decontam0.001.rds")

  Alpha(baselines_decB, "Bacteria", output_folder)
  Alpha(baselines_decA, "Archaea", output_folder)
  Alpha(baselines_decE, "Eukaryote", output_folder)
}

execute_alpha()
