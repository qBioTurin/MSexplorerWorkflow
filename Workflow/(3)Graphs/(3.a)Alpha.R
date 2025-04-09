source("utilities.R")
output_folder = "Output/NEW_PLOT_RDS/"
createFolder(output_folder)
baselines_decB = readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam.rds")
baselines_decA = readRDS(file = "Output/SUPERVISED_DEC/Archaea_Supervised_decontam.rds")
baselines_decE = readRDS(file = "Output/SUPERVISED_DEC/Eukaryota_Supervised_decontam.rds")

# Calculate alfa diversity
######################

Alpha <- function(baselines_dec,Domain,output_folder) {
    baselines_dec=baselines_decB
    
    baselines_dec_richness = estimate_richness(baselines_dec, split = TRUE, measures = c("Observed","Shannon", "Simpson", "Chao1"))
    baselines_dec_richness <- data.frame(id = row.names(baselines_dec_richness), baselines_dec_richness)

    baselines_dec_metadata = data.frame(baselines_dec@sam_data)
    baselines_dec_metadata = merge(baselines_dec_richness, baselines_dec_metadata, all.y=TRUE, by = "id")

# add column 
    baselines_dec_metadata$factor_days_gc_collection90 <- cut(baselines_dec_metadata$days_gc_collection,
                                      breaks = c(0, 90, Inf),
                                      labels = c("0-90", "Above 90"),
                                      right = FALSE)


  write.csv(baselines_dec_metadata,  gsub(" ","",paste("input/",Domain,"_alpha_metadata.csv")))
    baselines_dec_metadata1=baselines_dec_metadata%>%
    tidyr::gather(Observed,Simpson,Shannon,key="index", value="value" )


######healty vs ms

    custom_labels <- c("negative" = "Untreated", "positive" = "Treated")

    category = ggplot(baselines_dec_metadata1, aes(x = category, y = value, na.rm = TRUE)) +
    geom_boxplot(aes(fill = category)) +
    facet_wrap(~index, scales = "free_y", nrow = 1) +
    labs(x = "", y = "Diversity Indexes",fill="") +
    #geom_signif(comparisons = list(c("HEALTHY","MS")),map_signif_level = FALSE) +
    geom_signif(comparisons = list(c("HEALTHY", "MS")),
              map_signif_level = function(p) {
                if (p < 0.049) {
                  return(paste0("* (", signif(p, 2), ")"))
                } else if (p < 0.055) {
                    return(paste0("NS. (", signif(p, 2), ")"))
                }
                  else{ 
                    return(paste0("NS."))
                }
              })+
    theme_classic() +
    scale_x_discrete(labels = custom_labels) +
    scale_fill_manual(values=c("#6EE2FF99","#ff410D99"),labels=c("HD","MS")) +
    theme(
        axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15),   
        axis.text.y = element_text(size = 15), 
        plot.title = element_text(size = 20),    
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 20))

    saveRDS(category, gsub(" ","",paste(output_folder,Domain,"_alpha_category.rds")))


######Glucorticoid_treatment
    sum(is.na(baselines_dec_metadata1$gc_treatment))
    clean_baselines_dec_metadata1 <- baselines_dec_metadata1[!is.na(baselines_dec_metadata1$gc_treatment), ]

    MS1 = clean_baselines_dec_metadata1 %>%
    filter(category == "MS") %>%
    droplevels()

    gc_treatment = ggplot(MS1, aes(x = gc_treatment, y = value, na.rm = TRUE)) +
    geom_boxplot(aes(fill = gc_treatment)) +
    facet_wrap(~index, scales = "free_y", nrow = 1) +
    labs(x = "", y = "Diversity Indexes",fill="") +
    geom_signif(comparisons = list(c("negative", "positive")),
              map_signif_level = function(p) {
                if (p < 0.049) {
                  return(paste0("* (", signif(p, 2), ")"))
                } else if (p < 0.055) {
                    return(paste0("NS. (", signif(p, 2), ")"))
                }
                  else{ 
                    return(paste0("NS."))
                }
              })+
    theme_classic() +
    scale_x_discrete(labels = custom_labels) +
    scale_fill_manual(values=c("#6EE2FF99","#ff410D99")) +
    theme(
        axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15),   
        axis.text.y = element_text(size = 15), 
        plot.title = element_text(size = 20),    
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 20),
        legend.position = "none")

    saveRDS(gc_treatment, gsub(" ","",paste(output_folder,Domain,"_alpha_gc_treatment.rds")))


######Glucorticoid_days
    clean_MS1 <- MS1[!is.na(MS1$factor_days_gc_collection90), ]
    gc_days = ggplot(clean_MS1, aes(x = factor_days_gc_collection90, y = value, na.rm = TRUE)) +
    geom_boxplot(aes(fill = factor_days_gc_collection90)) +
    facet_wrap(~index, scales = "free_y", nrow = 1) +
    labs(x = "Days after GC collection", y = "Diversity Indexes",fill="") +
    geom_signif(comparisons = list(c("0-90", "Above 90")),
              map_signif_level = function(p) {
                if (p < 0.049) {
                  return(paste0("* (", signif(p, 2), ")"))
                } else if (p < 0.055) {
                    return(paste0("NS. (", signif(p, 2), ")"))
                }
                  else{ 
                    return(paste0("NS."))
                }
              })+
    theme_classic() +
    scale_x_discrete(labels = custom_labels) +
    scale_fill_manual(values=c("#6EE2FF99","#ff410D99")) +
    theme(
        axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15), 
        axis.text.x = element_text(size = 15),   
        axis.text.y = element_text(size = 15), 
        plot.title = element_text(size = 20),    
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 20),
        legend.position = "none")
    saveRDS(gc_days, gsub(" ","",paste(output_folder,Domain,"_alpha_gc_days.rds")))
}

Alpha(baselines_decB,"Bacteria",output_folder)
Alpha(baselines_decA,"Archaea",output_folder)
Alpha(baselines_decE,"Eukaryota",output_folder)


