source("Settings/utilities.R")
source("Settings/Packages.R")
B = readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
E = readRDS("Output/SUPERVISED_DEC/Eukaryote_Supervised_decontam0.001.rds")
A = readRDS("Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.001.rds")
baselines_dec = merge_phyloseq(A,B,E)

# Import metadata
baselines_metadata = read.csv("/Users/doratortarolo/Documents/Microbiome/MS/MS_HD/Data/Metadata/20241115/20241115_MetadataHSvsT0_modified.csv", 
                              header = TRUE, 
                              sep = ",",
                              na = c("", " ", "NA"), 
                              check.names = TRUE)
baselines_metadata = baselines_metadata[,-1]


baselines_metadata = baselines_metadata %>%
  mutate(
    across(id, .fns= as.character),
    across(.cols = c(category, naive, 
                     previous_therapy, sex, smoking_habit, physical_activity, edss,
                     antibiotic_use, spike_in, sample_type, therapy, sequencing_batch, 
                     glucocorticoid_treatment, age_terciles, pyr_median, simplified_batch,
                     bmi_classes, disease_actvity, edss_t0, edss_t12, edss_t24, 
                     gc_treatment, lesion_burden, bone_marrow_lesions, subtentorial_lesions, 
                     gadolinium_contrast, clinical_presentation, full_recovery, 
                     new_lesions_t12, relapse_t12, new_lesions_t24, relapse_t24, 
                     medas_percent, therapy_change, prognosis_at_onset, center),
           .fns = as.factor),
    across(.cols = c(age, bmi, pyr_mds, dna_quantification, days_gc_collection, medas_percent, 
                     percent_th17, treg_il10_pos, treg_cd39_pos), 
           .fns = as.numeric),
    across(c(sample_collection_date, onset_date, diagnosis_date, glc_date_before_t0), 
           as.POSIXct, format = "%Y-%m-%d %H:%M:%S"))

# convert names to new ids
new_ids = read_excel("/Users/doratortarolo/Documents/Microbiome/MS/MS_HD/Rebuttal/Correlation_lifestyle/conv_names.xlsx")
new_ids = new_ids %>%
  dplyr::rename(new_id = id, 
         id = old_id)

baselines_metadata = merge(baselines_metadata, new_ids, by = "id")

baselines_metadata = baselines_metadata %>%
  select(-id) %>%
  dplyr::rename(id = new_id) 


# Metadata analysis
# Remove non-naive patients from metadata
keep_samples = baselines_metadata$id[baselines_metadata$naive %in% c("yes", "healthy") & baselines_metadata$antibiotic_use == "no"]
baselines_metadata = baselines_metadata[baselines_metadata$id %in% keep_samples,]

# Remove HD
baselines_metadata = baselines_metadata %>%
  filter(category == "MS")

# pyrMDS
summary(baselines_metadata$pyr_mds)
skim(baselines_metadata$pyr_mds)

#physical activity
summary(baselines_metadata$physical_activity)
skim(baselines_metadata$physical_activity)


## Correlation
alpha = read.csv("Output/AlphaMetadata/Bacteria_alpha_metadata_obs.csv")
corr_alpha <- alpha %>%
  select(id, Observed, Shannon, Simpson) %>%
  merge(baselines_metadata, by = "id")



## PYR-MDS 
pyr_correlation <- corr_alpha %>%
  select(Observed, Shannon, Simpson, PyrMDS = pyr_mds) %>%
  na.omit() %>%
  as.data.frame()

# Helper to compute label text for Pearson r
pearson_label <- function(x, y) {
  ct <- cor.test(x, y, method = "pearson")
  r <- unname(ct$estimate)
  ci <- ct$conf.int
  sprintf("Pearson r = %.2f\n95%% CI [%.2f, %.2f]\np = %.3g", r, ci[1], ci[2], ct$p.value)
}

# Optional: Spearman label (uncomment if you want Spearman instead)
spearman_label <- function(x, y) {
  ct <- cor.test(x, y, method = "spearman", exact = FALSE)
  r <- unname(ct$estimate)
  sprintf("Spearman rho = %.2f\np = %.3g", r, ct$p.value)
}


base_font_size <- 12
theme_set(theme_minimal(base_size = base_font_size))

# (your existing plot code, unchanged except annotate size adjusted)
p1 <- ggplot(pyr_correlation, aes(x = Observed, y = PyrMDS)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, color = "#adbcc5", fill = "#bcc6d0") +
  labs(x = "Observed", y = "PyrMDS", title = "PyrMDS vs Observed") +
  theme(
    plot.title = element_text(size = base_font_size, face = "bold"),
    axis.title = element_text(size = base_font_size),
    axis.text = element_text(size = base_font_size)
  ) +
  annotate("text",
           x = Inf, y = -Inf,
           label = pearson_label(pyr_correlation$Observed, pyr_correlation$PyrMDS),
           hjust = 1.05, vjust = -0.5, size = 10 * 0.35) # Adjust size for better fit

p2 <- ggplot(pyr_correlation, aes(x = Shannon, y = PyrMDS)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, color = "#2c7fb8", fill = "#c6dbef") +
  labs(x = "Shannon index", y = "PyrMDS", title = "PyrMDS vs Shannon") +
  theme(
    plot.title = element_text(size = base_font_size, face = "bold"),
    axis.title = element_text(size = base_font_size),
    axis.text = element_text(size = base_font_size)
  ) +
  annotate("text",
           x = Inf, y = -Inf,
           label = pearson_label(pyr_correlation$Shannon, pyr_correlation$PyrMDS),
           hjust = 1.05, vjust = -0.5, size = 10 * 0.35) # Adjust size for better fit

p3 <- ggplot(pyr_correlation, aes(x = Simpson, y = PyrMDS)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, color = "#7fc97f", fill = "#d9f0d3") +
  labs(x = "Simpson index", y = "PyrMDS", title = "PyrMDS vs Simpson") +
  theme(
    plot.title = element_text(size = base_font_size, face = "bold"),
    axis.title = element_text(size = base_font_size),
    axis.text = element_text(size = base_font_size)
  ) +
  annotate("text",
           x = Inf, y = -Inf,
           label = pearson_label(pyr_correlation$Simpson, pyr_correlation$PyrMDS),
           hjust = 1.05, vjust = -0.5, size = 10 * 0.35) # Adjust size for better fit


# Combine side-by-side
p_pyr <- (p1 + p2 + p3) + plot_layout(ncol = 3, widths = c(1, 1))
ggsave("Output/Correlations_lifestyle/PyrMDS.png", plot = p_pyr, width = 7.27, height = 10.69, dpi = 300)

## PHYSICAL ACTIVITY
physact_correlation = corr_alpha %>%
  select(Observed, Shannon, Simpson, physical_activity) %>%
  na.omit(corr_alpha)

# normality
shapiro.test(physact_correlation$Observed) # not signif - normal
shapiro.test(physact_correlation$Shannon) # not signif - normal
shapiro.test(physact_correlation$Simpson) # signif - not normal

kruskal.test(Observed    ~ physical_activity, data = physact_correlation) # not signif
kruskal.test(Shannon    ~ physical_activity, data = physact_correlation) # not signif
kruskal.test(Simpson    ~ physical_activity, data = physact_correlation) # not signif

# Fit the model
anova_model <- aov(Observed ~ physical_activity, data = physact_correlation)
anova_model <- aov(Shannon ~ physical_activity, data = physact_correlation) # not signif
anova_model <- aov(Simpson ~ physical_activity, data = physact_correlation) # not signif

# View ANOVA table
summary(anova_model)

# visulization
phys_long <- physact_correlation %>%
  mutate(
    physical_activity = factor(
      physical_activity,
      levels = c("sedentary", "partially_active", "active"),
      labels = c("sedentary", "partially active", "active")
    )) %>%
  select(Observed, Shannon, Simpson, PhysicalActivity = physical_activity) %>%
  pivot_longer(
    cols      = c(Observed, Shannon, Simpson),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("Observed", "Shannon", "Simpson")
    )
  )

physact <- ggplot(phys_long, aes(x = PhysicalActivity, y = value, fill = PhysicalActivity)) +
  geom_boxplot() +
  facet_wrap(~ metric, scales = "free_y") +
  theme_classic() +
  labs(
    x    = "Physical Activity",
    y    = "Value",
    fill = "Activity Level"
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    strip.text   = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 12)
  )
ggsave("Output/Correlations_lifestyle/physact.png", plot = physact, width = 7.27, height = 10.69, dpi = 300)

