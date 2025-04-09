B = readRDS("~/Downloads/postProcess-main-5/Output/SUPERVISED_DEC/Bacteria_Supervised_decontam264.rds")
E = readRDS("~/Downloads/postProcess-main-5/Output/SUPERVISED_DEC/Eukaryota_Supervised_decontam.rds")
A = readRDS("~/Downloads/postProcess-main-5/Output/SUPERVISED_DEC/Archaea_Supervised_decontam.rds")
baselines_dec = merge_phyloseq(A,B,E)

baselines_dec_metadata = as.data.frame(baselines_dec@sam_data)

library(dplyr)

keep2 <- rownames(baselines_dec@sam_data)[grepl("^MS", rownames(baselines_dec@sam_data))]

baselines_dec = prune_samples(keep2, baselines_dec)

# Calculate alpha diversity
###########################
baselines_dec_richness = estimate_richness(baselines_dec, split = TRUE, measures = c("Observed","Shannon", "Simpson", "Chao1"))
baselines_dec_richness <- data.frame(id = row.names(baselines_dec_richness), baselines_dec_richness)

baselines_dec_metadata = data.frame(baselines_dec@sam_data)
#unisco metadati con le diversity calcolate in un unico dataframe
baselines_dec_metadata = merge(baselines_dec_richness, baselines_dec_metadata, all.y=TRUE, by = "id")


# Perform model regression analysis 
####################################
source("FunzioneLRM.R")

baselines_dec_metadata = baselines_dec_metadata %>%
  mutate(
    across(id, .fns= as.character),
    across(.cols = c(category, naive, 
                     sex, smoking_habit, physical_activity, edss,
                     antibiotic_use, spike_in, sample_type, therapy, sequencing_batch, age_terciles, pyr_median, simplified_batch,
                     bmi_classes, disease_actvity, edss_t0, edss_t12, edss_t24, 
                     gc_treatment, lesion_burden, bone_marrow_lesions, subtentorial_lesions, 
                     gadolinium_contrast, clinical_presentation, full_recovery, 
                     new_lesions_t12, relapse_t12, new_lesions_t24, relapse_t24, 
                     medas_percent, therapy_change, prognosis_at_onset, center),
           .fns = as.factor),
    across(.cols = c(Observed, Chao1, se.chao1, Shannon, Simpson, age, bmi, pyr_mds, 
                     dna_quantification, days_gc_collection, medas_percent, 
                     percent_th17, treg_il10_pos, treg_cd39_pos), 
           .fns = as.numeric),
    across(c(sample_collection_date, onset_date, diagnosis_date, glc_date_before_t0), 
           as.POSIXct, format = "%Y-%m-%d %H:%M:%S"))

# Define predictors and outcomes
################################
predittori = c("sex", "age", "bmi", "gc_treatment", "lesion_burden", "bone_marrow_lesions","subtentorial_lesions","gadolinium_contrast")
outcome = c("Observed", "Shannon", "Simpson")

# Run function
###############
res = LRM_microbiome(metadata = baselines_dec_metadata %>% as.data.frame(), predictors = predittori,outcomes = outcome)

# Save results
##############
filename = gsub("rds","txt", i)

sink(file = filename, type = c("output", "message"), split = TRUE)
print(paste("Observed", " ", filename))
print(res$Observed$Intercept)
print(paste("Shannon", " ", filename))
print(res$Shannon$Intercept)
print(paste("Simpson", " ", filename))
print(res$Simpson$Intercept)
sink()

}
