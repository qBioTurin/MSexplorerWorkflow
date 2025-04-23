source("Settings/utilities.R")
output_folder = "Output/Leaps/"
createFolder(output_folder)

B = readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
E = readRDS("Output/SUPERVISED_DEC/Eukaryote_Supervised_decontam0.001.rds")
A = readRDS("Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.001.rds")
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
source("Workflow/Step7_Leaps/FunctionLRM.R")

baselines_dec_metadata = baselines_dec_metadata %>%
       mutate(
              across(.cols = c(id), .fns = as.character),
    
              across(.cols = c(sex, category, clinical_presentation, gc_treatment,
                     subtentorial_lesions, spinal_cord_lesion, gadolinium_contrast,
                     sequencing_batch,WORSENING, EDSS_DIAGNOSI, EDSS_PROGRESSIONE, Event,
                     naive, previous_therapy, antibiotic_use, sample_type,lesion_burden), 
                     .fns = as.factor),
    
              across(.cols = c(age, bmi, EventTime), .fns = as.numeric),
              across(.cols = c(sample_collection_date), .fns = as.Date)
       )

# Define predictors and outcomes
################################
predictors = c("sex", "age", "bmi", "gc_treatment", "lesion_burden", "spinal_cord_lesion",
"subtentorial_lesions","gadolinium_contrast")
outcome = c("Observed", "Shannon", "Simpson")

# Run function
###############
res = LRM_microbiome(metadata = newbaselines_dec_metadata %>% as.data.frame(), predictors = predictors_ale,outcomes = outcome)

# Save results
##############
filename = gsub("rds","txt", i)
filename = paste0(output_folder, filename)

sink(file = filename, type = c("output", "message"), split = TRUE)
print(paste("Observed", " ", filename))
print(res$Observed$Intercept)
print(paste("Shannon", " ", filename))
print(res$Shannon$Intercept)
print(paste("Simpson", " ", filename))
print(res$Simpson$Intercept)

sink()


