source("utilities.R")
output_folder = "Output/RDS/"
createFolder(output_folder)
# Import Data
#############
# Import Biom
baselines <- import_biom("input/fin_hs_vs_MS_T0.biom")

# Import metadata
baselines_metadata = read.csv("input/20241205_MetadataHSvsT0_modified.csv", 
                              header = TRUE, 
                              sep = ",",
                              na = c("", " ", "NA"), 
                              check.names = TRUE)
baselines_metadata = baselines_metadata[,-1]


baselines_metadata = baselines_metadata %>%
  mutate(
    across(.cols=c(id,old_id), .fns= as.character),
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
  

#####################################
# BASELINE: Create Phyloseq Object #
####################################
# Prep taxa table
baselines@tax_table@.Data <- substring(baselines@tax_table@.Data, 4)
colnames(baselines@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
TAXA_baselines = tax_table(baselines@tax_table@.Data)

# Prep OTU table
OTU_baselines = otu_table(baselines@otu_table@.Data, taxa_are_rows = TRUE)
colnames(OTU_baselines) <- gsub("_bracken_species", "", colnames(OTU_baselines))
# Create phyloseq obj
ALL_baselines = phyloseq(OTU_baselines, TAXA_baselines)

otu = data.frame(otu_table(ALL_baselines))
taxa = data.frame(tax_table(ALL_baselines))

# Check that otu patients and metadata patients match
otu_patients = c(colnames(OTU_baselines))
metadata_patients = c(baselines_metadata$id)



# Create and add tree
TREE_baselines = rtree(ntaxa(ALL_baselines), rooted=TRUE, tip.label=taxa_names(ALL_baselines))

ALL_baselines = merge_phyloseq(ALL_baselines, TREE_baselines)

# Remove non-naive samples from biom + subjects with antibiotic use
keep_samples = baselines_metadata$id[baselines_metadata$naive %in% c("yes", "healthy") & baselines_metadata$antibiotic_use == "no"]
ALL_baselines = prune_samples(keep_samples, ALL_baselines)

# Create sample_data by associating metadata
##############################################
# Remove non-naive patients from metadata
baselines_metadata = baselines_metadata[baselines_metadata$id %in% keep_samples,]

# Remove previous_therapy column since we removed non-naive patients
baselines_metadata = subset(baselines_metadata, select = -c(`previous_therapy`))


sampledata = baselines_metadata[baselines_metadata$id %in% colnames(ALL_baselines@otu_table),] 
sampledata = sampledata %>% arrange(factor(id, levels = colnames(ALL_baselines@otu_table)))

baselines_sampledata = sample_data(data.frame(sampledata,
                                        row.names=sample_names(ALL_baselines@otu_table),
                                        stringsAsFactors=FALSE))

ALL_baselines = merge_phyloseq(ALL_baselines, baselines_sampledata)


# Quick overview of phyloseq object
microbiome::summarize_phyloseq(ALL_baselines)

# Save phyloseq with ALL kingdoms
saveRDS(ALL_baselines, file = gsub(" ","",paste(output_folder,"ALL_baseline_phylo.rds"))) 

# Subset by Kingdom
###################
# Subset Bacteria

bact_baselines = subset_taxa(ALL_baselines, Kingdom == "Bacteria")
saveRDS(bact_baselines, file =gsub(" ","",paste(output_folder,"Bact_baseline_phylo.rds")))

# Subset Eukaryotes
euk_baselines = subset_taxa(ALL_baselines, Kingdom == "Eukaryota")
saveRDS(euk_baselines, file =gsub(" ","",paste(output_folder,"Eukaryotes_baseline_phylo.rds"))) 


# Subset Archaea
archaea_baselines = subset_taxa(ALL_baselines, Kingdom == "Archaea")
saveRDS(archaea_baselines, file = gsub(" ","",paste(output_folder,"Archaea_baseline_phylo.rds")))

