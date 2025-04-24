source("Settings/utilities.R")
output_folder = "Output/RDS/"
createFolder(output_folder)
# Import Data
#############
# Import Biom
baselines <- import_biom("InputData/fin_hs_vs_MS_T0.biom")

# Import metadata
baselines_metadata = read.csv("InputData/metadataMS.csv", 
                              header = TRUE, 
                              sep = ",",
                              na = c("", " ", "NA"), 
                              check.names = TRUE)

baselines_metadata <- baselines_metadata %>%
  mutate(
    across(.cols = c(id), .fns = as.character),
    
    across(.cols = c(sex, category, clinical_presentation, gc_treatment,
                     subtentorial_lesions, spinal_cord_lesion, gadolinium_contrast,
                     sequencing_batch,WORSENING, EDSS_DIAGNOSI, EDSS_PROGRESSIONE, Event,
                     naive, previous_therapy, antibiotic_use, sample_type,lesion_burden), 
           .fns = as.factor),
    
    across(.cols = c(age, bmi ,EventTime
          ), .fns = as.numeric),
    
    across(.cols = c(sample_collection_date), .fns = as.Date)
  )

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
keep_samples = baselines_metadata$id[baselines_metadata$naive %in% c("yes", "healthy")  & baselines_metadata$antibiotic_use == "no"]
ALL_baselines = prune_samples(keep_samples, ALL_baselines)

# Create sample_data by associating metadata
##############################################
# Remove non-naive patients from metadata
baselines_metadata = baselines_metadata[baselines_metadata$id %in% keep_samples,]

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

