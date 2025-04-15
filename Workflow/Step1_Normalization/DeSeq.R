source("Settings/utilities.R")
output_folder = "Output/DESEQ_RDS/"
createFolder(output_folder)
# Input data
############
ALL_baselines = readRDS("Output/RDS/ALL_baseline_phylo.rds")

# Convert otu table to dataframe
#prova_df = as.data.frame(prova@otu_table)
ALL_baselines_df = as.data.frame(ALL_baselines@otu_table)

coldata_filt = as.data.frame(ALL_baselines@sam_data)
coldata_filt$sequencing_batch
# Convert into a DESeq Dataset
ALL_baselines_df <- DESeqDataSetFromMatrix(countData = round(ALL_baselines_df), 
                                         colData = coldata_filt, 
                                         design = ~ sequencing_batch)

# Normalize ALL_baselines_DDS
ALL_baselines_DDS = DESeq(ALL_baselines_df, test="Wald", fitType="parametric")

# Inspect Normalized Counts
normalized_counts <- counts(ALL_baselines_DDS, normalized = TRUE)

# Round to integer the results of DESeq2 variance-stabilization of counts
# and create OTU
OTU_ds <- round(otu_table(normalized_counts, taxa_are_rows = TRUE))

# Create new phyloseq
#####################
ALL_baselines_ds = phyloseq(OTU_ds, ALL_baselines@tax_table, ALL_baselines@sam_data)

# Save normalized data for ALL kingdoms phyloseq
################################################
saveRDS(ALL_baselines_ds, file = gsub(" ","",paste(output_folder,"ALL_DeSeq.rds"))) 


# Subset by Kingdom
###################
bact_baselines_ds = subset_taxa(ALL_baselines_ds, Kingdom == "Bacteria")
saveRDS(bact_baselines_ds, file = gsub(" ","",paste(output_folder,"Bact_DeSeq.rds")))

# Subset Archaea
archaea_baselines_ds = subset_taxa(ALL_baselines_ds, Kingdom == "Archaea")
saveRDS(archaea_baselines_ds, file = gsub(" ","",paste(output_folder,"Archaea_DeSeq.rds"))) 

# Subset Eukaryotes
euk_baselines_ds = subset_taxa(ALL_baselines_ds, Kingdom == "Eukaryota")
saveRDS(euk_baselines_ds, file = gsub(" ","",paste(output_folder,"Eukaryota_DeSeq.rds")))

