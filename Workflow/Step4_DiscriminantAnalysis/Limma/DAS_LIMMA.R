source("Settings/utilities.R")

output_folderMSHD="Output/LIMMA_score/MSHD/"
createFolder(output_folderMSHD)
output_folderGC="Output/LIMMA_score/GC/"
createFolder(output_folderGC)
output_folder = "Output/LIMMA_score/"
createFolder(output_folder)
output_folder_001 = "Output/LIMMA_score/001/"
createFolder(output_folder_001)
output_folder_01 = "Output/LIMMA_score/01/"
createFolder(output_folder_01)
output_folder_05 = "Output/LIMMA_score/05/"
createFolder(output_folder_05)
output_folder_GC_comp="Output/LIMMA_score/GC_comp/"
createFolder(output_folder_GC_comp)

das_limma <- function(baselines_dec,analisys,status,output_folder,Domain,info) {

# Keep only samples that have received gc treatment
samples = as.data.frame(sample_data(baselines_dec))

if(analisys!="category"){
  samples = samples[!samples$gc_treatment == "healthy",]#gc threatment only scheme not do for msvshd
  samples = samples[!is.na(samples$gc_treatment),]
}

baselines_dec = prune_samples(samples$id, baselines_dec)
# 264 taxa and 47 samples


################################# Prep matrix 
# Create merged taxa + otu table
# tax table
taxatables = as.data.frame(tax_table(baselines_dec))
taxatables = rownames_to_column(taxatables, var = "otu_id")

# otu table
otutables = data.frame(otu_table(baselines_dec))
otutables = as.data.frame(abundances(otutables, transform ="compositional"))
otutables = rownames_to_column(otutables, var = "otu_id")

# combine tables
combined_df = full_join(taxatables, otutables, by = "otu_id")

## from combined tables, keep only columns otu_id, Genus+species + otu table)
norm_data <- combined_df %>%
  mutate(Genus_species = paste(Genus, Species, sep = " "))
norm_data = norm_data[, -c(2:8)]
# 264 taxa x 47 samples

rownames(norm_data) = norm_data$Genus_species
norm_data = norm_data[,-c(1,length(norm_data))]
colnames(norm_data)
# 264 x 47 samples
  
  
## Select patients metadata 
metadata = read.csv("InputData/metadataMS.csv", 
                    header = TRUE, 
                    sep = ",",
                    na = c("", " ", "NA"), 
                    check.names = TRUE)
samples = c(colnames(norm_data))
metadata = metadata[metadata$id %in% samples,] # 71 samples
#metadata = metadata[,-1]
# 47 x 49
rownames(metadata) <- metadata$id


#####################
###for category #####
######################
# select metadata to visualize in heatmap
if(analisys=="category"){
  metadata_hm = metadata %>% 
    select(category) %>%
    mutate(category = as.factor(category))

  norm_data = norm_data %>%
    select(rownames(metadata_hm))

  log_data <- log2(norm_data + 1)

  design <- model.matrix(~ category, data = metadata_hm)
}

############################
###for analysis gluctrt#####
############################
# select metadata to visualize in heatmap
else if (analisys=="gc_treatment"){
  metadata_hm = metadata %>% 
    select(gc_treatment) %>%
    mutate(gc_treatment = as.factor(gc_treatment))

  norm_data = norm_data %>%
    select(rownames(metadata_hm))

  log_data <- log2(norm_data + 1)

  design <- model.matrix(~ gc_treatment, data = metadata_hm)
}

############################
###for analysis gluctrt pos/neg comparing lesion#####
############################
# select metadata to visualize in heatmapt
if(analisys=="lesion_burden"){
  if(status=="positive"){
        metadata_hm = metadata %>% 
            filter(gc_treatment =="positive") %>%
            select(lesion_burden) %>%
            mutate(lesion_burden = as.factor(lesion_burden))
  }else if (status == "negative")
  {
        metadata_hm = metadata %>% 
            filter(gc_treatment =="negative") %>%
            select(lesion_burden) %>%
            mutate(lesion_burden = as.factor(lesion_burden)) }
  else if (status == "both") 
      {
        metadata_hm = metadata %>% 
            filter(gc_treatment =="positive" | gc_treatment =="negative") %>%
            select(lesion_burden) %>%
            mutate(lesion_burden = as.factor(lesion_burden)) 
        }
  
  norm_data = norm_data %>%
  select(rownames(metadata_hm))

  log_data <- log2(norm_data + 1)

  design <- model.matrix(~ lesion_burden, data = metadata_hm)
}

############################


############################
###for analysis gluctrt pos/neg comparing spinal cord#####
############################
# select metadata to visualize in heatmapt
if(analisys=="spinal_cord_lesion"){
  if(status=="positive"){
        metadata_hm = metadata %>% 
            filter(gc_treatment =="positive") %>%
            select(spinal_cord_lesion) %>%
            mutate(spinal_cord_lesion = as.factor(spinal_cord_lesion))
  }else if (status == "negative"){
        metadata_hm = metadata %>% 
            filter(gc_treatment =="negative") %>%
            select(spinal_cord_lesion) %>%
            mutate(spinal_cord_lesion = as.factor(spinal_cord_lesion)) }
  else if (status == "both") 
      {
        metadata_hm = metadata %>% 
            filter(gc_treatment =="positive" | gc_treatment =="negative") %>%
            select(spinal_cord_lesion) %>%
            mutate(spinal_cord_lesion = as.factor(spinal_cord_lesion)) 
        }
  
  norm_data = norm_data %>%
  select(rownames(metadata_hm))

  log_data <- log2(norm_data + 1)

  design <- model.matrix(~ spinal_cord_lesion, data = metadata_hm)
}

############################
############################
###gaodolinum#####
############################
# select metadata to visualize in heatmapt

if(analisys=="gadolinium_contrast"){
  if(status=="positive"){
        metadata_hm = metadata %>% 
            filter(gc_treatment =="positive") %>%
            select(gadolinium_contrast) %>%
            mutate(gadolinium_contrast = as.factor(gadolinium_contrast))
  }else if (status == "negative"){
        metadata_hm = metadata %>% 
            filter(gc_treatment =="negative") %>%
            select(gadolinium_contrast) %>%
            mutate(gadolinium_contrast = as.factor(gadolinium_contrast)) }
  else if (status == "both") 
      {
        metadata_hm = metadata %>% 
            filter(gc_treatment =="positive" | gc_treatment =="negative") %>%
            select(gadolinium_contrast) %>%
            mutate(gadolinium_contrast = as.factor(gadolinium_contrast)) 
        }
  norm_data = norm_data %>%
  select(rownames(metadata_hm))

  log_data <- log2(norm_data + 1)

  design <- model.matrix(~ gadolinium_contrast, data = metadata_hm)
}

############################
############################
###subtentorial#####
############################
# select metadata to visualize in heatmapt
if(analisys=="subtentorial_lesions"){
  if(status=="positive"){
        metadata_hm = metadata %>% 
            filter(gc_treatment =="positive") %>%
            select(subtentorial_lesions) %>%
            mutate(subtentorial_lesions = as.factor(subtentorial_lesions))
  }else if (status == "negative"){
        metadata_hm = metadata %>% 
            filter(gc_treatment =="negative") %>%
            select(subtentorial_lesions) %>%
            mutate(subtentorial_lesions = as.factor(subtentorial_lesions)) }
  else if (status == "both") 
      {
        metadata_hm = metadata %>% 
            filter(gc_treatment =="positive" | gc_treatment =="negative") %>%
            select(subtentorial_lesions) %>%
            mutate(subtentorial_lesions = as.factor(subtentorial_lesions)) 
        }
  
  norm_data = norm_data %>%
  select(rownames(metadata_hm))

  log_data <- log2(norm_data + 1)

  design <- model.matrix(~ subtentorial_lesions, data = metadata_hm)
}
### LIMMA DAS analysis

# Fit linear models
fit <- lmFit(log_data, design)

# Apply empirical Bayes moderation
fit <- eBayes(fit)

top_table <- topTable(fit, coef=NULL, number = 100, adjust = "fdr")
top_table=top_table[top_table$P.Value < 0.05, ]
dim(top_table)
#print(top_table)
write.table(top_table, file=gsub(" ","",paste(output_folder,Domain,info,"_limma_",analisys,"_",status,"_top_table.tsv")), sep="\t", row.names=TRUE, col.names=NA)

if (info!="")
  info=paste0("_",info)
# Save DAS
#################
filterForHeatMap=norm_data[rownames((top_table)),]
if (status=="" | status=="both"){
    write.csv(filterForHeatMap, file=gsub(" ","",paste(output_folder,"/",Domain,"_",analisys,info,"_limma.csv")))}
else {    
    write.csv(filterForHeatMap, file=gsub(" ","",paste(output_folder,"/",Domain,"_",analisys,"_",status,info,"_limma.csv")))
  }
}


execute_limma <- function() {

  baselines_dec_001 = readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
  baselines_decA = readRDS(file = "Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.001.rds")
  baselines_decE = readRDS(file = "Output/SUPERVISED_DEC/Eukaryota_Supervised_decontam0.001.rds")

  baselines_dec_01 = readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds")
  baselines_dec_05 = readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds")
 
  analysis=c("lesion_burden","spinal_cord_lesion","gadolinium_contrast","subtentorial_lesions")

  das_limma(baselines_decA,"category","both",output_folderMSHD,"Archaea","")
  das_limma(baselines_decE,"category","both",output_folderMSHD,"Eukaryota","")
  das_limma(baselines_dec_001,"category","both",output_folderMSHD,"Bacteria","")
    
  das_limma(baselines_decA,"gc_treatment","both",output_folderGC,"Archaea","")
  das_limma(baselines_decE,"gc_treatment","both",output_folderGC,"Eukaryota","")
  das_limma(baselines_dec_001,"gc_treatment","both",output_folderGC,"Bacteria","")

  status=c("positive","negative")

  for(i in 1:length(analysis)){
    das_limma(baselines_dec_001,analysis[i],"both",output_folder_001,"Bacteria","001")
    das_limma(baselines_dec_01,analysis[i],"both",output_folder_01,"Bacteria","01")
    das_limma(baselines_dec_05,analysis[i],"both",output_folder_05,"Bacteria","05")
    for(j in 1:length(status)){
      das_limma(baselines_dec_001,analysis[i],status[j],output_folder_GC_comp,"Bacteria","")
      das_limma(baselines_decA,analysis[i],status[j],output_folder_GC_comp,"Archaea","")
      das_limma(baselines_decE,analysis[i],status[j],output_folder_GC_comp,"Eukaryota","")
    }
  }
  das_limma(baselines_dec_001,"gc_treatment","both",output_folder_001,"Bacteria","001")
  das_limma(baselines_dec_01,"gc_treatment","both",output_folder_01,"Bacteria","01")
  das_limma(baselines_dec_05,"gc_treatment","both",output_folder_05,"Bacteria","05")
}

execute_limma()


