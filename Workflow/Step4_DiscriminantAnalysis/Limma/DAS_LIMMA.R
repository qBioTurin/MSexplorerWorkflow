source("utilities.R")
output_folder = "Output/LIMMA_score/"
createFolder(output_folder)
output_folder264 = "Output/DAS_ONLY_LIMMA/264/"
createFolder(output_folder264)
output_folder167 = "Output/DAS_ONLY_LIMMA/167/"
createFolder(output_folder167)
output_folder116 = "Output/DAS_ONLY_LIMMA/116/"
createFolder(output_folder116)
output_folder39 = "Output/DAS_ONLY_LIMMA/39/"
createFolder(output_folder39)


output_folder = "Output/DAS_ONLY_LIMMA_HD/"
createFolder(output_folder)
output_folder_hd_264 = "Output/DAS_ONLY_LIMMA_HD/264/"
createFolder(output_folder264)
output_folder_hd_39 = "Output/DAS_ONLY_LIMMA_HD/39/"
createFolder(output_folder39)




baselines_decB = readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam.rds")
baselines_decA = readRDS(file = "Output/SUPERVISED_DEC/Archaea_Supervised_decontam.rds")
baselines_decE = readRDS(file = "Output/SUPERVISED_DEC/Eukaryota_Supervised_decontam.rds")

baselines_dec264 = readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontamOLD.rds")
print(baselines_dec264)
baselines_dec167 = readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontamMidWay.rds")
baselines_dec116 = readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam.rds")
baselines_dec39 = readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_dec_elabundance005.rds")
print(baselines_dec39)


das_limma <- function(baselines_dec,analisys,status,output_folder,Domain,info) {

baselines_dec<-baselines_dec264
analisys="gc_treatment"
status="healty"

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
metadata = read.csv("input/20241205_MetadataHSvsT0_modified.csv", 
                    header = TRUE, 
                    sep = ",",
                    na = c("", " ", "NA"), 
                    check.names = TRUE)
samples = c(colnames(norm_data))
metadata = metadata[metadata$id %in% samples,] # 71 samples
metadata = metadata[,-1]
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
###for analysis gluctrt pos/neg comparing bone marrow lesion#####
############################
# select metadata to visualize in heatmapt
if(analisys=="bone_marrow_lesions"){
  if(status=="positive"){
        metadata_hm = metadata %>% 
            filter(gc_treatment =="positive") %>%
            select(bone_marrow_lesions) %>%
            mutate(bone_marrow_lesions = as.factor(bone_marrow_lesions))
  }else if (status == "negative"){
        metadata_hm = metadata %>% 
            filter(gc_treatment =="negative") %>%
            select(bone_marrow_lesions) %>%
            mutate(bone_marrow_lesions = as.factor(bone_marrow_lesions)) }
  else if (status == "both") 
      {
        metadata_hm = metadata %>% 
            filter(gc_treatment =="positive" | gc_treatment =="negative") %>%
            select(bone_marrow_lesions) %>%
            mutate(bone_marrow_lesions = as.factor(bone_marrow_lesions)) 
        }
  
  norm_data = norm_data %>%
  select(rownames(metadata_hm))

  log_data <- log2(norm_data + 1)

  design <- model.matrix(~ bone_marrow_lesions, data = metadata_hm)
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

# Save DAS
#################
filterForHeatMap=norm_data[rownames((top_table)),]
if (status==""){
    write.csv(filterForHeatMap, file=gsub(" ","",paste(output_folder,Domain,info,"_limma_",analisys,".csv")))}
else {    
    write.csv(filterForHeatMap, file=gsub(" ","",paste(output_folder,Domain,info,"_limma_",analisys,"_",status,".csv")))
  }
}

array1=c("category","gc_treatment")
array2=c("lesion_burden","bone_marrow_lesions","gadolinium_contrast","subtentorial_lesions")
array3=c("positive","negative","both")
for (i in 1:length(array1)){
    das_limma(baselines_decB,array1[i],"",output_folder,"Bacteria")
 #   das_limma(baselines_decA,array1[i],"",output_folder,"Archaea")
 #   das_limma(baselines_decE,array1[i],"",output_folder,"Eukaryota")
}
for (j in 1:length(array2)){
        for (k in 1:length(array3)){
            das_limma(baselines_decB,array2[j],array3[k],output_folder,"Bacteria")
  #          das_limma(baselines_decA,array2[j],array3[k],output_folder,"Archaea")
  #          das_limma(baselines_decE,array2[j],array3[k],output_folder,"Eukaryota")
        }
}

for(i in 1:length(array2)){

    das_limma(baselines_dec264,array2[i],"both",output_folder264,"Bacteria","_264")
    das_limma(baselines_dec167,array2[i],"both",output_folder167,"Bacteria","_167")
    das_limma(baselines_dec116,array2[i],"both",output_folder116,"Bacteria","_116")
    das_limma(baselines_dec39,array2[i],"both",output_folder39,"Bacteria","_39")
}
das_limma(baselines_dec264,"gc_treatment","both",output_folder264,"Bacteria","_264")
das_limma(baselines_dec167,"gc_treatment","both",output_folder167,"Bacteria","_167")  
das_limma(baselines_dec116,"gc_treatment","both",output_folder116,"Bacteria","_116")
das_limma(baselines_dec39,"gc_treatment","both",output_folder39,"Bacteria","_39")


das_limma(baselines_dec264,"category","both",output_folder_hd_264,"Bacteria","_hd_264")
das_limma(baselines_dec39,"category","both",output_folder_hd_39,"Bacteria","_hd_39")
