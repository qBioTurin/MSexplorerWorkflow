# Carica i pacchetti necessari
library(readr)
library(maaslin3)
library(stringr)
# Percorsi ai file CSV
metadata_path <- "InputData/Metadata_MS.csv"
data_path <- "InputData/MS_pathabundance_cpm_230925.tsv"
rename_path <- "InputData/rename.csv"
# Leggi il file CSV per la metadata
metadata <- read_csv(metadata_path)
rename<-read_csv(rename_path)
rename<-as.data.frame(rename)
# Leggi il file CSV per i dati
data <- read_tsv(data_path)
data <- as.data.frame(data)
colnames(data) <- sub("_.+", "", colnames(data))
# Check the column names to find the correct one for pathway
print(colnames(data))
# Set rownames using the correct column name (adjust if needed)
rownames(data)<-data$Pathway
data<-data%>% select(-Pathway)

data1<-t(data)
data1 <- as.data.frame(data1)
data1$old_id <- rownames(data1)

data1 <- merge(data1, rename, by.x = "old_id", by.y = "old_id", all.x = TRUE)
data1 <- data1[!is.na(data1$id), ]
rownames(data1)<-data1$id
dataf <- data1 %>% select(-c(old_id,id))
dataf <- dataf[rownames(dataf) %in% metadata$id, ]
datf<-t(dataf)
datf <- as.data.frame(datf)
metadata<-as.data.frame(metadata)
rownames(metadata) <- metadata$id
head(data1)
if (!dir.exists("Output/Diff_humann/category/")) {
    dir.create("Output/Diff_humann/category/", recursive = TRUE)
}


datf <- datf[!rownames(datf) %in% c("UNMAPPED", "UNINTEGRATED"), ]


pathways <- datf[!grepl("UNINTEGRATED|UNMAPPED", rownames(datf)), ]

library(dplyr)
library(stringr)
library(tibble)

pathways2 <- pathways %>% tibble::rownames_to_column("Pathway")
pathw_filt <- pathways2 %>%
  filter(!str_detect(Pathway, "\\|g__")) %>%
  filter(!str_detect(Pathway, "unclassified"))
rownames(pathw_filt) <- pathw_filt$Pathway
pathw_filt <- pathw_filt %>% select(-Pathway)
res <- maaslin3(
  input_data = pathw_filt,       
  input_metadata = metadata,          
  output = "Output/Diff_humann/category/",
  formula = "~ category",
  normalization = "TSS",
  transform = "LOG",
  min_prevalence = 0.01,
  min_abundance = 0.0001
)


metadataMS <- metadata[metadata$category == "MS",]
metadataMSpos <- metadataMS[metadataMS$gc_treatment == "positive" & !is.na(metadataMS$gc_treatment), ]
metadataMSneg <- metadataMS[metadataMS$gc_treatment == "negative" & !is.na(metadataMS$gc_treatment), ]

dataMS<-t(pathw_filt)
dataMS <- dataMS[rownames(dataMS) %in% metadataMS$id, ]
dataMS <- t(dataMS)
dataMS <- as.data.frame(dataMS)

dataMSpos<-t(pathw_filt)
dataMSpos <- dataMSpos[rownames(dataMSpos) %in% metadataMSpos$id, ]
dataMSpos <- t(dataMSpos)
dataMSpos <- as.data.frame(dataMSpos)

dataMSneg<-t(pathw_filt)
dataMSneg <- dataMSneg[rownames(dataMSneg) %in% metadataMSneg$id, ]
dataMSneg <- t(dataMSneg)
dataMSneg <- as.data.frame(dataMSneg)


  metadataMSpos$gc_treatment=as.factor(metadataMSpos$gc_treatment)
  metadataMSpos$category=as.factor(metadataMSpos$category)
  metadataMSpos$lesion_burden=as.factor(metadataMSpos$lesion_burden)
  levels(metadataMSpos$lesion_burden) =c("low","high")
  metadataMSpos$spinal_cord_lesion=as.factor(metadataMSpos$spinal_cord_lesion)
  levels(metadataMSpos$spinal_cord_lesion) =c("BM_low","BM_high")
  metadataMSpos$gadolinium_contrast=as.factor(metadataMSpos$gadolinium_contrast)
  levels(metadataMSpos$gadolinium_contrast) =c("NoActive","Active")
  metadataMSpos$subtentorial_lesions=as.factor(metadataMSpos$subtentorial_lesions)
  levels(metadataMSpos$subtentorial_lesions) =c("No","Yes")

  metadataMSneg$gc_treatment=as.factor(metadataMSneg$gc_treatment)
  metadataMSneg$category=as.factor(metadataMSneg$category)
  metadataMSneg$lesion_burden=as.factor(metadataMSneg$lesion_burden)
  levels(metadataMSneg$lesion_burden) =c("low","high")
  metadataMSneg$spinal_cord_lesion=as.factor(metadataMSneg$spinal_cord_lesion)
  levels(metadataMSneg$spinal_cord_lesion) =c("BM_low","BM_high")
  metadataMSneg$gadolinium_contrast=as.factor(metadataMSneg$gadolinium_contrast)
  levels(metadataMSneg$gadolinium_contrast) =c("NoActive","Active")
  metadataMSneg$subtentorial_lesions=as.factor(metadataMSneg$subtentorial_lesions)
  levels(metadataMSneg$subtentorial_lesions) =c("No","Yes")

res <- maaslin3(
  input_data = dataMS,       
  input_metadata = metadataMS,          
  output = "Output/Diff_humann/GC_THREATMENT/",
  formula = "~ gc_treatment",
  normalization = "TSS",
  transform = "LOG",
  min_prevalence = 0.01,
  min_abundance = 0.0001
)

res <- maaslin3(
  input_data = dataMSpos,       
  input_metadata = metadataMSpos,          
  output = "Output/Diff_humann/GC_THREATMENT_POS_GADOLINIUM/",
  formula = "~ gadolinium_contrast",
  normalization = "TSS",
  transform = "LOG",
  min_prevalence = 0.01,
  min_abundance = 0.0001
)
res <- maaslin3(
  input_data = dataMSneg,       
  input_metadata = metadataMSneg,          
  output = "Output/Diff_humann/GC_THREATMENT_NEG_GADOLINIUM/",
  formula = "~ gadolinium_contrast",
  normalization = "TSS",
  transform = "LOG",
  min_prevalence = 0.01,
  min_abundance = 0.0001
)

res <- maaslin3(
  input_data = dataMSneg,       
  input_metadata = metadataMSneg,          
  output = "Output/Diff_humann/GC_THREATMENT_NEG_SUBTENTORIAL/",
  formula = "~ subtentorial_lesions",
  normalization = "TSS",
  transform = "LOG",
  min_prevalence = 0.01,
  min_abundance = 0.0001
)
res <- maaslin3(
  input_data = dataMSpos,       
  input_metadata = metadataMSpos,          
  output = "Output/Diff_humann/GC_THREATMENT_POS_SUBTENTORIAL/",
  formula = "~ subtentorial_lesions",
  normalization = "TSS",
  transform = "LOG",
  min_prevalence = 0.01,
  min_abundance = 0.0001
)
res <- maaslin3(
  input_data = dataMSpos,       
  input_metadata = metadataMSpos,          
  output = "Output/Diff_humann/GC_THREATMENT_SPINAL/",
  formula = "~ spinal_cord_lesion",
  normalization = "TSS",
  transform = "LOG",
  min_prevalence = 0.01,
  min_abundance = 0.0001
)

res <- maaslin3(
  input_data = dataMSneg,       
  input_metadata = metadataMSneg,          
  output = "Output/Diff_humann/GC_THREATMENT_NEG_SPINAL/",
  formula = "~ spinal_cord_lesion",
  normalization = "TSS",
  transform = "LOG",
  min_prevalence = 0.01,
  min_abundance = 0.0001
)

res <- maaslin3(
  input_data = dataMSneg,       
  input_metadata = metadataMSneg,          
  output = "Output/Diff_humann/GC_THREATMENT_NEG_LESION_BURDEN/",
  formula = "~ lesion_burden",
  normalization = "TSS",
  transform = "LOG",
  min_prevalence = 0.01,
  min_abundance = 0.0001
)

res <- maaslin3(
  input_data = dataMSpos,       
  input_metadata = metadataMSpos,          
  output = "Output/Diff_humann/GC_THREATMENT_POS_LESION_BURDEN/",
  formula = "~ lesion_burden",
  normalization = "TSS",
  transform = "LOG",
  min_prevalence = 0.01,
  min_abundance = 0.0001
)
head(metadata)
head(data)

# remove UN--pathways

# 12197 103
 

# Remove rows with species contribution or unclassified species
pathw_filt <- pathways %>%
  filter(!grepl("\\|g__", Pathway)) %>%
  filter(!str_detect(Pathway, "unclassified"))