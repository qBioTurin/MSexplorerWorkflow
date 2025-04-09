# Load necessary library
library(readr)

# Define the path to the TSV file
tab1<- read_tsv("Output/DAS_ONLY_LEFSE/ALL/BACT_BM_Lesion.res")
tab2<- read_tsv("Output/DAS_ONLY_LEFSE/ALL/BACT_GC.res")
tab3<-remove_common(tab1,tab2)
#remove_common<-
remove_common<- function(tab1,tab2){
tab1$Comparison <- paste(tab1$Genus, tab1$Species, sep = " ")
tab2$Comparison <- paste(tab2$Genus, tab2$Species, sep = " ")
# Subtract rows in tab1 that also appear in tab2
tab1 <- tab1[!tab1$Comparison %in% tab2$Comparison, ]
tab1$Comparison <- NULL
return(tab1)
}
# Load the RDS file
phyloseqAll<-readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontamOLD.rds")
load_taxaAll <- as.data.frame(phyloseqAll@tax_table)
load_otuAll <- as.data.frame(phyloseqAll@otu_table)
load_metadata <-as.data.frame(phyloseqAll@sam_data) 
patient <- load_metadata$id[load_metadata$gc_treatment %in% c("positive", "negative")]

load_otuAll <- load_otuAll[, colnames(load_otuAll) %in% patient]
