# Load RDS files
arch <- readRDS("Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.001.rds")
euk <- readRDS("Output/SUPERVISED_DEC/Eukaryote_Supervised_decontam0.001.rds")
bact <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
library(phyloseq)
taxa_arch<-tax_table(arch)
taxa_euk<-tax_table(euk)
taxa_bact<-tax_table(bact)
metadata<-sample_data(arch)

arch_prev <- readRDS("Output/RDS/Archaea_baseline_phylo.rds")
euk_prev <- readRDS("Output/RDS/Eukaryotes_baseline_phylo.rds")
bact_prev <- readRDS("Output/RDS/Bact_baseline_phylo.rds")
print(arch)
print(euk)
print(bact)
print(arch_prev)
print(euk_prev)
print(bact_prev)
euk_taxa<-tax_table(euk)
bact_taxa<-tax_table(bact)
arch_taxa<-tax_table(arch)

otu<- as.data.frame(otu_table(euk))
otu_bact<- as.data.frame(otu_table(bact))


bact_old<- readRDS("Output/RDS/Bact_baseline_phylo.rds")

# Append first 4 columns from merged_shannon_union.tsv to the bottom of merged_shannon_union_old.tsv
new_tsv <- read_tsv("Output/DAS_ALPHA_MAASLIN3/merged_shannon_union.tsv")
old_tsv <- read_tsv("Output/DAS_ALPHA_MAASLIN3/merged_shannon_union_old.tsv", col_names = FALSE)

rownames(old_tsv)<-old_tsv$X1

data1 <- read_csv("InputData/cluster_mapping_2026-03-20.csv")
data2<- read_csv("InputData/Metadata_MS.csv")

data1$id <- as.character(data1$id)
data2$id <- as.character(data2$id)
merged_data <- merge(data1, data2, by = "id", all = TRUE)
write.csv(merged_data, "InputData/merged_data_cluster.csv", row.names = FALSE)


tsv_file <- read_csv("InputData/merged_data_cluster.csv") %>%
    mutate(Cluster = ifelse(Cluster == 1, "bad_prognosis", ifelse(Cluster == 2, "good_prognosis", Cluster)))
write.csv(tsv_file, "InputData/merged_data_cluster_name.csv", row.names = FALSE)


gado<-rownames(tax_table(readRDS("Output/BOI/001/Bacteria_gadolinium_001.rds")))
lesion<-rownames(tax_table(readRDS("Output/BOI/001/Bacteria_lesion_001.rds")))
spinal_cord<-rownames(tax_table(readRDS("Output/BOI/001/Bacteria_spinal_cord_001.rds")))
subtentorial<-rownames(tax_table(readRDS("Output/BOI/001/Bacteria_subtentorial_lesions_001.rds")))
list_tax<-unique(c(gado, lesion, spinal_cord, subtentorial))
###46taxa