#BiocManager::install("biobakery/maaslin3")
library(maaslin3)
# Load an RDS file in R



Maaslin3 <- function(rds,output_folder,formulaInUse,categorySel, elements) {
# Create data for maaslin3
otu <- rds@otu_table
taxa <- rds@tax_table
metadata_df <- as(rds@sam_data, "data.frame")
ids_to_retain <- rownames(metadata_df[metadata_df[[categorySel]] %in% elements, ])
otu <- otu[, colnames(otu) %in% ids_to_retain]
metadata_df <- metadata_df[rownames(metadata_df) %in% ids_to_retain, ]
taxa <- as.data.frame(taxa)
taxa$Genus_Species <- paste(taxa$Genus, taxa$Species, sep = "_")
Tnew <- taxa %>% select(Genus_Species)
Tnew$taxid <- rownames(Tnew)

otu <- as.data.frame(otu)
otu$taxid <- rownames(otu)
# Merge the two data frames by taxid
merged_df <- merge(otu, Tnew, by = "taxid")
rownames(merged_df) <- merged_df$Genus_Species
merged_df <- merged_df %>% select(-c(taxid, Genus_Species))
merged_df <- t(merged_df)

# Create metadata for masslin3

class(metadata_df)
metadata_df <- metadata_df %>%
    select(category, gc_treatment, lesion_burden, spinal_cord_lesion, gadolinium_contrast, subtentorial_lesions)

metadata_df$gc_treatment <- as.factor(metadata_df$gc_treatment)
metadata_df$category <- as.factor(metadata_df$category)
metadata_df$lesion_burden <- as.factor(metadata_df$lesion_burden)
levels(metadata_df$lesion_burden) <- c("low", "high")
metadata_df$spinal_cord_lesion <- as.factor(metadata_df$spinal_cord_lesion)
levels(metadata_df$spinal_cord_lesion) <- c("BM_low", "BM_high")
metadata_df$gadolinium_contrast <- as.factor(metadata_df$gadolinium_contrast)
levels(metadata_df$gadolinium_contrast) <- c("NoActive", "Active")
metadata_df$subtentorial_lesions <- as.factor(metadata_df$subtentorial_lesions)
levels(metadata_df$subtentorial_lesions) <- c("No", "Yes")

if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
}
formulaInUse <- paste("~", formulaInUse, sep = "")
fit_out <- maaslin3(input_data = merged_df,
                    input_metadata = metadata_df,
                    output = output_folder,
                    formula = formulaInUse,
                    normalization = 'TSS',
                    transform = 'LOG',
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 250,
                    cores = 1)
}

rds_001 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")

Maaslin3(rds = rds_001, output_folder = "Output/MAASLIN3/001/001_lesion_burden_onlygc",
 formulaInUse = "lesion_burden", categorySel = "gc_treatment", elements = c("positive", "negative"))

Maaslin3(rds = rds_001, output_folder = "Output/MAASLIN3/001/001_spinal_cord_lesion_onlygc",
 formulaInUse = "spinal_cord_lesion", categorySel = "gc_treatment", elements = c("positive", "negative"))

Maaslin3(rds = rds_001, output_folder = "Output/MAASLIN3/001/001_gadolinium_contrast_onlygc",
 formulaInUse = "gadolinium_contrast", categorySel = "gc_treatment", elements = c("positive", "negative"))

Maaslin3(rds = rds_001, output_folder = "Output/MAASLIN3/001/001_subtentorial_lesions_onlygc",
 formulaInUse = "subtentorial_lesions", categorySel = "gc_treatment", elements = c("positive", "negative"))

rds_01 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds")

Maaslin3(rds = rds_01, output_folder = "Output/MAASLIN3/01/01_lesion_burden_onlygc",
 formulaInUse = "lesion_burden", categorySel = "gc_treatment", elements = c("positive", "negative"))

Maaslin3(rds = rds_01, output_folder = "Output/MAASLIN3/01/01_spinal_cord_lesion_onlygc",
 formulaInUse = "spinal_cord_lesion", categorySel = "gc_treatment", elements = c("positive", "negative"))

Maaslin3(rds = rds_01, output_folder = "Output/MAASLIN3/01/01_gadolinium_contrast_onlygc",
 formulaInUse = "gadolinium_contrast", categorySel = "gc_treatment", elements = c("positive", "negative"))

Maaslin3(rds = rds_01, output_folder = "Output/MAASLIN3/01/01_subtentorial_lesions_onlygc",
 formulaInUse = "subtentorial_lesions", categorySel = "gc_treatment", elements = c("positive", "negative"))


rds_05 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds")

Maaslin3(rds = rds_05, output_folder = "Output/MAASLIN3/05/05_lesion_burden_onlygc",
 formulaInUse = "lesion_burden", categorySel = "gc_treatment", elements = c("positive", "negative"))

Maaslin3(rds = rds_05, output_folder = "Output/MAASLIN3/05/05_spinal_cord_lesion_onlygc",
 formulaInUse = "spinal_cord_lesion", categorySel = "gc_treatment", elements = c("positive", "negative"))

Maaslin3(rds = rds_05, output_folder = "Output/MAASLIN3/05/05_gadolinium_contrast_onlygc",
 formulaInUse = "gadolinium_contrast", categorySel = "gc_treatment", elements = c("positive", "negative"))

Maaslin3(rds = rds_05, output_folder = "Output/MAASLIN3/05/05_subtentorial_lesions_onlygc",
 formulaInUse = "subtentorial_lesions", categorySel = "gc_treatment", elements = c("positive", "negative"))
