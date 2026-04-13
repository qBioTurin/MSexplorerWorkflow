
rds <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds")
otu <- rds@otu_table
taxa <- rds@tax_table

metadata_df <- as(rds@sam_data, "data.frame")

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

merged_df <- merged_df[rownames(merged_df) %in% c("Phocaeicola_vulgatus", "Bacteroides_ovatus", "Bacteroides_xylanisolvens", "Alistipes_onderdonkii", "Roseburia_intestinalis"), ]


t_abundance <- t(merged_df)

t_abundance_scaled <- scale(t_abundance)

pca_result <- prcomp(t_abundance_scaled, center = TRUE, scale. = TRUE)

pca_df <- as.data.frame(pca_result$x)
pca_df$SampleID <- rownames(pca_df)

metadata_df$SampleID <- rownames(metadata_df)  # Assicurati che SampleID sia colonna
pca_df <- merge(pca_df, metadata_df[, c("SampleID", "gc_treatment")], by = "SampleID")


ggplot(pca_df, aes(x = PC1, y = PC2, color = gc_treatment)) +
  geom_point(size = 3) +
  labs(title = "PCA of Abundance Data by GC Treatment",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)"),
       color = "GC Treatment") +
  theme_minimal() +
  theme(text = element_text(size = 14))
ggsave("Output/PCA/PCA_Abundance_GC_Treatment_5.jpg", width = 10, height = 8)
