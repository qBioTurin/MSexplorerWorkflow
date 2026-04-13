# Load three RDS files
otu1 <- readRDS("Output/merge_DAS/MSHD/Archaea_MsHd_merged.rds")
otu2 <- readRDS("Output/merge_DAS/MSHD/Bacteria_MsHd_merged.rds")
otu3 <- readRDS("Output/merge_DAS/MSHD/Eukaryote_MsHd_merged.rds")
taxa<-readRDS("Output/DESEQ_RDS/ALL_DeSeq.rds")
taxa_all <- as.data.frame(taxa@tax_table)
taxa_all$Genus_sp <- paste(taxa_all$Genus, taxa_all$Species, sep = " ")
taxa_all <- taxa_all["Genus_sp"]
taxa_all$taxid <- rownames(taxa_all)
otu1 <- as.data.frame(otu1@otu_table)
otu2 <- as.data.frame(otu2@otu_table)
otu3 <- as.data.frame(otu3@otu_table)
otu1$taxid <- rownames(otu1)
otu2$taxid <- rownames(otu2)
otu3$taxid <- rownames(otu3)
# Merge the three OTU tables by row names (assuming same samples/columns)
taxa1 <- merge(otu1, taxa_all, by = "taxid", all.x = TRUE)
taxa2 <- merge(otu2, taxa_all, by = "taxid", all.x = TRUE)
taxa3 <- merge(otu3, taxa_all, by = "taxid", all.x = TRUE)

# Move taxid column first and Genus_sp column second
move_genus_sp <- function(df) {
    cols <- colnames(df)
    taxid_idx <- which(cols == "taxid")
    genus_sp_idx <- which(cols == "Genus_sp")
    other_idx <- setdiff(seq_along(cols), c(taxid_idx, genus_sp_idx))
    new_order <- c(cols[taxid_idx], cols[genus_sp_idx], cols[other_idx])
    df <- df[, new_order]
    return(df)
}

taxa1 <- move_genus_sp(taxa1)
taxa2 <- move_genus_sp(taxa2)
taxa3 <- move_genus_sp(taxa3)

write.csv(taxa1, file = "Output/merge_DAS/MSHD/Archaea_MsHd_merged.csv", row.names = FALSE)
write.csv(taxa2, file = "Output/merge_DAS/MSHD/Bacteria_MsHd_merged.csv", row.names = FALSE)
write.csv(taxa3, file = "Output/merge_DAS/MSHD/Eukaryote_MsHd_merged.csv", row.names = FALSE)

# Example: merge all OTU tables (if needed)
otu_merged <- rbind(taxa1, taxa2, taxa3)
# Replace NA with 0 (optional, if you want to fill missing values)
otu_merged[is.na(otu_merged)] <- 0

# Assign merged table back to otu1 for further use
otu1 <- otu_merged
print(otu1)
# Print summaries to check
summary(otu1)
summary(otu2)
summary(otu3)