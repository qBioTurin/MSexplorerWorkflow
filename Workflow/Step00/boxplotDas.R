library(ggplot2)




create_boxplot<-function(rds1,name){
taxa1 = as.data.frame(rds1@tax_table)
taxa_names = rownames(taxa1)
rds_bact = readRDS("Output/RDS/merged_phylo_renamed.rds")
counts = as.data.frame(rds_bact@otu_table)
if (!phyloseq::taxa_are_rows(rds_bact)) {
    counts = t(counts)
}
sample_data <- as.data.frame(as.matrix(phyloseq::sample_data(rds_bact)))
group <- sample_data$category  # Modifica "category" se il nome è diverso
taxa<-as.data.frame(phyloseq::tax_table(rds_bact))

taxa$genusSpecies <- paste(taxa$Genus, taxa$Species, sep = " ")
# Filtra solo i campioni con gruppo "HEALTHY" o "MS"
valid_idx <- which(group %in% c("HEALTHY", "MS"))
counts <- counts[, valid_idx, drop = FALSE]
group <- group[valid_idx]
filename=paste0(name,".pdf")
pdf(file = file.path("Image/boxplot/", filename), width = 8, height = 6)
for (taxon in taxa_names) {
    if (taxon %in% rownames(counts)) {
        df <- data.frame(
            Count = as.numeric(counts[taxon, ]),
            Group = group
        )

        p <- ggplot(df, aes(x = Group, y = Count, fill = Group)) +
            geom_boxplot() +
            geom_signif(
                comparisons = list(c("HEALTHY", "MS")),
                map_signif_level = function(p) {
                    if (p < 0.001) {
                        return(paste0("*** (", signif(p, 2), ")"))
                    } else if (p < 0.01) {
                        return(paste0("** (", signif(p, 2), ")"))
                    } else if (p < 0.049) {
                        return(paste0("* (", signif(p, 2), ")"))
                    } else if (p < 0.055) {
                        return(paste0("NS. (", signif(p, 2), ")"))
                    } else {
                        return(paste0("NS."))
                    }
                }
            ) +
            theme_classic() +
            theme(
                axis.title.x = element_text(size = 15),
                axis.title.y = element_text(size = 15),
                axis.text.x = element_text(size = 15),
                axis.text.y = element_text(size = 15),
                plot.title = element_text(size = 20),
                legend.text = element_text(size = 15),
                legend.title = element_text(size = 15),
                strip.text = element_text(size = 20)
            )
        genus_species <- taxa[rownames(taxa) == taxon, "genusSpecies"]
        p <- p + ggtitle(genus_species)
        print(p)
    }
}
dev.off()
}

create_boxplot( readRDS("Output/merge_DAS/01/Bacteria_gadolinium_contrast_01_merged.rds"),"gadolinium_contrast_01")
create_boxplot(readRDS())


# Export the boxplot data to CSV for further analysis
export_boxplot_data <- function(rds1, name) {
    taxa1 = as.data.frame(rds1@tax_table)
    taxa_names = rownames(taxa1)
    rds_bact = readRDS("Output/RDS/merged_phylo_renamed.rds")
    counts = as.data.frame(rds_bact@otu_table)
    if (!phyloseq::taxa_are_rows(rds_bact)) {
        counts = t(counts)
    }
    sample_data <- as.data.frame(as.matrix(phyloseq::sample_data(rds_bact)))
    group <- sample_data$category
    taxa <- as.data.frame(phyloseq::tax_table(rds_bact))
    taxa$genusSpecies <- paste(taxa$Genus, taxa$Species, sep = " ")
    valid_idx <- which(group %in% c("HEALTHY", "MS"))
    counts <- counts[, valid_idx, drop = FALSE]
    group <- group[valid_idx]

    all_data <- data.frame()
    for (taxon in taxa_names) {
        if (taxon %in% rownames(counts)) {
            df <- data.frame(
                Taxon = taxon,
                GenusSpecies = taxa[rownames(taxa) == taxon, "genusSpecies"],
                Sample = colnames(counts),
                Count = as.numeric(counts[taxon, ]),
                Group = group
            )
            all_data <- rbind(all_data, df)
        }
    }
    write.csv(all_data, file = file.path("Image/boxplot/", paste0(name, "_boxplot_data.csv")), row.names = FALSE)
}

export_boxplot_data(readRDS("Output/merge_DAS/01/Bacteria_gadolinium_contrast_01_merged.rds"), "gadolinium_contrast_01")



p1<-read_csv("InputData/rename.csv")
name<-p1$old_id
name <- sort(name)
