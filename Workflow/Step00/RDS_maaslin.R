phylo<-readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
species<-read.csv("Output/MAASLIN3/Rebuttal/Bact_001/category_MAASLIN3_SR_features_K_T.csv")
taxa<-tax_table(phylo)
taxa<-as.data.frame(taxa)
taxa$genus_species<-paste(taxa$Genus,taxa$Species,sep="_")
matches <- rownames(taxa)[taxa$genus_species %in% species$Species]
phylo <- prune_taxa(matches, phylo)
print(phylo)
createFolder("Output/rebuttal")
saveRDS(phylo, "Output/rebuttal/Bacteria_category_heatmap.rds")




gado_old<-readRDS("InputData/old_boi/Bacteria_gadolinium_001.rds")
lesion_old<-readRDS("InputData/old_boi/Bacteria_lesion_001.rds")
subt_old<-readRDS("InputData/old_boi/Bacteria_subtentorial_lesions_001.rds")
spinal_old<-readRDS("InputData/old_boi/Bacteria_spinal_cord_001.rds")


gado<-readRDS("Output/BOI/001/Bacteria_gadolinium_001.rds")
lesion<-readRDS("Output/BOI/001/Bacteria_lesion_001.rds")
spinal_cord<-readRDS("Output/BOI/001/Bacteria_spinal_cord_001.rds")
subtentorial<-readRDS("Output/BOI/001/Bacteria_subtentorial_lesions_001.rds")

extract_species <- function(x) {
    if (is.null(x)) return(character(0))
    if (is.character(x)) return(x)
    if (is.factor(x)) return(as.character(x))

    build_genus_species <- function(tbl) {
        if (is.null(tbl)) return(NULL)
        df <- as.data.frame(tbl)
        cols <- colnames(df)
        if (is.null(cols)) return(NULL)

        genus_col <- cols[tolower(cols) == "genus"][1]
        species_col <- cols[tolower(cols) == "species"][1]
        if (is.na(genus_col) || is.na(species_col)) return(NULL)

        genus <- trimws(as.character(df[[genus_col]]))
        species <- trimws(as.character(df[[species_col]]))
        genus[is.na(genus)] <- ""
        species[is.na(species)] <- ""

        gs <- paste(genus, species, sep = " ")
        gs <- gsub("_+$", "", gs)
        gs <- gs[gs != ""]
        if (length(gs) == 0) return(NULL)
        unique(gs)
    }

    get_taxa_names <- function(obj) {
        tryCatch({
            if (exists("taxa_names", mode = "function")) {
                taxa_names(obj)
            } else if (requireNamespace("phyloseq", quietly = TRUE)) {
                phyloseq::taxa_names(obj)
            } else {
                NULL
            }
        }, error = function(e) NULL)
    }

    get_tax_table <- function(obj) {
        tryCatch({
            if (exists("tax_table", mode = "function")) {
                tax_table(obj)
            } else if (requireNamespace("phyloseq", quietly = TRUE)) {
                phyloseq::tax_table(obj)
            } else {
                NULL
            }
        }, error = function(e) NULL)
    }

    if (inherits(x, "phyloseq")) {
        gs <- build_genus_species(get_tax_table(x))
        if (!is.null(gs)) return(as.character(gs))
        return(as.character(get_taxa_names(x)))
    }
    if (isS4(x)) {
        gs <- build_genus_species(get_tax_table(x))
        if (!is.null(gs)) return(as.character(gs))
        tx <- get_taxa_names(x)
        if (!is.null(tx)) return(as.character(tx))
    }
    if (is.data.frame(x) || is.matrix(x)) {
        gs <- build_genus_species(x)
        if (!is.null(gs)) return(as.character(gs))
        rn <- rownames(x)
        if (!is.null(rn)) return(as.character(rn))
        return(as.character(x[, 1]))
    }
    if (is.atomic(x)) return(as.character(x))
    as.character(unlist(x, use.names = FALSE))
}

list_tax <- unique(c(
    extract_species(gado),
    extract_species(lesion),
    extract_species(spinal_cord),
    extract_species(subtentorial)
))
###46taxa

# Create comparison tables for each group
create_comparison <- function(old_data, new_data, group_name) {
    old_species <- unique(extract_species(old_data))
    new_species <- unique(extract_species(new_data))

    common <- intersect(old_species, new_species)
    only_old <- setdiff(old_species, new_species)
    only_new <- setdiff(new_species, old_species)
    
    result <- list(
        common = data.frame(Species = common, Type = "Common"),
        only_old = data.frame(Species = only_old, Type = paste("OLD")),
        only_new = data.frame(Species = only_new, Type = paste("NEW"))
    )
    return(result)
}

gado_comp <- create_comparison(gado_old, gado, "Gado")
print(gado_comp$common)
lesion_comp <- create_comparison(lesion_old, lesion, "Lesion")
spinal_comp <- create_comparison(spinal_old, spinal_cord, "Spinal")
subt_comp <- create_comparison(subt_old, subtentorial, "Subt")

# Create summary tables
gado_table <- rbind(gado_comp$common, gado_comp$only_old, gado_comp$only_new)
lesion_table <- rbind(lesion_comp$common, lesion_comp$only_old, lesion_comp$only_new)
spinal_table <- rbind(spinal_comp$common, spinal_comp$only_old, spinal_comp$only_new)
subt_table <- rbind(subt_comp$common, subt_comp$only_old, subt_comp$only_new)

# Union of all old and all new
all_old <- unique(c(
    extract_species(gado_old),
    extract_species(lesion_old),
    extract_species(spinal_old),
    extract_species(subt_old)
))
all_new <- unique(c(
    extract_species(gado),
    extract_species(lesion),
    extract_species(spinal_cord),
    extract_species(subtentorial)
))
common_all <- intersect(all_old, all_new)
only_old_all <- setdiff(all_old, all_new)
only_new_all <- setdiff(all_new, all_old)

union_table <- rbind(
    data.frame(Species = common_all, Type = "Common"),
    data.frame(Species = only_old_all, Type = "OLD"),
    data.frame(Species = only_new_all, Type = "NEW")
)

write.table(gado_table, "Output/rebuttal/gado_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lesion_table, "Output/rebuttal/lesion_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(spinal_table, "Output/rebuttal/spinal_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(subt_table, "Output/rebuttal/subt_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(union_table, "Output/rebuttal/union_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
