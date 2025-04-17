source("Settings/utilities.R")
output_folder = "Output/SUPERVISED_DEC/"
createFolder(output_folder)

decontam <- function(baselines, decontam_list, level, Domain, output_folder, dec_level) {
   
   # Remove controls
    baselines = subset_samples(baselines, sample_type == "sample")
    otu = data.frame(otu_table(baselines))
    tax = data.frame(tax_table(baselines))

    # Remove features with ambiguous phylum/genus annotation (NAs)
    if (level == "Phylum") {
        baselines_ds <- subset_taxa(baselines, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized", "NA"))
    } else if (level == "Genus") {
        baselines_ds <- subset_taxa(baselines, !is.na(Genus) & !Genus %in% c("", "uncharacterized", "NA"))
    }
    table(tax_table(baselines)[, level], exclude = NULL)

    # Compute prevalence of each feature (the number of samples in which a taxon appears at least once)
    Prev = apply(
        X = otu_table(baselines),
        MARGIN = ifelse(taxa_are_rows(baselines), yes = 1, no = 2),
        FUN = function(x) { sum(x > 0) }
    )
    # Adds taxonomy and total read counts, and stores as data frame
    Prev = data.frame(
        Prevalence = Prev,
        TotalAbundance = taxa_sums(baselines),
        tax_table(baselines)
    )

    # Remove microorganisms not from dataset
    if (level == "Phylum") {
        supervised_tax = tax[tax$Phylum %in% decontam_list, ]
    } else if (level == "Genus") {
        supervised_tax = tax[tax$Genus %in% decontam_list, ]
    }

    supervised = prune_taxa(rownames(supervised_tax), baselines)
    supervised_otu = data.frame(otu_table(supervised))

    # Relative Abundance filter (>0.001)
    ######################################
    # Convert to relative abundance
    rel = data.frame(abundances(supervised, transform = "compositional"))

    # Assign 0 to all OTUs that have relative abundance < 0.001
    threshold = dec_level
    rel[rel < threshold] <- 0

    # Remove all OTUs that have 0 counts in all columns
    rel <- rel %>%
        mutate(sum = rowSums(across(where(is.numeric)))) %>%
        filter(sum > 0)

    final = prune_taxa(rownames(rel), supervised)
    final_taxa = data.frame(tax_table(final))
    final_otu = data.frame(otu_table(final))
    table(final_taxa$Phylum)
    dec_level_ch <- as.character(dec_level)
    final_name = gsub(" ", "", paste(output_folder, Domain, "_Supervised_decontam",dec_level_ch,".rds"))
    saveRDS(final, file = final_name)
}

execute_sup_decontam <- function() {
    # Generale usando script prevalence
    ####################################

    # Manual selection Bacteria
    ##################
    # Manual selection using Encyclopedia of Life (https://eol.org/), Taxonomy Browser and Pubmed (Phylum name AND human)
    # Phyla associated with human:
    # Actinomycetota/actinobacteria
    # Bacillota/Firmicutes
    # Bacteroidota/Bacteroidetes
    # Campylobacterota
    # Chlamydiota
    # Chloroflexota
    # Lentisphaerota
    # Mycoplasmatota/Mycoplasma/Tenericutes
    # Pseudomonadota
    # Spirochaetota
    # Synergistota
    # Verrucomicrobiota

    # Phyla associated with human in only few ref:
    # Acidobacteriota
    # Candidatus Saccharibacteria
    # Elusimicriobiota
    # Gemmatimonadota
    # Ignavibacteriota
    # Myxococcota
    # Nitrospirota
    # Planctomycetota
    # Thermodesulfobacteriota

    human_phylaB = c(
        "Actinomycetota", "Bacillota", "Bacteroidota", "Campylobacterota",
        "Chlamydiota", "Chloroflexota", "Lentisphaerota", "Mycoplasmatota",
        "Pseudomonadota", "Spirochaetota", "Synergistota", "Verrucomicrobiota",
        "Acidobacteriota", "Candidatus Saccharibacteria", "Elusimicriobiota",
        "Gemmatimonadota", "Ignavibacteriota", "Myxococcota", "Nitrospirota",
        "Planctomycetota", "Thermodesulfobacteriota"
    )

    # Manual selection
    ##################
    # Manual selection using Encyclopedia of Life (https://eol.org/), Taxonomy Browser and Pubmed (Phylum name AND human)
    # Phyla associated with human:
    # Euryarchaeota
    # Thermoproteota

    # Phyla associated with human in only few ref:
    # Candidatus Korarchaeota
    # Nanoarchaeota

    human_phylaA = c("Euryarchaeota", "Thermoproteota", "Candidatus Korarchaeota", "Nanoarchaeota")

    # Manual selection
    # For this we consulted an expert of this field (there are very few info on them in the literature and for this we decided to use the genus and not phyla)
    human_genus = c(
        "Acremonium", "Agaricales", "Agaricomycetes", "Aspergillus",
        "Aureobasidium", "Byssochlamys", "Candida", "Cladosporium",
        "Claroideoglomus", "Clavispora", "Cyberlindnera", "Debaryomyces",
        "Geotrichum", "Hanseniaspora", "Hypocreales", "Itersonilia", "Kazachstania",
        "Kazchstania", "Kluyveromyces", "Kurtzmaniella", "Lachancea", "Lecanoromycetes",
        "Malassezia", "Malasseziales", "Metschnikowia", "Ogataea", "Penicillium", "Pichia",
        "Polyporales", "Rhizoplaca", "Rhodotorula", "Saccharomycetales", "Schizosaccharomyces",
        "Starmerella", "Torulaspora", "Trametes", "Vishniacozyma", "Wickerhamomyces",
        "Yamadazyma", "Yarrowia", "Saccharomyces"
    )

    bact_baselines_ds = readRDS(file = "Output/DESEQ_RDS/Bact_DeSeq.rds")
    archaea_baselines_ds = readRDS(file = "Output/DESEQ_RDS/Archaea_DeSeq.rds")
    euk_baselines_ds = readRDS(file = "Output/DESEQ_RDS/Eukaryote_DeSeq.rds")

    decontam(bact_baselines_ds, human_phylaB, "Phylum", "Bacteria", output_folder, 0.05)
    decontam(bact_baselines_ds, human_phylaB, "Phylum", "Bacteria", output_folder, 0.01)
    decontam(bact_baselines_ds, human_phylaB, "Phylum", "Bacteria", output_folder, 0.001)
    decontam(archaea_baselines_ds, human_phylaA, "Phylum", "Archaea", output_folder, 0.001)
    decontam(euk_baselines_ds, human_genus, "Genus", "Eukaryote", output_folder, 0.001)
}

execute_sup_decontam()
