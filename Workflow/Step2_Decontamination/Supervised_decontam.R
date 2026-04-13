source("Settings/utilities.R")
output_folder = "Output/SUPERVISED_DEC/"
output_folder_prev = "Output/PREVALENCE/"
createFolder(output_folder)
createFolder(output_folder_prev)

decontam <- function(baselines, Domain, output_folder, dec_level,prevalence = FALSE,contaminant=NULL) {
   
   # Remove controls
    baselines = subset_samples(baselines, sample_type == "sample")
    otu = data.frame(otu_table(baselines))
    tax = data.frame(tax_table(baselines))

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
    if(!is.null(contaminant)){

        baselines <- prune_taxa(!(taxa_names(baselines) %in% contaminant), baselines)
    }
    # # Remove microorganisms not from dataset
    # if (level == "Phylum") {
    #     supervised_tax = tax[tax$Phylum %in% decontam_list, ]
    # } else if (level == "Genus") {
    #     supervised_tax = tax[tax$Genus %in% decontam_list, ]
    # }

    #supervised = prune_taxa(rownames(supervised_tax), baselines)
    supervised <- baselines
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
    prova1<-final_taxa

    dec_level_ch <- as.character(dec_level)
    final_name = gsub(" ", "", paste(output_folder, Domain, "_Supervised_decontam",dec_level_ch,".rds"))
    
    saveRDS(final, file = final_name)
}

execute_sup_decontam <- function(){
    contE=c(460519)
    bact_baselines_ds = readRDS(file = "Output/DESEQ_RDS/Bact_DeSeq.rds")
    archaea_baselines_ds = readRDS(file = "Output/DESEQ_RDS/Archaea_DeSeq.rds")
    euk_baselines_ds = readRDS(file = "Output/DESEQ_RDS/Eukaryote_DeSeq.rds")

    decontam(baselines=bact_baselines_ds, Domain="Bacteria", output_folder=output_folder, dec_level=0.05)
    decontam(baselines=bact_baselines_ds, Domain="Bacteria", output_folder=output_folder, dec_level=0.01)
    decontam(baselines=bact_baselines_ds, Domain="Bacteria", output_folder=output_folder, dec_level=0.001)
    decontam(baselines=bact_baselines_ds, Domain="Bacteria", output_folder=output_folder, dec_level=0)
    decontam(baselines=archaea_baselines_ds, Domain="Archaea", output_folder=output_folder, dec_level=0)
    decontam(baselines=archaea_baselines_ds, Domain="Archaea", output_folder=output_folder, dec_level=0.001)
    decontam(baselines=euk_baselines_ds, Domain="Eukaryote", output_folder=output_folder, dec_level=0, contaminant = contE)
    decontam(baselines=euk_baselines_ds, Domain="Eukaryote", output_folder=output_folder, dec_level=0.001, contaminant = contE)
}

execute_sup_decontam()
prime<-readRDS("Output/DESEQ_RDS/Eukaryote_DeSeq.rds")
taxa<-as.data.frame(phyloseq::tax_table(prime))
