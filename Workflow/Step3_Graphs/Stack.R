source("utilities.R")
output_folder = "Output/NEW_PLOT_RDS/"
createFolder(output_folder)
baseline_decB = readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam.rds")
baseline_decA = readRDS(file = "Output/SUPERVISED_DEC/Archaea_Supervised_decontam.rds")
baseline_decE = readRDS(file = "Output/SUPERVISED_DEC/Eukaryota_Supervised_decontam.rds")

levelsB=c(  "Bifidobacterium adolescentis  ","Collinsella aerofaciens  ","Segatella copri  "
,"Phocaeicola dorei  ","Bifidobacterium longum  ","Faecalibacterium prausnitzii  ",
"Agathobacter rectalis  " ,"Bacteroides uniformis  ","Phocaeicola vulgatus  ", 
"Blautia wexlerae  ","Other" )
levelsA=c("Methanomicrobium antiquum  ", "Methanosarcina barkeri  ", 
"Methanococcus maripaludis  ", "Methanobrevibacter smithii  ", 
"Methanosphaera sp. BMS  ", "Methanimicrococcus sp. Hf6  ", "Methanobrevibacter sp. TLL-48-HuF1  ", 
"Methanolobus uncultured Methanolobus sp.  ", "Methanoregula uncultured Methanoregula sp.  ", 
"Haloferax volcanii  ","Other" )
levelsE=c( "Candida albicans  ", "Saccharomyces cerevisiae  ", "Aspergillus fumigatus  ", 
"Torulaspora globosa  ", "Malassezia japonica  ", "Saccharomyces mikatae  ", 
"Aspergillus ni
àdulans  ", "Penicillium oxalicum  ", 
"Schizosaccharomyces pombe  ", "Aspergillus puulaauensis  ", "Other")

#Domain="Bacteria"
stackbar <- function(baseline_dec, SP_levels,Domain, output_folder) {

  rel_abund = as.data.frame(abundances(baseline_dec, transform = "compositional"))
  baseline_dec = phyloseq(otu_table(rel_abund, taxa_are_rows = TRUE), tax_table(baseline_dec), sample_data(baseline_dec))

# Change category labels
  sample_data(baseline_dec)$category <- factor(sample_data(baseline_dec)$category, 
                                                        levels = c("HEALTHY", "MS"), 
                                                        labels = c("HD", "MS"))

  baselines_metadata = data.frame(sample_data(baseline_dec))


# Aggregate low represented phyla
  baseline_dc_species <- baseline_dec %>% microbiomeutilities::aggregate_top_taxa2(top=10, "Species")

# For species: merge genus and species names
  species = rownames(otu_table(baseline_dc_species))
  tax_df = as.data.frame(tax_table(baseline_dec))

  genus <- tax_df[tax_df$Species %in% species, c("Genus", "Species")]
  genus <- genus[match(species, genus$Species), ]

  genus_species = paste(genus[,"Genus"], species, " ")

  genus_species <- gsub("NA Other  ", "Other", genus_species)

# Create dataframe easily accessible by ggplot
  dat = as.data.frame(otu_table(baseline_dc_species))
  dat = dat %>%
    mutate(OTU_id = genus_species) %>%
    pivot_longer(-OTU_id, names_to = "id", values_to = "count")

# Add taxonomy
  tax = as.data.frame(tax_table(baseline_dc_species))
  tax = tax %>%
    mutate(OTU_id =genus_species)

# Merge counts with taxonomy information
  dat = dat %>%
    dplyr::left_join(tax, by = "OTU_id")

# Add metadata information
  dat = dat %>%
    dplyr::left_join(baselines_metadata, by = "id")


  dat$OTU_id=factor(dat$OTU_id, levels=SP_levels)
  palet_10=colors_10
  names(palet_10)=levels(dat$OTU_id)
  palet_10["Other"] = "#9b9a9a"
  pl=dat %>% 
    arrange(edss) %>%
    mutate(numeric_id= gsub("HD||MS", "", id)) %>%
    ggplot(aes(x= numeric_id , y = count)) + 
    facet_grid(~ category, scales = "free_x", space = "free_x") + 
    geom_bar(aes(fill = OTU_id), stat = "identity", position = "fill", width = 1) +
    xlab("Sample ID") + 
    scale_fill_manual(values=palet_10) +#inserire oter come unico colore
    scale_y_continuous(name = "Relative Abundance", 
                     labels = scales::percent) +
    theme_bw()+
    theme(legend.position="bottom",axis.text.x = element_text(angle = 0, size = 8),
        axis.text.y = element_text (color = "black",size=12),
        strip.text = element_text(face = "bold",size=20),
        strip.background = element_blank()) +
    guides(fill=guide_legend(title="Species"))
  pl
  saveRDS(pl, gsub(" ", "", paste(output_folder, Domain, "_StackedBar_Species.rds", sep = "")))
}

stackbar(baseline_decB, levelsB,"Bacteria", output_folder)
stackbar(baseline_decA, levelsA,"Archaea", output_folder)
stackbar(baseline_decE, levelsE,"Eukaryota", output_folder)
