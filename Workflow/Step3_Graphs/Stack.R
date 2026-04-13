source("Settings/utilities.R")
source("Settings/Packages.R")
source("Settings/Colorpalette.R")
output_folder <- "Output/NEW_PLOT_RDS/"
createFolder(output_folder)

stackbar <- function(baseline_dec, Domain, output_folder) {
  tax <- as.data.frame(phyloseq::tax_table(baseline_dec))
  baseline_dec <- phyloseq::prune_taxa(phyloseq::taxa_sums(baseline_dec) > 0, baseline_dec)
  # Save full sample list before pruning empty samples
  full_sample_data <- data.frame(phyloseq::sample_data(baseline_dec))
  full_sample_data$id <- rownames(full_sample_data)
  baseline_dec <- phyloseq::prune_samples(phyloseq::sample_sums(baseline_dec) > 0, baseline_dec)

  rel_abund <- as.data.frame(microbiome::abundances(baseline_dec, transform = "compositional"))
  baseline_dec <- phyloseq::phyloseq(
    phyloseq::otu_table(rel_abund, taxa_are_rows = TRUE),
    phyloseq::tax_table(baseline_dec),
    phyloseq::sample_data(baseline_dec)
  )
  # Change category labels
  phyloseq::sample_data(baseline_dec)$category <- factor(phyloseq::sample_data(baseline_dec)$category,
    levels = c("HEALTHY", "MS"),
    labels = c("HD", "MS")
  )

  baselines_metadata <- data.frame(phyloseq::sample_data(baseline_dec))

  # Update Species column in tax_table as "Genus Species"
  tax_upd <- as.data.frame(phyloseq::tax_table(baseline_dec))
  tax_upd$Species <- rownames(tax_upd)
  phylo2 <- baseline_dec
  phyloseq::tax_table(phylo2) <- phyloseq::tax_table(as.matrix(tax_upd))
  tax_upd2 <- as.data.frame(phyloseq::tax_table(baseline_dec))
  tax_upd2$Species <- paste(tax_upd2$Genus, tax_upd2$Species, sep = " ")
  phyloseq::tax_table(baseline_dec) <- phyloseq::tax_table(as.matrix(tax_upd2))

  # Aggregate low represented phyla
  baseline_dc_species <- phylo2 %>% microbiomeutilities::aggregate_top_taxa2(top = 10, "Species")

  # For species: merge genus and species names
  species <- rownames(phyloseq::otu_table(baseline_dc_species))
  tax_df <- as.data.frame(phyloseq::tax_table(baseline_dec))
  genus <- tax_df[rownames(tax_df) %in% species, c("Genus", "Species")]
  genus <- genus[match(species, rownames(genus)), ]


  genus_species <- trimws(genus$Species)

  genus_species[is.na(genus_species) | genus_species == ""] <- "Other"
  genus_species <- gsub("^NA Other$|^NA NA$", "Other", genus_species)

  save_array <- genus_species
  if ("Other" %in% genus_species) {
    save_array <- c(save_array[save_array != "Other"], "Other")
  } else {
    save_array <- save_array[save_array != "Other"]
  }

  # Create dataframe easily accessible by ggplot
  dat <- as.data.frame(phyloseq::otu_table(baseline_dc_species))
  dat <- dat %>%
    mutate(OTU_id = genus_species) %>%
    pivot_longer(-OTU_id, names_to = "id", values_to = "count")

  # Add taxonomy
  tax <- as.data.frame(phyloseq::tax_table(baseline_dc_species))
  tax <- tax %>%
    mutate(OTU_id = genus_species)

  # Merge counts with taxonomy information
  dat <- dat %>%
    dplyr::left_join(tax, by = "OTU_id")

  # Add metadata information
  dat <- dat %>%
    dplyr::left_join(baselines_metadata, by = "id")

  # dat <- dat %>%
  #  dplyr::filter(!is.na(count) & count > 0)


  dat$OTU_id <- factor(dat$OTU_id, levels = save_array)
  otu_levels <- levels(dat$OTU_id)
  palet_10 <- rep(colors_10, length.out = length(otu_levels))
  names(palet_10) <- otu_levels
  if ("Other" %in% names(palet_10)) palet_10["Other"] <- "#9b9a9a"
  dat_plot <- dat %>%
    dplyr::mutate(
      numeric_id = gsub("^(HD|MS)", "", id),
      numeric_id = ifelse(nchar(numeric_id) == 0, id, numeric_id),
      numeric_order = suppressWarnings(as.numeric(numeric_id))
    )

  # Add missing (empty) samples so they still appear on x-axis
  missing_ids <- setdiff(full_sample_data$id, unique(dat_plot$id))
  if (length(missing_ids) > 0) {
    full_sample_data$category <- factor(
      full_sample_data$category,
      levels = c("HEALTHY", "MS"), labels = c("HD", "MS")
    )
    missing_numeric_id <- gsub("^(HD|MS)", "", missing_ids)
    missing_numeric_id <- ifelse(nchar(missing_numeric_id) == 0, missing_ids, missing_numeric_id)
    missing_rows <- data.frame(
      id = missing_ids,
      OTU_id = "Other",
      count = 0,
      category = full_sample_data$category[match(missing_ids, full_sample_data$id)],
      numeric_id = missing_numeric_id,
      numeric_order = suppressWarnings(as.numeric(missing_numeric_id)),
      stringsAsFactors = FALSE
    )
    dat_plot <- dplyr::bind_rows(dat_plot, missing_rows)
  }

  ordered_levels <- dat_plot %>%
    dplyr::distinct(numeric_id, numeric_order) %>%
    dplyr::arrange(is.na(numeric_order), numeric_order, numeric_id) %>%
    dplyr::pull(numeric_id) %>%
    as.character()

  dat_plot$numeric_id <- factor(dat_plot$numeric_id, levels = ordered_levels)

  pl <- dat_plot %>%
    ggplot(aes(x = numeric_id, y = count)) +
    facet_grid(~category, scales = "free_x", space = "free_x") +
    geom_bar(aes(fill = OTU_id), stat = "identity", position = "fill", width = 1) +
    xlab("Sample ID") +
    scale_fill_manual(values = palet_10) +
    scale_y_continuous(
      name = "Relative Abundance",
      labels = scales::percent
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom", axis.text.x = element_text(angle = 45, size = 8),
      axis.text.y = element_text(color = "black", size = 12),
      strip.text = element_text(face = "bold", size = 20),
      strip.background = element_blank()
    ) +
    guides(fill = guide_legend(title = "Species"))
  pl
  saveRDS(pl, gsub(" ", "", paste(output_folder, Domain, "_StackedBar_Species.rds", sep = "")))
  return(pl)
}


execute_stackbar <- function() {
  baseline_decB <- readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.rds")
  baseline_decA <- readRDS(file = "Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.rds")
  baseline_decE <- readRDS(file = "Output/SUPERVISED_DEC/Eukaryote_Supervised_decontam0.rds")

  stackbar(baseline_dec = baseline_decB, Domain = "Bacteria", output_folder)
  stackbar(baseline_dec = baseline_decA, Domain = "Archaea", output_folder)
  stackbar(baseline_dec = baseline_decE, Domain = "Eukaryote", output_folder)
}

execute_stackbar()
