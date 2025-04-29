#####heatmaps bacteria score####
myfunction = function(input, filterhm) {
  
  baselines_dec = input
  
  taxatables = as.data.frame(tax_table(baselines_dec))
  taxatables = rownames_to_column(taxatables, var = "otu_id")
  
  # otu table
  otutables = data.frame(otu_table(baselines_dec))
  otutables = rownames_to_column(otutables, var = "otu_id")
  
  # combine tables
  combined_df = full_join(taxatables, otutables, by = "otu_id")
  
  ## from combined tables, keep only columns otu_id, Genus+species + otu table)
  norm_data = combined_df %>%
  mutate(Genus_species = paste(Genus, Species, sep = " "))
  metadata_kingdom = norm_data %>%
  select(Kingdom, Phylum, Genus, Genus_species)
  rownames(metadata_kingdom) = metadata_kingdom$Genus_species
  metadata_kingdom = metadata_kingdom %>% select(-Genus_species)
  norm_data = norm_data[, -c(2:8)]
  
  rownames(norm_data) = norm_data$Genus_species
  norm_data = norm_data[, -c(1, length(norm_data))]
  colnames(norm_data)
  
  metadata = read.csv("InputData/metadataMS.csv", 
            header = TRUE, 
            sep = ",",
            na = c("", " ", "NA"), 
            check.names = TRUE)
  baselines_metadata <- baselines_metadata %>%
  mutate(
    across(.cols = c(id), .fns = as.character),
    
    across(.cols = c(sex, category, clinical_presentation, gc_treatment,
             subtentorial_lesions, spinal_cord_lesion, gadolinium_contrast,
             sequencing_batch, WORSENING, EDSS_DIAGNOSI, EDSS_PROGRESSIONE, Event,
             naive, previous_therapy, antibiotic_use, sample_type, lesion_burden), 
       .fns = as.factor),
    
    across(.cols = c(age, bmi, EventTime), .fns = as.numeric),
    
    across(.cols = c(sample_collection_date), .fns = as.Date)
  )
  
  rownames(metadata) = metadata$id
  
  # select metadata to visualize in heatmap
  ###########################################
  # lesion burden
  metadata_hm = metadata %>%
  filter(category %in% c("HEALTHY", "MS")) %>%
  select(Sex = sex,
       Status = category,
       gc_treatment = gc_treatment,
       Lesion_Burden = lesion_burden,
       Spinal_cord_lesions = spinal_cord_lesion,
       gadolinium_contrast = gadolinium_contrast,
       subtentorial_lesions = subtentorial_lesions) %>%
  mutate(Sex = factor(Sex, levels = c("M", "F"), labels = c("M", "F")),
       Status = factor(Status, levels = c("HEALTHY", "MS"), labels = c("HV", "MS")),
       gc_treatment = factor(gc_treatment, levels = c("healthy", "positive", "negative"), labels = c("healthy", "positive", "negative")),
       Lesion_Burden = factor(Lesion_Burden, levels = c(0, 1), labels = c("Low", "High")),
       Spinal_cord_lesions = factor(Spinal_cord_lesions, levels = c(0, 1), labels = c("No", "Yes")),
       subtentorial_lesions = factor(subtentorial_lesions, levels = c(0, 1), labels = c("No", "Yes")),
       gadolinium_contrast = factor(gadolinium_contrast, levels = c(0, 1), labels = c("No", "Yes")))
  
  # #Spinal cord lesions
  norm_data_z = t(scale(t(norm_data)))
  
  ann_colors = list(
  Sex = setNames(c("blue", "pink"),
           levels(as.factor(metadata_hm$Sex))),
  Status = setNames(c("#bad7f2", "#ef6f6c"),
            levels(as.factor(metadata_hm$Status))),
  gc_treatment = setNames(c("#4D4D4D", "#D7D7D7"),
              c("positive", "negative")),
  Spinal_cord_lesions = setNames(c("brown", "#13913B"),
                   levels(as.factor(metadata_hm$Spinal_cord_lesions))),
  Lesion_Burden = setNames(c("orange", "purple"),
               levels(as.factor(metadata_hm$Lesion_Burden))),
  subtentorial_lesions = setNames(c("#233E6CFF", "#FFEA46FF"),
                  levels(as.factor(metadata_hm$subtentorial_lesions))),
  gadolinium_contrast = setNames(c("#E21277", "#0DC9A7"),
                   levels(as.factor(metadata_hm$gadolinium_contrast)))
  )
  
  paletteLength = 35
  
  myBreaks = c(seq(min(norm_data_z), 0, length.out = ceiling(paletteLength / 2) + 1),
         seq(-1 / paletteLength, max(norm_data_z), length.out = floor(paletteLength / 2)))
  
  myColors1 = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(paletteLength + 1))
  
  phylumPalette = viridisLite::viridis
  
  ann_colors$Phylum = setNames(phylumPalette(n = length(unique(metadata_kingdom$Phylum))),
                 levels(as.factor(metadata_kingdom$Phylum)))
  
  metadata_hm = filterhm(metadata_hm)
  
  metadata_kingdom <- metadata_kingdom %>%
  select(Phylum)
  
  heatmap <- pheatmap::pheatmap(norm_data_z, 
                annotation_col = metadata_hm, 
                annotation_colors = ann_colors,
                annotation_row = metadata_kingdom, 
                cluster_cols = T,
                cluster_rows = T,
                col = myColors1,
                border_color = NA,
                breaks = myBreaks, 
                silent = T)
  
  heatmap = ggplotify::as.ggplot(heatmap)
  
  return(heatmap)
}

Les = readRDS("Output/merge_DAS/001/Bacteria_lesion_burden_001_merged.rds")
Gad = readRDS("Output/merge_DAS/001/Bacteria_gadolinium_contrast_001_merged.rds")
SC = readRDS("Output/merge_DAS/001/Bacteria_spinal_cord_lesion_001_merged.rds")
ST = readRDS("Output/merge_DAS/001/Bacteria_subtentorial_lesions_001_merged.rds")

heatmap1 <- myfunction(
  Les, filterhm = function(df) df %>% select(Lesion_Burden))

heatmap2 <- myfunction(
  SC, filterhm = function(df) df %>% select(Spinal_cord_lesions))

heatmap3 <- myfunction(
  Gad, filterhm = function(df) df %>% select(gadolinium_contrast))

heatmap4 <- myfunction(
  ST, filterhm = function(df) df %>% select(subtentorial_lesions))

Les_ID = rownames(Les@tax_table)
Gad_ID = rownames(Gad@tax_table)
SC_ID = rownames(SC@tax_table)
ST_ID = row.names(ST@tax_table)

venn_list = list("Lesion" = Les_ID, "Gadolinium" = Gad_ID, "Spinal Cord" = SC_ID, "Subtentorial" = ST_ID)

venn <- ggVennDiagram(venn_list,
            set_color = c("purple", '#0DC9A7', '#13913B', '#aa9c2f'),
            show_intersect = F,
            label = "count",
            set_size = 5) 

overlaps <- VennDiagram::calculate.overlap(
  x = list(Les_ID, Gad_ID, SC_ID, ST_ID)
)

ggsave(plot = venn, filename = "Image/Heatmap/venn.pdf",
     height = 8, width = 12, limitsize = FALSE)
patch = heatmap1 | (heatmap2 / heatmap3 / heatmap4)

ggsave(plot = patch, filename = "Image/Heatmap/patchMod.pdf",
     height = 20, width = 24, limitsize = FALSE)
