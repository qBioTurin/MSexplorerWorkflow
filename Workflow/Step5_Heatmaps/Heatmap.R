source("Settings/utilities.R")
output_folder = "Image/Heatmap/"
createFolder(output_folder)
#### Functions ############
generate.heatmap = function(Bacteria,Archaea,Eukaryota,filename,filterRows = FALSE,output_folder)
{
  data_bact=data.generation(Bacteria)
  data_arch=data.generation(Archaea)
  data_euk=data.generation(Eukaryota)
  
  ### Kingdom together ###### 
  norm_data_z=rbind(rbind(data_bact[[1]],data_arch[[1]]),data_euk[[1]])
  metadata_hm=rbind(rbind(data_bact[[2]],data_arch[[2]]),data_euk[[2]])
  metadata_kingdom=rbind(rbind(data_bact[[3]],data_arch[[3]]),data_euk[[3]])
 
  ann_colors = list(
    Sex= setNames(c("blue", "pink"),
                  levels(as.factor(metadata_hm$Sex))),
    Status= setNames(c("#bad7f2","#ef6f6c"),
                     levels(as.factor(metadata_hm$Status))),
    gc_treatment = setNames(c("#4D4D4D","#D7D7D7"),
                                        c("positive" , "negative") ) ,
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
  
  ### Kingdom separated ###### 
  heatmap.kingdom = function(data,color,paletteLength, col_order = NULL,
                             phylumPalette =  viridisLite::inferno, USEann_colors = T,filterRows ){
    row_number=(nrow(data[[1]])>1)
    norm_data_z = data[[1]]
    metadata_hm= data[[2]]
    metadata_kingdom= data[[3]]
    
    norm_data_z[is.na(norm_data_z) | is.nan(norm_data_z)] = 0
    
    if(filterRows){
      norm_data_z = norm_data_z[, grep(x = colnames(norm_data_z),pattern = "MS")]
     }
    
    myBreaks = c(seq(min(norm_data_z), 0, length.out = ceiling(paletteLength/2) + 1),
                 seq(-1/paletteLength, max(norm_data_z), length.out = floor(paletteLength/2)))
    
    col_fun = circlize::colorRamp2(myBreaks,rev(color))
    
    if(!is.null(col_order) )  norm_data_z= norm_data_z[,col_order] 
    
    ann_colors$Phylum = setNames(phylumPalette(n = length(unique(metadata_kingdom$Phylum))),
                                 levels(as.factor(metadata_kingdom$Phylum)))
    
    metadata_hm = metadata_hm %>% select(subtentorial_lesions, gadolinium_contrast,
                                         Spinal_cord_lesions,Lesion_Burden,
                                         gc_treatment, Sex, Status)
    
    metadata_hm$gc_treatment[ metadata_hm$gc_treatment == "healthy"] = NA
    
    if(length(unique(metadata_hm$Status)) == 1)
      metadata_hm = metadata_hm %>% select(-Status)
    
    if(length(unique(metadata_hm$Lesion_Burden)) == 1)
      metadata_hm = metadata_hm %>% select(-Lesion_Burden)
    
    if(length(unique(metadata_hm$gc_treatment)) == 1)
      metadata_hm = metadata_hm %>% select(-gc_treatment)
    
    if(!USEann_colors){
      metadata_hm = NA
      ann_colors = ann_colors["Phylum"]
    }
    
    phil_name =  paste0("Phylum ", unique(metadata_kingdom$Kingdom))
    names(ann_colors)[names(ann_colors) == "Phylum"]  = phil_name
    metadata_kingdom = metadata_kingdom %>% select(Phylum)
    colnames(metadata_kingdom) = phil_name
    
    heatmap<-pheatmap::pheatmap(norm_data_z, 
                                annotation_col = metadata_hm, 
                                annotation_colors = ann_colors,
                                annotation_row = metadata_kingdom, 
                                cluster_cols = ifelse(is.null(col_order), T, F ), 
                                cluster_rows = row_number,
                                col=color,
                                border_color = NA,
                                breaks  = myBreaks,
                                silent = T
    )
    
    
    return(heatmap)
  }
  
  myColors1 = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(paletteLength+1))
  myColors2 = colorRampPalette(RColorBrewer::brewer.pal(9, "PuOr"))(paletteLength+1)
  myColors3 = colorRampPalette(RColorBrewer::brewer.pal(9, "BrBG"))(paletteLength+1)
  
  heatmap_bact = heatmap.kingdom(data_bact,myColors1, 35,
                                 phylumPalette = viridisLite::mako, filterRows = filterRows)
  

  heatmap_list = list()
  colOrder_bact = NULL
  
  if (!is.null(data_bact[[1]])) {
    heatmap_bact = heatmap.kingdom(data_bact, myColors1, 35, phylumPalette = viridisLite::mako, filterRows = filterRows)
    heatmap_list$Bacteria = heatmap_bact
    ggplotify::as.ggplot(heatmap_bact) -> heatmap_bact_ggplot
    colOrder_bact = heatmap_bact_ggplot$plot_env$plot$tree_col$labels[heatmap_bact_ggplot$plot_env$plot$tree_col$order]
  }
  
  if (!is.null(data_arch[[1]])) {
    heatmap_arch = heatmap.kingdom(data_arch,col_order = colOrder_bact, myColors3, 35,
                                   phylumPalette = viridisLite::cividis,
                                   USEann_colors = F, filterRows = filterRows)
    if (!is.null(data_bact[[1]]))  heatmap_arch$gtable$widths = heatmap_bact$gtable$widths 
    heatmap_list$Archaea = heatmap_arch
  }
  if (!is.null(data_euk[[1]])) {
    heatmap_euk = heatmap.kingdom(data_euk,
                                  col_order = colOrder_bact, myColors2, 35,
                                  phylumPalette = viridisLite::magma,
                                  USEann_colors = F, filterRows = filterRows)
    if (!is.null(data_bact[[1]])) heatmap_euk$gtable$widths <- heatmap_bact$gtable$widths 
    else if(!is.null(data_euk[[1]])) heatmap_euk$gtable$widths <- heatmap_arch$gtable$widths
    heatmap_list$Eukaryota = heatmap_euk
  }

  # row_counts <- c(nrow(data_bact[[1]]), nrow(data_arch[[1]]), nrow(data_euk[[1]]))
  # total_rows <- sum(row_counts)
  # relative_heights <- row_counts / total_rows + c(0.05,0,0.05)
  # 
  # tbac = textGrob("Bacteria", rot = 90, gp = gpar(fontsize=14))
  # teuk = textGrob("Eukaryota", rot = 90, gp = gpar(fontsize=14))
  # tarc = textGrob("Archaea", rot = 90, gp = gpar(fontsize=14))
  
  row_counts = sapply(heatmap_list, function(hm) if (!is.null(hm)) nrow(hm$gtable) else 0)
  total_rows = sum(row_counts)
  relative_heights = row_counts / total_rows #+ c(0.05, 0, 0.05)
  
  labels = names(heatmap_list)
  text_grobs = lapply(labels, function(lbl) textGrob(lbl, rot = 90, gp = gpar(fontsize = 14)))
  
  g = grid.arrange(grobs = c(rbind(text_grobs, lapply(heatmap_list, function(hm) hm$gtable))),
                   ncol = 2, heights = relative_heights, widths = c(0.05, 1.5))
  
  ggsave(plot = g, filename = gsub(" ", "", paste0(output_folder, filename, ".pdf")),
         height = 22, width = 18, limitsize = FALSE)
  
  return(NULL)
}
data.generation=function(baselines_dec){
  if(is.null(baselines_dec))
    return(list(NULL,NULL,NULL))
    
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

  metadata_kingdom = norm_data%>%
    select(Kingdom,Phylum,Genus,Genus_species)
  rownames(metadata_kingdom) = metadata_kingdom$Genus_species
  metadata_kingdom = metadata_kingdom%>%select(-Genus_species)
  norm_data = norm_data[, -c(2:8)]

  
  rownames(norm_data) = norm_data$Genus_species
  norm_data = norm_data[,-c(1,length(norm_data))]
  colnames(norm_data)

   metadata = read.csv("InputData/metadataMS.csv", 
                      header = TRUE, 
                      sep = ",",
                      na = c("", " ", "NA"), 
                      check.names = TRUE)
 metadata <- metadata %>%
  mutate(
    across(.cols = c(id), .fns = as.character),
    
    across(.cols = c(sex, category, clinical_presentation, gc_treatment,
                     subtentorial_lesions, spinal_cord_lesion, gadolinium_contrast,
                     sequencing_batch,WORSENING, EDSS_DIAGNOSI, EDSS_PROGRESSIONE, Event,
                     naive, previous_therapy, antibiotic_use, sample_type,lesion_burden), 
           .fns = as.factor),
    across(.cols = c(age, bmi, EventTime
          ), .fns = as.numeric),
    across(.cols = c(sample_collection_date), .fns = as.Date)
  )
  
  samples = c(colnames(norm_data))
  metadata = metadata[metadata$id %in% samples,] 
  rownames(metadata) = metadata$id
  
  # select metadata to visualize in heatmap
  ###########################################
  # lesion burden
  metadata_hm = metadata %>%
    filter(category %in% c("HEALTHY","MS") ) %>%
    select(Sex = sex,
           Status=category,
           gc_treatment = gc_treatment,
           Lesion_Burden = lesion_burden,
           Spinal_cord_lesions = spinal_cord_lesion,
           gadolinium_contrast = gadolinium_contrast,
           subtentorial_lesions = subtentorial_lesions) %>%
    mutate(Sex= factor(Sex, levels = c("M", "F"), labels = c("M", "F")),
           Status= factor(Status, levels=c("HEALTHY","MS"),labels=c("HV","MS")),
           gc_treatment= factor(gc_treatment, levels = c("healthy","positive", "negative"), labels = c("healthy","positive", "negative")),
           Lesion_Burden = factor(Lesion_Burden, levels = c(0, 1), labels = c("Low", "High")),
           Spinal_cord_lesions = factor(Spinal_cord_lesions, levels = c(0, 1), labels = c("No", "Yes")),
           subtentorial_lesions = factor(subtentorial_lesions, levels = c(0, 1), labels = c("No", "Yes")),
           gadolinium_contrast = factor(gadolinium_contrast, levels = c(0, 1), labels = c("No", "Yes")))
  # #bone marrow lesions
  norm_data_z = t(scale(t(norm_data)))
  return(list(norm_data_z,metadata_hm,metadata_kingdom))
}

######################
##########ALL###########
generate.heatmap(
  Bacteria = readRDS(file = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds"),
  Archaea = readRDS(file = "Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.001.rds"),
  Eukaryota = readRDS(file = "Output/SUPERVISED_DEC/Eukaryote_Supervised_decontam0.001.rds"),,
  filename = "ALL", output_folder)

########################
######### MSHD ############
generate.heatmap(
  Bacteria = readRDS(file = "Output/merge_DAS/MSHD/Bacteria_MsHd_merged.rds"),
  Archaea = readRDS(file = "Output/merge_DAS/MSHD/Archaea_MsHd_merged.rds"),
  Eukaryota = readRDS(file = "Output/merge_DAS/MSHD/Eukaryote_MsHd_merged.rds"),,
  filename = "category", output_folder)
#########GC_LESION############
generate.heatmap(
  Bacteria = readRDS("Output/merge_DAS/GC/Bacteria_GC_merged.rds"),
  Archaea = readRDS("Output/merge_DAS/GC/Archaea_GC_merged.rds"),
  Eukaryota = readRDS("Output/merge_DAS/GC/Eukaryote_GC_merged.rds"),,
  filename = "GC_LESION", output_folder)


#########lesion_burden############
 generate.heatmap(
        Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_lesion_burden_positive_merged.rds")) ,
        Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_lesion_burden_positive_merged.rds")),
        Eukaryota = readRDS(paste0("Output/merge_DAS/GC_comp/Eukaryote_lesion_burden_positive_merged.rds")),,
        filename = paste0("lesion_burden_positive"), output_folder )

 generate.heatmap(
        Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_lesion_burden_negative_merged.rds")) ,
        Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_lesion_burden_negative_merged.rds")),
        Eukaryota = readRDS(paste0("Output/merge_DAS/GC_comp/Eukaryote_lesion_burden_negative_merged.rds")),,
        filename = paste0("lesion_burden_negative"), output_folder )


#########spinal_cord_lesion############
generate.heatmap(
        Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_spinal_cord_lesion_positive_merged.rds")) ,
        Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_spinal_cord_lesion_positive_merged.rds")),
        Eukaryota =NULL,,
        filename = paste0("spinal_cord_lesion_positive"), output_folder )########lol


generate.heatmap(
        Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_spinal_cord_lesion_negative_merged.rds")) ,
        Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_spinal_cord_lesion_negative_merged.rds")),
        Eukaryota = readRDS(paste0("Output/merge_DAS/GC_comp/Eukaryote_spinal_cord_lesion_negative_merged.rds")),,
        filename = paste0("spinal_cord_lesion_negative"), output_folder )

#########subtentorial_lesions############
generate.heatmap(
        Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_subtentorial_lesions_positive_merged.rds")) ,
        Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_subtentorial_lesions_positive_merged.rds")),
        Eukaryota = readRDS(paste0("Output/merge_DAS/GC_comp/Eukaryote_subtentorial_lesions_positive_merged.rds")),,
        filename = paste0("subtentorial_lesions_positive"), output_folder )

generate.heatmap(
        Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_subtentorial_lesions_negative_merged.rds")) ,
        Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_subtentorial_lesions_negative_merged.rds")),
        Eukaryota = NULL,,
        filename = paste0("subtentorial_lesions_negative"), output_folder )#######LOL
#########gadolinium_contrast############
generate.heatmap(
        Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_gadolinium_contrast_positive_merged.rds")) ,
        Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_gadolinium_contrast_positive_merged.rds")),
        Eukaryota = readRDS(paste0("Output/merge_DAS/GC_comp/Eukaryote_gadolinium_contrast_positive_merged.rds")),,
        filename = paste0("gadolinium_contrast_positive"), output_folder )
generate.heatmap(
        Bacteria = readRDS(paste0("Output/merge_DAS/GC_comp/Bacteria_gadolinium_contrast_negative_merged.rds")) ,
        Archaea = readRDS(paste0("Output/merge_DAS/GC_comp/Archaea_gadolinium_contrast_negative_merged.rds")),
        Eukaryota = readRDS(paste0("Output/merge_DAS/GC_comp/Eukaryote_gadolinium_contrast_negative_merged.rds")),,
        filename = paste0("gadolinium_contrast_negative"), output_folder )
