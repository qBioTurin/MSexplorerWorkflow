source("utilities.R")
source("Colorpalette.R")
output_folder = "Image/AlphaDasMod/"####output folder
createFolder(output_folder)
library(ggrepel)
metadata <- read_delim("input/Bacteria_alpha_metadata.csv", delim = ",", escape_double = FALSE,  trim_ws = TRUE)###metadata folder

  metadata$gc_treatment=as.factor(metadata$gc_treatment)
  metadata$category=as.factor(metadata$category)
  metadata$lesion_burden=as.factor(metadata$lesion_burden)
  levels(metadata$lesion_burden) =c("Lesion Low","Lesion High")
  metadata$bone_marrow_lesions=as.factor(metadata$bone_marrow_lesions)
  levels(metadata$bone_marrow_lesions) =c("Spinal cord Low","Spinal cord High")
  metadata$gadolinium_contrast=as.factor(metadata$gadolinium_contrast)
  levels(metadata$gadolinium_contrast) =c("Gadolinum No","Gadolinum Yes")
  metadata$subtentorial_lesions=as.factor(metadata$subtentorial_lesions)
  levels(metadata$subtentorial_lesions) =c("Subtentorial No","Subtenorial Yes")  
  metadata$fused_col <- paste(metadata$id,metadata$lesion_burden, metadata$bone_marrow_lesions, metadata$gadolinium_contrast,metadata$subtentorial_lesions, sep = "_")
  metaSel <- metadata %>%
    select(id, fused_col) %>%
    as.data.frame()
 
create_alpa_cluster <- function(file,metaData,output_folder) {
    data_new <- read.csv(file)
    rownames(data_new) <- data_new$X 
    data_new <- data_new[, -1]
    set.seed(1234)
    data_new <- t(data_new)
    data_new <- as.data.frame(data_new)
    data_new2 <- data_new %>% mutate(id = rownames(.))
    merged_table <- merge(data_new2, metadata, by = "id")%>%
        select(id, lesion_burden, bone_marrow_lesions, gadolinium_contrast, subtentorial_lesions)
    rownames(merged_table) <- merged_table$id
    merged_table <- merged_table[, -1]
    write.csv(merged_table, file = paste0(output_folder, "metadata_table.csv"), row.names = TRUE)
    data_new <- cbind(data_new, id = rownames(data_new))
    #data_new <- as.data.frame(data_new)
    merged_table <- data_new %>% 
    select(Lesion, Spinal_Cord, Gadolinium, Subtentorial)
    toSaveTable<-merged_table
    merged_table <- merged_table %>%
    mutate(across(everything(), as.numeric))

    sil_values <- fviz_nbclust(merged_table, kmeans, method = "silhouette")$data
    optimal_k_sil <- sil_values$clusters[which.max(sil_values$y)]
    pca_res<-prcomp(merged_table, scale = TRUE)
    pca_res <- prcomp(merged_table, scale = TRUE)

    plot <- autoplot(kmeans(merged_table, as.numeric(optimal_k_sil)), 
            data = merged_table, label = FALSE, frame = TRUE, 
            frame.border.color = "black",loadings.label = TRUE, loadings.label.size  = 3,frame.type = 'norm') +
            geom_point(color = "black") +  # Ensures all points are black
            geom_text_repel(aes(label = rownames(merged_table)), color = "black", 
                    size = 4, max.overlaps = Inf) +
            theme_classic() +
            theme(legend.position = "none") +
            ggtitle("K-means clustering") +
            scale_color_manual(values = rep("black", as.numeric(optimal_k_sil)))  
plot
    file_name <- tools::file_path_sans_ext(basename(file))
    ggsave(paste0(output_folder, file_name, ".jpeg"), plot = plot, width = 10, height = 10, dpi = 300)
}


files <- list.files(path = "Output/modDAS_ALPHA", pattern = "\\.csv$", full.names = TRUE)####input folder

for (file in files) {
    create_alpa_cluster(file,metaData,output_folder)
}


