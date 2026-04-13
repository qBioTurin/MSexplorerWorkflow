# Load data
data <- read.csv("InputData/Humann_20path_HDvsMS_norm.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(data) <- data$X
data <- data %>% select(-X)
sample_pres <- colnames(data)

metadata <- read.csv("InputData/Metadata_MS.csv", header = TRUE, stringsAsFactors = FALSE)
metadata <- metadata %>% filter(id %in% sample_pres & !is.na(category))
metadata$category[metadata$category == "HEALTHY"] <- "HV"
data <- data[, metadata$id]

# Convert data to matrix
data_matrix <- as.matrix(data)
data_matrix[is.na(data_matrix) | is.nan(data_matrix)] <- 0

# Reorder columns: first HV, then MS alphabetically
metadata <- metadata %>% arrange(category, id)
colOrder <- match(metadata$id, colnames(data_matrix))
data_matrix <- data_matrix[, colOrder]

# Prepare column annotation
col_anno <- data.frame(category = metadata$category)
rownames(col_anno) <- metadata$id
ann_colors_col <- list(category = c(HV = "#bad7f2", MS = "#ef6f6c"))

# Calcolo range simmetrico
abs_max <- max(abs(data_matrix), na.rm = TRUE)

# Numero di step complessivi
n_total <- 35  

# Breaks simmetrici: da -abs_max a +abs_max
myBreaks <- seq(-abs_max, abs_max, length.out = n_total + 1)

# Palette continua (blu → bianco → rosso) con 35 gradienti
myColors1 <- colorRampPalette(c("blue", "white", "red"))(n_total)

# Creazione heatmap
heatmap <- pheatmap(
  mat = data_matrix,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  annotation_col = col_anno,
  annotation_colors = ann_colors_col,
  col = myColors1,
  breaks = myBreaks,
  fontsize = 14,
  fontsize_row = 18,
  fontsize_col = 12,
  fontsize.legend = 14,
  border_color = NA,
  silent = TRUE
)

if(!dir.exists("Image/Humann")) dir.create("Image/Humann", recursive = TRUE)
pdf("Image/Humann/humann_category_norm.pdf", width = 35, height = 15)
grid::grid.draw(heatmap$gtable)
dev.off()

########################################################HEALTHY
# Load data
data <- read.csv("InputData/Humann_20path_HDvsMS.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(data) <- data$X
data <- data %>% select(-X)
sample_pres <- colnames(data)

metadata <- read.csv("InputData/Metadata_MS.csv", header = TRUE, stringsAsFactors = FALSE)
metadata <- metadata %>% filter(id %in% sample_pres & !is.na(category) & category != "HEALTHY")
metadata$category[metadata$category == "HEALTHY"] <- "HV"
data <- data[, metadata$id]

# Convert data to matrix
data_matrix <- as.matrix(data)
data_matrix[is.na(data_matrix) | is.nan(data_matrix)] <- 0

# Reorder columns: first HV, then MS alphabetically
metadata <- metadata %>% arrange(category, id)
colOrder <- match(metadata$id, colnames(data_matrix))
data_matrix <- data_matrix[, colOrder]

# Prepare column annotation
col_anno <- data.frame(category = metadata$category)
rownames(col_anno) <- metadata$id
#ann_colors_col <- list(category = c(HV = "#bad7f2", MS = "#ef6f6c"))

# Calcolo range simmetrico
abs_max <- max(abs(data_matrix), na.rm = TRUE)

# Numero di step complessivi
n_total <- 35  

# Breaks simmetrici: da -abs_max a +abs_max
myBreaks <- seq(-abs_max, abs_max, length.out = n_total + 1)

# Palette continua (blu → bianco → rosso) con 35 gradienti
myColors1 <- colorRampPalette(c("blue", "white", "red"))(n_total)

# Creazione heatmap
heatmap <- pheatmap(
  mat = data_matrix,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
 # annotation_col = col_anno,
 # annotation_colors = ann_colors_col,
  col = myColors1,
  breaks = myBreaks,
  fontsize = 14,
  fontsize_row = 18,
  fontsize_col = 12,
  fontsize.legend = 14,
  border_color = NA,
  silent = TRUE
)

if(!dir.exists("Image/Humann")) dir.create("Image/Humann", recursive = TRUE)
pdf("Image/Humann/humann_category_MS.pdf", width = 35, height = 20)
grid::grid.draw(heatmap$gtable)
dev.off()

###########################################################MS

data <- read.csv("InputData/Humann_20path_HDvsMS.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(data) <- data$X
data <- data %>% select(-X)
sample_pres <- colnames(data)

metadata <- read.csv("InputData/Metadata_MS.csv", header = TRUE, stringsAsFactors = FALSE)
metadata <- metadata %>% filter(id %in% sample_pres & !is.na(category) & category != "MS")
#metadata$category[metadata$category == "HEALTHY"] <- "HV"
data <- data[, metadata$id]

# Convert data to matrix
data_matrix <- as.matrix(data)
data_matrix[is.na(data_matrix) | is.nan(data_matrix)] <- 0

# Reorder columns: first HV, then MS alphabetically
metadata <- metadata %>% arrange(category, id)
colOrder <- match(metadata$id, colnames(data_matrix))
data_matrix <- data_matrix[, colOrder]

# Prepare column annotation
col_anno <- data.frame(category = metadata$category)
rownames(col_anno) <- metadata$id
#ann_colors_col <- list(category = c(HV = "#bad7f2", MS = "#ef6f6c"))

# Calcolo range simmetrico
abs_max <- max(abs(data_matrix), na.rm = TRUE)

# Numero di step complessivi
n_total <- 35  

# Breaks simmetrici: da -abs_max a +abs_max
myBreaks <- seq(-abs_max, abs_max, length.out = n_total + 1)

# Palette continua (blu → bianco → rosso) con 35 gradienti
myColors1 <- colorRampPalette(c("blue", "white", "red"))(n_total)

# Creazione heatmap
heatmap <- pheatmap(
  mat = data_matrix,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
 # annotation_col = col_anno,
 # annotation_colors = ann_colors_col,
  col = myColors1,
  breaks = myBreaks,
  fontsize = 14,
  fontsize_row = 18,
  fontsize_col = 12,
  fontsize.legend = 14,
  border_color = NA,
  silent = TRUE
)

if(!dir.exists("Image/Humann")) dir.create("Image/Humann", recursive = TRUE)
pdf("Image/Humann/humann_category_healthy.pdf", width = 35, height = 20)
grid::grid.draw(heatmap$gtable)
dev.off()


##########################GC treatment##########################################################
# Load data
data <- read.csv("InputData/Humann_20path_MS.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(data) <- data$X
data <- data %>% select(-X)
sample_pres <- colnames(data)

metadata <- read.csv("InputData/Metadata_MS.csv", header = TRUE, stringsAsFactors = FALSE)
metadata <- metadata %>%
  filter(id %in% sample_pres & !is.na(gc_treatment) & gc_treatment != "HEALTHY")
data <- data[, metadata$id]

# Convert data to matrix
data_matrix <- as.matrix(data)
data_matrix[is.na(data_matrix) | is.nan(data_matrix)] <- 0

# Reorder columns: by gc_treatment, then alphabetically
metadata <- metadata %>% arrange(gc_treatment, id)
colOrder <- match(metadata$id, colnames(data_matrix))
data_matrix <- data_matrix[, colOrder]

# Calcolo range simmetrico
abs_max <- max(abs(data_matrix), na.rm = TRUE)

# Numero di step totali (50 negativi + 50 positivi)
n_total <- 35  

# Breaks simmetrici
myBreaks <- seq(-abs_max, abs_max, length.out = n_total + 1)

# Palette continua simmetrica (blu → bianco → rosso)
myColors1 <- colorRampPalette(c("blue", "white", "red"))(n_total)

# Prepare column annotation
col_anno <- data.frame(gc_treatment = metadata$gc_treatment)
rownames(col_anno) <- metadata$id
ann_colors_col <- list(gc_treatment = setNames(c("#4D4D4D", "#D7D7D7"),
                                              c("positive", "negative")))

# Creazione heatmap
heatmap <- pheatmap(
  mat = data_matrix,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  annotation_col = col_anno,
  annotation_colors = ann_colors_col,
  col = myColors1,
  breaks = myBreaks,
  fontsize = 14,
  fontsize_row = 18,
  fontsize_col = 12,
  fontsize.legend = 14,
  border_color = NA,
  silent = TRUE
)

if(!dir.exists("Image/Humann")) dir.create("Image/Humann", recursive = TRUE)
pdf("Image/Humann/humann_gc_treatment.pdf", width = 35, height = 15)
grid::grid.draw(heatmap$gtable)
dev.off()
##############################positive

data <- read.csv("InputData/Humann_20path_MS.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(data) <- data$X
data <- data %>% select(-X)
sample_pres <- colnames(data)

metadata <- read.csv("InputData/Metadata_MS.csv", header = TRUE, stringsAsFactors = FALSE)
metadata <- metadata %>%
  filter(id %in% sample_pres & !is.na(gc_treatment) & gc_treatment != "HEALTHY" & gc_treatment != "negative")
data <- data[, metadata$id]

# Convert data to matrix
data_matrix <- as.matrix(data)
data_matrix[is.na(data_matrix) | is.nan(data_matrix)] <- 0

# Reorder columns: by gc_treatment, then alphabetically
metadata <- metadata %>% arrange(gc_treatment, id)
colOrder <- match(metadata$id, colnames(data_matrix))
data_matrix <- data_matrix[, colOrder]

# Calcolo range simmetrico
abs_max <- max(abs(data_matrix), na.rm = TRUE)

# Numero di step totali (50 negativi + 50 positivi)
n_total <- 35  

# Breaks simmetrici
myBreaks <- seq(-abs_max, abs_max, length.out = n_total + 1)

# Palette continua simmetrica (blu → bianco → rosso)
myColors1 <- colorRampPalette(c("blue", "white", "red"))(n_total)

# Prepare column annotation
col_anno <- data.frame(gc_treatment = metadata$gc_treatment)
rownames(col_anno) <- metadata$id
#ann_colors_col <- list(gc_treatment = setNames(c("#4D4D4D", "#D7D7D7"),
#                                              c("positive", "negative")))

# Creazione heatmap
heatmap <- pheatmap(
  mat = data_matrix,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
 # annotation_col = col_anno,
 # annotation_colors = ann_colors_col,
  col = myColors1,
  breaks = myBreaks,
  fontsize = 14,
  fontsize_row = 18,
  fontsize_col = 12,
  fontsize.legend = 14,
  border_color = NA,
  silent = TRUE
)

if(!dir.exists("Image/Humann")) dir.create("Image/Humann", recursive = TRUE)
pdf("Image/Humann/humann_gc_treatment_positive.pdf", width = 35, height = 15)
grid::grid.draw(heatmap$gtable)
dev.off()


##############################negative
data <- read.csv("InputData/Humann_20path_MS.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(data) <- data$X
data <- data %>% select(-X)
sample_pres <- colnames(data)

metadata <- read.csv("InputData/Metadata_MS.csv", header = TRUE, stringsAsFactors = FALSE)
metadata <- metadata %>%
  filter(id %in% sample_pres & !is.na(gc_treatment) & gc_treatment != "HEALTHY" & gc_treatment != "positive")
data <- data[, metadata$id]

# Convert data to matrix
data_matrix <- as.matrix(data)
data_matrix[is.na(data_matrix) | is.nan(data_matrix)] <- 0

# Reorder columns: by gc_treatment, then alphabetically
metadata <- metadata %>% arrange(gc_treatment, id)
colOrder <- match(metadata$id, colnames(data_matrix))
data_matrix <- data_matrix[, colOrder]

# Calcolo range simmetrico
abs_max <- max(abs(data_matrix), na.rm = TRUE)

# Numero di step totali (50 negativi + 50 positivi)
n_total <- 35  

# Breaks simmetrici
myBreaks <- seq(-abs_max, abs_max, length.out = n_total + 1)

# Palette continua simmetrica (blu → bianco → rosso)
myColors1 <- colorRampPalette(c("blue", "white", "red"))(n_total)

# Prepare column annotation
col_anno <- data.frame(gc_treatment = metadata$gc_treatment)
rownames(col_anno) <- metadata$id
#ann_colors_col <- list(gc_treatment = setNames(c("#4D4D4D", "#D7D7D7"),
   #                                            c("positive", "negative")))

# Creazione heatmap
heatmap <- pheatmap(
  mat = data_matrix,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  #annotation_col = col_anno,
  #annotation_colors = ann_colors_col,
  col = myColors1,
  breaks = myBreaks,
  fontsize = 14,
  fontsize_row = 18,
  fontsize_col = 12,
  fontsize.legend = 14,
  border_color = NA,
  silent = TRUE
)

if(!dir.exists("Image/Humann")) dir.create("Image/Humann", recursive = TRUE)
pdf("Image/Humann/humann_gc_treatment_negative.pdf", width = 35, height = 15)
grid::grid.draw(heatmap$gtable)
dev.off()




library(gridExtra)

# Leggi il file TSV dei risultati significativi
sig_results <- read.table("Output/Diff_humann/category/significant_results.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sig_results <- sig_results[, c("feature", "metadata", "pval_joint", "qval_joint")]
# Crea una tabella per il PDF

table_grob <- tableGrob(sig_results, rows = NULL)


# Salva la tabella in PDF
pdf("Image/Humann/significant_results_category_table.pdf", width = 12, height = 35)
grid::grid.draw(table_grob)
dev.off()
