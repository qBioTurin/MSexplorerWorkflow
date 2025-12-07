DASLimma <- function(baselines_dec, metadata, analisys,
                     subcategory, sub_name, output_folder, discriminant) {
    if (!dir.exists(output_folder)) {
        dir.create(output_folder, recursive = TRUE)
    }
    cat("TYPE baselines_dec:", class(baselines_dec), "\n")
    cat("TYPE metadata:", class(metadata), "\n")
    subname_paste <- sub_name_creator(sub_name)
    domain <- find_domain(baselines_dec)
    baselines_dec <- readRDS(baselines_dec)

    filename <- gsub(" ", "", paste0(
        "limma_", domain, "_",
        analisys, ".csv"
    ))
    fullname <- gsub(" ", "", paste0(output_folder, filename))
    # Estrai i dati di sample_data e converti in un tibble
    samples <- as.data.frame(as_tibble(sample_data(baselines_dec)))

    # Filtra solo se subcategory è diverso da "none"
    if (subcategory != "none" | subcategory != "") {
        samples <- samples %>% filter(.data[[subcategory]] %in% sub_name)
    }
    baselines_dec <- prune_samples(samples$id, baselines_dec)
    print(baselines_dec)

    # Prepare taxonomic and OTU tables
    taxatables <- as.data.frame(tax_table(baselines_dec)) %>%
        rownames_to_column(var = "otu_id")

    otutables <- data.frame(otu_table(baselines_dec))
    otutables <- as.data.frame(abundances(otutables, transform = "compositional"))
    otutables <- rownames_to_column(otutables, var = "otu_id")

    combined_df <- full_join(taxatables, otutables, by = "otu_id") %>%
        mutate(Genus_species = paste(Genus, Species, sep = " ")) %>%
        select(Genus_species, everything(), -otu_id, -c(2:8))

    norm_data <- combined_df
    colnames(norm_data) <- gsub("[^[:alnum:]_]+", "_", colnames(norm_data))

    rownames(norm_data) <- norm_data$Genus_species
    norm_data <- norm_data[, -c(1, length(norm_data))]

    # Read metadata
    metadata <- read.csv(metadata,
        header = TRUE,
        sep = ",",
        na = c("", " ", "NA"),
        check.names = TRUE
    )
    samples <- colnames(norm_data)
    metadata$id <- gsub("[^[:alnum:]_]+", "_", metadata$id)
    # replace any non-alphanumeric or underscore characters with underscore
    samples <- gsub("[^[:alnum:]_]+", "_", samples)
    metadata <- metadata[metadata$id %in% samples, ]

    if (length(sub_name) > 0) {
        metadata_hm <- metadata
        rownames(metadata_hm) <- metadata_hm$id
        metadata_hm <- metadata_hm %>% filter(.data[[subcategory]] %in% sub_name)
        metadata_hm <- metadata_hm %>% mutate(across(all_of(analisys), as.factor))
    } else {
        metadata_hm <- metadata %>%
            filter(id %in% colnames(norm_data))

        rownames(metadata_hm) <- metadata_hm$id

        metadata_hm <- metadata_hm %>%
            select(all_of(c(analisys, discriminant))) %>%
            mutate(across(all_of(analisys), as.factor))
    }

    # Ensure data consistency
    norm_data <- norm_data %>%
        select(rownames(metadata_hm))

    log_data <- log2(norm_data + 1)

    formula <- paste0("~ ", analisys)
    for (el in discriminant) {
        formula <- paste0(formula, " + ", el)
    }

    design <- model.matrix(as.formula(formula), data = metadata_hm)


    cat("Dimension log_data:", dim(log_data), "\n")
    cat("Dimension design:", dim(design), "\n")
    ncol(log_data)
    array1 <- colnames(log_data)
    array2 <- rownames(design)
    diff <- setdiff(array1, array2)
    if (length(diff) > 0) {
        cat("Samples in log_data not in design:", diff, "\n")
        log_data <- log_data[, !colnames(log_data) %in% diff]
    } else {
        cat("All samples in log_data are present in design.\n")
    }

    ### LIMMA DAS Analysis
    fit <- lmFit(log_data, design)
    fit <- eBayes(fit)

    # remove NA coefficients
    estimable_coefs <- colnames(fit$coefficients)[!is.na(fit$coefficients[1, ])]

    if (length(estimable_coefs) == 0) {
        stop("Nessun coefficiente stimabile nel modello LIMMA.")
    }

    top_table <- topTable(fit,
        coef = estimable_coefs,
        number = Inf,
        adjust = "fdr"
    )
    write.csv(top_table, file = gsub(".csv", "_toptable.csv", fullname))
    top_table_p <- top_table %>%
        filter(P.Value < 0.05)    
    top_table_adj<- top_table %>%
        filter(adj.P.Val < 0.05)
    # Save results
    filter_for_heatmap_adj <- norm_data[rownames(top_table_adj), ]
    filter_for_heatmap <- norm_data[rownames(top_table_p), ]
    write.csv(filter_for_heatmap, file = gsub(" ", "", fullname))
    write.csv(filter_for_heatmap_adj, file = gsub(" ", "", gsub(".csv", "_adj.csv", fullname)))
    return(fullname)
}
source("Settings/Packages.R")
source("Settings/utilities.R")
source("Workflow/Step4_DiscriminantAnalysis/Limma/function.R")

output_folderMSHD <- "Output/LIMMA_score_corr/MSHD/"
createFolder(output_folderMSHD)
output_folder_001 <- "Output/LIMMA_score_corr/001/"
createFolder(output_folder_001)
output_folder_01 <- "Output/LIMMA_score_corr/01/"
createFolder(output_folder_01)
output_folder_05 <- "Output/LIMMA_score_corr/05/"
createFolder(output_folder_05)
library(tibble)

DASLimma(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "lesion_burden",
    subcategory = "gc_treatment",
    sub_name = c("positive", "negative"),
    output_folder = output_folder_001,
    discriminant = c("age", "sex", "bmi")
)
DASLimma(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "lesion_burden",
    subcategory = "gc_treatment",
    sub_name = c("positive", "negative"),
    output_folder = output_folder_01,
    discriminant = c("age", "sex", "bmi")
)
DASLimma(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "lesion_burden",
    subcategory = "gc_treatment",
    sub_name = c("positive", "negative"),
    output_folder = output_folder_05,
    discriminant = c("age", "sex", "bmi")
)
DASLimma(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "gadolinium_contrast",
    subcategory = "gc_treatment",
    sub_name = c("positive", "negative"),
    output_folder = output_folder_001,
    discriminant = c("age", "sex", "bmi")
)
DASLimma(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "gadolinium_contrast",
    subcategory = "gc_treatment",
    sub_name = c("positive", "negative"),
    output_folder = output_folder_01,
    discriminant = c("age", "sex", "bmi")
)
DASLimma(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "gadolinium_contrast",
    subcategory = "gc_treatment",
    sub_name = c("positive", "negative"),
    output_folder =    output_folder_05,
    discriminant = c("age", "sex", "bmi")
)
DASLimma(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "spinal_cord_lesion",
    subcategory = "gc_treatment",
    sub_name = c("positive", "negative"),
    output_folder = output_folder_001,
    discriminant = c("age", "sex", "bmi")
)
DASLimma(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "spinal_cord_lesion",
    subcategory = "gc_treatment",
    sub_name = c("positive", "negative"),
    output_folder = output_folder_01,
    discriminant = c("age", "sex", "bmi")
)
DASLimma(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "spinal_cord_lesion",
    subcategory = "gc_treatment",
    sub_name = c("positive", "negative"),
    output_folder = output_folder_05,
    discriminant = c("age", "sex", "bmi")
)
DASLimma(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "subtentorial_lesions",
    subcategory = "gc_treatment",
    sub_name = c("positive", "negative"),
    output_folder = output_folder_001,
    discriminant = c("age", "sex", "bmi")
)

DASLimma(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "category",
    subcategory = "category",
    sub_name = c("MS", "HEALTHY"),
    output_folder = output_folderMSHD,
    discriminant = c("age", "sex", "bmi")
)

DASLimma(
    baselines_dec = "Output/SUPERVISED_DEC/Eukaryote_Supervised_decontam0.001.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "category",
    subcategory = "category",
    sub_name = c("MS", "HEALTHY"),
    output_folder = output_folderMSHD,
    discriminant = c("age", "sex", "bmi ")
)

DASLimma(
    baselines_dec = "Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.001.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "category",
    subcategory = "category",
    sub_name = c("MS", "HEALTHY"),
    output_folder = output_folderMSHD,
    discriminant = c("age", "sex", "bmi ")
)

