DASLimmaVoom <- function(baselines_dec,
                     metadata,
                     analisys,
                     subcategory = "none",
                     sub_name = NULL,
                     output_folder = ".",
                     discriminant = character(0)) {
    # Caricamento librerie richieste
    requireNamespace("phyloseq", quietly = TRUE)
    requireNamespace("edgeR", quietly = TRUE)
    requireNamespace("limma", quietly = TRUE)
    requireNamespace("dplyr", quietly = TRUE)
    requireNamespace("tibble", quietly = TRUE)
    requireNamespace("rlang", quietly = TRUE)

    # --- 1) Preparazione cartella output ---
    if (!dir.exists(output_folder)) {
        dir.create(output_folder, recursive = TRUE)
    }

    # --- 2) Lettura dell'oggetto phyloseq ---
    # baselines_dec può essere un path a RDS o un oggetto phyloseq già caricato
    if (is.character(baselines_dec) && file.exists(baselines_dec)) {
        ps <- readRDS(baselines_dec)
    } else if (inherits(baselines_dec, "phyloseq")) {
        ps <- baselines_dec
    } else {
        stop("baselines_dec deve essere un path a un .rds contenente un oggetto phyloseq o un oggetto phyloseq.")
    }

    # Optional: trova domain per nome file di output (se hai la funzione find_domain)
    domain <- tryCatch({
        if (exists("find_domain")) find_domain(baselines_dec) else "domain"
    }, error = function(e) "domain")

    filename <- gsub(" ", "", paste0("limma_", domain, "_", analisys, ".csv"))
    fullname <- file.path(output_folder, filename)

    # --- 3) Estrai sample_data e applica filtro subcategory (se richiesto) ---
    samples_df <- as.data.frame(phyloseq::sample_data(ps))


    if (!is.null(subcategory) && subcategory != "none" && subcategory != "" && !is.null(sub_name)) {
        samples_df <- samples_df[samples_df[[subcategory]] %in% sub_name, , drop = FALSE]
    }

    # prune i sample nel phyloseq in modo coerente
    ps <- phyloseq::prune_samples(samples_df$id, ps)

    # --- 4) Prepara tax_table (per nomi specie) e conti grezzi ---
    tax_tab <- NULL
    if (!is.null(phyloseq::tax_table(ps, errorIfNULL = FALSE))) {
        tax_tab <- as.data.frame(phyloseq::tax_table(ps))
        tax_tab <- tibble::rownames_to_column(tax_tab, var = "otu_id")
    }

    # Preleva la otu_table e assicurati che le colonne siano campioni (samples in colonne)
    otu_mat <- as.matrix(phyloseq::otu_table(ps))
    if (phyloseq::taxa_are_rows(ps)) {
        # taxa x samples -> trasponi per avere samples in colonne
        otu_mat <- t(otu_mat)
    }
    # ora otu_mat: rows = samples, cols = taxa? Wait: we want counts matrix taxa x samples (rows = taxa)
    # For edgeR we need matrix with rows = features (genes/taxa) and cols = samples:
    # after as.matrix(otu_table(ps)) when taxa_are_rows(ps) TRUE -> rows = taxa, cols = samples
    # after transpose above, we have rows = samples, cols = taxa -> so let's ensure rows = taxa
    # Better: get otu2 as as.matrix(otu_table(ps)); if taxa_are_rows TRUE keep it; else transpose
    otu2 <- as.matrix(phyloseq::otu_table(ps))
    if (!phyloseq::taxa_are_rows(ps)) {
        otu2 <- t(otu2)
    }
    counts <- otu2  # rows = taxa, cols = samples

    # Sanitize sample names to match metadata
    colnames(counts) <- gsub("[^[:alnum:]_]+", "_", colnames(counts))

    # --- 5) Read metadata (file path) or use provided dataframe ---
    if (is.character(metadata) && file.exists(metadata)) {
        meta <- read.csv(metadata,
                         header = TRUE,
                         sep = ",",
                         na.strings = c("", " ", "NA"),
                         check.names = TRUE,
                         stringsAsFactors = FALSE)
    } else if (is.data.frame(metadata)) {
        meta <- metadata
    } else {
        stop("metadata deve essere un path a csv o un data.frame.")
    }

    # sanitize metadata id and keep only samples present in counts
    if (!"id" %in% colnames(meta)) {
        stop("metadata deve contenere la colonna 'id' con i sample IDs.")
    }
    meta$id <- gsub("[^[:alnum:]_]+", "_", meta$id)
    samples_names <- colnames(counts)
    meta_sub <- meta[meta$id %in% samples_names, , drop = FALSE]

    if (nrow(meta_sub) == 0) {
        stop("Nessun sample del metadata corrisponde ai nomi colonne dei counts.")
    }

    # Imposta rownames del metadata per matching con design e ordina coerentemente
    rownames(meta_sub) <- meta_sub$id

    # Se era stato richiesto un subset tramite sub_name e subcategory, assicurati che metadata contenga solo quelli
    if (!is.null(subcategory) && subcategory != "none" && subcategory != "" && !is.null(sub_name)) {
        meta_sub <- meta_sub %>% dplyr::filter(.data[[subcategory]] %in% sub_name)
    }

    # --- 6) Costruzione metadata_hm usata per design (in base a analisys e discriminant) ---
    # Verifica che le colonne richieste esistano
    cols_needed <- unique(c(analisys, discriminant))
    missing_cols <- setdiff(cols_needed, colnames(meta_sub))
    if (length(missing_cols) > 0) {
        stop(paste("Le seguenti colonne metadata sono mancanti:", paste(missing_cols, collapse = ", ")))
    }

    # Scegli i campioni coerenti e ordina counts per metadata
    samples_keep <- rownames(meta_sub)
    counts <- counts[, samples_keep, drop = FALSE]

    # Converti la variabile di interesse in factor se serve
    meta_sub <- meta_sub %>% dplyr::mutate(across(all_of(analisys), as.factor))

    # Crea metadata per il design: prendi almeno analisys + discriminant
    design_meta <- meta_sub %>% dplyr::select(all_of(cols_needed))
    # Assicurati che rownames siano gli id
    rownames(design_meta) <- rownames(meta_sub)

    # --- 7) Costruzione formula e design matrix ---
    formula_str <- paste0("~ ", analisys)
    if (length(discriminant) > 0) {
        for (el in discriminant) {
            formula_str <- paste0(formula_str, " + ", el)
        }
    }
    design <- model.matrix(as.formula(formula_str), data = design_meta)

    # Controllo di coerenza design vs counts
    samples_design <- rownames(design)
    samples_counts <- colnames(counts)
    if (!all(samples_counts %in% samples_design) || !all(samples_design %in% samples_counts)) {
        # mantieni intersezione e riordina
        common <- intersect(samples_counts, samples_design)
        if (length(common) == 0) stop("Nessuna corrispondenza tra samples in counts e metadata/design.")
        counts <- counts[, common, drop = FALSE]
        design <- design[common, , drop = FALSE]
    }

    # --- 8) Esegui limma-voom ---

    library(edgeR)
    library(limma)
    counts<-as.data.frame(counts)

    dge <- edgeR::DGEList(counts = counts)
    # Calcola fattori di normalizzazione TMM

    dge <- edgeR::calcNormFactors(dge)

    # Voom (ottiene log-CPM + pesi)
    v <- limma::voom(dge, design = design, plot = FALSE)

    limma_fit <- limma::lmFit(v, design)
    limma_fit <- limma::eBayes(limma_fit)

    # --- 9) Estrai coefficienti stimabili ---
    if (is.null(limma_fit$coefficients)) stop("limma_fit$coefficients è NULL dopo lmFit.")
    estimable_coefs <- colnames(limma_fit$coefficients)[!is.na(limma_fit$coefficients[1, ])]
    # rimuovi intercept se presente dalle testate
    estimable_coefs_no_intercept <- estimable_coefs[estimable_coefs != "(Intercept)"]

    if (length(estimable_coefs_no_intercept) == 0) {
        stop("Nessun coefficiente stimabile (diverso dall'intercetta) nel modello LIMMA.")
    }

    # topTable: se più coefficienti testati -> usa test F per tutte le colonne,
    # altrimenti testa il singolo coefficiente
    TopTableClean <- function(fit,
                              coef = NULL,
                              number = Inf,
                              adjust.method = "BH",
                              sort.by = "P",
                              p.value = 1) {
        if (!inherits(fit, "MArrayLM")) stop("fit must be an MArrayLM object")
        if (is.null(fit$coefficients)) stop("No coefficients found in fit object")

        if (is.null(coef)) {
            # default: test F su tutti i coefficienti tranne intercept
            tt <- limma::topTable(fit, coef = NULL, number = number,
                                  adjust.method = adjust.method, sort.by = sort.by, p.value = p.value)
            tested <- colnames(fit$coefficients)[colnames(fit$coefficients) != "(Intercept)"]
        } else if (length(coef) == 1) {
            tt <- limma::topTable(fit, coef = coef, number = number,
                                  adjust.method = adjust.method, sort.by = sort.by, p.value = p.value)
            tested <- coef
        } else {
            # se ci sono piu' coef, esegui un test F generale (coef = NULL) ma annota quali sono
            tt <- limma::topTable(fit, coef = NULL, number = number,
                                  adjust.method = adjust.method, sort.by = sort.by, p.value = p.value)
            tested <- coef
        }

        tt$Tested_Coefficients <- paste(tested, collapse = "; ")
        return(tt)
    }

    # scegli cosa passare: se hai più coefficienti stimabili diversi da intercept -> usa NULL (F-test)
    coef_for_test <- if (length(estimable_coefs_no_intercept) == 1) estimable_coefs_no_intercept else NULL

    top_table <- TopTableClean(limma_fit,
                               coef = coef_for_test,
                               number = Inf,
                               adjust.method = "fdr",
                               sort.by = "F")

    # --- 10) Filtri e salvataggio risultati ---
    # Forma nomi output
    write.csv(top_table, file = gsub(".csv", "_toptable.csv", fullname), row.names = TRUE)

    # Filtra per p-value e adjusted p-value
    # top_table ha colonne: "P.Value" e "adj.P.Val" quando topTable è chiamata
    top_table_p <- dplyr::filter(top_table, .data$P.Value < 0.05)

    taxa<-as.data.frame(phyloseq::tax_table(ps))
    taxa$Genus_Species <- paste(taxa$Genus, taxa$Species, sep = " ")
    taxa<-taxa %>%
        select(Genus_Species)
    taxa$taxid <- rownames(taxa)
    top_table_p$taxid <- rownames(top_table_p)
    table<-merge(top_table_p, taxa, by="taxid")
    table_fin <- table %>%
        select(Genus_Species, taxid, P.Value)
    # salva file finali
    write.csv(table_fin, file = gsub(".csv", "_adj.csv", fullname), row.names = TRUE)
    return(invisible(NULL))
}


source("Settings/Packages.R")
source("Settings/utilities.R")
source("Workflow/Step4_DiscriminantAnalysis/Limma/function.R")

output_folderMSHD <- "Output/LIMMA_voom/MSHD/"
createFolder(output_folderMSHD)
output_folder_001 <- "Output/LIMMA_voom/001/"
createFolder(output_folder_001)
output_folder_01 <- "Output/LIMMA_voom/01/"
createFolder(output_folder_01)
output_folder_05 <- "Output/LIMMA_voom/05/"
createFolder(output_folder_05)
library(tibble)

DASLimmaVoom(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "lesion_burden",
    subcategory = "gc_treatment",
    sub_name = c("positive", "negative"),
    output_folder = output_folder_001,
    discriminant = c("gc_treatment", "age", "sex", "bmi")
)

DASLimmaVoom(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "lesion_burden",
    subcategory = "gc_treatment",
    sub_name = c("positive", "negative"),
    output_folder = output_folder_05,
    discriminant = c("gc_treatment", "age", "sex", "bmi")
)

DASLimmaVoom(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "lesion_burden",
    subcategory = "gc_treatment",
    sub_name = c("positive", "negative"),
    output_folder = output_folder_01,
    discriminant = c("gc_treatment", "age", "sex", "bmi")
)

DASLimmaVoom(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "gadolinium_contrast",
    subcategory = "gc_treatment",
    sub_name = c("positive", "negative"),
    output_folder = output_folder_001,
    discriminant = c("gc_treatment", "age", "sex", "bmi")
)
DASLimmaVoom(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "gadolinium_contrast",
    subcategory = "gc_treatment",
    sub_name = c("positive", "negative"),
    output_folder = output_folder_05,
    discriminant = c("gc_treatment", "age", "sex", "bmi")
)
DASLimmaVoom(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "gadolinium_contrast",
    subcategory = "gc_treatment",
    sub_name = c("positive", "negative"),
    output_folder = output_folder_01,
    discriminant = c("gc_treatment", "age", "sex", "bmi")

)
DASLimmaVoom(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "spinal_cord_lesion",
    subcategory = "gc_treatment",
    sub_name = c("positive", "negative"),
    output_folder = output_folder_001,
    discriminant = c("gc_treatment", "age", "sex", "bmi")
)
DASLimmaVoom(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "spinal_cord_lesion",
    subcategory = "gc_treatment",
    sub_name = c("positive", "negative"),
    output_folder = output_folder_05,
    discriminant = c("gc_treatment", "age", "sex", "bmi")
)
DASLimmaVoom(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "spinal_cord_lesion",
    subcategory = "gc_treatment",
    sub_name = c("positive", "negative"),
    output_folder = output_folder_01,
    discriminant = c("gc_treatment", "age", "sex", "bmi")
)
DASLimmaVoom(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "subtentorial_lesions",  
    subcategory = "gc_treatment",
    sub_name = c("positive", "negative"),
    output_folder = output_folder_001,
    discriminant = c("gc_treatment", "age", "sex", "bmi")
)
DASLimmaVoom(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "subtentorial_lesions",  
    subcategory = "gc_treatment",
    sub_name = c("positive", "negative"),
    output_folder = output_folder_05,
    discriminant = c("gc_treatment", "age", "sex", "bmi")
)
DASLimmaVoom(
    baselines_dec = "Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds",
    metadata = "InputData/NewMetadataMS_fin.csv",
    analisys = "subtentorial_lesions",  
    subcategory = "gc_treatment",
    sub_name = c("positive", "negative"),   
    output_folder = output_folder_01,
    discriminant = c("gc_treatment", "age", "sex", "bmi")
)

# fine script
