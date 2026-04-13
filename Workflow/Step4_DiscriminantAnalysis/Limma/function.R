###########################################################################################################################
### MAIN FUNCTION
###########################################################################################################################
#' Process BIOM and Metadata Files with Decontamination and Phyloseq Conversion
#'
#' This function processes microbiome data using a Kraken pipeline. It imports
#' a BIOM file and associated metadata, performs DESeq2-based supervised decontamination,
#' and saves the results as `.rds` files.
#'
#' @param biom_file File input object (from `shiny::fileInput`) representing the BIOM file.
#' @param metadata_file File input object representing the metadata CSV file.
#' @param Bacteria_level Integer. Taxonomic level to retain for Bacteria (e.g., 6 for Genus).
#' @param Archaea_level Integer. Taxonomic level to retain for Archaea.
#' @param Eukaryota_level Integer. Taxonomic level to retain for Eukaryota.
#' @param selectDesign Character. The column in the metadata used as experimental design.
#' @param out_dir Character. Base output directory path where results will be saved.
#'
#' @return A list containing:
#' \describe{
#'   \item{list_of_decontam}{A list of filtered phyloseq objects after supervised DESeq2 decontamination.}
#'   \item{tab_list}{The imported phyloseq object list created from the raw data.}
#' }
#'
#' @seealso \code{\link{KrakenDataImport}}, \code{\link{KrakenDeSeq}}, \code{\link{createFolder}}

process_file <- function(biom_file, metadata_file, Bacteria_level = NULL, Archaea_level = NULL, Eukaryota_level = NULL, selectDesign = NULL, out_dir, do_normalization = FALSE, tab_list) {
  output_folder <- paste0(out_dir, "RDS/")
  output_folderD <- paste0(out_dir, "DESEQ_RDS/")
  createFolder(output_folderD)
  allPath <- paste0(output_folder, "/ALL_baseline_phylo.rds")
  list_of_decontam <- KrakenDeSeq(allPath, output_folderD, Bacteria_level, Archaea_level, Eukaryota_level, selectDesign, do_normalization)


  return(list(list_of_decontam = list_of_decontam, tab_list = tab_list))
}

#' Create a Dynamic Checkbox Group Input
#'
#' This function creates a `checkboxGroupInput` UI element only if the choices are valid.
#' It pre-selects values that match the globally defined `selected_genera` vector.
#'
#' @param id Character. The input ID to assign to the checkbox group.
#' @param label Character. The label to display above the checkbox group.
#' @param choices Character vector. The available choices to display as checkboxes.
#'
#' @return A `shiny` UI object (checkbox group) or `NULL` if `choices` is empty.
#'
#' @note This function depends on a global variable `selected_genera` being defined.
#'
#' @seealso \code{\link[shiny]{checkboxGroupInput}}

createCheckboxGroup <- function(id, label, choices) {
  if (is.null(choices) || length(choices) == 0) {
    return(NULL)
  }

  checkboxGroupInput(
    id, label,
    choices = choices, # Mostra sempre tutti i generi
    selected = choices[choices %in% selected_genera] # Seleziona solo quelli presenti in selected_genera
  )
}

#' Create Differential Abundance Results using LEfSe and/or LIMMA
#'
#' This function performs differential abundance analysis using LEfSe and/or LIMMA
#' on decontaminated microbiome data, based on the specified pipeline (`kraken` or `metaphlan`).
#' It merges the results into a single file for downstream analysis.
#'
#' @param analisys Character. Type of analysis performed (e.g., "disease_vs_control").
#' @param sub_category Character. Name of the metadata column used for subsetting (e.g., "Treatment").
#' @param sub_name Character vector. Values within `sub_category` to subset (e.g., c("Treated", "Control")).
#' @param program Character. Which method(s) to use: `"lefse"`, `"limma"`, or `"both"`.
#' @param kingdom Integer. Kingdom index: 1 = Archaea, 2 = Bacteria, 3 = Eukaryota.
#' @param pipeline Character. Pipeline used: `"kraken"` or `"metaphlan"`.
#' @param output_folder Character. Path to the output root folder.
#'
#' @return The path to the merged result file as returned by `DasMerge()`.
#'
#' @details The function:
#' \itemize{
#'   \item Loads supervised decontaminated abundance data and matching alpha diversity metadata.
#'   \item Performs differential abundance analysis using LEfSe and/or LIMMA depending on `program`.
#'   \item Creates an empty table for the method not selected (to ensure consistent merge).
#'   \item Calls `DasMerge()` to create a unified table and saves it to a `.rds` file.
#' }
#'
#' @seealso \code{\link{DASLefse}}, \code{\link{DASLimma}}, \code{\link{DasMerge}}


create_das <- function(analisys, sub_category, sub_name, program, kingdom, pipeline, output_folder, metadata, discriminants = NULL) {
  # Percorsi di output
  output_folder_lefse <- paste0(output_folder, "/LEFSE/step0/")
  output_folder_limma <- paste0(output_folder, "/LIMMA/")
  output_folder_maaslin3 <- paste0(output_folder, "/MAASLIN3/")
  output_folder_merge <- paste0(output_folder, "/MERGE_DAS/")

  kingdom <- as.integer(kingdom)
  if (!dir.exists(output_folder_merge)) dir.create(output_folder_merge, recursive = TRUE)

  # File path per pipeline
  if (pipeline == "kraken") {
    baselines_dec <- c(
      paste0(output_folder, "/SUPERVISED_DEC/Archaea_Supervised_decontam.rds"),
      paste0(output_folder, "/SUPERVISED_DEC/Bacteria_Supervised_decontam.rds"),
      paste0(output_folder, "/SUPERVISED_DEC/Eukaryota_Supervised_decontam.rds")
    )
    metadata_path <- c(
      paste0(output_folder, "/ALPHA/Archaea_alpha_metadata.csv"),
      paste0(output_folder, "/ALPHA/Bacteria_alpha_metadata.csv"),
      paste0(output_folder, "/ALPHA/Eukaryota_alpha_metadata.csv")
    )
  } else if (pipeline == "metaphlan") {
    baselines_dec <- c(
      paste0(output_folder, "/META_SUPERVISED_DEC/Archaea_Supervised_meta_no_strain.rds"),
      paste0(output_folder, "/META_SUPERVISED_DEC/Bacteria_Supervised_meta_no_strain.rds"),
      paste0(output_folder, "/META_SUPERVISED_DEC/Eukaryota_Supervised_meta_no_strain.rds")
    )
    metadata_path <- c(
      paste0(output_folder, "/ALPHA/Archaea_meta_alpha_metadata.csv"),
      paste0(output_folder, "/ALPHA/Bacteria_meta_alpha_metadata.csv"),
      paste0(output_folder, "/ALPHA/Eukaryota_meta_alpha_metadata.csv")
    )
  } else {
    stop("Pipeline not recognized. Use 'kraken' or 'metaphlan'.")
  }

  metadata_hm <- read.csv(metadata_path[kingdom],
    header = TRUE,
    sep = ",",
    na = c("", " ", "NA"),
    check.names = TRUE
  )
  if (!is.null(discriminants)) {
    # trim whitespace and remove empty strings
    discriminants <- trimws(discriminants)
    discriminants <- discriminants[discriminants != ""]
    if (length(discriminants) == 0) discriminants <- NULL
  }
  if (is.null(discriminants)) {
    print("No discriminants provided")
  } else {
    print(paste("Discriminants:", paste(discriminants, collapse = ", ")))
  }
  domain <- find_domain(baselines_dec[kingdom])

  # --- LEfSe ---
  lefse_df <- tryCatch(
    {
      if (grepl("lefse", program, ignore.case = TRUE)) {
        lefse_path <- DASLefse(
          baselines_dec[kingdom], metadata_path[kingdom],
          analisys, sub_category, sub_name,
          output_folder_lefse, discriminants
        )
        read_tsv(lefse_path)
      } else {
        stop("Skipping LEfSe")
      }
    },
    error = function(e) {
      message("Warning: LEfSe failed - ", e$message)
      names <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
      as.data.frame(matrix(ncol = length(names), nrow = 0, dimnames = list(NULL, names)))
    }
  )

  # --- LIMMA ---
  limma_df <- tryCatch(
    {
      if (grepl("limma", program, ignore.case = TRUE)) {
        limma_path <- DASLimma(
          baselines_dec[kingdom], metadata_path[kingdom],
          analisys, sub_category, sub_name,
          output_folder_limma, discriminants
        )
        read.csv(limma_path)
      } else {
        stop("Skipping LIMMA")
      }
    },
    error = function(e) {
      message("Warning: LIMMA failed - ", e$message)
      ids <- c(metadata_hm$id)
      names <- c("X", ids)
      as.data.frame(matrix(ncol = length(names), nrow = 0, dimnames = list(NULL, names)))
    }
  )

  # --- Maaslin3 ---
  maaslin3_df <- tryCatch(
    {
      if (grepl("maaslin3", program, ignore.case = TRUE)) {
        DASMaaslin3(
          rds_path = baselines_dec[kingdom],
          metadata_path = metadata_path[kingdom],
          main_el = analisys,
          discr_els = discriminants,
          categorySel = sub_category,
          elements = sub_name,
          output_folder = output_folder_maaslin3
        )
      } else {
        stop("Skipping Maaslin3")
      }
    },
    error = function(e) {
      message("Warning: Maaslin3 failed - ", e$message)
      names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
      as.data.frame(matrix(ncol = length(names), nrow = 0, dimnames = list(NULL, names)))
    }
  )

  # --- Merge finale ---
  final_output <- gsub(" ", "", paste0(
    output_folder_merge, "fuse_", domain, "_",
    sub_name_creator(sub_name), "_", sub_category, "_",
    analisys, ".rds"
  ))
  
  # Limita la lunghezza del nome file (mantenendo path e estensione)
  max_filename_length <- 122
  dir_path <- dirname(final_output)
  filename <- basename(final_output)
  
  if (nchar(filename) > max_filename_length) {
    extension <- ".rds"
    max_name_length <- max_filename_length - nchar(extension)
    filename <- paste0(substr(filename, 1, max_name_length), extension)
    final_output <- file.path(dir_path, filename)
  }

  return(DasMerge(
    baselines_dec[kingdom], lefse_df, limma_df, maaslin3_df,
    final_output, sub_category, sub_name
  ))
}

find_domain <- function(path) {
  domain <- "Unknown"
  if (grepl("bact", path, ignore.case = TRUE)) {
    domain <- "bacteria"
  } else if (grepl("euk", path, ignore.case = TRUE)) {
    domain <- "eukaryote"
  } else if (grepl("arch", path, ignore.case = TRUE)) {
    domain <- "archaea"
  }
  return(domain)
}

sub_name_creator <- function(sub_name) {
  if (length(sub_name) > 1) {
    subname_paste <- paste(sub_name, collapse = "_")
  } else {
    subname_paste <- sub_name
  }
  return(subname_paste)
}
change_type <- function(variable, df, levels, colTypes, ord_factor) {
  column_class <- colTypes[variable]
  col_data <- df[[variable]]

  if (column_class == "character") {
    new_col <- as.character(col_data)
  } else if (column_class == "date") {
    new_col <- suppressWarnings(as.Date(col_data, format = "%Y-%m-%d"))
    if (all(is.na(new_col))) {
      new_col <- suppressWarnings(as.Date(col_data, format = "%d/%m/%Y"))
    }
    if (all(is.na(new_col))) {
      new_col <- suppressWarnings(as.Date(col_data, format = "%m/%d/%Y"))
    }
  } else if (column_class == "numeric") {
    new_col <- suppressWarnings(as.numeric(col_data))
  } else if (column_class == "factor") {
    if (ord_factor[variable]) {
      new_col <- factor(col_data,
        ordered = TRUE,
        levels = levels[[variable]]
      )
    } else {
      new_col <- as.factor(col_data)
    }
  } else if (column_class == "logic") {
    if (any(is.na(suppressWarnings(as.logical(levels[[variable]]))))) {
      new_col <- as.factor(col_data)
    } else {
      new_col <- suppressWarnings(as.logical(col_data))
    }
  } else {
    # fallback di sicurezza se qualcosa va storto
    new_col <- as.character(col_data)
  }

  return(new_col)
}

callback <- c(
  "var colnames = table.columns().header().to$().map(function() {",
  "return this.innerHTML;",
  "}).get();",
  "Shiny.onInputChange('colnames', colnames);",
  "table.on('dblclick.dt', 'thead th', function(e) {",
  "var $th = $(this);",
  "var index = $th.index();",
  "var colname = $th.text(), newcolname = colname;",
  "var $input = $('<input type=\"text\">')",
  "$input.val(colname);",
  "$th.empty().append($input);",
  "$input.on('change', function() {",
  "newcolname = $input.val();",
  "if(newcolname != colname){",
  "$(table.column(index).header()).text(newcolname);",
  "colnames[index] = newcolname;",
  "Shiny.onInputChange('colnames', colnames);",
  "}",
  "$input.remove();",
  "}).on('blur', function() {",
  "$(table.column(index).header()).text(newcolname);",
  "$input.remove();",
  "});",
  "})"
)

golden <- (1 + sqrt(5)) / 2
