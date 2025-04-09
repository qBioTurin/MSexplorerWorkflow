library("tools")
source("utilities.R")
limma_folder = "Output/modLimma/"
lefse_folder = "Output/modLefse/"
createFolder(limma_folder)
createFolder(lefse_folder)


Limma_GC <- "Output/LIMMA/Bacteria_limma_gc_treatment.csv"

limma_Lesion<- "Output/LIMMA/Bacteria_limma_lesion_burden_both.csv"
limma_BM_Lesion <- "Output/LIMMA/Bacteria_limma_bone_marrow_lesions_both.csv"
limma_Gadolinium <- "Output/LIMMA/Bacteria_limma_gadolinium_contrast_both.csv"
limma_Subtentorial <- "Output/LIMMA/Bacteria_limma_subtentorial_lesions_both.csv"

lefse_GC <- "Output/LEFSE/final_output/BACT_GC.res"

lefse_Lesion <- "Output/LEFSE/final_outputv2/BACT_Lesion.res"
lefse_BM_Lesion <- "Output/LEFSE/final_outputv2/BACT_BM_Lesion.res"
lefse_Gadolinium <- "Output/LEFSE/final_outputv2/BACT_Gadolinium.res"
lefse_Subtentorial <- "Output/LEFSE/final_outputv2/BACT_Subtentorial.res"

remove_common_species_limma <- function(file1, file2, outputFolder){
    #Read csv file
    read1 <- read.csv(file1)
    read2 <- read.csv(file2)
    pre <- nrow(read1)
    # Find all common rows
    common_elements <- intersect(read1$X, read2$X)
    read1 <- read1[!read1$X %in% common_elements, ]

    # count nrows after sostitution
    post <- nrow(read1)

    # Print info
    print(paste("Removed Rows:", length(common_elements)))
    print("Common species removed:")
    print(common_elements)  # Stampa i nomi delle specie rimosse

    # Check length of removed rows
    if ((pre - post) != length(common_elements)) {
        print("Error: mismatch after row removal!")
    } else {
        file1 <- basename(file1) 
        filename_no_ext <- file_path_sans_ext(file1)

        path <- gsub(" ", "", paste0(outputFolder, filename_no_ext, "_mod.csv"))

        write.csv(read1, path, row.names = FALSE)
        print("Execution completed!")
    }
}

remove_common_species_lefse <- function(file1, file2, output_folder) {
    #Read tsv file
    
    read1 <- read.csv(file1)
    read2 <- read.csv(file2)
    
    # starting number of rows
    pre <- nrow(read1)

    # Find all common rows
    common_elements <- intersect(read1$X, read2$X)
    read1 <- read1[!read1$Species %in% common_elements, ]

    # count nrows after sostitution
    post <- nrow(read1)

    # Print info
    print(paste("Removed Rows:", length(common_elements)))
    print("Common species removed:")
    print(common_elements)  # Stampa i nomi delle specie rimosse

    # Check length of removed rows
    if ((pre - post) != length(common_elements)) {
        print("Error: mismatch after row removal!")
    } else {
        filename_no_ext <- file_path_sans_ext(basename(file1))

        # Output folder
        path <- file.path(output_folder, paste0(filename_no_ext, "_mod.rds"))

        write.table(read1, path, sep = "\t", row.names = FALSE, quote = FALSE)
        print("Execution completed!")
    }
}



remove_common_species_limma(limma_Lesion, limma_GC, limma_folder)
remove_common_species_limma(limma_BM_Lesion,limma_GC ,limma_folder)
remove_common_species_limma( limma_Gadolinium,limma_GC, limma_folder)
remove_common_species_limma(limma_Subtentorial,limma_GC , limma_folder)

remove_common_species_lefse(lefse_Lesion, lefse_GC, lefse_folder)
remove_common_species_lefse(lefse_BM_Lesion, lefse_GC, lefse_folder)
remove_common_species_lefse(lefse_Gadolinium, lefse_GC, lefse_folder)
remove_common_species_lefse(lefse_Subtentorial, lefse_GC, lefse_folder)
