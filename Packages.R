############################################################
#######AUTOMATIZZARE SCRIPT################################
############################################################

source("Colorpalette.R")
library(dplyr)
library(ggplot2)
library(phyloseq)
library(ggsignif)
library(ggsignif)
library(patchwork)
library(ape)
library(DESeq2)
library(microbiome)
library(tidyr)
library(MetBrewer)
library(ggbreak)
library(tibble)
library(limma)
library(microbiome)
library(readr)
library(pheatmap)
library(microbiomeutilities)
library(ggpubr)
library(grid)
library(leaps)
library(UpSetR)
library(NbClust)
library(ggfortify)
library(factoextra)
library(cluster)
library(ggalt)
library(gridExtra)
library(FDRestimation)


# Lista delle librerie richieste
packages <- c(
  "dplyr", "ggplot2", "phyloseq", "ggsignif", "patchwork", "ape", "DESeq2", 
  "microbiome", "tidyr", "MetBrewer", "ggbreak", "tibble", "limma", 
  "readr", "pheatmap", "microbiomeutilities", "ggpubr", "grid", 
  "leaps", "UpSetR"
)

# Funzione per installare e caricare i pacchetti
install_and_load <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installando il pacchetto:", pkg))
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# Itera sulla lista e verifica/installazione/caricamento delle librerie
for (pkg in packages) {
  install_and_load(pkg)
}
###############controllare lista
# Lista dei pacchetti necessari
packages <- c(
  "dplyr", "ggplot2", "phyloseq", "ggsignif", "patchwork",
  "ape", "DESeq2", "microbiome", "tidyr", "MetBrewer", 
  "ggbreak", "tibble"
)

# Installazione dei pacchetti mancanti
install_missing <- function(packages) {
  installed <- installed.packages()[, "Package"]
  to_install <- packages[!packages %in% installed]
  if (length(to_install) > 0) {
    install.packages(to_install, dependencies = TRUE)
  } else {
    message("Tutti i pacchetti sono già installati.")
  }
}

# Installare Bioconductor per phyloseq e DESeq2, se necessario
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("phyloseq", "DESeq2"), update = FALSE,force=TRUE)

# Installazione dei pacchetti restanti
install_missing(packages)

# Verifica di eventuali problemi
suppressMessages(lapply(packages, library, character.only = TRUE))
# Forzare la reinstallazione di DESeq2
BiocManager::install("DESeq2", force = TRUE)

BiocManager::install("microbiome")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma")
install.packages("readr")
install.packages("pheatmap")

install.packages(c("phyloseq", "ggplot2", "dplyr", "tibble"))
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("joey711/phyloseq")
install.packages(c("nloptr", "lme4", "pbkrtest", "car", "rstatix", "ggpubr"))
remotes::install_github("microsud/microbiomeutilities")

install.packages(rstatix)


install.packages("ggpubr")
install.packages("rstatix")  # Install rstatix if not already installed
library(rstatix)
install.packages("factoextra", repos="http://cran.rstudio.com/")