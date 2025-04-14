############################################################
#######AUTOMATIZZARE SCRIPT################################
############################################################

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


# CRAN packages
cran_packages <- c(
  "dplyr", "ggplot2", "tidyr", "tibble", "patchwork", "gridExtra", 
  "ggsignif", "ggbreak", "MetBrewer", "ggalt", "cluster", "NbClust",
  "ggfortify", "factoextra", "UpSetR", "leaps", "readr", "grid", 
  "ggpubr", "rstatix", "pheatmap", "nloptr", "lme4", "pbkrtest", "car"
)

# Bioconductor Packages
bioc_packages <- c("phyloseq", "DESeq2", "microbiome", "limma")

# GitHub Packages
github_packages <- list(
  phyloseq = "joey711/phyloseq",
  microbiomeutilities = "microsud/microbiomeutilities"
)

#  CRAN installation
install_if_missing_cran <- function(pkgs) {
  installed <- rownames(installed.packages())
  to_install <- pkgs[!pkgs %in% installed]
  if (length(to_install) > 0) {
    install.packages(to_install, dependencies = TRUE)
  }
}

# Bioconductor installation
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
install_if_missing_bioc <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, update = FALSE, force = TRUE)
    }
  }
}

# GitHub installation
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
install_github_packages <- function(pkgs) {
  for (pkg in names(pkgs)) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      remotes::install_github(pkgs[[pkg]])
    }
  }
}

install_if_missing_cran(cran_packages)
install_if_missing_bioc(bioc_packages)
install_github_packages(github_packages)

# Check installation
main_packages <- c(cran_packages, bioc_packages, names(github_packages))
suppressMessages(invisible(lapply(main_packages, function(p) {
  try(library(p, character.only = TRUE), silent = TRUE)
})))