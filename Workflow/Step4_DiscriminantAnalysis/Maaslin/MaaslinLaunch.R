source("Workflow/Step4_DiscriminantAnalysis/Maaslin/Maaslin22.R")
source("Settings/Packages.R")
#source("Workflow/Step4_DiscriminantAnalysis/Maaslin/Maaslin3_funct.R")

rds_0 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.rds")
library(dplyr)

Maaslin3(
    rds = rds_0, output_folder = "Output/MAASLIN3/0_disc/0_lesion_burden_onlygc_plus_gc",
    main_el = "lesion_burden", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "0"
)

Maaslin3(
    rds = rds_0, output_folder = "Output/MAASLIN3/0_disc/0_spinal_cord_lesion_onlygc_plus_gc",
    main_el = "spinal_cord_lesion", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "0"
)

Maaslin3(
    rds = rds_0, output_folder = "Output/MAASLIN3/0_disc/0_gadolinium_contrast_onlygc_plus_gc",
    main_el = "gadolinium_contrast", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "0"
)

Maaslin3(
    rds = rds_0, output_folder = "Output/MAASLIN3/0_disc/0_subtentorial_lesions_onlygc_plus_gc",
    main_el = "subtentorial_lesions", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "0"
)












rds_001 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")

Maaslin3(
    rds = rds_001, output_folder = "Output/MAASLIN3/001_disc/001_lesion_burden_onlygc_plus_gc",
    main_el = "lesion_burden", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "001"
)

Maaslin3(
    rds = rds_001, output_folder = "Output/MAASLIN3/001_disc/001_spinal_cord_lesion_onlygc_plus_gc",
    main_el = "spinal_cord_lesion", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "001"
)

Maaslin3(
    rds = rds_001, output_folder = "Output/MAASLIN3/001_disc/001_gadolinium_contrast_onlygc_plus_gc",
    main_el = "gadolinium_contrast", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "001"
)

Maaslin3(
    rds = rds_001, output_folder = "Output/MAASLIN3/001_disc/001_subtentorial_lesions_onlygc_plus_gc",
    main_el = "subtentorial_lesions", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "001"
)







rds_01 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds")

Maaslin3(
    rds = rds_01, output_folder = "Output/MAASLIN3/01_disc/01_lesion_burden_onlygc_plus_gc",
    main_el = "lesion_burden", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "01"
)

Maaslin3(
    rds = rds_01, output_folder = "Output/MAASLIN3/01_disc/01_spinal_cord_lesion_onlygc_plus_gc",
    main_el = "spinal_cord_lesion", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "01"
)

Maaslin3(
    rds = rds_01, output_folder = "Output/MAASLIN3/01_disc/01_gadolinium_contrast_onlygc_plus_gc",
    main_el = "gadolinium_contrast", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "01"
)

Maaslin3(
    rds = rds_01, output_folder = "Output/MAASLIN3/01_disc/01_subtentorial_lesions_onlygc_plus_gc",
    main_el = "subtentorial_lesions", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "01"
)

rds_05 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds")

Maaslin3(
    rds = rds_05, output_folder = "Output/MAASLIN3/05_disc/05_lesion_burden_onlygc_plus_gc",
    main_el = "lesion_burden", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
    perc = "05"
)

Maaslin3(
    rds = rds_05, output_folder = "Output/MAASLIN3/05_disc/05_spinal_cord_lesion_onlygc_plus_gc",
    main_el = "spinal_cord_lesion", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
    perc = "05"
)

Maaslin3(
    rds = rds_05, output_folder = "Output/MAASLIN3/05_disc/05_gadolinium_contrast_onlygc_plus_gc",
    main_el = "gadolinium_contrast", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
    perc = "05"
)

Maaslin3(
    rds = rds_05, output_folder = "Output/MAASLIN3/05_disc/05_subtentorial_lesions_onlygc_plus_gc",
    main_el = "subtentorial_lesions", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
    perc = "05"
)


rds_Euk <- readRDS("Output/SUPERVISED_DEC/Eukaryote_Supervised_decontam0.rds")
Maaslin3(
    rds = rds_Euk, output_folder = "Output/MAASLIN3/Eukaryote_MSHD_disc",
    main_el = "category", categorySel = "category", elements = c("HEALTHY", "MS"), discr_els = c("sex", "age", "bmi")
)
rds_Arch <- readRDS("Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.rds")
Maaslin3(
    rds = rds_Arch, output_folder = "Output/MAASLIN3/Archaea_MSHD_disc",
    main_el = "category", categorySel = "category", elements = c("HEALTHY", "MS"), discr_els = c("sex", "age", "bmi")
)
Maaslin3(
    rds = rds_0, output_folder = "Output/MAASLIN3/Bacteria_MSHD_disc",
    main_el = "category", categorySel = "category", elements = c("HEALTHY", "MS"), discr_els = c("sex", "age", "bmi")
)
Maaslin3(
    rds = rds_001, output_folder = "Output/MAASLIN3/Bacteria_MSHD_disc_001",
    main_el = "category", categorySel = "category", elements = c("HEALTHY", "MS"), discr_els = c("sex", "age", "bmi")
)

########################### false


#rds_001 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
library(dplyr)

Maaslin3(
    rds = rds_0, output_folder = "Output/MAASLIN3/0_disc/0_lesion_burden_onlygc_F",
    main_el = "lesion_burden", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
    perc = "0", median_comparison_abundance_in = FALSE
)

Maaslin3(
    rds = rds_0, output_folder = "Output/MAASLIN3/0_disc/0_spinal_cord_lesion_onlygc_F",
    main_el = "spinal_cord_lesion", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
    perc = "0", median_comparison_abundance_in = FALSE
)

Maaslin3(
    rds = rds_0, output_folder = "Output/MAASLIN3/0_disc/0_gadolinium_contrast_onlygc_F",
    main_el = "gadolinium_contrast", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
    perc = "0", median_comparison_abundance_in = FALSE
)
Maaslin3(
    rds = rds_0, output_folder = "Output/MAASLIN3/0_disc/001_subtentorial_lesions_onlygc_F",
    main_el = "subtentorial_lesions", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
    perc = "0", median_comparison_abundance_in = FALSE
)




Maaslin3(
    rds = rds_001, output_folder = "Output/MAASLIN3/001_disc/001_lesion_burden_onlygc_F",
    main_el = "lesion_burden", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
    perc = "001", median_comparison_abundance_in = FALSE
)

Maaslin3(
    rds = rds_001, output_folder = "Output/MAASLIN3/001_disc/001_spinal_cord_lesion_onlygc_F",
    main_el = "spinal_cord_lesion", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
    perc = "001", median_comparison_abundance_in = FALSE
)

Maaslin3(
    rds = rds_001, output_folder = "Output/MAASLIN3/001_disc/001_gadolinium_contrast_onlygc_F",
    main_el = "gadolinium_contrast", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
    perc = "001", median_comparison_abundance_in = FALSE
)
Maaslin3(
    rds = rds_001, output_folder = "Output/MAASLIN3/001_disc/001_subtentorial_lesions_onlygc_F",
    main_el = "subtentorial_lesions", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
    perc = "001", median_comparison_abundance_in = FALSE
)
#rds_01 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.01.rds")

Maaslin3(
    rds = rds_01, output_folder = "Output/MAASLIN3/01_disc/01_lesion_burden_onlygc_F",
    main_el = "lesion_burden", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
    perc = "01", median_comparison_abundance_in = FALSE
)

Maaslin3(
    rds = rds_01, output_folder = "Output/MAASLIN3/01_disc/01_spinal_cord_lesion_onlygc_F",
    main_el = "spinal_cord_lesion", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
    perc = "01", median_comparison_abundance_in = FALSE
)
Maaslin3(
    rds = rds_01, output_folder = "Output/MAASLIN3/01_disc/01_gadolinium_contrast_onlygc_F",
    main_el = "gadolinium_contrast", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
    perc = "01", median_comparison_abundance_in = FALSE
)
Maaslin3(
    rds = rds_01, output_folder = "Output/MAASLIN3/01_disc/01_subtentorial_lesions_onlygc_F",
    main_el = "subtentorial_lesions", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
    perc = "01", median_comparison_abundance_in = FALSE
)
#rds_05 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.05.rds")

Maaslin3(
    rds = rds_05, output_folder = "Output/MAASLIN3/05_disc/05_lesion_burden_onlygc_F",
    main_el = "lesion_burden", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
    perc = "05", median_comparison_abundance_in = FALSE
)

Maaslin3(
    rds = rds_05, output_folder = "Output/MAASLIN3/05_disc/05_spinal_cord_lesion_onlygc_F",
    main_el = "spinal_cord_lesion", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
    perc = "05", median_comparison_abundance_in = FALSE
)
Maaslin3(
    rds = rds_05, output_folder = "Output/MAASLIN3/05_disc/05_gadolinium_contrast_onlygc_F",
    main_el = "gadolinium_contrast", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
    perc = "05", median_comparison_abundance_in = FALSE
)
Maaslin3(
    rds = rds_05, output_folder = "Output/MAASLIN3/05_disc/05_subtentorial_lesions_onlygc_F",
    main_el = "subtentorial_lesions", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("age", "sex", "bmi"),
    perc = "05", median_comparison_abundance_in = FALSE
)

rds_Euk <- readRDS("Output/SUPERVISED_DEC/Eukaryote_Supervised_decontam0.001.rds")
Maaslin3(
    rds = rds_Euk, output_folder = "Output/MAASLIN3/Rebuttal/Euk_001",
    main_el = "category", categorySel = "category", elements = c("HEALTHY", "MS"), discr_els = c("sex", "age", "bmi"),median_comparison_abundance_in = TRUE
)

rds_Arch <- readRDS("Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.001.rds")

Maaslin3(
    rds = rds_Arch, output_folder = "Output/MAASLIN3/Rebuttal/Arch_001",
    main_el = "category", categorySel = "category", elements = c("HEALTHY", "MS"), discr_els = c("sex", "age", "bmi")
)
Maaslin3(
    rds = rds_001, output_folder = "Output/MAASLIN3/Rebuttal/Bact_001",
    main_el = "category", categorySel = "category", elements = c("HEALTHY", "MS"), discr_els = c("sex", "age", "bmi")
)






rds_001 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
metadata<- data.frame(sample_data(rds_001))
print(unique(metadata$Cluster))
colnames(rds_001$metadata)
Maaslin3Fin(
    rds = rds_001, output_folder = "Output/MAASLIN3/001_disc_cluster/001_cluster_gc",
    main_el = "Cluster", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("sex", "age", "bmi"),
    perc = "001"
)
library(maaslin3)
