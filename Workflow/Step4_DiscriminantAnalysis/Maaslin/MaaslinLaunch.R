source("Workflow/Step4_DiscriminantAnalysis/Maaslin/Maaslin22.R")
source("Settings/Packages.R")
source("Workflow/Step4_DiscriminantAnalysis/Maaslin/Maaslin3.R")
rds_Euk <- readRDS("Output/SUPERVISED_DEC/Eukaryote_Supervised_decontam0.001.rds")
rds_001 <- readRDS("Output/SUPERVISED_DEC/Bacteria_Supervised_decontam0.001.rds")
rds_Arch <- readRDS("Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.001.rds")
library(maaslin3)


Maaslin3(
    rds = rds_001, output_folder = "Output/MAASLIN3/001_disc/001_lesion_burden_onlygc_plus_gc",
    main_el = "lesion_burden", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("gc_treatment", "age", "sex", "bmi"),
    perc = "001"
)
print(readRDS("Output/MAASLIN3/001_disc/001_lesion_burden_onlygc_plus_gc/lesion_burden_MAASLIN3_SR_features_001_T.rds"))

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


Maaslin3(
    rds = rds_Euk, output_folder = "Output/MAASLIN3/Eukaryote_MSHD_disc",
    main_el = "category", categorySel = "category", elements = c("HEALTHY", "MS"), discr_els = c("sex", "age", "bmi")
)
rds_Arch <- readRDS("Output/SUPERVISED_DEC/Archaea_Supervised_decontam0.001.rds")
Maaslin3(
    rds = rds_Arch, output_folder = "Output/MAASLIN3/Archaea_MSHD_disc",
    main_el = "category", categorySel = "category", elements = c("HEALTHY", "MS"), discr_els = c("sex", "age", "bmi")
)
Maaslin3(
    rds = rds_001, output_folder = "Output/MAASLIN3/Bacteria_MSHD_disc_001",
    main_el = "category", categorySel = "category", elements = c("HEALTHY", "MS"), discr_els = c("sex", "age", "bmi")
)

Maaslin3(
    rds = rds_001, output_folder = "Output/MAASLIN3/Bacteria_gc_treatment_disc",
    main_el = "gc_treatment", categorySel = "gc_treatment", elements = c("positive","negative"), discr_els = c("sex", "age", "bmi")
)
Maaslin3(
    rds = rds_Euk, output_folder = "Output/MAASLIN3/Eukaryote_gc_treatment_disc",
    main_el = "category", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("sex", "age", "bmi")
)
Maaslin3(
    rds = rds_Arch, output_folder = "Output/MAASLIN3/Archaea_gc_treatment_disc",
    main_el = "category", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("sex", "age", "bmi")
)


metadata<- data.frame(sample_data(rds_001))
print(unique(metadata$Cluster))

Maaslin3(
    rds = rds_001, output_folder = "Output/MAASLIN3/001_disc_cluster/001_cluster_gc",
    main_el = "Cluster", categorySel = "gc_treatment", elements = c("positive", "negative"), discr_els = c("sex", "age", "bmi"),
    perc = "001"
)

