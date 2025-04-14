
# Loading all the necessary 
source("Settings/utilities.R")
source("Settings/Packages.R")
source("Settings/Colorpalette.R")


# Step 0: Loading the data and metadata
source("Workflow/Step0_DataImport/DataImport.R")

# Step 1: use Deseq to do unsupervised decontamination
source("Workflow/Step1_Normalization/DeSeq.R")

# Step 2: use supervised decontamination
source("Workflow/Step2_Decontamination/Supervised_decontam.R")

# Step 3: create graphs of Alpha, Beta diversity and Stackbars
#Alpha
source("Workflow/Step3_Graphs/Alpha.R") 
#Beta
source("Workflow/Step3_Graphs/Beta.R")
#Stackbars
source("Workflow/Step3_Graphs/Stack.R")
#patchwork
source("Workflow/Step3_Graphs/Patchwork.R")

# Step 4: creation of DAS using limma and lefse
#limma
source("Workflow/Step4_DiscriminantAnalysis/Limma/DAS_LIMMA.R")
#lefse
source("Workflow/Step4_DiscriminantAnalysis/Lefse/PreLefse.R")
#furthermore for lefse you need other steps see README
#fuse DAS
source("Workflow/Step4_DiscriminantAnalysis/Merge_DAS/Merge_DAS.R")

# Step 5: create of heatmaps
source("Workflow/Step5_Heatmaps/HeatmapNew.R")

# Step 6 calculate Alpha diversity DAS
source("Workflow/Step6_DAS_Alpha/DAS_ALPHA.R")

# Step 7: LEAPS?
source("Workflow/Step7_Leaps/script_LEAPS_nuovo_perpaper.R")
