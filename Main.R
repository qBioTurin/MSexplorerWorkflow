
# Loading all the necessary
source("Settings/utilities.R")
source("Settings/Packages.R")
# add install.packages('VennDiagram')

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
source("Workflow/Step3_Graphs/BetaNOhv.R")
#Stackbars
source("Workflow/Step3_Graphs/Stack.R")
#patchwork
source("Workflow/Step3_Graphs/Patchwork.R")

# Step 4: creation of DAS using Maaslin3
source("Workflow/Step4_DiscriminantAnalysis/Maaslin/Maaslin3_funct.R")
source("Workflow/Step4_DiscriminantAnalysis/Maaslin/MaaslinLaunch.R")

# Step 5: create of heatmaps
source("Workflow/Step5_Heatmaps/Heatmap.R")
source("Workflow/Step5_Heatmaps/heatmap_Venn_patchmod.R")

# Step 6 calculate Alpha diversity DAS
source("Workflow/Step6_DAS_Alpha/DAS_ALPHA.R")

# Step 7: LEAPS
source("Workflow/Step7_Leaps/LEAPS.R")
