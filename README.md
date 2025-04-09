# MSexplorerWorkflow

This is the post-processing workflow we created and used for our pipeline, which includes Kraken2, Bracken, and KrakenBio. It consists of an initial decontamination step, followed by various plots that help provide a better understanding of the results obtained from the taxonomic classifier.

# Steps

In the following steps, I will do my best to explain all the processes required to recreate our post-processing workflow, from the initial requirements to the final stages. If anything doesn’t work for you, please feel free to open an issue on GitHub.

## Pre-condition

This is the essential framework needed to begin the post-processing step. Everything included here is strictly necessary and mandatory for the successful execution of the next phase, so please be methodical and precise.

### [R](https://www.r-project.org/)

To use this pipeline, you need a BIOM file created from the output of the Bracken program* and a TSV file containing the corresponding metadata.

To properly use the programs included in your R environment, you must ensure that all the libraries listed in the Packages.R file are installed. To achieve this, simply run the file; it will check if all dependencies are present and, if not, install them automatically (NOTE: automatic installation is not yet implemented, so you need to install them manually). This process is fairly time- and resource-intensive, so make sure you have a stable internet connection. If you are working on a laptop, ensure the charger is connected.

Once you have all the necessary files, you can begin with the DataImport step. In this step, you must set the path to the BIOM file and the metadata file. Be careful to assign the correct type to each variable. Lastly, if needed, remove excess text in the Bracken file column names—for example, in our case, we removed "_bracken_species" to ensure the names match those in the metadata file.

Next, you can remove a list of samples that are not relevant to your research. This step is necessary for our study but is completely optional. Finally, an RDS file is created for each Kingdom, which will be required in subsequent steps.

*The BIOM file must result from the merging of multiple Bracken files, and the column names must correspond to the names of the individual Bracken files.


### [Lefse](https://huttenhower.sph.harvard.edu/lefse/)

In order to install LEfSe, you need Python 2, and depending on your operating system, this could be a bit of a challenge since it has been deprecated for over half a decade now. Anyway, in the following steps, I will guide you through its installation on Debian 12 (this should be quite similar for other Unix-like systems):

First and foremost, you need Conda installed on your machine. To do this, simply follow their [Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

Now, install python2.7 and create a symbolic link called python2(not necesssary but very handy) using the following two commands:
```bash
sudo apt-get install python2.7
sudo ln -s /usr/bin/python2.7 /usr/bin/python2
```
Create the environment for lefse using **explicitly** python 2.7, then activate it 
```bash
conda create -n lefse_env python=2.7 -y
conda activate lefse_env
```
Now, add the default, conda-forge, and bioconda channels, and then install LEfSe with the last command:
```bash
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda install -c bioconda lefse -y
```
You are now done! The following lines are an example of how to use LEfSe:
```bash
conda activate lefse_env #only if you havent already
lefse-format_input.py file_input file_input.in -c 1 -s -1 -u 2 -o 1000000
run_lefse.py -l 2  file_input.in file_input.res
```
For any additional info or issues, please refer to the SegataLab [github page](https://github.com/SegataLab/lefse) or [biobackery tutorial](https://github.com/biobakery/biobakery/wiki/lefse)

**IMPORTANT**The next few steps will have a folder for each one containing all the necessary files. To avoid running into errors, it is recommended to load all the libraries listed in the [Package.R](https://github.com/franky2204/postProcess/blob/main/Packages.R). One last thing: every path is relative, so the paths do not need to be modified unless the user wants to use a custom input. **IMPORTANT**


**TODO** update the Package file to automatize it
## [Step 0 Data import](https://github.com/franky2204/postProcess/blob/main/(0)DataImport/DataImport.R)

This script takes as input the biom file created by [bracken](https://github.com/jenniferlu717/Bracken) and the metadata of your files in CSV format. Once uploaded, the user must define a type for each column. This step is only necessary if the metadata differs from the one provided in this script. Additionally, if a custom input is used, it may be necessary to change the OTU table column names.

Lastly, for this dataset, we remove the naive patients, who will not appear in the future scripts.

Once completed, 4 RDS files will be created:

- Bacteria
- Archaea
- Eukaryota
- All 3 together

## [Step 1 Normalization + Unsupervised decontamination(Deseq)](https://github.com/franky2204/postProcess/blob/main/(1)Normalization/DeSeq.R)

In this step, we take as input the combined RDS generated by the previous step and pass it through[DEseq](https://github.com/thelovelab/DESeq2), a statistical tool used for analyzing DNA sequencing data. It identifies differentially expressed genes between conditions by modeling the counts of sequencing reads and accounting for variability in the data. Finally, we separate the 3 kingdoms to create one RDS file for each of them.

## [Step 2 Decontamination(Supervised)](https://github.com/franky2204/postProcess/blob/main/(2)Decontamination/Supervised_decontam.R)

This step takes as input the three RDS files created in the previous step and three lists of microorganisms associated with humans, for which we identified sources using the [Encyclopedia of Life](https://eol.org/). Additionally, we annotated each of the two lists to indicate whether a microorganism was found in only one source or in more than one. For archaea, the process is slightly more complex due to the limited studies available. The few sources, combined with decontamination at the phylum level, leave too many contaminants. To address this, we consulted an expert to obtain a more specific list of genera.

This script takes the aforementioned RDS files, removes controls and rows with ambiguous phylum/genus annotations (such as NULL values), and computes the prevalence for each feature, adding taxonomy and row counts. Lastly, it removes all phyla/genera not present in the tables and applies a relative abundance filter (>0.001) to eliminate additional potential contaminants. The program then saves this updated table as an RDS file (one for each domain).

## [Step 3 Graphs(Alfa&Beta+Stackbar)](https://github.com/franky2204/postProcess/blob/main/(3)Graphs)
In these four scripts, we focused on identifying and generating graphs that would not only be explanatory but also easy to interpret. The primary objective of this analysis was to provide a comprehensive description of the population we identified, highlighting any significant differences between patients undergoing glucocorticoid treatment and those who were not. By leveraging visual tools, we aimed to make complex data more accessible, allowing for a clearer understanding of the variations within the population.

Our approach was to select graph types that could effectively capture the nuances of the data while maintaining clarity. In particular, we were interested in how glucocorticoid treatment might influence the observed patterns, both at the individual and group levels. This allowed us to assess the treatment’s impact on various biological and clinical factors in a manner that was both statistically rigorous and visually intuitive.

Through the use of these graphs, we sought to make the results not only informative but also actionable, enabling further hypothesis testing or clinical applications. Ultimately, the goal was to ensure that the findings could be easily understood and communicated to a broad audience, whether researchers, clinicians, or stakeholders involved in the study.
### [Alpha diversity](https://github.com/franky2204/postProcess/blob/main/(3)Graphs/(3.a)Alpha.R)
This script takes as input the manually decontaminated RDS file and utilizes ggplot2 to create three distinct graphs. These graphs will be saved as RDS files and serve to visually demonstrate the differences in alpha diversity across various conditions. Specifically, the graphs will illustrate:

- The comparison in alpha diversity between healthy individuals and those with multiple sclerosis (MS).
- The impact of glucocorticoid (GC) treatment on the microbial diversity.
- The correlation between alpha diversity and the number of days since glucocorticoid treatment was administered.

By using ggplot2, the script ensures that these plots are not only informative but also aesthetically clear, allowing for easy interpretation of the relationships between the biological variables and treatment conditions. The resulting RDS files can be used for further analysis or incorporated into reports to communicate the findings effectively.
### [Beta diversity](https://github.com/franky2204/postProcess/blob/main/(3)Graphs/(3.a)Beta.R)
Similar to the alpha diversity analysis, this script takes the manually decontaminated RDS file as input and uses ggplot2 to generate a beta diversity graph. The graph will be created based on two key factors:

- The comparison between healthy individuals and those diagnosed with multiple sclerosis (MS).
- The effect of glucocorticoid (GC) treatment on the microbial community.

The resulting beta diversity plot will visually represent the dissimilarities between samples, providing insights into how these two conditions (health status and GC treatment) influence the microbiome structure. The use of ggplot2 ensures the graph is not only informative but also clear and easy to interpret for further analysis.
The final outputs will be two RDS file.
### [StackBar](https://github.com/franky2204/postProcess/blob/main/(3)Graphs/(3.a)Stack.R)

Stackbars are designed to visually represent the top 10 most abundant species found within each domain, making it easier to compare species distributions across different datasets. To construct these stackbars in a clear and comprehensible manner, we created an array that includes each of the top 10 species, along with an additional category labeled 'Other.' This category accounts for the sum of all species that do not appear in the top 10, ensuring that no relevant data is omitted while maintaining visual clarity.

To enhance readability and consistency, we standardized the color scheme by assigning the same shade of grey to the 'Other' category across all stackbars. This allows viewers to quickly differentiate between the most abundant species and the aggregated remaining species. Additionally, our array structure enables us to assign a distinct, predetermined color to each species. By ensuring that each species is consistently represented by a unique color, we improve both the interpretability and overall comprehensibility of the visualization. This approach makes it easier to compare results across different domains, identify patterns, and extract meaningful insights from the data.

### [Patchwork fuse graph](https://github.com/franky2204/postProcess/blob/main/(3)Graphs/(3.b)Patchwork_ggplot.R)

This step creates a visually harmonious and easy-to-read composition by integrating the previously generated graphs into a unified structure. It achieves this result by merging identical labels and efficiently managing space to arrange all elements into a structured table of parallel components.

To provide a brief recap, the composition includes a total of nine graphs, categorized as follows:

Bacteria:

- P1: Stack bar chart for bacterial species
- P2: Alpha diversity graphs for the bacterial category
- P3: Beta diversity graphs for the bacterial category

Archaea:

- P4: Stack bar chart for archaeal species
- P5: Alpha diversity graphs for the archaeal category
- P6: Beta diversity graphs for the archaeal category

Eukaryota:

- P7: Stack bar chart for eukaryotic species
- P8: Alpha diversity graphs for the eukaryotic category
- P9: Beta diversity graphs for the eukaryotic category

This structured approach enhances clarity, allowing for a more intuitive comparison across different domains while ensuring a cohesive and well-organized visual representation.



## [Step 4.1 Lefse](https://github.com/franky2204/postProcess/blob/main/(4.1)Lefse)

LEfSe is the acronime od Linear discrciminant ananlysis effect size and its used to determines the features most likely to explain difference between classes. The way it was implemented is:

### [Step 1 (4.1.a) - PreLefSe](https://github.com/franky2204/postProcess/blob/main/(4.1)Lefse/(4.1.a)PreLefse.R)

The PreLefSe script is an R-based preprocessing tool designed to prepare data for LEfSe (Linear Discriminant Analysis Effect Size) analysis. It takes as input:

Supervised decontaminated data (processed in Step 2)
Metadata (created in Step 3)
The script extracts the taxa table from the decontaminated file and assigns levels to a key categorical field, which will be used in the LEfSe pipeline. Based on this selected field, the script generates a structured table where:

The first row contains the chosen metadata field.
The second row lists the patient IDs.
The following rows contain phylogenetic classifications, followed by the relative abundance of each species.
Generated Tables and Field Selection
For this study, six different metadata fields are used, resulting in a total of 10 tables per domain, categorized as follows:

Primary Comparisons (Across All Patients)
- category – Comparison between Healthy vs. MS patients.
- gc_treatment – Comparison between GC-positive vs. GC-negative patients.

Subset Analyses (Separate Tables for GC-Positive and GC-Negative 
Patients)These four additional metadata fields are analyzed separately for GC-positive and GC-negative groups, leading to eight tables in total:

- lesion_burden – Analysis based on lesion burden classification.
- bone_marrow_lesion – Presence or absence of bone marrow lesions.
- gadolinium_contrast – Gadolinium contrast enhancement patterns.
- subtenorial_lesion – Lesions located in the subtenorial region.

This structured approach ensures that each comparison is appropriately categorized, facilitating meaningful insights from the LEfSe analysis.

### [Step 2 - lefseEX](https://github.com/franky2204/postProcess/blob/main/(4.1)Lefse/(4.1.b)LefseElaborator/lefseEx.sh)

The lefseEX script is a Bash-based automation tool designed to streamline the execution of LEfSe without manual intervention. It automatically runs the necessary commands to process the input data, perform statistical analysis, and generate results.

Requirements

The LEfSe Conda environment must be activated before running the script. If not activated, the script will not function correctly.
Customization Options
Users can easily modify input and output locations by adjusting the following script variables:

- DIRECTORY – Path to the input data folder.
- OUT1 – Path for intermediate output files.
- OUT2 – Path for final LEfSe results.

By automating the LEfSe pipeline, this script enhances efficiency, reducing manual errors and ensuring reproducibility across analyses.

In order to execute this file you must be in the postProcess folder(the base of the project) and the command must be : 
```bash
conda activate lefse_env
bash \(4.1\)Lefse/\(4.1.b\)LefseElaborator/lefseEx.sh 
```
### [Step 3 - main.py](https://github.com/franky2204/postProcess/blob/main/(4.1)Lefse/(4.1.b)LefseElaborator/main.py)

This step involves the initial manipulation of the LEfSe output, specifically focusing on two tasks:

**Trimming unimportant species.**

Filtering out species that do not meet the significance criteria.
Converting the phylogenetic tree format – Changing the format of the tree from pipe-separated (|) to tab-separated (\t) in order to create a TSV file that can later be merged with limma results.
Customization
To modify the input and output folder paths, simply adjust the variables input_folder and output_folder inside the script.
Like for the previous step you must be in the postProcess folder and execute:
```bash
/bin/python3 "(4.1)Lefse/(4.1.b)LefseElaborator/main.py"
```
**Expected Output**

Once the script has been run, it will output a series of files containing the phylogenetic paths of the relevant species. These files are essential for the subsequent analyses and will be used for further processing, including integration with the limma results.

**Important Considerations**

**Special Character Handling in LEfSe Output:**

LEfSe modifies all special characters in species names by replacing them with underscores (_). This could lead to discrepancies when matching species names across different datasets.
Therefore, it is imperative that the user manually review the output and replace any underscores with the correct special characters (e.g., spaces, hyphens, etc.).
To assist with this, you can refer to one of the RDS files containing the taxa tables to view the original species names and ensure proper character restoration.

NOTE: This step is crucial for maintaining data integrity and preventing mismatches in downstream analyses.

## [Step 4.2 Limma](https://github.com/franky2204/postProcess/blob/main/(4.2)Limma/DAS_LIMMA.R)

Limma is a powerful tool for data analysis, particularly suited for linear models and differential expression analysis in omics data. It leverages a systematic approach to identify significant changes between conditions by applying statistical models.

**Input and Processing**
Limma takes as input the supervised decontaminated data (from previous steps), performing the following steps:

- Normalization – The data is normalized to adjust for technical variation and ensure comparability across samples.
- Log Transformation – The data is then log-transformed (base 2) to stabilize variance and make the data more suitable for linear modeling.
Metadata Integration
The same metadata used in the LEfSe analysis is applied here, enabling the creation of a set of 10 tables per domain. These tables are generated for each of the following categories, mirroring the structure of LEfSe:

Primary comparisons (e.g., Healthy vs. MS, GC-positive vs. GC-negative)
Subset analyses (e.g., lesion burden, bone marrow lesion, gadolinium contrast, subtenorial lesion)

**Output**
The result of the Limma analysis will be a set of tables that highlight differential expression across domains, helping to identify key species or features that exhibit significant differences under various conditions. These tables will serve as the foundation for further interpretation and integration with other tools and analysis pipelines.

## [Step 5 Merge DAS](https://github.com/franky2204/postProcess/blob/main/(5)Merge_DAS/Merge_DAS.R)

This script takes the output from LEfSe and Limma and merges them into 10 distinct lists, making a merge of the two. The process works as follows:

The script retrieves all the species identified by LEfSe and Limma for each category (e.g., gc_treatment, lesion_burden, etc.), without filtering out species that appear in only one of the two results.
For each of the 10 categories (which correspond to the primary comparisons and subset analyses), it merges the species lists from both LEfSe and Limma into a single comprehensive list.
The script then searches for the abundance values of all these species in the supervised decontaminated data based on the coort(eg. in the GC will be selected only patience with gc_treatment postive)for the specific domain, assigning the appropriate abundance values to each species.
The final output consists of 10 RDS files, each containing the full list of species for a given category, along with their corresponding abundance values.
This process ensures that all species are considered, providing a complete view of the species distribution across various conditions, regardless of whether they appear in both the LEfSe and Limma results.

## [Step 6 Heatmaps](https://github.com/franky2204/postProcess/blob/main/(6)Heatmaps/HeatmapNew.R)
This script processes the DAS generated during the merge DAS step and creates a series of heatmaps designed to identify significant clusters. These clusters play a crucial role in refining the search for key species across different comparative analyses. By visually highlighting areas of interest, the heatmaps help researchers better interpret patterns and distributions within the dataset.

To run this script, the DAS files produced in Step 5 are required as input. The output consists of a series of PDF document that presents important elements, such as gadolinium and subtentorial regions, each represented with a distinct color. The color scale is carefully structured to enhance interpretability: white indicates zero presence, while progressively more saturated hues correspond to increasing concentrations. This visual representation allows for a more intuitive understanding of the data, making it easier to detect significant variations and trends.

## [Step 7 Leaps](https://github.com/franky2204/postProcess/blob/main/(7)Leaps/DAS_Leaps_ALE_alltaxatogetherNEW.R)


## [Step 8 DAS Alpha diversity](https://github.com/franky2204/postProcess/blob/main/(8)DAS_Alpha)
WIP

## [Step 9 Score ratio](https://github.com/franky2204/postProcess/blob/main/(9)Score_Ratio/F:B ratio.R)
WIP
