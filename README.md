# WGCNexplore

### Hello and welcome!
#### This app was designed to facilitate user exploration and comparison of gene or transcript modules generated by Weighted Gene Co-Expression Analysis (WGCNA)", [(Langfelder & Horvath, BMC Bioinformatics 2008)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559). Note: the app does NOT run WGCNA for you; users are recommended to use the [WGCNA R package](https://cran.r-project.org/web/packages/WGCNA/index.html) tutorials for which are provided here: [WGCNA tutorials](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/). Example WGCNA scripts are also available in the scripts directory.

### Overview 
#### Weighted Gene Co-expression Network Analysis (WGCNA) is a widely used bioinformatics method that defines gene clusters (modules) based on expression levels, allowing the user to study the relationships between these modules and compare network topologies across expression datasets. This approach overcomes challenges of gene-based studies, but it is difficult to compare WGCNA results within & across datasets, as generated modules are highly dependent on user parameter inputs. As such, validation of modules based on module sizes (including grey module), number of modules in set, and biological data including cell-type and functional enrichment is an essential step in WGCNA analyses. We developed this application to facilitate and encourage user exploration of modules generated by their own WGCNA analyses of human transcriptomic data.

### Usage instructions
#### In this application, users may upload a CSV file delineating gene/transcript-module assignemnts as generated by WGCNA. Genes and transcripts should be indicated by ensembl ID. Data may be uploaded in wide format or long format:

![alt text](https://github.com/rlsmith1/WGCNexplore/blob/main/www/example_wide.png?raw=true)
![alt text](https://github.com/rlsmith1/WGCNexplore/blob/main/www/example_long.png?raw=true)

#### (These are just example data, not the result of any analysis). Note the ensembl ID column should be named either ensembl\_gene\_id or ensembl\_transcript\_id, depending on if the data are genes or transcripts, respectively. Module sets and modules may be named as the user desires.

#### This application compares module sets in terms of module number and size as well as hypergeometric overlap between each module of each set. Z-scored cell-type enrichment for each module is calculated using the [Lake et al (2018)](https://www.nature.com/articles/nbt.4038) cell-type dataset. The user also has the option to provide their own cell-type data, in the following long format: ![alt text](https://github.com/rlsmith1/WGCNexplore/blob/main/www/df_cell_type.png?raw=true) 

#### GO Functional enrichment using the [topGO R package](https://bioconductor.org/packages/release/bioc/html/topGO.html) (version 2.48.0.) is also performed on each module of interest, as well as genes that intersect between modules of choice
#### Once plots have been generated and functional enrichment has been run for a given module/module set combo, these plots & results are stored and do not need to be rerun (within a session). All plots are downloadable as PDFs and all data tables (hypergeometric test and functional enrichment results) are downloadable as CSVs. Please contact me with any comments, issues, suggestions, or feedback. Hope you find this useful!!

App initially generated with: R version 4.2.0, RStudio 2022.07.1+554 Release (7872775ebddc40635780ca1ed238934c3345c5de, 2022-07-22) for macOS, shiny version 1.7.2, shinydashboard version 0.7.2

