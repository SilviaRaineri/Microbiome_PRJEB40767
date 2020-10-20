# Microbiome_PRJEB40767
Analysis from a shotgun metagenomics experiment of gut microbiome samples from obese female rats treated with a panel of weightloss drugs for 42 Days. Samples were collected at Day -3 (before the start of the treatment) and Day 42 (end of the treatment). Results, methods and experimental design will be published in a dedicated paper (manuscript in preparation). 

##  Dependencies
The analysis was carried out using the R packages listed below:
- R = 3.6.1 
- Deseq2 = 1.24.0
- Phyloseq = 1.28.0 
- fgsea = 1.10.1
- tidyverse = 1.3.0 
- reshape2 = 1.4.4
- car = 3.0-9
At the end of each rmd script, a sessionInfo output can be found, with further specifications.

## Scripts 
The repository contains five different scripts, each focusing on a different portion of the analysis. 
- Weight_diversity =  weight loss and gut microbiome diversity analyses 
- OGTT = Oral Glucose Tolerance Test results
- Food_intake 
- Diff_species = differentially abundant species and gut microbiome composition
- Diff_genes = differentially abundant genes and fgsea (functional analysis)
NB: Before running each script, make sure to change the working directory path in the first block of code of each rmd script.

## Running the analysis
To run the analysis you should follow the steps below:

### Clone this repository 
`git clone https://github.com/SilviaRaineri/Microbiome_PRJEB40767.git` 

### Create a working directory
`mkdir gut_microbiome_analysis` 

### Copy relevant files to working directory
```
cp <path-to-Microbiome_PRJEB40767>/*rmd # insert specific name of the script you want to run
cp <path-to-Microbiome_PRJEB40767>/config.R # contains links to data tables stored on figshare
cp <path-to-Microbiome_PRJEB40767>/R/* # contains script with helper functions
cp <path-to-Microbiome_PRJEB40767>/Data/* # Contains files needed in Diff_genes (Species_colors) or Weight_diversity (july19_map and RATFAT04_weights)
```

### Render markdown
Once all the necessary files are downloaded in the appropriate folder(s) and the proper working directory is indicated in the .Rmd script of interest, one should be able to reproduce the figures and results detailed in the corresponding html report (see Microbiome_PRJEB40767/reports). 

`rmarkdown::render("Diff_species.Rmd", output_format="html_document")` 
