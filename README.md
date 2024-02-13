# CMDMultiOmics

The `CMDMultiOmics` package provides an expandable framework for the assembly, storage, integration, and visualization of transcriptomic, proteomic, methylome, and metabolome datasets. Currently the `CMDMultiOmics` app houses custom selection of data that are related to a group of cardiometabolic diseases from trascriptomics, proteomics, methylation or metabolomic data types. This shiny application is paired with transcript expression and proteomics processing pipelines for the addition of user selected data sets. 

# Installation

## 1) Install dependencies
install.packages(c("plotly", "shinyWidgets", "data.table", "feather", "textclean", "dplyr", "ggplot2", "RColorBrewer", "stringr", "readxl", "R.utils", "umap", "tidyr", "tidyverse", "makeunique", "radiant.data", "pander")) <br> 
 <br> 
BiocManager::install(c("Biobase", "limma", "DESeq2", "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "GEOquery", "mia", "clariomdhumantranscriptcluster.db", "hugene10sttranscriptcluster.db", "hugene20sttranscriptcluster.db", "hugene11sttranscriptcluster.db", "hgu133plus2.db", "hwgcod.db", "illuminaHumanv4.db", "hta20transcriptcluster.db", "hgu133a.db", "AnnotationDbi", "rpx", "RforProteomics", "BiocFileCache", "DEP", "NormalyzerDE", "SummarizedExperiment", "rtracklayer")) <br> 

## 2) Install app
devtools::install_github("brgongol/CMDMultiOmics") <br> 

## 3) Load vignette
library(pander) <br>
openFileInOS(system.file("extdata", "Vignette/CMDMultiOmicsVignette.html", package="MultiOmics")) <br>

## 4) Execute app deployment
library(CMDMultiOmics) <br> 
homedir <- "path to database/CMDMultiOmicsAnalysis" <br> 
makeDirectory(homedir) <br> 
setwd(homedir) <br> 
Open and execute the "WorkflowScript.R" file in the "Scripts" directory. <br> 




![alt text](https://github.com/brgongol/CMDMultiOmics/images/Figure2.png?raw=true)










