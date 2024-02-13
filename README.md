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


# Overview


![Overview](/images/Figure1.png)

**Overview:** The figure above illustrates the data types and potential applications of the CMDMultiOmics app.  Data types illustrated in the shaded areas to the left and right are contained in separate data repositories and loadable across tabular pages to explore individual data types between datasets as illustrated in the center section. Once loaded, exploratory analysis can be performed on selected datasets including visualization of raw expression profiles or differential expression results across a custom selection of datasets and genes. Alternatively, individual datasets and genes can be visualized with volcano plots.  Additional applications of the loaded data include applications to AI driven analytics after data download.  




![Features](/images/Figure2.png)

App Features: The figure above illustrates the vistalization tools featured in the app. (A-B) Violin plots of the expression levels of a single gene across a single (A) or several datasets (B). (C) Volcano plot of the deferential abundance levels of a gene, metabolite, or protein. (D) A table illustrating the deferential abundance levels of a gene, metabolite, or protein. (E) heatmap indicating the log2 fold change levels of selected genes across selected datasets. (F) Bar plot ordered from left to right first by the number of times a gene or protein was identified in a selected dataset and second by the average expression level. (G) Table indicating the position and score of a differentially methylated gene.  Color coding at the top of the figure indicates which app exploration modalities, here illustrated column wise across panels A-F, are available for each data type.  While all analytic tools re available for transcriptomics data, volcano plots (C) and differential abundance tables (D) are available for metabolomic data, and volcano plots (C) and differential abundance tables (D) heatmaps (E) and bar plots (F) are available for proteomics data.










