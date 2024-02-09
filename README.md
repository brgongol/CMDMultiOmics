# MultiOmics

Step by step instructions for app installation. <br> 
The following three steps should be executed in a R terminal: <br> 

#### 1) Install dependencies
install.packages(c("plotly", "shinyWidgets", "data.table", "feather", "textclean", "dplyr", "ggplot2", "RColorBrewer", "stringr", "readxl", "R.utils", "umap", "tidyr", "tidyverse", "makeunique", "radiant.data")) <br> 
 <br> 
BiocManager::install(c("Biobase", "limma", "DESeq2", "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "GEOquery", "mia", "clariomdhumantranscriptcluster.db", "hugene10sttranscriptcluster.db", "hugene20sttranscriptcluster.db", "hugene11sttranscriptcluster.db", "hgu133plus2.db", "hwgcod.db", "illuminaHumanv4.db", "hta20transcriptcluster.db", "hgu133a.db", "AnnotationDbi", "rpx", "RforProteomics", "BiocFileCache", "DEP", "NormalyzerDE", "SummarizedExperiment", "rtracklayer")) <br> 

#### 2) Install app
devtools::install_github("brgongol/MultiOmics") <br> 

#### 3) Execute app deployment
library(MultiOmics) <br> 
homedir <- "path to database/MultiOmicsAnalysis" <br> 
makeDirectory(homedir) <br> 
setwd(homedir) <br> 
Open and execute the "WorkflowScript.R" file in the "Scripts" directory. <br> 


