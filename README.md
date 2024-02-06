# MultiOmics

Step by step instructions for app installation: <br> 

#### 1) Install dependencies
install.packages(c("plotly", "shinyWidgets", "data.table", "feather", "textclean", "dplyr", "ggplot2", "RColorBrewer", "stringr" <br> 
	"readxl", "R.utils", "umap", "tidyr", "tidyverse", "makeunique", "radiant.data")) <br> 
 <br> 
BiocManager::install(c("Biobase", "limma", "DESeq2", "org.Hs.eg.db", "org.Mm.eg.db", "GEOquery", "mia", "clariomdhumantranscriptcluster.db", <br> 
	"hugene10sttranscriptcluster.db", "hugene20sttranscriptcluster.db", "hugene11sttranscriptcluster.db", "hgu133plus2.db", <br> 
	"hwgcod.db", "illuminaHumanv4.db", "hta20transcriptcluster.db", "hgu133a.db", "AnnotationDbi", "rpx", "RforProteomics", "BiocFileCache", <br> 
	"DEP", "NormalyzerDE", "SummarizedExperiment")) <br> 

#### 2) Install app
devtools::install_github("brgongol/MultiOmics") <br> 

#### 3) Execute app deployment
library(MultiOmics) <br> 
homedir <- "path to database/MultiOmicsAnalysis" <br> 
makeDirectory(homedir) <br> 
setwd(homedir) <br> 
Open and execute the "WorkflowScript.R" file in the "Scripts" directory. <br> 


