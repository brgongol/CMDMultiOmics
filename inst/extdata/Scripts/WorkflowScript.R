
############################
#### set up environment ####
############################
install.packages("path to tar file/MultiOmics_0.1.0.tar.gz", repos = NULL, type="source")
library(MultiOmics); library(data.table); library(shiny)
homedir <- "path to database/MultiOmicsAnalysis"
makeDirectory(homedir)
setwd(homedir)

#####################################
#### Obtain platform information ####
#####################################
PlatAnnotInfo <- PlatformAnnotationLoad(PlatInfo=fread(file.path(homedir, "OverviewFiles", "GEOPlatformInfo.csv"), header = TRUE))
fwrite(PlatAnnotInfo, file.path(homedir, "OverviewFiles", "GEOPlatformAnnotation.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
PlatAnnotInfo <- fread(file.path(homedir, "OverviewFiles", "GEOPlatformAnnotation.txt"))

############################
#### Load overview file ####
############################
overview <- fread(file.path(homedir, "OverviewFiles", "GEODataOverview3.csv"), header = TRUE)
overview <- overview[GEO2R == TRUE,]
#### check for duplicates assay names ####
overview$AssayNames <-gsub(" ", "", paste(overview$ID, overview$FCColumnNames, overview$Tissue, overview$Disease, sep = "_"))
overview[duplicated(AssayNames),][!AssayNames == "__",] # check should return nothing

################################
#### Download data from GEO ####
################################
Compiled <- GEOCompile(DS=overview$ID,
                       gpl=overview$GPLNumber,
                       gsm=overview$ComparisonVector,
                       namestr=overview$FCColumnNames,
                       nameraw=overview$RawColumnNames,
                       PlatAnnotInfo = PlatAnnotInfo,
                       destdir = file.path(homedir, "GEOcache"),
                       filename = NULL,
                       writeDB=TRUE,
                       writeRaw=overview$DownloadRawData,
                       GenerateMetaData=overview$DownloadMetaData,
                       MetaDataPath = file.path(homedir, "ArrayMetaData"),
                       writeMetaData=TRUE,
                       DBPath=file.path(homedir, "AppData"),
                       Technology = overview$Technology,
                       renameRaw = FALSE,
                       subsetRaw = FALSE)

#################################################################################
#### double check that the fold change calculations were performed correctly ####
#################################################################################
checkFile <- GEO2RDirectionCheck(DBPath=file.path(homedir, "AppData"),
                                 DS=overview$ID,
                                 namestr=overview$FCColumnNames,
                                 gsm=overview$ComparisonVector,
                                 Technology = overview$Technology,
                                 GraphPath = file.path(homedir, "DirectionCheck"),
                                 subsetRaw = FALSE)
checkFile[[1]][correlation < 0,]
fwrite(checkFile[[1]], file.path("./OverviewFiles/ArrayCheck.xls"), row.names = FALSE, quote = FALSE, sep = "\t")

####################################################################################################
#### Create or update Summarized experiment objects that integrate all data sets: DEG data only #### Checks for comparisons that are already in the database
#################################################################################################### Duplicated comparisons are not deleted from "5_AppData"
(DESE <- DESEGenerate(DEGDatapath=file.path(homedir, "AppData"), SEPath = file.path(homedir, "ProcessFiles", "SumarizedExp_DB.rds") ))
(DESE <- readRDS(file.path(homedir, "ProcessFiles", "SumarizedExp_DB.rds")))

#########################################
#### Create or update Raw data files #### Expression levels of duplicated genes are averaged. (Usually occurs with datasets generated using array technology if multiple probes measure the expression level of a gene.)
######################################### Note: sleep is only required if using a computer that synchronizes the compiled table to a cloud server.
#### Add Array data ####
RawDataCompile(Fpath = file.path(homedir, "AppData"), outPath = file.path("./ProcessFiles/RawData.txt"), CheckIdenticalData = FALSE, StartAt = 1, sleep = 60, DataType = "Array")
RawArrayComplete <- fread(file.path("./ProcessFiles/RawData.txt"))

#############################################
#### remove files from staging directory #### # Optional
#############################################
# files <- list.files(file.path(homedir, "AppData"))
# # files <- files[grepl(".rds", files) & grepl("_Raw", files)] # un-comment to remove only raw files
# for(b in 1:length(files)){fn <- file.path(homedir, "AppData", files[b]); if (file.exists(fn)){ file.remove(fn)}}

############################################
#### Create meta data for Raw data file #### # note RNAseq RunInfo files need to be downloaded manually and stored in the
############################################ # RunInfoFiles directory with the naming convention: "GSEXXXXX_SraRunTable.txt".
metaDF <- MetDataCompile(RNAseqFilePath= file.path(homedir, "RunInfo"), ArrayFilePath= file.path(homedir, "ArrayMetaData"), overview = overview)

#####################################################
#### Create raw data summarizedExperiment object ####
#####################################################
library(SummarizedExperiment)
RawSE <- GenerateRawSE(df = metaDF, ArrayDT=RawArrayComplete)
saveRDS(RawSE, file=file.path(homedir, "ProcessFiles", "SumarizedExp_RawDB.rds"))
RawSE <- readRDS(file.path(homedir, "ProcessFiles", "SumarizedExp_RawDB.rds"))
assay(RawSE)[1:5,1:5]
rowData(RawSE)
colData(RawSE)

################################################################################
################################################################################
#### Proteomic data workflow ###################################################
################################################################################
################################################################################
library(DEP); library(readxl); library(readr); library(stringr)
library(rpx); library(RforProteomics); library(BiocFileCache); library(data.table); library(readxl); library(textclean)

############################
#### Load overview file ####
############################
overview <- fread(file.path(homedir, "OverviewFiles", "GEODataOverview3.csv"), header = TRUE)
overviewpProteomics <- overview[DataType == "Proteomic",]

##################################
#### Download proteomics data ####
##################################
# devtools::install_version("dbplyr", version = "2.3.4")
# ProteomicsDataDownload(path = file.path(homedir, "Proteomic_1"), DS = overviewpProteomics$ID)

############################################
#### Loop through files and format data ####
############################################
library(makeunique)
# MZList <- FormatMaxQuant(path = file.path(homedir, "Proteomic_1"))
# names(MZList)

########################################
#### Save formatted proteomics data ####
########################################
# proteomicMZSave(MZList=MZList, path = file.path(homedir, "Proteomic_2"))

###########################################
#### Create experimental design object ####
###########################################
DesignDT <- DesignMatrixFromNames(Fpath = file.path(homedir, "Proteomic_2"))
fwrite(DesignDT, file.path(homedir, "OverviewFiles", "DesignMatrix.xls"), row.names = FALSE, quote = FALSE, sep = "\t")

################################################
#### load into summarized experiment object ####
################################################
library(SummarizedExperiment)
tissueSplitList <- ProtSELoad(DesignDT=DesignDT, Fpath=file.path(homedir, "Proteomic_2"))
assay(tissueSplitList[[2]])
rowData(tissueSplitList[[2]])

##################################################################
#### Perform overall protein abundance assessments ###############
##################################################################
#### compile data together to assess protein expression level ####
TotMel <- RowDataCompile(tissueSplitList=tissueSplitList)
TotMel

#############################################################################################
#### Filter for proteins that are identified in all replicates of at least one condition ####
#############################################################################################
library(radiant.data); library(tidyr); library(dplyr); library(tidyverse)
dataFiltList <- DataFilter(dataSeList=tissueSplitList, thr = 0)
assay(dataFiltList[[3]])
rowData(dataFiltList[[3]])

###############################################
#### Perform N-peptides per protein cutoff ####
###############################################
PepCutOffList <- NPeptideThreshold(dataFiltList=dataFiltList, Npeptides = 2)
PepCutOffList[[2]]                  # return the percent of records remaining
dataPepCutOff <- PepCutOffList[[1]] # return the data
colnames(assay(dataPepCutOff[[1]]))

###############################
#### Perform normalization ####
###############################
library(NormalyzerDE)
MultiNormalizeList <- MultiNormalization(dataPepCutOff=dataPepCutOff)
names(MultiNormalizeList)

#### create normalization density plots ####
############################################
densityPlotList <- densityPlotFromList(MultiNormalizeList)
names(densityPlotList)
densityPlotList[["mean"]]
densityPlotList[["median"]]
densityPlotList[["vsn"]]
densityPlotList[["DEPvsn"]]
densityPlotList[["loess"]]
densityPlotList[["rlr"]]
densityPlotList[["smad"]]

###############################################################
#### Select normalization method ##############################
###############################################################
NormalizedSE <- MultiNormalizeList[["median"]]
assay(NormalizedSE[[3]])
rowData(NormalizedSE[[3]])

#############################
#### Impute missing data ####
#############################
#### Determine if any data sets have missing values that need to be imputed ####
(Missing <- DetermineMising(data=NormalizedSE))
#### Impute missing values ####
set.seed(1)
NormImpAll <- DataImpute(dataFiltList = NormalizedSE, type = "MinProb")

##################################################################
#### remove samples if there are less than 50 samples in them ####
##################################################################
VarRMThreshList <- LowSampleCountRmove(VarRMList=NormImpAll, cut = 50)

#############################################################################
#### calculate Fold changes on imputed, thresholded, and normalized data ####
#############################################################################
FCcut <- log2(1); Pcut <- 0.05; sigCol = "BHCorrection" # "p.val" "p.adj" "BHCorrection"
Comparisons <- fread(file.path(homedir, "OverviewFiles", "ProteomicComparisons.xls"), header = TRUE)
compList2 <- Comparisons$Comparison
names(compList2) <- Comparisons$dataset
compList2 <- compList2[names(compList2) %in% names(VarRMThreshList)]
dataDiffNorm <- DEAnalysis(DataList = VarRMThreshList, type = "manual", ComparisonList=compList2)
dataDiffNorm

###############################
#### Save data to database ####
###############################
SaveToProteomicDB(SEList=dataDiffNorm, Path=file.path(homedir, "Proteomic_3"))

#################################################
#### Load selected Sumarized experiment file ####
#################################################
files <- list.files(file.path(homedir, "Proteomic_3"))
SelectedList <- ProteomicSELoad(Path=file.path(homedir, "Proteomic_3"), Fname=files[1])

###################################################################
#### Denote significant proteins based on user defined cutoffs ####
###################################################################
DEPList <- SigDEAnnotate(dataDiff=dataDiffNorm, alpha = Pcut, lfc = FCcut, sigCol= sigCol)
DEPList[[1]]
colData(DEPList[[1]])
rowData(DEPList[[1]])
assay(DEPList[[1]])

#######################
#### volcano plots ####
#######################
if(sigCol == "p.adj"){adj <- TRUE; BHadj <- FALSE
} else if(sigCol == "p.val"){ adj <- FALSE; BHadj <- FALSE
} else if(sigCol == "BHCorrection"){ adj <- FALSE; BHadj <- TRUE}
PlotList <- VolcanoPlot(DEPList = DEPList, ComparisonList=compList2, ymin = 0, ymax = 15, xmin = -7.5, xmax = 7.5, addNames = FALSE, Adjusted = adj, BHadjust = BHadj)
PlotList[[1]]

################################################################################
################################################################################
#### Load app ##################################################################
################################################################################
################################################################################
Appetup(homedir)
runApp("./Scripts")







