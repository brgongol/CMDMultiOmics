
############################
#### set up environment ####
############################
library(MultiOmics); library(data.table)
homedir <- "path to database/MultiOmicsAnalysis"
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
                                 subsetRaw = FALSE,
                                 writeRaw=overview$DownloadRawData,
                                 RawQCPath = file.path(homedir, "RawQC"))
checkFile[[1]][correlation < 0,]
fwrite(checkFile[[1]], file.path("./OverviewFiles/ArrayCheck.xls"), row.names = FALSE, quote = FALSE, sep = "\t")

####################################################################################################
#### Create or update Summarized experiment objects that integrate all data sets: DEG data only #### Checks for comparisons that are already in the database
#################################################################################################### Duplicated comparisons are not deleted from "AppData"
(DESE <- DESEGenerate(DEGDatapath=file.path(homedir, "AppData"), SEPath = file.path(homedir, "ProcessFiles", "SumarizedExp_DB.rds") ))
(DESE <- readRDS(file.path(homedir, "ProcessFiles", "SumarizedExp_DB.rds")))

#########################################
#### Create or update Raw data files ####
#########################################

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

################################################################################
################################################################################
#### Load app ##################################################################
################################################################################
################################################################################
library(shiny)
AppSetup(homedir)
runApp("./Scripts")








