

#' Create directories
#'
#' @param homedir The directory where the database will be stored
#' @import R.utils
#' @export
makeDirectory <- function(homedir){
  dirs <- c(homedir,
            file.path(homedir, "ArraymetaData"),
            file.path(homedir, "AppData"),
            file.path(homedir, "GEOcache"),
            file.path(homedir, "DirectionCheck"),
            file.path(homedir, "Proteomic_1"),
            file.path(homedir, "Proteomic_3") )
  sapply(dirs, dir.create)
  SysFpath <- system.file("extdata/", package = "MultiOmics")
  from <- file.path(SysFpath, list.files(SysFpath))
  to <- file.path(homedir, list.files(SysFpath))
  sapply(seq_along(from), function(x){
    R.utils::copyDirectory(from=from[x], to=to[x], recursive=TRUE)
  }) }

#' Compile Array annotation information
#'
#' @param PlatInfo A csv file containing a table of the Array information
#' @import data.table
#' @import clariomdhumantranscriptcluster.db
#' @import hugene10sttranscriptcluster.db
#' @import hugene20sttranscriptcluster.db
#' @import hugene11sttranscriptcluster.db
#' @import hgu133plus2.db
#' @import hwgcod.db
#' @import illuminaHumanv4.db
#' @import hta20transcriptcluster.db
#' @import hgu133a.db
#' @import org.Hs.eg.db
#' @import GEOquery
#' @import dplyr
#' @export
PlatformAnnotationLoad <- function(PlatInfo){
  AnnotationDT <- data.table()
  for(i in 1:nrow(PlatInfo)){
    if(PlatInfo$GPLNumber[i] == "GPL23126"){
      # BiocManager::install("clariomdhumantranscriptcluster.db")
      library(clariomdhumantranscriptcluster.db)
      x<-clariomdhumantranscriptclusterGENENAME
      mapped_probes<-mappedkeys(x)
      GENENAME <- unlist(as.list(x[mapped_probes]))
      DTNames <- data.table(ID=names(GENENAME), GENENAME)
      x<-clariomdhumantranscriptclusterENTREZID
      mapped_probes<-mappedkeys(x)
      ENTREZ<-  unlist(as.list(x[mapped_probes]))
      DTENTREZ <- data.table(ID=names(ENTREZ), ENTREZID = ENTREZ)
      mer <- merge(DTNames, DTENTREZ, by = "ID")
      symbols <- as.character(mer$ENTREZ)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"),"ENTREZID" ) )
      GPL23126 <- merge(mer, goHuman, by = "ENTREZID")
      GPL23126 <- GPL23126[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL23126$GPLID <- "GPL23126"
      AnnotationDT <- rbind(AnnotationDT, GPL23126)
    } else if(PlatInfo$GPLNumber[i] == "GPL30511"){
      z <- getGEO("GPL30511")
      temp <- as.data.table(z@dataTable@table)
      setnames(temp, "Symbol", "SYMBOL")
      symbols <- temp$SYMBOL
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("ENTREZID"),"SYMBOL"))
      GPL30511 <- merge(temp, goHuman, by = "SYMBOL")
      GPL30511 <- GPL30511[,c("SYMBOL", "ID", "ENTREZID"), with = FALSE]
      GPL30511 <- GPL30511[!is.na(ENTREZID),]
      GPL30511$ENTREZID <- as.character(GPL30511$ENTREZID)
      Values <- GPL30511$ENTREZID
      mappings <- AnnotationDbi::select(org.Hs.eg.db, keys=Values, columns=c("GENENAME"),keytype="ENTREZID")
      GPL30511 <- merge(GPL30511, mappings, by = "ENTREZID")
      GPL30511 <- GPL30511[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL30511$GPLID <- "GPL30511"
      AnnotationDT <- rbind(AnnotationDT, GPL30511)
    } else if(PlatInfo$GPLNumber[i] == "GPL28577"){
      z <- getGEO("GPL28577")
      temp <- as.data.table(z@dataTable@table)
      symbols <- temp$ID
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("ENTREZID"),"SYMBOL"))
      GPL28577 <- cbind(temp, goHuman)
      GPL28577 <- GPL28577[!is.na(ENTREZID),]
      GPL28577$ENTREZID <- as.character(GPL28577$ENTREZID)
      Values <- GPL28577$ENTREZID
      mappings <- AnnotationDbi::select(org.Hs.eg.db, keys=Values, columns=c("GENENAME"),keytype="ENTREZID")
      GPL28577 <- merge(GPL28577, mappings, by = "ENTREZID")
      GPL28577 <- GPL28577[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL28577$GPLID <- "GPL28577"
      AnnotationDT <- rbind(AnnotationDT, GPL28577)
    } else if(PlatInfo$GPLNumber[i] == "GPL29503"){
      z <- getGEO("GPL29503")
      temp <- as.data.table(z@dataTable@table)
      symbols <- temp$ID
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("ENTREZID", "GENENAME"),"SYMBOL"))
      GPL29503 <- cbind(temp, goHuman)
      GPL29503 <- GPL29503[,c("ID", "ENTREZID", "SYMBOL", "GENENAME"), with = FALSE]
      GPL29503$GPLID <- "GPL29503"
      AnnotationDT <- rbind(AnnotationDT, GPL29503)
    } else if(PlatInfo$GPLNumber[i] == "GPL14951"){
      z <- getGEO("GPL14951")
      temp <- as.data.table(z@dataTable@table)
      setnames(temp, c("Symbol", "Entrez_Gene_ID"), c("SYMBOL", "ENTREZID"))
      GPL14951 <- temp
      GPL14951 <- GPL14951[,c("ID", "SYMBOL", "ENTREZID"), with = FALSE]
      GPL14951 <- unique(GPL14951[!is.na(ENTREZID),])
      GPL14951 <- GPL14951[!is.na(ENTREZID),]
      GPL14951$ENTREZID <- as.character(GPL14951$ENTREZID)
      Values <- GPL14951$ENTREZID
      mappings <- AnnotationDbi::select(org.Hs.eg.db, keys=Values, columns=c("GENENAME"),keytype="ENTREZID")
      GPL14951 <- merge(GPL14951, mappings, by = "ENTREZID")
      GPL14951 <- GPL14951[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL14951$GPLID <- "GPL14951"
      AnnotationDT <- rbind(AnnotationDT, GPL14951)
    } else if(PlatInfo$GPLNumber[i] == "GPL6244"){
      # BiocManager::install("hugene10sttranscriptcluster.db")
      library(hugene10sttranscriptcluster.db)
      x<-hugene10sttranscriptclusterGENENAME
      mapped_probes<-mappedkeys(x)
      GENENAME <- unlist(as.list(x[mapped_probes]))
      DTNames <- data.table(ID=names(GENENAME), GENENAME)
      x<-hugene10sttranscriptclusterENTREZID
      mapped_probes<-mappedkeys(x)
      ENTREZ<-  unlist(as.list(x[mapped_probes]))
      DTENTREZ <- data.table(ID=names(ENTREZ), ENTREZID = ENTREZ)
      mer <- merge(DTNames, DTENTREZ, by = "ID")
      symbols <- as.character(mer$ENTREZ)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"),"ENTREZID" ) )
      GPL6244 <- merge(mer, goHuman, by = "ENTREZID")
      GPL6244 <- GPL6244[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL6244$GPLID <- "GPL6244"
      AnnotationDT <- rbind(AnnotationDT, GPL6244)
    } else if(PlatInfo$GPLNumber[i] == "GPL16686"){
      # BiocManager::install("hugene20sttranscriptcluster.db")
      library(hugene20sttranscriptcluster.db)
      x<-hugene20sttranscriptclusterGENENAME
      mapped_probes<-mappedkeys(x)
      GENENAME <- unlist(as.list(x[mapped_probes]))
      DTNames <- data.table(ID=names(GENENAME), GENENAME)
      x<-hugene20sttranscriptclusterENTREZID
      mapped_probes<-mappedkeys(x)
      ENTREZ<-  unlist(as.list(x[mapped_probes]))
      DTENTREZ <- data.table(ID=names(ENTREZ), ENTREZID = ENTREZ)
      mer <- merge(DTNames, DTENTREZ, by = "ID")
      symbols <- as.character(mer$ENTREZ)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"),"ENTREZID" ) )
      GPL16686 <- merge(mer, goHuman, by = "ENTREZID")
      GPL16686 <- GPL16686[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL16686$GPLID <- "GPL16686"
      AnnotationDT <- rbind(AnnotationDT, GPL16686)
    } else if(PlatInfo$GPLNumber[i] == "GPL570"){
      z <- getGEO("GPL570")
      temp <- as.data.table(z@dataTable@table)
      temp <- temp[,c("ID", "Gene Symbol", "ENTREZ_GENE_ID"), with = FALSE]
      temp$`Gene Symbol` <- gsub(" ///.+", "", temp$`Gene Symbol`)
      temp$ENTREZ_GENE_ID <- gsub(" ///.+", "", temp$ENTREZ_GENE_ID)
      temp <- temp[!temp$`Gene Symbol` == "",]
      setnames(temp, c("Gene Symbol", "ENTREZ_GENE_ID"), c("SYMBOL", "ENTREZID"))
      GPL570 <- temp
      GPL570 <- GPL570[!is.na(ENTREZID),]
      GPL570 <- GPL570[!ENTREZID == "",]
      GPL570$ENTREZID <- as.character(GPL570$ENTREZID)
      Values <- GPL570$ENTREZID
      mappings <- AnnotationDbi::select(org.Hs.eg.db, keys=Values, columns=c("GENENAME"),keytype="ENTREZID")
      GPL570 <- unique(merge(GPL570, mappings, by = "ENTREZID", allow.cartesian = TRUE))
      GPL570 <- GPL570[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL570$GPLID <- "GPL570"
      AnnotationDT <- rbind(AnnotationDT, GPL570)
    } else if(PlatInfo$GPLNumber[i] == "GPL11532"){
      # BiocManager::install("hugene11sttranscriptcluster.db")
      library(hugene11sttranscriptcluster.db)
      x<-hugene11sttranscriptclusterGENENAME
      mapped_probes<-mappedkeys(x)
      GENENAME <- unlist(as.list(x[mapped_probes]))
      DTNames <- data.table(ID=names(GENENAME), GENENAME)
      x<-hugene11sttranscriptclusterENTREZID
      mapped_probes<-mappedkeys(x)
      ENTREZ<-  unlist(as.list(x[mapped_probes]))
      DTENTREZ <- data.table(ID=names(ENTREZ), ENTREZID = ENTREZ)
      mer <- merge(DTNames, DTENTREZ, by = "ID")
      symbols <- as.character(mer$ENTREZ)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"),"ENTREZID" ) )
      GPL11532 <- merge(mer, goHuman, by = "ENTREZID")
      GPL11532 <- GPL11532[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL11532$GPLID <- "GPL11532"
      AnnotationDT <- rbind(AnnotationDT, GPL11532)
    } else if(PlatInfo$GPLNumber[i] == "GPL14877"){
      # BiocManager::install("hgu133plus2.db")
      library(hgu133plus2.db)
      x<-hgu133plus2GENENAME
      mapped_probes<-mappedkeys(x)
      GENENAME <- unlist(as.list(x[mapped_probes]))
      DTNames <- data.table(ID=names(GENENAME), GENENAME)
      x<-hgu133plus2ENTREZID
      mapped_probes<-mappedkeys(x)
      ENTREZ<-  unlist(as.list(x[mapped_probes]))
      DTENTREZ <- data.table(ID=names(ENTREZ), ENTREZID = ENTREZ)
      mer <- merge(DTNames, DTENTREZ, by = "ID")
      symbols <- as.character(mer$ENTREZ)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"),"ENTREZID" ) )
      goHuman <- unique(goHuman)
      GPL14877 <- merge(mer, goHuman, by = "ENTREZID")
      GPL14877 <- GPL14877[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL14877$GPLID <- "GPL14877"
      AnnotationDT <- rbind(AnnotationDT, GPL14877)
    } else if(PlatInfo$GPLNumber[i] == "GPL2895"){
      # BiocManager::install("hwgcod.db")
      library(hwgcod.db)
      x<-hwgcodGENENAME
      mapped_probes<-mappedkeys(x)
      GENENAME <- unlist(as.list(x[mapped_probes]))
      DTNames <- data.table(ID=names(GENENAME), GENENAME)
      x<-hwgcodENTREZID
      mapped_probes<-mappedkeys(x)
      ENTREZ<-  unlist(as.list(x[mapped_probes]))
      DTENTREZ <- data.table(ID=names(ENTREZ), ENTREZID = ENTREZ)
      mer <- merge(DTNames, DTENTREZ, by = "ID")
      symbols <- as.character(mer$ENTREZ)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"),"ENTREZID" ) )
      goHuman <- unique(goHuman)
      GPL2895 <- merge(mer, goHuman, by = "ENTREZID")
      GPL2895 <- GPL2895[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL2895$ID <- as.character(GPL2895$ID)
      GPL2895$ID <- gsub("GE", "", GPL2895$ID)
      GPL2895$GPLID <- "GPL2895"
      AnnotationDT <- rbind(AnnotationDT, GPL2895)
    } else if(PlatInfo$GPLNumber[i] == "GPL10558"){
      z <- getGEO("GPL10558")
      temp <- as.data.table(z@dataTable@table)
      setnames(temp, c("ID", "Entrez_Gene_ID"), c("ID", "ENTREZID"))
      symbols <- as.character(temp$ENTREZID)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL", "GENENAME"),"ENTREZID" ) )
      goHuman <- unique(goHuman)
      goHuman$ENTREZID <- as.character(goHuman$ENTREZID)
      temp$ENTREZID <- as.character(temp$ENTREZID)
      GPL10558 <- merge(temp, goHuman, by = "ENTREZID")
      GPL10558 <- GPL10558[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL10558 <- GPL10558[!is.na(ENTREZID),]
      GPL10558$GPLID <- "GPL10558"
      AnnotationDT <- rbind(AnnotationDT, GPL10558)
    } else if(PlatInfo$GPLNumber[i] == "GPL14550"){
      z <- getGEO("GPL14550")
      z <- z@dataTable@table[, c("ID", "GENE", "GENE_SYMBOL", "GENE_NAME")]
      z <- z[!z$GENE_NAME == "",]
      setnames(z, c("ID", "GENE", "GENE_SYMBOL", "GENE_NAME"), c("ID", "ENTREZID", "SYMBOL", "GENENAME"))
      GPL14550 <- as.data.table(z)
      GPL14550 <- GPL14550[,c("ID", "ENTREZID", "SYMBOL", "GENENAME"), with = FALSE]
      GPL14550$GPLID <- "GPL14550"
      AnnotationDT <- rbind(AnnotationDT, GPL14550)
    } else if(PlatInfo$GPLNumber[i] == "GPL17586"){
      # BiocManager::install("hta20transcriptcluster.db")
      library(hta20transcriptcluster.db)
      x<-hta20transcriptclusterGENENAME
      mapped_probes<-mappedkeys(x)
      GENENAME <- unlist(as.list(x[mapped_probes]))
      DTNames <- data.table(ID=names(GENENAME), GENENAME)
      x<-hta20transcriptclusterENTREZID
      mapped_probes<-mappedkeys(x)
      ENTREZ<-  unlist(as.list(x[mapped_probes]))
      DTENTREZ <- data.table(ID=names(ENTREZ), ENTREZID = ENTREZ)
      mer <- merge(DTNames, DTENTREZ, by = "ID")
      symbols <- as.character(mer$ENTREZ)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"),"ENTREZID" ) )
      goHuman <- unique(goHuman)
      GPL17586 <- merge(mer, goHuman, by = "ENTREZID")
      GPL17586 <- GPL17586[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL17586$GPLID <- "GPL17586"
      AnnotationDT <- rbind(AnnotationDT, GPL17586)
    } else if(PlatInfo$GPLNumber[i] == "GPL8910"){
      z <- getGEO("GPL8910")
      temp <- as.data.table(z@dataTable@table)
      Values <- unique(temp$`Reference Accession`)
      mappings <- AnnotationDbi::select(org.Hs.eg.db, keys=Values, columns=c("ENTREZID","SYMBOL","GENENAME"),keytype="REFSEQ")
      t2 <- temp[,c("ID", "Reference Accession"), with = FALSE] %>% setnames("Reference Accession", "REFSEQ")
      GPL8910 <- merge(t2, mappings, by = "REFSEQ")
      GPL8910 <- GPL8910[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL8910$GPLID <- "GPL8910"
      AnnotationDT <- rbind(AnnotationDT, GPL8910)
    } else if(PlatInfo$GPLNumber[i] == "GPL127"){
      z <- getGEO("GPL127")
      temp <- as.data.table(z@dataTable@table)
      setnames(temp, c("ID", "GENE"), c("ID", "ENTREZID"))
      symbols <- as.character(temp$ENTREZID)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"),"ENTREZID" ) )
      goHuman <- unique(goHuman)
      goHuman$ENTREZID <- as.character(goHuman$ENTREZID)
      temp$ENTREZID <- as.character(temp$ENTREZID)
      GPL127 <- merge(temp, goHuman, by = "ENTREZID", all.x = TRUE)
      GPL127 <- GPL127[,c("ENTREZID", "ID", "GENE_NAME", "SYMBOL"), with = FALSE]
      setnames(GPL127, "GENE_NAME", "GENENAME")
      GPL127 <- GPL127[!ENTREZID == "",]
      GPL127 <- GPL127[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL127$GPLID <- "GPL127"
      AnnotationDT <- rbind(AnnotationDT, GPL127)
    } else if(PlatInfo$GPLNumber[i] == "GPL128"){
      z <- getGEO("GPL128")
      temp <- as.data.table(z@dataTable@table)
      setnames(temp, c("ID", "GENE"), c("ID", "ENTREZID"))
      symbols <- as.character(temp$ENTREZID)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"),"ENTREZID" ) )
      goHuman <- unique(goHuman)
      goHuman$ENTREZID <- as.character(goHuman$ENTREZID)
      temp$ENTREZID <- as.character(temp$ENTREZID)
      GPL128 <- merge(temp, goHuman, by = "ENTREZID", all.x = TRUE)
      GPL128 <- GPL128[,c("ENTREZID", "ID", "GENE_NAME", "SYMBOL"), with = FALSE]
      setnames(GPL128, "GENE_NAME", "GENENAME")
      GPL128 <- GPL128[!ENTREZID == "",]
      GPL128 <- GPL128[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL128$GPLID <- "GPL128"
      AnnotationDT <- rbind(AnnotationDT, GPL128)
    } else if(PlatInfo$GPLNumber[i] == "GPL131"){
      z <- getGEO("GPL131")
      temp <- as.data.table(z@dataTable@table)
      setnames(temp, c("ID", "GENE"), c("ID", "ENTREZID"))
      symbols <- as.character(temp$ENTREZID)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"),"ENTREZID" ) )
      goHuman <- unique(goHuman)
      goHuman$ENTREZID <- as.character(goHuman$ENTREZID)
      temp$ENTREZID <- as.character(temp$ENTREZID)
      GPL131 <- merge(temp, goHuman, by = "ENTREZID", all.x = TRUE)
      GPL131 <- GPL131[,c("ENTREZID", "ID", "GENE_NAME", "SYMBOL"), with = FALSE]
      setnames(GPL131, "GENE_NAME", "GENENAME")
      GPL131 <- GPL131[!ENTREZID == "",]
      GPL131 <- GPL131[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL131$GPLID <- "GPL131"
      AnnotationDT <- rbind(AnnotationDT, GPL131)
    } else if(PlatInfo$GPLNumber[i] == "GPL549"){
      z <- getGEO("GPL549")
      temp <- as.data.table(z@dataTable@table)
      setnames(temp, c("ID", "GENE"), c("ID", "ENTREZID"))
      symbols <- as.character(temp$ENTREZID)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"),"ENTREZID" ) )
      goHuman <- unique(goHuman)
      goHuman$ENTREZID <- as.character(goHuman$ENTREZID)
      temp$ENTREZID <- as.character(temp$ENTREZID)
      GPL549 <- merge(temp, goHuman, by = "ENTREZID", all.x = TRUE)
      GPL549 <- GPL549[,c("ENTREZID", "ID", "Gene Name", "SYMBOL"), with = FALSE]
      setnames(GPL549, "Gene Name", "GENENAME")
      GPL549 <- GPL549[!ENTREZID == "",]
      GPL549 <- GPL549[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL549$GPLID <- "GPL549"
      AnnotationDT <- rbind(AnnotationDT, GPL549)
    } else if(PlatInfo$GPLNumber[i] == "GPL96"){
      z <- getGEO("GPL96")
      temp <- as.data.table(z@dataTable@table)
      GPL96 <- temp[,c("ID", "Gene Symbol", "Gene Title", "ENTREZ_GENE_ID"), with = FALSE]
      setnames(GPL96, c("Gene Symbol", "Gene Title", "ENTREZ_GENE_ID"), c("SYMBOL", "GENENAME", "ENTREZID"))
      GPL96$ENTREZID <- gsub(" /// .+", "", GPL96$ENTREZID)
      GPL96$SYMBOL <- gsub(" /// .+", "", GPL96$SYMBOL)
      GPL96 <- GPL96[!ENTREZID == "",]
      GPL96 <- GPL96[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL96$GPLID <- "GPL96"
      AnnotationDT <- rbind(AnnotationDT, GPL96)
    } else if(PlatInfo$GPLNumber[i] == "GPL3050"){
      z <- getGEO("GPL3050")
      temp <- as.data.table(z@dataTable@table)
      temp <- temp[,c("ID", "ENSEMBL_ID"), with = FALSE]
      setnames(temp, c("ID", "ENSEMBL_ID"), c("ID", "ENSEMBL"))
      symbols <- as.character(temp$ENSEMBL)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL", "ENTREZID"), "ENSEMBL") )
      goHuman <- unique(goHuman)
      GPL3050 <- merge(temp, goHuman, by = "ENSEMBL")
      GPL3050 <- unique(GPL3050[,!c("ENSEMBL"), with = FALSE])
      GPL3050 <- GPL3050[!is.na(ENTREZID),]
      Values <- unique(GPL3050$ENTREZID)
      mappings <- AnnotationDbi::select(org.Hs.eg.db, keys=Values, columns=c("GENENAME"),keytype="ENTREZID")
      GPL3050 <- merge(GPL3050, mappings, by = "ENTREZID")
      GPL3050 <- GPL3050[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL3050$GPLID <- "GPL3050"
      AnnotationDT <- rbind(AnnotationDT, GPL3050)
    } else if(PlatInfo$GPLNumber[i] == "GPL20115"){
      z <- getGEO("GPL20115")
      temp <- as.data.table(z@dataTable@table)
      temp <- temp[, c("ID", "EntrezGeneID", "GeneName"), with = FALSE]
      setnames(temp, c("EntrezGeneID", "GeneName"), c("ENTREZID", "GENENAME"))
      temp$ENTREZID <- as.character(temp$ENTREZID)
      temp <- temp[!is.na(ENTREZID),]
      symbols <- as.character(temp$ENTREZID)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"), "ENTREZID") )
      goHuman <- unique(goHuman)
      GPL20115 <- unique(merge(temp, goHuman, by = "ENTREZID"))
      GPL20115 <- GPL20115[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL20115$GPLID <- "GPL20115"
      AnnotationDT <- rbind(AnnotationDT, GPL20115)
    } else if(PlatInfo$GPLNumber[i] == "GPL18056"){
      z <- getGEO("GPL18056")
      temp <- as.data.table(z@dataTable@table)
      temp <- temp[,c("ID", "Gene Symbol", "Gene Name"), with = FALSE]
      setnames(temp, c("ID", "Gene Symbol", "Gene Name"), c("ID", "SYMBOL", "GENENAME"))
      symbols <- as.character(temp$SYMBOL)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("ENTREZID"), "SYMBOL") )
      goHuman <- unique(goHuman)
      GPL18056 <- unique(merge(temp, goHuman, by = "SYMBOL"))
      GPL18056 <- GPL18056[!is.na(ENTREZID),]
      GPL18056 <- GPL18056[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL18056$GPLID <- "GPL18056"
      AnnotationDT <- rbind(AnnotationDT, GPL18056)
    } else if(PlatInfo$GPLNumber[i] == "GPL19886"){
      z <- getGEO("GPL19886")
      temp <- as.data.table(z@dataTable@table)
      temp <- temp[, c("ID", "Gene symbol", "EntrezGeneID"), with = FALSE]
      setnames(temp, c("ID", "Gene symbol", "EntrezGeneID"), c("ID", "SYMBOL", "ENTREZID"))
      temp$ENTREZID <- as.character(temp$ENTREZID)
      Values <- temp$ENTREZID
      mappings <- AnnotationDbi::select(org.Hs.eg.db, keys=Values, columns=c("GENENAME"),keytype="ENTREZID")
      GPL19886 <- merge(temp, mappings, by = "ENTREZID")
      GPL19886 <- GPL19886[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL19886$GPLID <- "GPL19886"
      AnnotationDT <- rbind(AnnotationDT, GPL19886)
    } else if (PlatInfo$GPLNumber[i] == "GPL23159"){
      z <- getGEO("GPL23159")
      z <- as.data.table(z@dataTable@table)
      REFSEQ <- gsub(" ", "", gsub("//.+", "" , z[[10]]))
      z$REFSEQ <- REFSEQ
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, REFSEQ, c("REFSEQ", "ENTREZID", "SYMBOL", "GENENAME"),"REFSEQ"))
      goHuman <- unique(goHuman[!is.na(ENTREZID),])
      z <- z[, c("ID", "REFSEQ"), with = FALSE]
      z <- z[!grepl("--CONTROL|--normgene", REFSEQ, ignore.case = TRUE),]
      GPL23159 <- merge(z, goHuman, by = "REFSEQ")
      GPL23159 <- GPL23159[,c("ID", "ENTREZID", "SYMBOL", "GENENAME"), with = FALSE]
      GPL23159$GPLID <- "GPL23159"
      AnnotationDT <- rbind(AnnotationDT, GPL23159)
    } else if (PlatInfo$GPLNumber[i] == "GPL19109"){
      z <- getGEO("GPL19109")
      z <- as.data.table(z@dataTable@table)
      setnames(z, "ENTREZ_GENE_ID", "ENTREZID")
      z <- z[!is.na(ENTREZID),]
      z$ENTREZID <- as.character(z$ENTREZID)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, z$ENTREZID, c("SYMBOL", "GENENAME"),"ENTREZID"))
      goHuman <- unique(goHuman[!is.na(ENTREZID),])
      GPL19109 <- unique(merge(z, goHuman, by = "ENTREZID"))
      GPL19109 <- GPL19109[,c("ID", "ENTREZID", "SYMBOL", "GENENAME"), with = FALSE]
      GPL19109$GPLID <- "GPL19109"
      AnnotationDT <- rbind(AnnotationDT, GPL19109)
    } else if (PlatInfo$GPLNumber[i] == "GPL17692"){
      z <- getGEO("GPL17692")
      z <- as.data.table(z@dataTable@table)
      z <- z[,c("ID", "GB_ACC"), with = FALSE][grepl("NR_|NM_", GB_ACC),]
      setnames(z, "GB_ACC", "REFSEQ")
      z <- z[!grepl("--CONTROL|--normgene", REFSEQ, ignore.case = TRUE),]
      REFSEQ <- z$REFSEQ
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, REFSEQ, c("REFSEQ", "ENTREZID", "SYMBOL", "GENENAME"),"REFSEQ"))
      goHuman <- unique(goHuman[!is.na(ENTREZID),])
      GPL17692 <- merge(z, goHuman, by = "REFSEQ")
      GPL17692 <- GPL17692[,c("ID", "ENTREZID", "SYMBOL", "GENENAME"), with = FALSE]
      GPL17692$GPLID <- "GPL17692"
      AnnotationDT <- rbind(AnnotationDT, GPL17692)
    } else if (PlatInfo$GPLNumber[i] == "GPL13667"){
      z <- getGEO("GPL13667")
      z <- as.data.table(z@dataTable@table)
      z <- z[,c("ID", "Entrez Gene"), with = FALSE]
      z <- z[!grepl("---", `Entrez Gene`, ignore.case = TRUE),]
      z <- z[!duplicated(`Entrez Gene`),]
      z$`Entrez Gene` <- as.character(z$`Entrez Gene`)
      setnames(z, "Entrez Gene", "ENTREZID")
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, z$ENTREZID, c("ENTREZID", "SYMBOL", "GENENAME"),"ENTREZID"))
      goHuman <- goHuman[!is.na(SYMBOL),]
      GPL13667 <- merge(z, goHuman, by = "ENTREZID")
      GPL13667 <- GPL13667[,c("ID", "ENTREZID", "SYMBOL", "GENENAME"), with = FALSE]
      GPL13667$GPLID <- "GPL13667"
      AnnotationDT <- rbind(AnnotationDT, GPL13667)
    } else if (PlatInfo$GPLNumber[i] == "GPL10335"){
      z <- getGEO("GPL10335")
      temp <- as.data.table(z@dataTable@table)
      setnames(temp, c("ID", "GeneID"), c("ID", "ENTREZID"))
      temp <- temp[, c("ID", "ENTREZID"), with = FALSE]
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, temp$ENTREZID, c("ENTREZID", "SYMBOL", "GENENAME"),"ENTREZID"))
      goHuman <- goHuman[!is.na(SYMBOL),]
      GPL10335 <- merge(temp, goHuman, by = "ENTREZID")
      GPL10335 <- GPL10335[,c("ID", "ENTREZID", "SYMBOL", "GENENAME"), with = FALSE]
      GPL10335$GPLID <- "GPL10335"
      AnnotationDT <- rbind(AnnotationDT, GPL10335)
      } else {
      print(paste("There is no annotation information available for:", PlatInfo$GPLNumber[i])) } }
  return(AnnotationDT)
}

#' GEO data download
#'
#' @param DS A vector containing the GEO GSE numbers of the desired datasets.
#' @param gpl A vector containing the GEO GPL numbers of the desired datasets.
#' @param gsm A vector containing the Comparison vector used in GEO2R describing how to parse the data and calculate fold changes.
#' @param namestr A vector containing the column names of the fold change calculations that are exported.
#' @param nameraw An optional vector containing the column names of the raw data columns that are exported.
#' @param PlatAnnotInfo The name of the object housing the plate annotation information exported from the `PlatformAnnotationLoad` function.
#' @param destdir The path to the directory to cache GEO data downloads.
#' @param filename description
#' @param writeDB Either TRUE or FALSE indicating if the fold change data should be saved to drive.
#' @param writeRaw A logical vector indicating if the raw data should be saved to drive.
#' @param GenerateMetaData A logical vector indicating if the meta data should be saved to drive.
#' @param MetaDataPath The path to the directory where the meta data should be saved.
#' @param writeMetaData Either TRUE or FALSE indicating if the meta data should be saved to drive.
#' @param DBPath The path to the directory where the raw and fold change data should be written.
#' @param Technology A character vector indicating either "Array" or "RNAseq" indicating if the data was generated using a minroArray or RNA sequencing technology.
#' @param renameRaw Either TRUE or FALSE indicating if the columns of the raw data should be renamed.
#' @param subsetRaw Either TRUE or FALSE indicating if the raw data should be subsetted to only the columns used in fold change calculations.
#' @param sleep The number of seconds each iteration should wait before starting the next.
#' @import Biobase
#' @import limma
#' @import DESeq2
#' @import GEOquery
#' @export
GEOCompile <- function(DS, gpl, gsm, namestr, nameraw, PlatAnnotInfo, destdir, filename=NULL,
                       writeDB=TRUE, writeRaw=TRUE, GenerateMetaData, MetaDataPath, writeMetaData=TRUE,
                       DBPath, Technology, renameRaw=FALSE, subsetRaw=FALSE, sleep = 30){
  DataListTemp <- list()
  RawListTemp <- list()
  MetaDataListTemp <- list()
  DataList <- list()
  RawList <- list()
  MetaDataList <- list()
  for(a in 1:length(DS)){

    if(Technology[a] == "Array"){
      nam <- namestr[a]
      nam <- gsub(" ", "", nam)
      nam <- strsplit(nam, ",")[[1]]
      gset <- suppressMessages(getGEO(DS[a], filename = filename[a], destdir = destdir))
      if(length(gset) > 1) idx <- grep(gpl[a], attr(gset, "names")) else idx <- 1
      gset <- gset[[idx]]
      fvarLabels(gset) <- make.names(fvarLabels(gset))
      comp <- gsub(" ", "", gsm[a])
      comp <- gsub(",", "", comp)
      gsms <- paste0(comp)
      #### Set up raw names ####
      if(renameRaw){
        namraw1 <- nameraw[a]
        namraw1 <- gsub(" ", "", namraw1)
        namraw1 <- make.unique(strsplit(namraw1, ",|;")[[1]])
        namraw1 <- gsub("\n", "", namraw1)
        namraw1 <- namraw1[!strsplit(gsms, split="")[[1]] == "X"] }
      sml <- c()
      for(i in 1:nchar(gsms)){ sml[i] <- substr(gsms,i,i)}
      ex <- exprs(gset)
      qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
      LogC <- (qx[5] > 100) ||
        (qx[6]-qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] < 2)
      if(LogC){ ex[which(ex <= 0)] <- NaN
      exprs(gset) <- log2(ex) }
      sml <- paste("G", sml, sep="")
      f1 <- as.factor(sml)
      gset$description2 <- f1
      design <- model.matrix(~description2 + 0, gset)
      colnames(design) <- levels(f1)
      fit <- lmFit(gset, design)
      cont.matrix <- makeContrasts(G1-G0, levels = design)
      fit2 <- contrasts.fit(fit, cont.matrix)
      fit2 <- eBayes(fit2, 0.01)
      tT <- topTable(fit2, adjust="fdr", sort.by = "B", number = 25000000000)
      #### subset ####
      ex2 <- data.table(subset(tT, select=c("ID", "logFC", "P.Value", "adj.P.Val")))
      ex2$ID <- as.character(ex2$ID)
      #### annotate with gene names ####
      plat <- PlatAnnotInfo[GPLID == gpl[a],][,!"GPLID", with = FALSE]
      if(nrow(plat) == 0){ print(paste("There is no annotation information available for", gpl[a])) }
      plat$ID <- as.character(plat$ID)
      plat <- plat[!duplicated(plat$ID),]
      ex2 <- merge(plat, ex2, by = "ID")
      ex2 <- ex2[,.SD[which.min(P.Value)], by = "SYMBOL"] # duplicated gene names are dealt with by taking the most significant record
      setnames(ex2, c("logFC", "P.Value", "adj.P.Val"), nam)
      ex2$ID <- as.character(ex2$ID)
      exraw <- data.table(ex)
      if(subsetRaw){ exraw <- exraw[,grepl("[0-9]", strsplit(as.character(gsms), split = "")[[1]]), with = FALSE] }
      if(renameRaw){ setnames(exraw, colnames(exraw), namraw1)
      } else { setnames(exraw, colnames(exraw), paste(DS[a], colnames(exraw), sep = "_")) }
      exraw$ID <- as.character(rownames(ex))
      #### annotate raw data with gene names ####
      exraw <- merge(plat, exraw, by = "ID")
      if(writeDB){
        fpath <- file.path(DBPath, paste(paste(DS[a], gsub("-", "_", gsub("_.+", "", nam[1])), sep = "_"), ".txt", sep = ""))
        saveRDS(ex2, file = gsub(".txt", ".rds", fpath) ) }
      if(writeRaw[a]){
        if(subsetRaw){
       fpath <- file.path(DBPath, paste(paste(DS[a], gsub("-", "_", gsub("_.+", "", nam[1])), "Raw", sep = "_"), ".txt", sep = ""))
        } else {
          fpath <- file.path(DBPath, paste(paste(DS[a], "Raw", sep = "_"), ".txt", sep = ""))
        }
       saveRDS(exraw, file = gsub(".txt", ".rds", fpath) ) }
      if(GenerateMetaData[a]){
        Pdat <- pData(gset)
        #### Add Meta data ####
        colMet <- as.data.table(Pdat)
        if(writeMetaData){
          fpath <- file.path(MetaDataPath, paste(paste(DS[a], gsub("-", "_", gsub("_.+", "", nam[1])), "MetaData", sep = "_"), ".txt", sep = ""))
          saveRDS(colMet, file = gsub(".txt", ".rds", fpath) ) }
        MetaDataListTemp[[1]] <- colMet
        MetaDataList <- c(MetaDataList, MetaDataListTemp)  }
      DataListTemp[[1]] <- ex2
      RawListTemp[[1]] <- exraw
      DataList <- c(DataList, DataListTemp)
      RawList <- c(RawList, RawListTemp) }

    if(Technology[a] == "RNAseq"){
      ACC <- paste("acc=", DS[a], sep = "")
      file <- paste("file=", DS[a], "_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep = "")
      comp <- gsub(" ", "", gsm[a])
      comp <- gsub(",", "", comp)
      gsms <- paste0(comp)
      #### Set up DEG names ####
      nam <- namestr[a]
      nam <- gsub(" ", "", nam)
      nam <- strsplit(nam, ",")[[1]]
      #### Set up raw names ####
      if(renameRaw){
        namraw1 <- nameraw[a]
        namraw1 <- gsub(" ", "", namraw1)
        namraw1 <- make.unique(strsplit(namraw1, ",|;")[[1]])
        namraw1 <- gsub("\n", "", namraw1) }
      # load counts table from GEO
      urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
      path <- paste(urld, ACC, file, sep="&");
      tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")
      apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
      annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
      rownames(annot) <- annot$GeneID
      sml <- strsplit(gsms, split="")[[1]]
      sel <- which(sml != "X")
      sml <- sml[sel]
      if(!subsetRaw){ exraw <- tbl }
      tbl <- tbl[ ,sel]
      gs <- factor(sml)
      groups <- make.names(c("Ctrl", "Tx"))
      levels(gs) <- groups
      sample_info <- data.frame(Group = gs, row.names = colnames(tbl))
      keep <- rowSums( tbl >= 10 ) >= min(table(gs))
      tbl <- tbl[keep, ]
      ds <- DESeqDataSetFromMatrix(countData=tbl, colData=sample_info, design= ~Group)
      ds <- DESeq(ds, test="Wald", sfType="poscount")
      r <- results (ds, contrast=c("Group", groups[2], groups[1]), alpha=0.05, pAdjustMethod ="fdr")
      tT <- r[order(r$padj)[1:length(r$padj)],]
      tT <- merge(as.data.frame(tT), annot, by=0, sort=F)
      tT <- subset(tT, select=c("GeneID","padj","pvalue","lfcSE","stat","log2FoldChange","baseMean","Symbol","Description"))
      #### subset ####
      ex2 <- data.table(subset(tT, select=c("GeneID", "Symbol", "Description", "log2FoldChange", "pvalue", "padj")))
      #### Adjust columnnames and order ####
      setnames(ex2, c("log2FoldChange", "pvalue", "padj"), nam)
      setnames(ex2, c("GeneID", "Symbol", "Description"), c("ENTREZID", "SYMBOL", "GENENAME"))
      ex2$ID <- ex2$ENTREZID
      ex2 <- ex2[,c("SYMBOL", "ID", "ENTREZID", "GENENAME", nam), with = FALSE]
      #### Get Raw data ####
      if(subsetRaw){ exraw <- counts(ds, normalized = F)  }
      GeneID <- as.integer(rownames(exraw))
      exraw <- as.data.table(exraw)
      #### Update column names ####
      if(renameRaw){
        namraw1 <- namraw1[gsub(".+_", "", namraw1) %in% colnames(exraw)]
        if(identical(colnames(exraw),gsub(".+_", "", namraw1))){
          setnames(exraw, colnames(exraw), namraw1)
        } else { print(paste("new colnames are not equivalent for", DS[a], nam[1])) }
      } else { setnames(exraw, colnames(exraw), paste(DS[a], colnames(exraw), sep = "_")) }
      exraw$ENTREZID <- GeneID
      #### merge FC and raw data together ####
      mer <- merge(ex2, exraw, by = "ENTREZID")
      if(writeDB){
        fpath <- file.path(DBPath, paste(paste(DS[a], gsub("-", "_", gsub("_.+", "", nam[1])), sep = "_"), ".txt", sep = ""))
        saveRDS(ex2, file = gsub(".txt", ".rds", fpath) ) }
      if(writeRaw[a]){
        merRaw <- mer[,!grepl("logFC|Pvalue|PValue", colnames(mer)), with = FALSE]
        if(subsetRaw){
        fpath <- file.path(DBPath, paste(paste(DS[a], gsub("-", "_", gsub("_.+", "", nam[1])), "Raw", sep = "_"), ".txt", sep = ""))
        } else {
          fpath <- file.path(DBPath, paste(paste(DS[a], "Raw", sep = "_"), ".txt", sep = ""))
        }
        saveRDS(merRaw, file = gsub(".txt", ".rds", fpath) ) }
      DataListTemp[[1]] <- mer
      RawListTemp[[1]] <- merRaw
      DataList <- c(DataList, DataListTemp)
      RawList <- c(RawList, RawListTemp) }
    print(paste("Completed", a, "of", length(DS)))
    Sys.sleep(sleep) }
  names(DataList) <- DS
  names(RawList) <- DS
 # if(length(MetaDataList) > 0){ names(MetaDataList) <- DS}
  return(list(FCData = DataList, RawData = RawList)) }

#' double check that the fold change calculations were performed correctly
#'
#' @param DBPath The path to the directory where the raw and fold change data are stored.
#' @param DS A vector containing the GEO GSE numbers of the desired datasets.
#' @param gsm A vector containing the Comparison vector used in GEO2R describing how to parse the data and calculate fold changes.
#' @param namestr A vector containing the column names of the fold change calculations that are exported.
#' @param Technology A vector containing the technologies used for each dataset.
#' @param GraphPath The directory to save the plots that are generated.
#' @import stringr
#' @import ggplot2
#' @export
GEO2RDirectionCheck <- function(DBPath, DS, namestr, gsm, Technology, GraphPath, subsetRaw = FALSE){
  plotList <- list()
  plotnames <- NULL
  checkFile <- data.table()
  for(a in 1:length(DS)){
    #### load data ####
    nam <- namestr[a]
    nam <- gsub(" ", "", nam)
    nam <- strsplit(nam, ",")[[1]]
    if(subsetRaw){
      rawfpath <- file.path(DBPath, paste(paste(DS[a], gsub("-", "_", gsub("_.+", "", nam[1])), "Raw", sep = "_"), ".rds", sep = ""))
    } else {
      rawfpath <- file.path(DBPath, paste(paste(DS[a], "Raw", sep = "_"), ".rds", sep = ""))
    }
    fcfpath <- file.path(DBPath, paste(paste(DS[a], gsub("-", "_", gsub("_.+", "", nam[1])), sep = "_"), ".rds", sep = ""))
    raw <- readRDS(rawfpath)
    FC <- readRDS(fcfpath)
    #### contain col subsetting vector ####
    comp <- gsub(" ", "", gsm[a])
    comp <- gsub(",", "", comp)
    gsms <- paste0(comp)
    spl <- str_split(gsub("X", "", gsms), "")[[1]]
    #### take fold change ####
    if(sum(raw[complete.cases(raw),5:ncol(raw)] < 0) > 0){
      raw <- cbind(raw[,c("ID", "SYMBOL", "ENTREZID", "GENENAME")], as.data.frame(as.matrix(raw[,!c("ID", "SYMBOL", "ENTREZID", "GENENAME")]) + 1)) }
    raw2 <- raw[,!c("ID", "SYMBOL", "ENTREZID", "GENENAME"), with = FALSE]
    comp <- gsub("_.+", "", colnames(FC)[ncol(FC)])
    spl <- gsub("[0-9]+$", "", str_split(comp, "-")[[1]])
    if(Technology[a] == "RNAseq"){
      #### perform normalization for RNAseq samples ####
      sml <- strsplit(gsms, split="")[[1]]
      sel <- which(sml != "X")
      sml <- sml[sel]
      tbl <- as.matrix(raw2[ , ..sel])
      gs <- factor(sml)
      groups <- make.names(c("Ctrl", "Tx"))
      levels(gs) <- groups
      sample_info <- data.frame(Group = gs, row.names = colnames(tbl))
      keep <- rowSums( tbl >= 10 ) >= min(table(gs))
      tbl <- tbl[keep, ]
      ds <- DESeqDataSetFromMatrix(countData=tbl, colData=sample_info, design= ~Group)
      ds <- DESeq(ds, test="Wald", sfType="poscount")
      #### Obtain normalized count values ####
      est <- estimateSizeFactors(ds)
      raw2 <- as.data.table(counts(est, normalized = TRUE))
      #### perform differential expression from normalized counts ####
      onetx <- rowMeans(raw2[,sml == 1, with = FALSE])
      zeroCtrl <- rowMeans(raw2[,sml == 0, with = FALSE])
      ratio <- onetx/zeroCtrl }
    if(Technology[a] == "Array"){
      #### perform normalization for Array samples ####
      sml <- c()
      for(i in 1:nchar(gsms)){ sml[i] <- substr(gsms,i,i)}
      ex <- raw2#exprs(gset)
      qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
      LogC <- (qx[5] > 100) ||
        (qx[6]-qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] < 2)
      if(LogC){ ex[which(ex <= 0)] <- NaN
      raw2 <- log2(ex) }
      onetx <- rowMeans(raw2[,str_split(gsms, "")[[1]] == 1, with = FALSE])
      zeroCtrl <- rowMeans(raw2[,str_split(gsms, "")[[1]] == 0, with = FALSE])
      ratio <- onetx/zeroCtrl }
    #### log2 ####
    log <- NULL
    if(sum(ratio[!is.na(ratio)] < 0) == 0){
      for(i in 1:length(ratio)){
        if(!is.na(ratio[i]) & !ratio[i] == "Inf"){
          if(ratio[i] > 1){
            log[i] <- log2(ratio[i])
          }  else if(ratio[i] > 0){
            log[i] <- -log2(1/ratio[i])
          } else if(ratio[i] == 0){
            log[i] <- 0 }
        } else{
          log[i] <- NA
        } } } else {
          log <- ratio }
    manFC <- data.table(SYMBOL = raw$SYMBOL, log = log)
    mer <- merge(FC, manFC, by = "SYMBOL")
    mer <- mer[!is.na(mer[[colnames(mer)[grepl("_logFC", colnames(mer))]]]),]
    mer <- mer[!is.na(mer$log),]
    FCcol <- colnames(mer)[grepl("_logFC", colnames(mer))]
    p <- print(ggplot(mer, aes(x=mer[[colnames(mer)[grepl("_logFC", colnames(mer))]]], y=log)) + geom_point()+ xlab(FCcol))
    #### write graph to file ####
    graphpath <- file.path(GraphPath, paste(a, "_", paste(DS[a], gsub("-", "_", gsub("_.+", "", nam[1])), sep = "_"), ".png", sep = "")) # ".xls"
    png(filename = graphpath, width = 2000, height = 1500, units = "px", res = 300); print(p); dev.off()
    plotList <- c(plotList, list(p))
    p <- NULL
    plotnames[a] <- gsub(".rds", "", gsub(".+/", "", fcfpath))
    #### fit a linear model ####
    linear <- lm(as.formula(paste( paste("`", FCcol, "`", sep = ""), "~", "log")), mer)
    LMIntercept <- linear$coefficients[2]
    #### quantify differences ####
    corelation <- stats::cor(mer[[FCcol]], mer$log)
    temp <- data.table(dataset = paste(paste(DS[a], gsub("-", "_", gsub("_.+", "", nam[1])), sep = "_"), sep = ""),
                       correct = sum(nrow(mer[mer[[FCcol]] > 0 & log > 0,]), nrow(mer[mer[[FCcol]] < 0 & log < 0,])),
                       incorrect = sum(nrow(mer[mer[[FCcol]] > 0 & log < 0,]), nrow(mer[mer[[FCcol]] > 0 & log < 0,])),
                       correlation= corelation,
                       LMEstimate = LMIntercept,
                       iteration = a)
    checkFile <- rbind(checkFile, temp )
    print(paste(a, "of", length(DS), "Completed")) }
  names(plotList) <- plotnames
  return(list(checkFile=checkFile, plotList = plotList)) }

#' Create a summarized experiment object that integrates all DEG data.
#'
#' @param DEGDatapath The path to the directory where the fold change data are stored.
#' @param SEPath The path to the directory where the compiled SummarizedExperiment object should be stored.
#' @import data.table
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @import SummarizedExperiment
#' @import AnnotationDbi
#' @export
DESEGenerate <- function(DEGDatapath, SEPath){
  #### Check if SE object exists ####
  tryCatch({ SE1 <- readRDS(SEPath)
  }, warning = function(w) {Scratch <<- TRUE; print("SE file not found. Compiling from scratch")})
  if(exists("SE1")){SE <- SE1; Scratch <- FALSE; print("Updating SE object") }
  #### obtain fold change file names ####
  files <- list.files(DEGDatapath)
  files <- files[grepl(".rds", files)]
  files <- files[!grepl("_Raw|_Spotfire|_RPKMRaw|_CountRaw", files)]
  if(Scratch){
    #### map Human annotation information to gene names ####
    x<-org.Hs.egSYMBOL; symbols<-mappedkeys(x)
    goHuman <- as.data.table( suppressMessages(AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"), "ENTREZID" ) ) %>% setnames(c("ENTREZID"), paste(c("ENTREZID"), "_Human", sep = ""))) # , "ENSEMBL"
    goHuman <- goHuman[!is.na(goHuman$SYMBOL),]
    #### map Mouse annotation information to gene names ####
    x<-org.Mm.egSYMBOL; symbols<-mappedkeys(x)
    goMouse <- as.data.table( suppressMessages(AnnotationDbi::select(org.Mm.eg.db, symbols, c("SYMBOL") , "ENTREZID") ) %>% setnames(c("ENTREZID"), paste(c("ENTREZID"), "_Mouse", sep = ""))) #, "ENSEMBL"
    goMouse <- goMouse[!is.na(goMouse$SYMBOL),]
    # #### map Rat annotation information to gene names ####
    # x<-org.Rn.egSYMBOL; symbols<-mappedkeys(x)
    # goRat <- as.data.table( select(org.Rn.eg.db, symbols, c("SYMBOL"), "ENTREZID") ) %>% setnames(c("ENTREZID"), paste(c("ENTREZID"), "_Rat", sep = "")))
    # goRat <- goRat[!is.na(goRat$ENTREZID),]
    #### set up rowData information ####
    ####################################
    print("formatting rowData information.")
    goMouse$SYMBOL <- toupper(goMouse$SYMBOL)
    goHuman$SYMBOL <- toupper(goHuman$SYMBOL)
    mer <- merge(goHuman, goMouse, by = "SYMBOL", all = TRUE)
    mer <- mer[!grepl("RIK$", mer$SYMBOL),]
    mer <- mer[!grepl("RIK[0-9]$", mer$SYMBOL),]
    mer <- mer[!grepl("---", mer$SYMBOL),]
    mer <- mer[!grepl("1-DEC", mer$SYMBOL),]
    mer <- mer[!grepl("1-MAR", mer$SYMBOL),]
    #### remove duplicated records ####
    dups <- unique(mer[duplicated(mer$SYMBOL),]$SYMBOL)
    dupRMDT <- data.table()
    for(i in 1:length(dups)){
      temp <- mer[mer$SYMBOL == dups[i],]
      if(nrow(temp[complete.cases(temp)]) == 0){
        df <- as.data.frame(is.na(as.matrix(temp[,2:3, with = FALSE])))
        FAL <- apply(df, 1, sum)
        if(sum(FAL == 2) == 4){ temp <- temp[1,]
        } else { temp <- unique(temp[FAL ==1,]) }
        dupRMDT <- rbind(dupRMDT, temp)
      } else { dupRMDT <- rbind(dupRMDT, temp[complete.cases(temp),]) } }
    dupRMDT
    #### Combine duplicated records with non-duplicated records ####
    dupRMDT <- dupRMDT[!duplicated(dupRMDT$SYMBOL),]
    mer2 <- mer[!(mer$SYMBOL %in% dupRMDT$SYMBOL),]
    mer2 <- rbind(mer2, dupRMDT)
    #### format SummarizedExperiment rowData data frame ####
    rowData <- as.data.frame(mer2)
    row.names(rowData) <- rowData$SYMBOL
    rowData <- rowData[order(rownames(rowData), decreasing = FALSE),]
    print("Compiling list of FC data.")
    checkIdent <- NULL
    SEList <- list()
    Names <- NULL
    for(i in 1:length(files)){
      #### load and format file ####
      temp <- as.data.table(readRDS(file.path(DEGDatapath, files[i])))
      temp$SYMBOL <- toupper(temp$SYMBOL)
      temp <- temp[toupper(temp$SYMBOL) %in% mer2$SYMBOL,]
      colnames(temp) <- gsub(".+_", "", colnames(temp))
      temp2 <- temp[,c("SYMBOL", "logFC", "Pvalue", "AdjPValue"), with = FALSE]
      #### if there are gene duplications, keep the more significant of the two ####
      temp2 <- temp2[,.SD[which.min(Pvalue)], by = "SYMBOL"]
      #### normalize row lengths across all names ####
      tempMer <- merge(mer2, temp2, by = "SYMBOL", all.x = TRUE)
      tempMer <- tempMer[,c("SYMBOL", "logFC", "Pvalue", "AdjPValue"), with = FALSE]
      #### format data frame ####
      tempMer <- as.data.frame(tempMer)
      rownames(tempMer) <- tempMer$SYMBOL
      tempMer$SYMBOL <- NULL
      tempMer <- tempMer[order(rownames(tempMer), decreasing = FALSE),]
      #### Save to list ####
      SEList[[i]] <- tempMer
      names(SEList)[i] <- gsub(".rds", "", files[i])
      print(paste("completed", i, "of", length(files))) }
    #### Create SummarizedExperiment object ####
    SE <- SummarizedExperiment(assays = SimpleList(SEList), rowData = rowData)
  } else {
    #### Check for duplicated names ####
    dups <- files[gsub(".rds", "", files) %in% names(assays(SE))]
    if(length(dups) > 0){
      print(paste("Ignoring", dups, "which alerady exist in the database."))
      files <- files[!(files %in% dups)] }
    if(length(files) > 0){
      rowData = as.data.table(rowData(SE))
      print("Compiling list of FC data.")
      SEList <- list()
      Names <- NULL
      for(i in 1:length(files)){
        #### load and format file ####
        temp <- as.data.table(readRDS(file.path(DEGDatapath, files[i])))
        temp$SYMBOL <- toupper(temp$SYMBOL)
        temp <- temp[toupper(temp$SYMBOL) %in% rowData$SYMBOL,]
        colnames(temp) <- gsub(".+_", "", colnames(temp))
        temp2 <- temp[,c("SYMBOL", "logFC", "Pvalue", "AdjPValue"), with = FALSE]
        #### if there are gene duplications, keep the more significant of the two ####
        temp2 <- temp2[,.SD[which.min(Pvalue)], by = "SYMBOL"]
        #### normalize row lengths across all names ####
        tempMer <- merge(rowData, temp2, by = "SYMBOL", all.x = TRUE)
        tempMer <- tempMer[,c("SYMBOL", "logFC", "Pvalue", "AdjPValue"), with = FALSE]
        #### format data frame ####
        tempMer <- as.data.frame(tempMer)
        rownames(tempMer) <- tempMer$SYMBOL
        tempMer$SYMBOL <- NULL
        tempMer <- tempMer[order(rownames(tempMer), decreasing = FALSE),]
        #### Save to list ####
        SEList[[i]] <- tempMer
        names(SEList)[i] <- gsub(".rds", "", files[i])
        print(paste("completed", i, "of", length(files))) }
      #### Create SummarizedExperiment object ####
      SEList1 <- list()
      SEass <- assays(SE)
      for(b in 1:length(SEass)){ SEList1[[b]] <- SEass[[b]];
      names(SEList1)[b] <- names(SEass)[b] }
      SE <- SummarizedExperiment(assays = SimpleList(c(SEList1, SEList)), rowData = rowData) } }
  saveRDS(SE, file=SEPath)
  return(SE) }

#' Combine all raw data into one data table.
#'
#' @param Fpath The file path to the directory where the raw data is stored.
#' @param outPath The path to the assembled raw data file.
#' @param StartAt A numeric value indicating the position of the raw datasets to start compiling raw data from.
#' @param sleep The number of seconds each iteration should wait before starting the next.
#' @param DataType Either "Array", "RNAseqRPKM", or "RNAseqCount" indicating the type of raw data to assemble into a data frame.  Only "Array" is used for data generated using GEO2R.
#' @param SampleSelect A logical vector indicating which raw datasets to incorporate into the app.
#' @import plyr
#' @import data.table
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @import AnnotationDbi
#' @export
RawDataCompile <- function(Fpath = file.path(homedir, "5_AppData"), outPath = file.path("./6_ProcessFiles/RawData.txt"), StartAt = 1, sleep = 20, DataType = "Array", CheckIdenticalData = FALSE){ #SampleSelect = overview$DownloadRawData,
  if(DataType == "Array"){
    tryCatch({ RawArrayComplete <- as.data.frame(fread(outPath)) #RawArrayComplete <- read_feather(outPath)# as.data.table()
    }, warning = function(w) {print("Raw Database file not found. Compiling from scratch"); Scratch <<- TRUE
    }, error = function(e) {    print("Raw Database file not found. Compiling from scratch"); Scratch <<- TRUE  })
    if(exists("RawArrayComplete")){ print("Updating Raw Database. Data will be ignored for duplicated column names with identical data."); Scratch <- FALSE }
    if(Scratch){
      #### map Human annotation information to gene names ####
      x<-org.Hs.egSYMBOL; symbols<-mappedkeys(x)
      goHuman <- as.data.table( suppressMessages(AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"), "ENTREZID" ) ) %>% setnames(c("ENTREZID"), paste(c("ENTREZID"), "_Human", sep = "")))
      goHuman <- goHuman[!is.na(goHuman$SYMBOL),]
      #### map Mouse annotation information to gene names ####
      x<-org.Mm.egSYMBOL; symbols<-mappedkeys(x)
      goMouse <- as.data.table( suppressMessages(AnnotationDbi::select(org.Mm.eg.db, symbols, c("SYMBOL") , "ENTREZID") ) %>% setnames(c("ENTREZID"), paste(c("ENTREZID"), "_Mouse", sep = "")))
      goMouse <- goMouse[!is.na(goMouse$SYMBOL),]
      # #### map Rat annotation information to gene names ####
      # x<-org.Rn.egSYMBOL; symbols<-mappedkeys(x)
      # goRat <- as.data.table( select(org.Rn.eg.db, symbols, c("SYMBOL"), "ENTREZID") ) %>% setnames(c("ENTREZID"), paste(c("ENTREZID"), "_Rat", sep = "")))
      # goRat <- goRat[!is.na(goRat$ENTREZID),]
      #### set up rowData information ####
      ####################################
      print("formatting rowData information.")
      goMouse$SYMBOL <- toupper(goMouse$SYMBOL)
      goHuman$SYMBOL <- toupper(goHuman$SYMBOL)
      mer <- merge(goHuman, goMouse, by = "SYMBOL", all = TRUE)
      mer <- mer[!grepl("RIK$", mer$SYMBOL),]
      mer <- mer[!grepl("RIK[0-9]$", mer$SYMBOL),]
      mer <- mer[!grepl("---", mer$SYMBOL),]
      mer <- mer[!grepl("1-DEC", mer$SYMBOL),]
      mer <- mer[!grepl("1-MAR", mer$SYMBOL),]
      #### remove duplicated records ####
      dups <- unique(mer[duplicated(mer$SYMBOL),]$SYMBOL)
      dupRMDT <- data.table()
      for(i in 1:length(dups)){
        temp <- mer[mer$SYMBOL == dups[i],]
        if(nrow(temp[complete.cases(temp)]) == 0){
          df <- as.data.frame(is.na(as.matrix(temp[,2:3, with = FALSE])))
          FAL <- apply(df, 1, sum)
          if(sum(FAL == 2) == 4){ temp <- temp[1,]
          } else { temp <- unique(temp[FAL ==1,]) }
          dupRMDT <- rbind(dupRMDT, temp)
        } else { dupRMDT <- rbind(dupRMDT, temp[complete.cases(temp),]) } }
      #### Combine duplicated records with non-duplicated records ####
      dupRMDT <- dupRMDT[!duplicated(dupRMDT$SYMBOL),]
      mer2 <- mer[!(mer$SYMBOL %in% dupRMDT$SYMBOL),]
      mer2 <- rbind(mer2, dupRMDT)
      #### format SummarizedExperiment rowData data frame ####
      rowData <- as.data.frame(mer2)
      row.names(rowData) <- rowData$SYMBOL
      rowData <- rowData[order(rownames(rowData), decreasing = FALSE),] }
    #### Loop through files and incorporate them into the raw data master file ####
    files <- list.files(Fpath)
    files <- files[grepl(".rds", files) & grepl("_Raw", files)]
    # files <- files[SampleSelect]
    files <- files[StartAt:length(files)]
    for(i in 1:length(files)){
      print("Loading experiment raw table")
      one <- as.data.table(readRDS(file.path(Fpath, files[i])))
      one <- one[,!c("ID", "GENENAME"), with = FALSE]
      one <- one[!is.na(SYMBOL),][!is.na(ENTREZID),]
      one$SYMBOL <- toupper(one$SYMBOL)
      ### remove duplicated genes ####
      print("Averaging duplicated genes")
      annot <- one[,(colnames(one) %in% c("ENTREZID","SYMBOL")), with = FALSE]
      annot <- annot[!duplicated(annot$ENTREZID),]
      one <- one[,!(colnames(one) %in% c("SYMBOL")), with = FALSE]
      one <- one[, lapply(.SD, mean), by = ENTREZID]
      one <- merge(annot, test, by = "ENTREZID")
      if(!exists("RawArrayComplete")){
        if(sum(rowData$ENTREZID_Human %in% one$ENTREZID) > 150){
          RawArrayComplete <- merge(rowData, one[,!c("SYMBOL"), with = FALSE], by.x = "ENTREZID_Human", by.y = "ENTREZID", all.x = TRUE) }
        if(sum(rowData$ENTREZID_Mouse %in% one$ENTREZID) > 150){
          RawArrayComplete <- merge(rowData, one[,!c("SYMBOL"), with = FALSE], by.x = "ENTREZID_Mouse", by.y = "ENTREZID", all.x = TRUE) }
      } else {
        ##### systematically check to see if data in duplicated columns are identical ####
        print("Checking for and resolving duplicated columns")
        if(CheckIdenticalData){
        for(c in 1:20){
          if(ncol(one) > 1){
            intColNames <- intersect(colnames(one), colnames(RawArrayComplete)); intColNames <- intColNames[!(intColNames %in% c("ENTREZID", "SYMBOL"))]
            if(length(intColNames) > 0){
              for(b in 1:length(intColNames)){
                t1 <- one[,c(intColNames[b], "SYMBOL"), with = FALSE]
                t2 <- RawArrayComplete[,c(intColNames[b], "SYMBOL")]
                t2 <- t2[!is.na(t2[[intColNames[b]]]),]
                mer <- merge(t1, t2, by = "SYMBOL")
                if(cor(mer[[2]], mer[[3]], use="complete.obs") < 0.95 ){
                  print(paste("Data for column name:", intColNames[b], "is not identical for", files[i], i, sep = " "))
                  setnames(one, intColNames[b], paste(intColNames[b], ".99", sep = "")) } }
              one <- one[,!(colnames(one) %in% intColNames), with = FALSE] } }
          }
        } else {
          intColNames <- intersect(colnames(one), colnames(RawArrayComplete)); intColNames <- intColNames[!(intColNames %in% c("ENTREZID", "SYMBOL"))]
          one <- one[,!(colnames(one) %in% intColNames), with = FALSE]
        }
        if(!ncol(one) > 2){ one <- data.table() }
        if(ncol(one) > 0){
          if(sum(RawArrayComplete$ENTREZID_Human %in% one$ENTREZID) > 150){
            RawArrayComplete <- merge(RawArrayComplete, one[,!c("SYMBOL"), with = FALSE], by.x = "ENTREZID_Human", by.y = "ENTREZID", all.x = TRUE) }
          if(sum(RawArrayComplete$ENTREZID_Mouse %in% one$ENTREZID) > 150){
            RawArrayComplete <- merge(RawArrayComplete, one[,!c("SYMBOL"), with = FALSE], by.x = "ENTREZID_Mouse", by.y = "ENTREZID", all.x = TRUE)
          } } }
      print(paste("incorporated", i, "of", length(files), sep = " "))
      fwrite(RawArrayComplete, gsub("feather", "txt", outPath), row.names = FALSE, quote = FALSE, sep = "\t")
      print(paste("System sleeping for", sleep, "seconds"))
      Sys.sleep(sleep) }
    print("Raw data updated") }
  if(DataType == "RNAseqRPKM"){
    tryCatch({ RawRPKMComplete <- as.data.frame(fread(outPath))#<- read_feather(outPath)
    }, warning = function(w) {print("Raw Database file not found. Compiling from scratch"); ScratchRPKM <<- TRUE
    }, error = function(e) {    print("Raw Database file not found. Compiling from scratch"); ScratchRPKM <<- TRUE  })
    if(exists("RawRPKMComplete")){ print("Updating Raw Database. Data will be ignored for duplicated column names with identical data."); ScratchRPKM <- FALSE }
    if(ScratchRPKM){
      #### map Human annotation information to gene names ####
      x<-org.Hs.egSYMBOL; symbols<-mappedkeys(x)
      goHuman <- as.data.table( suppressMessages(AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"), "ENTREZID" ) ) %>% setnames(c("ENTREZID"), paste(c("ENTREZID"), "_Human", sep = ""))) # , "ENSEMBL"
      goHuman <- goHuman[!is.na(goHuman$SYMBOL),]
      #### map Mouse annotation information to gene names ####
      x<-org.Mm.egSYMBOL; symbols<-mappedkeys(x)
      goMouse <- as.data.table( suppressMessages(AnnotationDbi::select(org.Mm.eg.db, symbols, c("SYMBOL") , "ENTREZID") ) %>% setnames(c("ENTREZID"), paste(c("ENTREZID"), "_Mouse", sep = ""))) #, "ENSEMBL"
      goMouse <- goMouse[!is.na(goMouse$SYMBOL),]
      # #### map Rat annotation information to gene names ####
      # x<-org.Rn.egSYMBOL; symbols<-mappedkeys(x)
      # goRat <- as.data.table( select(org.Rn.eg.db, symbols, c("SYMBOL"), "ENTREZID") ) %>% setnames(c("ENTREZID"), paste(c("ENTREZID"), "_Rat", sep = "")))
      # goRat <- goRat[!is.na(goRat$ENTREZID),]
      #### set up rowData information ####
      ####################################
      print("formatting rowData information.")
      goMouse$SYMBOL <- toupper(goMouse$SYMBOL)
      goHuman$SYMBOL <- toupper(goHuman$SYMBOL)
      mer <- merge(goHuman, goMouse, by = "SYMBOL", all = TRUE)
      mer <- mer[!grepl("RIK$", mer$SYMBOL),]
      mer <- mer[!grepl("RIK[0-9]$", mer$SYMBOL),]
      mer <- mer[!grepl("---", mer$SYMBOL),]
      mer <- mer[!grepl("1-DEC", mer$SYMBOL),]
      mer <- mer[!grepl("1-MAR", mer$SYMBOL),]
      #### remove duplicated records ####
      dups <- unique(mer[duplicated(mer$SYMBOL),]$SYMBOL)
      dupRMDT <- data.table()
      for(i in 1:length(dups)){
        temp <- mer[mer$SYMBOL == dups[i],]
        if(nrow(temp[complete.cases(temp)]) == 0){
          df <- as.data.frame(is.na(as.matrix(temp[,2:3, with = FALSE])))
          FAL <- apply(df, 1, sum)
          if(sum(FAL == 2) == 4){ temp <- temp[1,]
          } else { temp <- unique(temp[FAL ==1,]) }
          dupRMDT <- rbind(dupRMDT, temp)
        } else { dupRMDT <- rbind(dupRMDT, temp[complete.cases(temp),]) } }
      #### Combine duplicated records with non-duplicated records ####
      dupRMDT <- dupRMDT[!duplicated(dupRMDT$SYMBOL),]
      mer2 <- mer[!(mer$SYMBOL %in% dupRMDT$SYMBOL),]
      mer2 <- rbind(mer2, dupRMDT)
      #### format SummarizedExperiment rowData data frame ####
      rowData <- as.data.frame(mer2)
      row.names(rowData) <- rowData$SYMBOL
      rowData <- rowData[order(rownames(rowData), decreasing = FALSE),] }
    #### Loop through files and incorporate them into the raw data master file ####
    files <- list.files(Fpath)
    files <- files[grepl(".rds", files) & grepl("_RPKMRaw", files)]
    files <- files[StartAt:length(files)]
    for(i in 1:length(files)){
      print("Loading experiment raw table")
      one <- as.data.table(readRDS(file.path(Fpath, files[i])))
      one <- one[,!c("ID", "GENENAME"), with = FALSE]
      one <- one[!is.na(SYMBOL),][!is.na(ENTREZID),]
      one$SYMBOL <- toupper(one$SYMBOL)
      setnames(one, colnames(one)[!colnames(one) %in% c("ENTREZID", "SYMBOL")], paste(gsub("_.+", "", files[i]), colnames(one)[!colnames(one) %in% c("ENTREZID", "SYMBOL")], sep = "_"))
      #### remove duplicated genes ####
      print("Averaging duplicated genes")
      annot <- one[,(colnames(one) %in% c("ENTREZID","SYMBOL")), with = FALSE]
      annot <- annot[!duplicated(annot$ENTREZID),]
      one <- one[,!(colnames(one) %in% c("SYMBOL")), with = FALSE]
      one <- one[, lapply(.SD, mean), by = ENTREZID]
      one <- merge(annot, test, by = "ENTREZID")
      if(!exists("RawRPKMComplete")){
        if(sum(rowData$ENTREZID_Human %in% one$ENTREZID) > 150){
          RawRPKMComplete <- merge(rowData, one[,!c("SYMBOL"), with = FALSE], by.x = "ENTREZID_Human", by.y = "ENTREZID", all.x = TRUE) }
        if(sum(rowData$ENTREZID_Mouse %in% one$ENTREZID) > 150){
          RawRPKMComplete <- merge(rowData, one[,!c("SYMBOL"), with = FALSE], by.x = "ENTREZID_Mouse", by.y = "ENTREZID", all.x = TRUE) }
      } else {
        ##### systematically check to see if data in duplicated columns are identical ####
        print("Checking for and resolving duplicated columns")
        if(CheckIdenticalData){
        for(c in 1:20){
          if(ncol(one) > 1){
            intColNames <- intersect(colnames(one), colnames(RawRPKMComplete)); intColNames <- intColNames[!(intColNames %in% c("ENTREZID", "SYMBOL"))]
            if(length(intColNames) > 0){
              for(b in 1:length(intColNames)){
                t1 <- one[,c(intColNames[b], "SYMBOL"), with = FALSE]
                t2 <- RawRPKMComplete[,c(intColNames[b], "SYMBOL")]
                t2 <- t2[!is.na(t2[[intColNames[b]]]),]
                mer <- merge(t1, t2, by = "SYMBOL")
                if(cor(mer[[2]], mer[[3]], use="complete.obs") < 0.95 ){
                  print(paste("Data for column name:", intColNames[b], "is not identical for", files[i], i, sep = " "))
                  setnames(one, intColNames[b], paste(intColNames[b], ".99", sep = ""))  } }
              one <- one[,!(colnames(one) %in% intColNames), with = FALSE] } } }
        } else {
          intColNames <- intersect(colnames(one), colnames(RawArrayComplete)); intColNames <- intColNames[!(intColNames %in% c("ENTREZID", "SYMBOL"))]
          one <- one[,!(colnames(one) %in% intColNames), with = FALSE]
        }
        if(!ncol(one) > 2){ one <- data.table() }
        if(ncol(one) > 0){
          if(sum(RawRPKMComplete$ENTREZID_Human %in% one$ENTREZID) > 150){
            RawRPKMComplete <- merge(RawRPKMComplete, one[,!c("SYMBOL"), with = FALSE], by.x = "ENTREZID_Human", by.y = "ENTREZID", all.x = TRUE) }
          if(sum(RawRPKMComplete$ENTREZID_Mouse %in% one$ENTREZID) > 150){
            RawRPKMComplete <- merge(RawRPKMComplete, one[,!c("SYMBOL"), with = FALSE], by.x = "ENTREZID_Mouse", by.y = "ENTREZID", all.x = TRUE) } } }
      print(paste("incorporated", i, "of", length(files), sep = " "))
      fwrite(RawRPKMComplete, gsub("feather", "txt", outPath), row.names = FALSE, quote = FALSE, sep = "\t")
      print(paste("System sleeping for", sleep, "seconds"))
      Sys.sleep(sleep) }
    print("RPKM data updated")
  }
  if(DataType == "RNAseqCount"){
    tryCatch({ RawCountComplete <- as.data.frame(fread(outPath))#<- read_feather(outPath)
    }, warning = function(w) {print("Raw Database file not found. Compiling from scratch"); ScratchCount <<- TRUE
    }, error = function(e) {    print("Raw Database file not found. Compiling from scratch"); ScratchCount <<- TRUE  })
    if(exists("RawCountComplete")){ print("Updating Raw Database. Data will be ignored for duplicated column names with identical data."); ScratchCount <- FALSE }
    if(ScratchCount){
      #### map Human annotation information to gene names ####
      x<-org.Hs.egSYMBOL; symbols<-mappedkeys(x)
      goHuman <- as.data.table( suppressMessages(AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"), "ENTREZID" ) ) %>% setnames(c("ENTREZID"), paste(c("ENTREZID"), "_Human", sep = ""))) # , "ENSEMBL"
      goHuman <- goHuman[!is.na(goHuman$SYMBOL),]
      #### map Mouse annotation information to gene names ####
      x<-org.Mm.egSYMBOL; symbols<-mappedkeys(x)
      goMouse <- as.data.table( suppressMessages(AnnotationDbi::select(org.Mm.eg.db, symbols, c("SYMBOL") , "ENTREZID") ) %>% setnames(c("ENTREZID"), paste(c("ENTREZID"), "_Mouse", sep = ""))) #, "ENSEMBL"
      goMouse <- goMouse[!is.na(goMouse$SYMBOL),]
      # #### map Rat annotation information to gene names ####
      # x<-org.Rn.egSYMBOL; symbols<-mappedkeys(x)
      # goRat <- as.data.table( select(org.Rn.eg.db, symbols, c("SYMBOL"), "ENTREZID") ) %>% setnames(c("ENTREZID"), paste(c("ENTREZID"), "_Rat", sep = "")))
      # goRat <- goRat[!is.na(goRat$ENTREZID),]
      #### set up rowData information ####
      ####################################
      print("formatting rowData information.")
      goMouse$SYMBOL <- toupper(goMouse$SYMBOL)
      goHuman$SYMBOL <- toupper(goHuman$SYMBOL)
      mer <- merge(goHuman, goMouse, by = "SYMBOL", all = TRUE)
      mer <- mer[!grepl("RIK$", mer$SYMBOL),]
      mer <- mer[!grepl("RIK[0-9]$", mer$SYMBOL),]
      mer <- mer[!grepl("---", mer$SYMBOL),]
      mer <- mer[!grepl("1-DEC", mer$SYMBOL),]
      mer <- mer[!grepl("1-MAR", mer$SYMBOL),]
      #### remove duplicated records ####
      dups <- unique(mer[duplicated(mer$SYMBOL),]$SYMBOL)
      dupRMDT <- data.table()
      for(i in 1:length(dups)){
        temp <- mer[mer$SYMBOL == dups[i],]
        if(nrow(temp[complete.cases(temp)]) == 0){
          df <- as.data.frame(is.na(as.matrix(temp[,2:3, with = FALSE])))
          FAL <- apply(df, 1, sum)
          if(sum(FAL == 2) == 4){ temp <- temp[1,]
          } else { temp <- unique(temp[FAL ==1,]) }
          dupRMDT <- rbind(dupRMDT, temp)
        } else { dupRMDT <- rbind(dupRMDT, temp[complete.cases(temp),]) } }
      #### Combine duplicated records with non-duplicated records ####
      dupRMDT <- dupRMDT[!duplicated(dupRMDT$SYMBOL),]
      mer2 <- mer[!(mer$SYMBOL %in% dupRMDT$SYMBOL),]
      mer2 <- rbind(mer2, dupRMDT)
      #### format SummarizedExperiment rowData data frame ####
      rowData <- as.data.frame(mer2)
      row.names(rowData) <- rowData$SYMBOL
      rowData <- rowData[order(rownames(rowData), decreasing = FALSE),] }
    #### Loop through files and incorporate them into the raw data master file ####
    files <- list.files(Fpath)
    files <- files[grepl(".rds", files) & grepl("_CountRaw", files)]
    files <- files[StartAt:length(files)]
    for(i in 1:length(files)){
      print("Loading experiment raw table")
      one <- as.data.table(readRDS(file.path(Fpath, files[i])))
      one <- one[,!c("ID", "GENENAME"), with = FALSE]
      one <- one[!is.na(SYMBOL),][!is.na(ENTREZID),]
      one$SYMBOL <- toupper(one$SYMBOL)
      setnames(one, colnames(one)[!colnames(one) %in% c("ENTREZID", "SYMBOL")], paste(gsub("_.+", "", files[i]), colnames(one)[!colnames(one) %in% c("ENTREZID", "SYMBOL")], sep = "_"))
      #### remove duplicated genes ####
      print("Averaging duplicated genes")
      annot <- one[,(colnames(one) %in% c("ENTREZID","SYMBOL")), with = FALSE]
      annot <- annot[!duplicated(annot$ENTREZID),]
      one <- one[,!(colnames(one) %in% c("SYMBOL")), with = FALSE]
      one <- one[, lapply(.SD, mean), by = ENTREZID]
      one <- merge(annot, test, by = "ENTREZID")
      if(!exists("RawCountComplete")){
        if(sum(rowData$ENTREZID_Human %in% one$ENTREZID) > 150){
          RawCountComplete <- merge(rowData, one[,!c("SYMBOL"), with = FALSE], by.x = "ENTREZID_Human", by.y = "ENTREZID", all.x = TRUE) }
        if(sum(rowData$ENTREZID_Mouse %in% one$ENTREZID) > 150){
          RawCountComplete <- merge(rowData, one[,!c("SYMBOL"), with = FALSE], by.x = "ENTREZID_Mouse", by.y = "ENTREZID", all.x = TRUE) }
      } else {
        ##### systematically check to see if data in duplicated columns are identical ####
        print("Checking for and resolving duplicated columns")
        if(CheckIdenticalData){
        for(c in 1:20){
          if(ncol(one) > 1){
            intColNames <- intersect(colnames(one), colnames(RawCountComplete)); intColNames <- intColNames[!(intColNames %in% c("ENTREZID", "SYMBOL"))]
            if(length(intColNames) > 0){
              for(b in 1:length(intColNames)){
                t1 <- one[,c(intColNames[b], "SYMBOL"), with = FALSE]
                t2 <- RawCountComplete[,c(intColNames[b], "SYMBOL")]
                t2 <- t2[!is.na(t2[[intColNames[b]]]),]
                mer <- merge(t1, t2, by = "SYMBOL")
                if(cor(mer[[2]], mer[[3]], use="complete.obs") < 0.95 ){
                  print(paste("Data for column name:", intColNames[b], "is not identical for", files[i], i, sep = " "))
                  setnames(one, intColNames[b], paste(intColNames[b], ".99", sep = ""))  } }
              one <- one[,!(colnames(one) %in% intColNames), with = FALSE] } } }
        } else {
          intColNames <- intersect(colnames(one), colnames(RawArrayComplete)); intColNames <- intColNames[!(intColNames %in% c("ENTREZID", "SYMBOL"))]
          one <- one[,!(colnames(one) %in% intColNames), with = FALSE]
        }
        if(!ncol(one) > 2){ one <- data.table() }
        if(ncol(one) > 0){
          if(sum(RawCountComplete$ENTREZID_Human %in% one$ENTREZID) > 150){
            RawCountComplete <- merge(RawCountComplete, one[,!c("SYMBOL"), with = FALSE], by.x = "ENTREZID_Human", by.y = "ENTREZID", all.x = TRUE) }
          if(sum(RawCountComplete$ENTREZID_Mouse %in% one$ENTREZID) > 150){
            RawCountComplete <- merge(RawCountComplete, one[,!c("SYMBOL"), with = FALSE], by.x = "ENTREZID_Mouse", by.y = "ENTREZID", all.x = TRUE) } } }
      print(paste("incorporated", i, "of", length(files), sep = " "))
      fwrite(RawCountComplete, gsub("feather", "txt", outPath), row.names = FALSE, quote = FALSE, sep = "\t")
      print(paste("System sleeping for", sleep, "seconds"))
      Sys.sleep(sleep) }
    print("Count data updated")  } }

#' Compile all meta data into one table
#'
#' @param RNAseqFilePath The file path where the RNAseq meta data is stored.  Files should be named according to the following naming convention: "GSEXXXXX_SraRunTable.txt"
#' @param ArrayFilePath The file path where the Array meta data is stored.
#' @import data.table
#' @export
MetDataCompile <- function(RNAseqFilePath, ArrayFilePath, overview = overview){
  #### RNAseq data ####
  files <- list.files(RNAseqFilePath)
  files <- file.path(RNAseqFilePath, files)
  dt1 <- data.table()
  for(i in 1:length(files)){
    t <- fread(files[i]);
    t$Dataset <- gsub(".+/|_.+", "", files[i])
    t$RawColumnNames <- paste(t$Dataset, t$`Sample Name`, sep = "_")
    dt1 <- rbind(dt1, t, fill = TRUE)}
  dt1 <- dt1[`Assay Type` == "RNA-Seq",]
  #### Array data ####
  files <- list.files(ArrayFilePath)
  files <- file.path(ArrayFilePath, files)
  dt2 <- data.table()
  for(i in 1:length(files)){
    # t <- fread(files[i]);
    t <- readRDS(files[i]);
    t$Dataset <- gsub(".+/|_.+", "", files[i])
    t$RawColumnNames <- paste(t$Dataset, t$geo_accession, sep = "_")
    dt2 <- rbind(dt2, t, fill = TRUE)}
  dt2 <- as.data.table(dt2)
  dt2$`Assay Type` <- "Array"
  #### Combine together ####
  dt <- rbind(dt1, dt2, fill = TRUE)
  #### clean up column names ####
  cname <- colnames(dt)
  cname <- gsub(":.+", "", cname)
  cname <- gsub("erythrocyte  eicosapentaenoic acid (% of total lipids)", "erythrocyte  eicosapentaenoic acid % of total lipids", cname, fixed = TRUE)
  cname <- gsub("erythrocyte arachidonic acid (% of total lipids)", "erythrocyte arachidonic acid % of total lipids", cname, fixed = TRUE)
  cname <- gsub("erythrocyte docosahexaenoic acid (% of total lipids)", "erythrocyte docosahexaenoic acid % of total lipids", cname, fixed = TRUE)
  cname <- gsub("liver arachidonic acid (% of diacylglycerols)", "liver arachidonic acid % of diacylglycerols", cname, fixed = TRUE)
  cname <- gsub("liver arachidonic acid (% of phospholipids)", "liver arachidonic acid % of phospholipids", cname, fixed = TRUE)
  cname <- gsub("liver arachidonic acid (% of triacylglycerols)", "liver arachidonic acid % of triacylglycerols", cname, fixed = TRUE)
  cname <- gsub("liver docosahexaenoic acid (% of diacylglycerols)", "liver docosahexaenoic acid % of diacylglycerols", cname, fixed = TRUE)
  cname <- gsub("liver docosahexaenoic acid (% of phospholipids)", "liver docosahexaenoic acid % of phospholipids", cname, fixed = TRUE)
  cname <- gsub("liver docosahexaenoic acid (% of triacylglycerols)", "liver docosahexaenoic acid % of triacylglycerols", cname, fixed = TRUE)
  cname <- gsub("liver eicosapentaenoic acid (% of diacylglycerols)", "liver eicosapentaenoic acid % of diacylglycerols", cname, fixed = TRUE)
  cname <- gsub("liver eicosapentaenoic acid (% of phospholipids)", "liver eicosapentaenoic acid % of phospholipids", cname, fixed = TRUE)
  cname <- gsub("liver eicosapentaenoic acid (% of triacylglycerols)", "liver eicosapentaenoic acid % of triacylglycerols", cname, fixed = TRUE)
  cname <- gsub(" \\(.+\\)", "", cname)
  cname <- gsub("\\..+", "", cname)
  cname <- tolower(cname)
  setnames(dt, colnames(dt), cname)
  dt <- dt[,!colnames(dt) %in% c("characteristics_ch1"), with = FALSE,]
  #### merge duplicated column names into one column ####
  cdups <- unique(colnames(dt)[duplicated(colnames(dt))])
  for(a in 1:length(cdups)){
    #### remove duplicated column names ####
    MerAtrribute <- NULL
    temp <- dt[,colnames(dt) %in% cdups[a], with = FALSE]
    for(i in 1:nrow(temp)){
      rec <- NULL
      for(b in 1:ncol(temp[i,])){ rec <- c(rec, temp[i,][[b]]) }
      if(sum(is.na(rec)) == length(rec)){ MerAtrribute[i] <- NA
      } else { MerAtrribute[i] <- paste(rec[!is.na(rec)], collapse = " / ")  } }
    dt <- dt[,!colnames(dt) %in% cdups[a], with = FALSE]
    dt[[cdups[a]]] <- MerAtrribute
    setnames(dt, colnames(dt), gsub("\\..+", "", colnames(dt))) }
  dt <- dt[!duplicated(rawcolumnnames),]
  df <- as.data.frame(dt)
  rownames(df) <- df$rawcolumnnames
  #### remove columns containing only missing values ####
  df <- df[,!apply(df, 2, function(x){sum(is.na(x))}) == nrow(df)]
  #### Add Tissue and disease annotations ####
  df <- df[df$dataset %in% unique(overview$ID),]
  df$disease <- NULL; df$tissue <- NULL
  df$disease <- NA; df$tissue <- NA
  O2 <- overview[,c("ID", "Tissue", "Disease"), with = FALSE]
  #### add tissue and disease columns ####
  tempID <- unique(O2$ID)
  for(i in 1:length(tempID)){
    tempDF <- df[df$dataset == tempID[i],]
    if(nrow(tempDF) > 0){
      tempAn <- unique(O2[O2$ID == tempID[i],])
      tempDF$disease <- tempAn$Disease
      tempDF$tissue <- tempAn$Tissue
      df <- rbind(df[!df$dataset == tempID[i],], tempDF) } }
  return(df) }

#' Generate raw data summarized experiment object.
#'
#' @param df The meta data data frame created using the `MetDataCompile` function.
#' @param ArrayDT The Array raw data frame created using the `RawDataCompile` function
#' @import SummarizedExperiment
#' @import data.table
#' @export
GenerateRawSE <- function(df, ArrayDT){
  #### remove records from annotation table not in data array ####
  remove <- setdiff(rownames(df), colnames(ArrayDT))
  dfsub <- df[!(rownames(df) %in% remove),]
  #### add records to annotation table that are in data array ####
  add <- setdiff(colnames(ArrayDT), rownames(df))
  add <- add[!(add %in% c("ENTREZID_Human", "SYMBOL", "ENTREZID_Mouse"))]
  addDT <- as.data.frame(matrix(data=NA, nrow = length(add), ncol = length(colnames(df))))
  colnames(addDT) <- colnames(df)
  rownames(addDT) <- add
  dfsub <- rbind(dfsub, addDT)
  #### Set up summarized experiment rowData, columnData, and final data frame ####
  Rdata <- as.data.frame(ArrayDT[,c("ENTREZID_Human", "SYMBOL", "ENTREZID_Mouse")])
  row.names(Rdata) <- ArrayDT$SYMBOL
  FDat <- ArrayDT[,!c("ENTREZID_Human", "SYMBOL", "ENTREZID_Mouse")]
  row.names(FDat) <- ArrayDT$SYMBOL
  colDat <- dfsub[match(colnames(FDat), row.names(dfsub)),]
  #### create SummarizedExperiment object ####
  RawSE <- SummarizedExperiment(assays=list(RawData=FDat), colData=colDat, rowData = Rdata)
  return(RawSE) }

#' Set up app files
#'
#' @param homedir the path to the database location.
#' @export
Appetup <- function(homedir){
  t <- readLines("./Scripts/server.R")
  t <- gsub("<path to data base>", file.path(homedir), t)
  writeLines(t, "./Scripts/server.R")
  t <- readLines("./Scripts/UI.R")
  t <- gsub("<path to data base>", file.path(homedir), t)
  writeLines(t, "./Scripts/UI.R")
  }

#' Create an index file to reference meta data information
#'
#' @param metaData A data table containing meta data information.
#' @import data.table
#' @export
CreateMetaIndex <- function(metaData){
  #### remove columns that are not meta data columns ####
  Nmeta <- c("run","assay type","avgspotlen","bases","bioproject","biosample","bytes","center name","consent",
             "datastore filetype","datastore provider","datastore region","experiment","instrument","librarylayout","libraryselection",
             "librarysource","organism","platform","releasedate","create_date","version","source_name","sra study","dataset","rawcolumnnames",
             "geo_accession", "relation diagnosis body mass index sample name")
  NmetaDF <- metaData[,!(colnames(metaData) %in% Nmeta)]
  #### create index Table ####
  Cnames <- colnames(NmetaDF)
  metaIndexDT <- data.table()
  for(i in 1:length(Cnames)){
    temvec <- NmetaDF[,i]
    names(temvec) <- rownames(NmetaDF)
    temvec <- temvec[!(temvec == "" | is.na(temvec))]
    tempDT <- data.table(rowName = names(temvec), variable = temvec, Colname = Cnames[i])
    if(i == 1){ metaIndexDT <- tempDT } else {  metaIndexDT <- rbind(metaIndexDT, tempDT)  } }
  metaIndexDT$Dataset <- sapply(metaIndexDT$rowName, function(x){str_split(x, "_")[[1]][1]})
  metaIndexDT$Sample <- sapply(metaIndexDT$rowName, function(x){str_split(x, "_")[[1]][2]})
  return(metaIndexDT) }

#' Select genes across features and datasets
#'
#' @param SE A sumarized experiment object containing expression level data.
#' @param feature The feature interested in exploring
#' @param datasets A vector of data sets to query.
#' @param genes A vector of genes to explore.
#' @param export Either "SE", "DF", or "plot" indicating if a summarizedExperiment object, data frame, or Plot should be returned.
#' @import data.table
#' @import SummarizedExperiment
#' @import ggplot2
#' @import RColorBrewer
#' @export
SEFeatureSelect <- function(SE, feature, datasets, genes = NA, export = "SE"){
  #### select feature ####
  featSelect <- colData(SE)[[feature]]
  SEtemp <- SE[,!(featSelect == "" | is.na(featSelect))]
  #### select datasets ####
  SEtemp <- SEtemp[,colData(SEtemp)[["dataset"]] %in% datasets]
  #### update colData ####
  colData(SEtemp) <- colData(SEtemp)[,c("rawcolumnnames", "dataset", feature)]
  #### update assay data ####
  SEtemp <- SEtemp[!(apply(assay(SEtemp), 1, function(x){sum(is.na(x))}) == ncol(assay(SEtemp))),]
  #### perform gene selection ####
  if(!is.na(genes)[1]){ SEtemp <- SEtemp[rownames(rowData(SEtemp)) %in% genes,] }
  if(nrow(assay(SEtemp)) > 0){
    #### melt data for graphing ####
    if(export == "DF" | export == "plot"){
      dat <- assay(SEtemp)
      dat$SYMBOL <- rownames(dat)
      dat <- reshape2::melt(dat, id.vars = "SYMBOL")
      #### Annotate data frame with annotation information
      annot <- colData(SEtemp)
      colnames(annot) <- c("variable", "dataset", "feature")
      datAnnot <- merge(dat, annot, by = "variable") }
    if(export == "plot"){
      #### remove missing values ####
      datAnnotClean <- as.data.frame(datAnnot[!is.na(datAnnot$value),])
      box_count = length(unique(datAnnotClean$feature))
      coul = brewer.pal(9, "Set3")
      coul = colorRampPalette(coul)(box_count)
      p=ggplot(data = datAnnotClean, aes(x=feature, y=value, fill=feature))+
        geom_violin(scale = "width", trim = FALSE, alpha = 0.7)+
        geom_boxplot(outlier.shape = NA,coef = 0, fill="gray", width=0.3)+
        coord_flip()+
        ylab("expression")+
        xlab("")+
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text.y = element_text(colour = 'black', size = 10),
              axis.text.x = element_text(colour = 'black', size = 10, vjust = 1, hjust=1),
              legend.position = "none",
              strip.text.y = element_text(angle = 0, size = 10))+
        stat_boxplot(geom='errorbar', linetype=1, width=0.2, position = "dodge2")+
        scale_fill_manual(values= coul) +
        facet_grid(dataset ~ ., scales = "free", space = "free") }
    if(export == "SE"){ return(SEtemp) }
    if(export == "DF"){ return(datAnnot) }
    if(export == "plot"){ return(p) }
  } else { print("None of the selected genes are in the dataset of interest. Please select a different gene")  } }

#' Select genes across features and datasets
#'
#' @param SE A sumarized experiment object containing expression level data.
#' @param feature The feature interested in exploring
#' @param datasets A vector of data sets to query.
#' @param metaIndexDT A data table of indexed meta data information.
#' @import data.table
#' @import SummarizedExperiment
#' @import ggplot2
#' @import RColorBrewer
#' @export
AvailableGenes <- function(SE, metaIndexDT, feature, datasets){
  (select1 <- metaIndexDT[metaIndexDT$Colname == feature,]) #### return data sets with selected features #### Shiny App menu that updates based on prior selection
  # (datasets <- unique(select1$Dataset))
  #### Add gene selection ####
  featSelect <- colData(SE)[[feature]]
  SEtemp <- SE[,!(featSelect == "" | is.na(featSelect))]
  #### select datasets ####
  SEtemp <- SEtemp[,colData(SEtemp)[["dataset"]] %in% datasets]
  #### update colData ####
  colData(SEtemp) <- colData(SEtemp)[,c("rawcolumnnames", "dataset", feature)]
  #### update assay data ####
  SEtemp <- SEtemp[!(apply(assay(SEtemp), 1, function(x){sum(is.na(x))}) == ncol(assay(SEtemp))),]
  all <- cbind(rowData(SEtemp), assay(SEtemp))
  genes <- all[apply(assay(SEtemp), 1, function(x){!sum(is.na(x)) == ncol(assay(SEtemp))}),]$SYMBOL
  return(genes) }

#' Return available genes for specific data sets.
#'
#' @param SE A sumarized experiment object containing expression level data.
#' @param feature The feature interested in exploring
#' @param datasets A vector of data sets to query.
#' @param metaIndexDT A data table of indexed meta data information.
#' @import data.table
#' @import SummarizedExperiment
#' @import ggplot2
#' @import RColorBrewer
#' @export
AvailableGenes <- function(SE, metaIndexDT, feature, datasets){
  (select1 <- metaIndexDT[metaIndexDT$Colname == feature,])
  #### Add gene selection ####
  featSelect <- colData(SE)[[feature]]
  SEtemp <- SE[,!(featSelect == "" | is.na(featSelect))]
  #### select datasets ####
  SEtemp <- SEtemp[,colData(SEtemp)[["dataset"]] %in% datasets]
  #### update colData ####
  colData(SEtemp) <- colData(SEtemp)[,c("rawcolumnnames", "dataset", feature)]
  #### update assay data ####
  SEtemp <- SEtemp[!(apply(assay(SEtemp), 1, function(x){sum(is.na(x))}) == ncol(assay(SEtemp))),]
  all <- cbind(rowData(SEtemp), assay(SEtemp))
  genes <- all[apply(assay(SEtemp), 1, function(x){!sum(is.na(x)) == ncol(assay(SEtemp))}),]$SYMBOL
  return(genes) }

#' Return available genes for specific data sets.
#'
#' @param path The directory where the downloaded data is stored.
#' @param DS The ProteomeExchange "PDX" numbers of the data to download.
#' @param DownloadRawData A logical value indicating if the raw data shoud be downloaded.
#' @import rpx
#' @import RforProteomics
#' @import BiocFileCache
#' @import data.table
#' @export
ProteomicsDataDownload <- function(path, DS, DownloadRawData = FALSE){
  bfc <- BiocFileCache(path, ask = FALSE)
  for(i in 1:length(DS)){
    tryCatch( suppressWarnings(px1 <- PXDataset(DS[i])),
              error = function(e){ print(paste(DS[i], "could not be downloaded. Check for errors in the PXD number"))} )
    if(exists("px1")){
      if(sum(list.files(file.path(path)) %in% DS[i]) ==0){ dir.create(file.path(path, DS[i])) }
      files <- suppressMessages(pxfiles(px1))
      if(!DownloadRawData){ files <- files[!grepl(".raw$|raw|.baf|RAW|.xml|.sne|.mgf|.wiff|.group", files, ignore.case = TRUE)] }
      print(files)
      mztab <- pxget(px1, files, cache = bfc)
      #### move files into approperiate file ####
      print(paste("moving files to", DS[i], "directory"))
      file.copy(from = mztab, to   = file.path(path, DS[i], files))
      print("removing files from old directory")
      file.remove(mztab)
      rm(px1)
    }
    print(paste("Completed", i, "of", length(DS)))
  }
  file.remove(file.path(path, list.files(path)[grepl("BiocFileCache.sqlite", list.files(path))]))
}

#' Format proteomics data. Note: due to the lack of consistency between data sets, each proteomics data set must be formatted individually. This function provides examples of how to format the data.
#'
#' @param path The directory where the downloaded data is stored.
#' @import data.table
#' @export
FormatMaxQuant <- function(path){
  files <- gsub(".+/", "", list.dirs(path, recursive = FALSE) )
  MZList <- list()
  for(i in 1:length(files)){
    if(files[i] == "PXD005847"){
      fpath <- file.path(path, files[i])
      fpath <- file.path(fpath, list.files(fpath))
      nam <- fread(fpath[grepl("experimentalDesignTemplate.txt", fpath)])
      data <- fread(fpath[grepl("proteinGroups.txt", fpath)])
      data <- data[,c("Majority protein IDs","Peptide counts (all)", colnames(data)[grep("Intensity", colnames(data), ignore.case = TRUE)],
                      "Protein names", "Gene names"), with = FALSE]
      data$Sample <- "HDL"
      data <- cbind(data[,!grepl("intensity", colnames(data), ignore.case = TRUE), with = FALSE],
                    data[,grepl("LFQ", colnames(data), ignore.case = TRUE), with = FALSE])
      data <- data[,c("Sample", "Majority protein IDs", "Gene names", "Protein names", "Peptide counts (all)",
                      colnames(data)[grepl("LFQ", colnames(data))]), with = FALSE]
      setnames(data, colnames(data)[1:5], c("Sample", "ProteinID", "GeneSymbol", "Description", "Numberofpeptides"))
      data$Numberofpeptides <- sapply(data$Numberofpeptides, function(x){ sum(as.numeric(strsplit(x, ";")[[1]])) })
      Cnames <- gsub("_.+", "", colnames(data) )
      Cnames <- makeunique::make_unique(Cnames, "_", wrap_in_brackets = FALSE)
      setnames(data, colnames(data), Cnames)
      MZList[[i]] <- data
      names(MZList)[i] <- files[i]
    }
    if(files[i] == "PXD008934"){
      fpath <- file.path(path, files[i])
      fpath <- file.path(fpath, list.files(fpath))
      data <- fread(fpath[grepl("Prosser_HumanHearts_MassSpec.txt", fpath)])#fpath[2])
      data$Sample <- "Heart"
      data <- cbind(data[,!grepl("intensity", colnames(data), ignore.case = TRUE), with = FALSE],
                    data[,grepl("LFQ", colnames(data), ignore.case = TRUE), with = FALSE])
      data <- data[,c("Sample", "Majority protein IDs",  "Gene names", "Protein names","Peptides",
                      colnames(data)[grepl("LFQ", colnames(data))]), with = FALSE]
      setnames(data, colnames(data)[1:5], c("Sample", "ProteinID", "GeneSymbol", "Description", "Numberofpeptides"))
      meta <- fread(file.path(path, files[i], "sample_details.txt"), skip = 2)
      meta <- meta[!meta$`File #` == "",]
      meta <- meta[1:34,1:7, with = FALSE]
      Cnames <- gsub(" S", " ", colnames(data))
      Cnames <- mgsub(Cnames, meta$`File #`, meta$Etiology)
      Cnames <- makeunique::make_unique(Cnames, "_", wrap_in_brackets = FALSE)
      setnames(data, colnames(data), Cnames)
      MZList[[i]] <- data
      names(MZList)[i] <- files[i]
    }
    if(files[i] == "PXD011283"){
      fpath <- file.path(path, files[i])
      fpath <- file.path(fpath, list.files(fpath))
      # unzip(fpath[1], exdir = file.path(path, files[i]))
      fpath <- file.path(path, files[i])
      fpath <- file.path(fpath, list.files(fpath))
      fpath <- file.path(fpath[1], list.files(fpath[1]))
      # meta <- fread(fpath[grepl("experimentalDesignTemplate.txt", fpath)])
      data <- fread(fpath[grepl("peptides.txt", fpath)])
      data$Sample <- "MV_urine"
      data <- cbind(data[,!grepl("Experiment|Identification ", colnames(data), ignore.case = TRUE), with = FALSE],
                    data[,grepl("Intensity MV", colnames(data), ignore.case = TRUE), with = FALSE])
      data <- data[,c("Sample", "Proteins",  "Gene names", "Protein names","Length",
                      colnames(data)[grepl("Intensity MV", colnames(data))]), with = FALSE]
      setnames(data, colnames(data)[1:5], c("Sample", "ProteinID", "GeneSymbol", "Description", "Numberofpeptides"))
      Cnames <- colnames(data)
      Cnames <- gsub("MV_urine_", "",  Cnames)
      Cnames <- gsub("_.+", "",  Cnames)
      colnames(data) <- makeunique::make_unique(Cnames, "_", wrap_in_brackets = FALSE)
      MZList[[i]] <- data
      names(MZList)[i] <- files[i]
    }
    if(files[i] == "PXD011839"){
      fpath <- file.path(path, files[i])
      fpath <- file.path(fpath, list.files(fpath))
      # nam <- fread(fpath[grepl("sdrf.tsv", fpath)])
      data <- fread(fpath[grepl("proteinGroups.txt", fpath)])
      data <- data[,c("Majority protein IDs","Peptide counts (all)", colnames(data)[grep("Intensity", colnames(data), ignore.case = TRUE)],
                      "Protein names", "Gene names"), with = FALSE]
      data$Sample <- "Plasma"
      data <- cbind(data[,!grepl("intensity", colnames(data), ignore.case = TRUE), with = FALSE],
                    data[,grepl("LFQ", colnames(data), ignore.case = TRUE), with = FALSE])
      data <- data[,c("Sample", "Majority protein IDs", "Gene names", "Protein names", "Peptide counts (all)",
                      colnames(data)[grepl("LFQ", colnames(data))]), with = FALSE]
      setnames(data, colnames(data)[1:5], c("Sample", "ProteinID", "GeneSymbol", "Description", "Numberofpeptides"))
      data$Numberofpeptides <- sapply(data$Numberofpeptides, function(x){ sum(as.numeric(strsplit(x, ";")[[1]])) })
      data <- data[,!grepl("Liver_Rep", colnames(data)), with = FALSE]
      data <- data[,!grepl("_Top14Top6Depl_", colnames(data)), with = FALSE]
      data <- data[,!grepl("_Top6Top14Depl_", colnames(data)), with = FALSE]
      Cnames <- gsub("20170405.+_Human_", "", colnames(data))
      Cnames <- gsub("_F.+", "", Cnames)
      colnames(data) <- makeunique::make_unique(Cnames, "_", wrap_in_brackets = FALSE)
      data$ProteinID <- sapply(data$ProteinID, function(x){ strsplit(x, ";")[[1]][1] })
      data$GeneSymbol <- sapply(data$GeneSymbol, function(x){ strsplit(x, ";")[[1]][1] })
      data$Description <- sapply(data$Description, function(x){ strsplit(x, ";")[[1]][1] })
      MZList[[i]] <- data
      names(MZList)[i] <- files[i]
    }
    if(files[i] == "PXD022545"){
      fpath <- file.path(path, files[i])
      fpath <- file.path(fpath, list.files(fpath))
      # system(paste("7za x", fpath[2]))
      fpath <- file.path(path, files[i])
      fpath <- file.path(fpath, list.files(fpath))
      fpath <- file.path(fpath[4], list.files(fpath[4]))
      data <- fread(fpath[grepl("proteinGroups.txt", fpath)])
      data <- data[,c("Majority protein IDs","Peptide counts (all)", colnames(data)[grep("Intensity", colnames(data), ignore.case = TRUE)], "Fasta headers"), with = FALSE]
      data$GeneSymbol <- sapply(data$`Fasta headers`, function(x){test <- gsub("sp.+\\|| OS.+|HUMAN ", "", x); spl <- str_split(test, "_")[[1]]; spl[1] })
      data$Description <- sapply(data$`Fasta headers`, function(x){test <- gsub("sp.+\\|| OS.+|HUMAN ", "", x); spl <- str_split(test, "_")[[1]]; spl[2] })
      data$`Fasta headers` <- NULL
      data$Sample <- "retinum"
      data <- cbind(data[,!grepl("intensity", colnames(data), ignore.case = TRUE), with = FALSE],
                    data[,grepl("LFQ", colnames(data), ignore.case = TRUE), with = FALSE])
      data <- data[,c("Sample", "Majority protein IDs", "GeneSymbol", "Description", "Peptide counts (all)",
                      colnames(data)[grepl("LFQ", colnames(data))]), with = FALSE]
      setnames(data, colnames(data)[1:5], c("Sample", "ProteinID", "GeneSymbol", "Description", "Numberofpeptides"))
      data$Numberofpeptides <- sapply(data$Numberofpeptides, function(x){ sum(as.numeric(strsplit(x, ";")[[1]])) })
      data$ProteinID <- sapply(data$ProteinID, function(x){ strsplit(x, ";")[[1]][1] })
      data$GeneSymbol <- sapply(data$GeneSymbol, function(x){ strsplit(x, ";")[[1]][1] })
      data$Description <- sapply(data$Description, function(x){ strsplit(x, ";")[[1]][1] })
      colnames(data) <- gsub("_0", "_", colnames(data))
      MZList[[i]] <- data
      names(MZList)[i] <- files[i]
    }
    if(files[i] == "PXD029200"){
      fpath <- file.path(path, files[i])
      fpath <- file.path(fpath, list.files(fpath))
      meta <- as.data.table(read_excel(fpath[grepl("Pride_explanatory_sample_file.xlsx", fpath)]))
      data <- suppressWarnings(as.data.table(read_xlsx(fpath[grepl("Data_from_maxquant_used_for_filtering.xlsx", fpath)])))
      data$Sample <- "Hepatocyte"
      data <- cbind(data[,!grepl("intensity", colnames(data), ignore.case = TRUE), with = FALSE],
                    data[,grepl("LFQ", colnames(data), ignore.case = TRUE), with = FALSE])
      data <- data[,c("Sample", "Protein IDs", "Gene names", "Protein names", "Peptide counts (all)",
                      colnames(data)[grepl("LFQ", colnames(data))]), with = FALSE]
      setnames(data, colnames(data)[1:5], c("Sample", "ProteinID", "GeneSymbol", "Description", "Numberofpeptides"))
      metasub <- meta[,c("Name in data analysis", "Day", "Treatment"), with = FALSE]
      metasub$Treatment <- gsub(" ", "", gsub("Low glucose, ", "", metasub$Treatment))
      metasub$Day <- paste("Day", metasub$Day, sep = "")
      metasub$SampleTx <- paste(metasub$Treatment, metasub$Day, sep = "_")
      metasub$colname <- paste("LFQ intensity", metasub$`Name in data analysis`)
      metasub$colname[1:3] <- c("LFQ intensity 011", "LFQ intensity 012", "LFQ intensity 013")
      setnames(data, metasub$colname, make.unique(metasub$SampleTx, "_"))
      data$Numberofpeptides <- sapply(data$Numberofpeptides, function(x){ sum(as.numeric(strsplit(x, ";")[[1]])) })
      Cnames <- colnames(data)
      colnames(data) <- gsub("_Day", "Day",  Cnames)
      colnames(data) <- makeunique::make_unique(gsub("_[0-9]+", "", colnames(data)), "_", wrap_in_brackets = FALSE)
      colnames(data)[6:length(colnames(data))] <- paste("LFQ intensity", colnames(data)[6:length(colnames(data))])
      MZList[[i]] <- data
      names(MZList)[i] <- files[i]
    }
    if(files[i] == "PXD030764"){
      fpath <- file.path(path, files[i])
      fpath <- file.path(fpath, list.files(fpath))
      data <- fread(fpath[grepl("proteinGroups.txt", fpath)])
      data$Sample <- "Hepatocyte"
      data <- cbind(data[,!grepl("intensity", colnames(data), ignore.case = TRUE), with = FALSE],
                    data[,grepl("LFQ", colnames(data), ignore.case = TRUE), with = FALSE])
      data <- data[,c("Sample", "Protein IDs", "Gene names", "Protein names", "Peptide counts (all)",
                      colnames(data)[grepl("LFQ", colnames(data))]), with = FALSE]
      setnames(data, colnames(data)[1:5], c("Sample", "ProteinID", "GeneSymbol", "Description", "Numberofpeptides"))
      data$Numberofpeptides <- sapply(data$Numberofpeptides, function(x){ sum(as.numeric(strsplit(x, ";")[[1]])) })
      colnames(data)[6:length(colnames(data))] <- mgsub(colnames(data)[6:length(colnames(data))], c("M", "P", "O"), c("myristic", "palmitic", "oleic"))
      Cnames <- makeunique::make_unique(gsub("_[0-9]+", "", colnames(data)), "_", wrap_in_brackets = FALSE)
      Cnames <- gsub("LFQ intensity palmiticoleic12", "LFQ intensity palmiticoleic12_1", Cnames)
      Cnames <- gsub("LFQ intensity palmiticmyristic11", "LFQ intensity palmiticmyristic11_1", Cnames)
      setnames(data, colnames(data), Cnames)
      data <- data[!GeneSymbol == "",]
      data$ProteinID <- sapply(data$ProteinID, function(x){ strsplit(x, ";")[[1]][1] })
      data$GeneSymbol <- sapply(data$GeneSymbol, function(x){ strsplit(x, ";")[[1]][1] })
      data$Description <- sapply(data$Description, function(x){ strsplit(x, ";")[[1]][1] })
      MZList[[i]] <- data
      names(MZList)[i] <- files[i]
    }
  }
  MZList <- MZList[!names(MZList) == ""]
  return(MZList)
}

#' Save formatted proteomics data
#'
#' @param path The directory where the downloaded data is stored.
#' @param MZList Formatted proteomics data stored as a list.
#' @import data.table
#' @export
proteomicMZSave <- function(MZList, path){
  for(i in 1:length(MZList)){ fwrite(MZList[[i]], file.path(path, paste(names(MZList)[i], ".xls", sep = "")), row.names=FALSE, quote=FALSE, sep="\t") } }

#' Design Matrix from file names
#'
#' @param Fpath The directory where the data is stored.
#' @import data.table
#' @export
DesignMatrixFromNames <- function(Fpath){
  files <- list.files(Fpath)
  DesignDT <- data.table()
  for(i in 1:length(files)){
    Cnames <- colnames(fread(file.path(Fpath, files[i])))
    Cnames <- Cnames[!(Cnames %in% c("Sample", "ProteinID", "GeneSymbol", "Description", "Numberofpeptides"))]
    DT <- data.table(dataset = gsub(".xls", "", files[i]),
                     label = Cnames,
                     condition = gsub("_[0-9]+$", "", gsub(".+ ", "", Cnames)),
                     replicate = gsub(".+_", "", Cnames) )
    DesignDT <- rbind(DesignDT, DT)
  }
  replicate <- gsub(".+_", "", makeunique::make_unique(DesignDT$condition, "_", wrap_in_brackets = FALSE))
  replicate[grepl("[a-z]", replicate, ignore.case = TRUE)] <- "1"
  DesignDT$replicate <- replicate
  return(DesignDT)}

#' Create summarized experiment objects
#'
#' @param DesignDT A design table containing the dataset, column label, and the treatment condition for samples.
#' @param Fpath The directory where the data is stored.
#' @import data.table
#' @import DEP
#' @export
ProtSELoad <- function(DesignDT, Fpath){
  SEList <- list()
  files <- list.files(Fpath)
  for(i in 1:length(files)){
    data <- fread(file.path(Fpath, files[i]))
    data_unique <- DEP::make_unique(data, "GeneSymbol", "ProteinID", delim = ";")
    LFQ_columns <- grep("LFQ.", colnames(data_unique))
    if(length(LFQ_columns) == 0){
      LFQ_columns <- grep("Intensity.", colnames(data_unique), ignore.case = TRUE)
    }
    experimental_design <- DesignDT[DesignDT$dataset == gsub(".xls", "", files[i]),]
    experimental_design$dataset <- NULL
    data_se <- make_se(as.data.frame(data_unique), LFQ_columns, as.data.frame(experimental_design))
    SEList[[i]] <- data_se
    names(SEList)[i] <- gsub(".xls", "", files[i])
  }
  return(SEList) }

#' Obtain row bound data
#'
#' @param tissueSplitList A list of summarizedExperiment objects.
#' @import data.table
#' @import DEP
#' @export
RowDataCompile <- function(tissueSplitList){
  TotMel <- data.table()
  for(i in 1:length(tissueSplitList)){
    temp <- as.data.frame(assay(tissueSplitList[[i]]))
    temp$GeneSymbol <- rownames(temp)
    temp <- as.data.table(temp)
    mel <- suppressWarnings(suppressMessages(reshape2::melt(temp)))
    mel$Comb <- paste(mel$variable, gsub("_.+", "", names(tissueSplitList)[i]), sep = ":")
    mel$Cell <- gsub("_.+", "", names(tissueSplitList)[i])
    TotMel <- rbind(TotMel, mel)
  }
  TotMel$CellTissue <- paste(TotMel$Cell, gsub("_.+", "", TotMel$variable), sep = ":")
  TotMel$Tissue <- gsub("_.+", "", TotMel$variable)
  TotMel$Diet <- gsub("Plasma_|Liver_", "", TotMel$variable)
  return(TotMel) }

#' Data filtering functions
#'
#' @param dataSeList A list of summarizedExperiment objects.
#' @param thr The number of missing values to filter by.
#' @import data.table
#' @import DEP
#' @import tidyr
#' @import dplyr
#' @import tidyverse
#' @export
DataFilter <- function(dataSeList, thr = 0){
  dataFiltList <- list()
  for(i in 1:length(dataSeList)){
    se = dataSeList[[i]]
    if (is.integer(thr)){
      thr <- as.numeric(thr)}
    assertthat::assert_that(inherits(se, "SummarizedExperiment"), is.numeric(thr), length(thr) == 1)
    if (any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
      stop("'name' and/or 'ID' columns are not present in '", deparse(substitute(se)), "'\nRun make_unique() and make_se() to obtain the required columns", call. = FALSE)
    }
    if (any(!c("label", "condition", "replicate") %in% colnames(colData(se)))) {
      stop("'label', 'condition' and/or 'replicate' columns are not present in '", deparse(substitute(se)), "'\nRun make_se() or make_se_parse() to obtain the required columns", call. = FALSE)
    }
    max_repl <- max(colData(se)$replicate)
    if (thr < 0 | thr > max_repl) { stop("invalid filter threshold applied", "\nRun filter_missval() with a threshold ranging from 0 to ", max_repl)
    }
    bin_data <- assay(se)
    idx <- is.na(assay(se))
    bin_data[!idx] <- 1
    bin_data[idx] <- 0
    rownames(bin_data) <- 1:nrow(bin_data)
    keep <- bin_data %>% data.frame() %>% radiant.data::rownames_to_column() %>%
      gather(ID, value, -rowname) %>% left_join(., data.frame(colData(se)), by = "ID") %>%
      group_by(rowname, condition) %>% summarize(miss_val = n() - sum(value)) %>%
      filter(miss_val <= thr) %>% spread(condition, miss_val)
    dataFiltList[[i]] <- se[as.numeric(keep$rowname), ]
    }
  names(dataFiltList) <- names(dataSeList)
  return(dataFiltList) }

#' N peptide filtering
#'
#' @param dataFiltList A list of summarizedExperiment objects.
#' @param Npeptides The number of peptides to filter by.
#' @import data.table
#' @import DEP
#' @export
NPeptideThreshold <- function(dataFiltList, Npeptides = 2){
  PepCutOff <- dataFiltList
  PercentRemain <- NULL
  for(i in 1:length(PepCutOff)){
    PepCutOff[[i]] <- PepCutOff[[i]][rowData(PepCutOff[[i]])$Numberofpeptides >= Npeptides,]
    PercentRemain[i] <- nrow(assay(PepCutOff[[i]]))/nrow(assay(dataFiltList[[i]]))
  }
  remaining <- data.table(dataset = names(PepCutOff), percentRemaining = PercentRemain)
  return(list(data = PepCutOff, remaining = remaining)) }

#' Multiple normalization
#'
#' @param dataPepCutOff A list of summarizedExperiment objects.
#' @import data.table
#' @import DEP
#' @import NormalyzerDE
#' @export
MultiNormalization <- function(dataPepCutOff){
  #### MEAN ####
  MeanNormList <- dataPepCutOff
  for(i in 1:length(MeanNormList)){
    assayTemp <- 2^assay(MeanNormList[[i]])
    norm <- meanNormalization(assayTemp)
    rownames(norm) <- rownames(assayTemp)
    assay(MeanNormList[[i]]) <- norm }
  #### MEDIAN ####
  MedianNormList <- dataPepCutOff
  for(i in 1:length(MedianNormList)){
    assayTemp <- 2^assay(MedianNormList[[i]])
    norm <- medianNormalization(assayTemp)
    rownames(norm) <- rownames(assayTemp)
    assay(MedianNormList[[i]]) <- norm }
  #### VSN ####
  VSNNormList <- dataPepCutOff
  for(i in 1:length(VSNNormList)){
    assayTemp <- 2^assay(VSNNormList[[i]])
    norm <- performVSNNormalization(assayTemp)
    rownames(norm) <- rownames(assayTemp)
    assay(VSNNormList[[i]]) <- norm }
  ### DEP package VSN ####
  DEPVSNNormList <- dataPepCutOff
  for(i in 1:length(DEPVSNNormList)){
    DEPVSNNormList[[i]] <- normalize_vsn(DEPVSNNormList[[i]]) }
  #### Loess ####
  LoessNormList <- dataPepCutOff
  for(i in 1:length(LoessNormList)){
    assayTemp <- assay(LoessNormList[[i]])
    norm <- performCyclicLoessNormalization(assayTemp, noLogTransform = TRUE)
    rownames(norm) <- rownames(assayTemp)
    assay(LoessNormList[[i]]) <- norm }
  #### RLR ####
  RLRNormList <- dataPepCutOff
  for(i in 1:length(RLRNormList)){
    assayTemp <- assay(RLRNormList[[i]])
    norm <- performGlobalRLRNormalization(assayTemp, noLogTransform = TRUE)
    rownames(norm) <- rownames(assayTemp)
    assay(RLRNormList[[i]]) <- norm }
  #### SMAD ####
  SMADNormList <- dataPepCutOff
  for(i in 1:length(SMADNormList)){
    assayTemp <- assay(SMADNormList[[i]])
    norm <- performSMADNormalization(assayTemp, noLogTransform = TRUE)
    rownames(norm) <- rownames(assayTemp)
    assay(SMADNormList[[i]]) <- norm }
  return(list(mean=MeanNormList, median=MedianNormList, vsn=VSNNormList, DEPvsn=DEPVSNNormList, loess=LoessNormList, rlr=RLRNormList, smad=SMADNormList))
}

#' Density plots from a list
#'
#' @param MultiNormalizeList A list of summarizedExperiment objects.
#' @import data.table
#' @import DEP
#' @import ggplot2
#' @export
densityPlotFromList <- function(MultiNormalizeList){
  densityPlotList <- list()
  for(i in 1:length(MultiNormalizeList)){
    temp <- RowDataCompile(tissueSplitList=MultiNormalizeList[[i]])
    densityPlotList[[i]] <- ggplot(temp, aes(x=value, color=variable))+
      geom_density()+
      theme(legend.position="none",
            panel.spacing = unit(0.1, "lines"),
            strip.text.x = element_text(size = 8) )+
      facet_wrap(~CellTissue)
  }
  names(densityPlotList) <- names(MultiNormalizeList)
  return(densityPlotList) }

#' Count missing values
#'
#' @param data A list of summarizedExperiment objects.
#' @import DEP
#' @export
DetermineMising <- function(data){
  Missing <- NULL
  for(i in 1:length(data)){
    if(any(is.na(assay(data[[i]])))){ Missing[i] <- i } }
  if(is.null(Missing)){ print("There are no datasets with missing values") }
  Missing <- Missing[!is.na(Missing)] }

#' Impute data
#'
#' @param dataFiltList A list of summarizedExperiment objects.
#' @param type Either "MinProb", "man", or "knn" indicating the type of imputation to perform.
#' @import DEP
#' @export
DataImpute <- function(dataFiltList, type = "MinProb"){
  ImputeList <- list()
  for(i in 1:length(dataFiltList)){
    if(type == "MinProb"){
      # Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
      ImputeList[[i]] <- impute2(dataFiltList[[i]], fun = "MinProb")}
    if(type == "man"){
      # Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
      ImputeList[[i]] <- impute2(dataFiltList[[i]], fun = "man", shift = 1.8, scale = 0.3) }
    if(type == "knn"){
      # Impute missing data using the k-nearest neighbour approach (for MAR)
      ImputeList[[i]] <- impute2(dataFiltList[[i]], fun = "knn", rowmax = 0.9) }
  }
  names(ImputeList) <- names(dataFiltList)
  return(ImputeList) }

#' Impute data
#'
#' @param se A list of summarizedExperiment objects.
#' @param fun Either "MinProb", "man", or "knn" indicating the type of imputation to perform.
#' @import DEP
#' @export
impute2 <- function (se, fun = c("bpca", "knn", "QRILC", "MLE", "MinDet","MinProb", "man", "min", "zero", "mixed", "nbavg"), ...)  {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"), is.character(fun))
  fun <- match.arg(fun)
  if (any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '", deparse(substitute(se)), "'\nRun make_unique() and make_se() to obtain the required columns", call. = FALSE)
  }
  if (!any(is.na(assay(se)))) {
    warning("No missing values in '", deparse(substitute(se)), "'. ", "Returning the unchanged object.", call. = FALSE)
    return(se)
  }
  rowData(se)$imputed <- apply(is.na(assay(se)), 1, any)
  rowData(se)$num_NAs <- rowSums(is.na(assay(se)))
  if (fun == "man") { se <- manual_impute(se)
  } else {
    MSnSet_data <- as(se, "MSnSet")
    MSnSet_imputed <- MSnbase::impute(MSnSet_data, method = fun, ...)
    assay(se, withDimnames=FALSE) <- MSnbase::exprs(MSnSet_imputed)
  }
  return(se)}

#' Remove samples with a low count
#'
#' @param VarRMList A list of summarizedExperiment objects.
#' @param cut The minimum cutoff for samples.
#' @import DEP
#' @export
LowSampleCountRmove <- function(VarRMList, cut){
  Keep <- NULL
  for(i in 1:length(VarRMList)){
    if(nrow(assay(VarRMList[[i]])) > cut){ Keep[i] <- TRUE
    } else { Keep[i] <- FALSE } }
  print(paste(names(VarRMList)[!Keep], "samples removed"))
  return(VarRMList[Keep])}

#' Calculate differential expression
#'
#' @param DataList A list of summarizedExperiment objects.
#' @param type Either "control", "all" or "manual" indicating the differential expression analysis algorithm to use.
#' @import DEP
#' @export
DEAnalysis <- function(DataList = NormImpAll, type = "manual", ComparisonList=ComparisonList){
  data_diff <- list()
  for(i in 1:length(DataList)){
    if(type == "control"){
      data_diff[[i]] <- test_diff(DataList[[i]], type = "control", control = "Plasma_ND") }
    if(type == "all"){
      data_diff[[i]] <- test_diff(DataList[[i]], type = "all") }
    if(type == "manual"){
      comps <- ComparisonList[names(ComparisonList) ==  names(DataList)[i]]#[[1]] #gsub("_normalized", "", gsub("_Raw", "", gsub("scaled_", "", names(DataList)[i] ) )) ][[1]]
      for(b in 1:length(comps)){
        data_diff[[length(data_diff)+1]] <- suppressWarnings(test_diff(DataList[[i]], type = "manual", test = comps[b]))
        names(data_diff)[length(data_diff)] <- paste(names(DataList)[i], comps[b], sep = "-")
      } } }
  #### perform manual P-value correction ####
  for(i in 1:length(data_diff)){
    dat <- rowData(data_diff[[i]])
    if(sum(grepl("p.val", colnames(dat))) == 1){
      datsub <-   as.data.frame(dat[,grepl("p.val", colnames(dat))])
      colnames(datsub) <- colnames(dat)[grepl("p.val", colnames(dat))]
    } else { datsub <- dat[,grepl("p.val", colnames(dat))] }
    for(a in 1:ncol(datsub)){
      adjusted <- p.adjust(datsub[[a]], method = "BH", n = length(datsub[[a]]))
      if(a == 1){
        adjustPvals <- as.data.frame(adjusted)
      } else { adjustPvals <- cbind(adjustPvals, as.data.frame(adjusted)) }  }
    colnames(adjustPvals) <- gsub("p.val", "BHCorrection", colnames(datsub))
    rowData(data_diff[[i]]) <- cbind(dat, adjustPvals)   }
  return(data_diff) }

#' Save DEP expression summarizedExperiment objects
#'
#' @param SEList A list of summarizedExperiment objects.
#' @param Path The path to where the data should be saved.
#' @import DEP
#' @export
SaveToProteomicDB <- function(SEList, Path){
  for(i in 1:length(SEList)){ saveRDS(SEList[[i]],  file.path(Path, paste(names(SEList[i]), ".rds", sep = ""))) } }

#' Load proteomics summarized experiment object from DB
#'
#' @param Fname The file name of the dataset to load.
#' @param Path The path to where the data should be saved.
#' @import DEP
#' @export
ProteomicSELoad <- function(Path=file.path(homedir, "Proteomic_3"), Fname=files[1]){
  SelectedList <- list()
  SelectedList[[1]] <- readRDS(file.path(Path, Fname))
  names(SelectedList) <- gsub(".rds", "", files[1])
  return(SelectedList) }

#' Denote DE samples
#'
#' @param dataDiff A list of summarizedExperiment objects
#' @param alpha The significance cutoff level to use.
#' @param lfc The log2 fold change cutoff to use.
#' @param sigCol Either "p.val", "p.adj", or "BHCorrection" indicating the significance metric to use.
#' @import DEP
#' @export
SigDEAnnotate <- function(dataDiff, alpha = 0.05, lfc = log2(2), sigCol="p.adj"){
  dep <- list()
  for(i in 1:length(dataDiff)){ dep[[i]] <- add_rejections2(dataDiff[[i]], alpha = alpha, lfc = lfc, sigCol=sigCol) }
  names(dep) <- names(dataDiff)
  return(dep) }

#' Denote DE samples
#'
#' @param diff A list of summarizedExperiment objects
#' @param alpha The significance cutoff level to use.
#' @param lfc The log2 fold change cutoff to use.
#' @param sigCol Either "p.val", "p.adj", or "BHCorrection" indicating the significance metric to use.
#' @import DEP
#' @export
add_rejections2 <- function (diff=dataDiff[[i]], alpha = 0.05, lfc = 1, sigCol="p.adj")  { # "p.val"
  if (is.integer(alpha)){alpha <- as.numeric(alpha)}
  if (is.integer(lfc)){lfc <- as.numeric(lfc)}
  assertthat::assert_that(inherits(diff, "SummarizedExperiment"), is.numeric(alpha), length(alpha) == 1, is.numeric(lfc), length(lfc) == 1)
  row_data <- rowData(diff, use.names = FALSE) %>% as.data.frame()
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present in '", deparse(substitute(diff)), "'\nRun make_unique() and make_se() to obtain the required columns", call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '", deparse(substitute(diff)), "'\nRun test_diff() to obtain the required columns", call. = FALSE)
  }
  cols_p <- grep(sigCol, colnames(row_data))
  cols_diff <- grep("_diff", colnames(row_data))
  if (length(cols_p) == 1) {
    rowData(diff)$significant <- row_data[, cols_p] <= alpha & abs(row_data[, cols_diff]) >= lfc
    rowData(diff)$contrast_significant <- rowData(diff, use.names = FALSE)$significant
    colnames(rowData(diff))[ncol(rowData(diff, use.names = FALSE))] <- gsub(sigCol, "significant", colnames(row_data)[cols_p])
  }
  if (length(cols_p) > 1) {
    p_reject <- row_data[, cols_p] <= alpha
    p_reject[is.na(p_reject)] <- FALSE
    diff_reject <- abs(row_data[, cols_diff]) >= lfc
    diff_reject[is.na(diff_reject)] <- FALSE
    sign_df <- p_reject & diff_reject
    sign_df <- cbind(sign_df, significant = apply(sign_df, 1, function(x) any(x)))
    colnames(sign_df) <- gsub(paste("_", sigCol, sep = ""), "_significant", colnames(sign_df))
    sign_df <- cbind(name = row_data$name, as.data.frame(sign_df))
    rowData(diff) <- merge(rowData(diff, use.names = FALSE),
                           sign_df, by = "name")
  }
  return(diff) }

#' Volcano plot
#'
#' @param DEPList A list of summarizedExperiment objects
#' @param ComparisonList A named vector of comparisons that were made.
#' @param ymin The minimum y axis value.
#' @param ymax The maximum y axis value.
#' @param xmin The minimum x axis value.
#' @param xmax The maximum x axis value.
#' @param addNames A logical value indicating if data point names should be added to the plot.
#' @param Adjusted A logical value indicating if adjusted significnce values should be plotted.
#' @param BHadjust A logical value indicating if the Benjamini–Hochberg method of adjusting p-values should be plotted.
#' @import DEP
#' @export
VolcanoPlot <- function(DEPList = DEPList, ComparisonList=ComparisonList, ymin = 0, ymax = 15, xmin = -3, xmax = 3, addNames = TRUE, Adjusted = TRUE, BHadjust){
  PlotList <- list()
  comps <- names(DEPList)# for(i in 1:length(DEPList)){ # comps <- ComparisonList[names(ComparisonList) == names(DEPList)[i]][[1]]#  gsub("_normalized", "", gsub("_Raw", "", gsub("scaled_", "", names(DEPList)[i] ) )) ][[1]]
  for(b in 1:length(comps)){
    tempLit <- list()
    nam <- names(DEPList)[b]#paste(names(DEPList)[i], comps[b], sep = "-")
    tempLit[[1]] <- plot_volcano2(dep = DEPList[[b]], contrast = gsub(".+-", "", comps[b]), label_size = 3, add_names = addNames, adjusted = Adjusted, BHadjusted = BHadjust) + ggtitle(nam) +
      ylim(ymin,ymax) + xlim(xmin,xmax)
    names(tempLit) <- nam
    PlotList <- c(PlotList, tempLit)
  } # }
  return(PlotList) }

#' Volcano plot
#'
#' @param DEP A summarizedExperiment object
#' @param contrast A named vector of comparisons that were made.
#' @param label_size The size of the labels on the plot.
#' @param add_names A logical value indicating if data point names should be added to the plot.
#' @param adjusted A logical value indicating if adjusted significnce values should be plotted.
#' @param BHadjusted A logical value indicating if the Benjamini–Hochberg method of adjusting p-values should be plotted.
#' @param plot A logical value indicating if a plot should be returned.
#' @import DEP
#' @export
plot_volcano2 <- function (dep, contrast, label_size = 3, add_names = TRUE, adjusted = FALSE, plot = TRUE, BHadjusted) {
  if (is.integer(label_size))
    label_size <- as.numeric(label_size)
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.character(contrast), length(contrast) == 1, is.numeric(label_size),
                          length(label_size) == 1, is.logical(add_names), length(add_names) ==
                            1, is.logical(adjusted), length(adjusted) == 1, is.logical(plot),
                          length(plot) == 1)
  row_data <- rowData(dep, use.names = FALSE)
  if(any(!c("name", "ID") %in% colnames(row_data))){
    stop(paste0("'name' and/or 'ID' columns are not present in '", deparse(substitute(dep)), "'.\nRun make_unique() to obtain required columns."), call. = FALSE)
  }
  if(length(grep("_p.adj|_diff", colnames(row_data))) < 1){
    stop(paste0("'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '", deparse(substitute(dep)), "'.\nRun test_diff() to obtain the required columns."), call. = FALSE)
  }
  if(length(grep("_significant", colnames(row_data))) < 1){
    stop(paste0("'[contrast]_significant' columns are not present in '", deparse(substitute(dep)), "'.\nRun add_rejections() to obtain the required columns."), call. = FALSE)
  }
  if(length(grep(paste(contrast, "_diff", sep = ""), colnames(row_data))) ==  0) {
    valid_cntrsts <- row_data %>% data.frame() %>% select(ends_with("_diff")) %>% colnames(.) %>% gsub("_diff", "", .)
    valid_cntrsts_msg <- paste0("Valid contrasts are: '", paste0(valid_cntrsts, collapse = "', '"), "'")
    stop("Not a valid contrast, please run `plot_volcano()` with a valid contrast as argument\n",
         valid_cntrsts_msg, call. = FALSE) }
  diff <- grep(paste(contrast, "_diff", sep = ""), colnames(row_data))
  if (adjusted) { p_values <- grep(paste(contrast, "_p.adj", sep = ""), colnames(row_data))
  } else if(BHadjusted){ p_values <- grep(paste(contrast, "_BHCorrection", sep = ""), colnames(row_data))
  } else {
    p_values <- grep(paste(contrast, "_p.val", sep = ""), colnames(row_data)) }
  signif <- grep(paste(contrast, "_significant", sep = ""), colnames(row_data))
  df <- data.frame(x = row_data[, diff], y = -log10(row_data[, p_values]), significant = row_data[, signif], name = row_data$name) %>%
    filter(!is.na(significant)) %>% arrange(significant)
  name1 <- gsub("_vs_.*", "", contrast)
  name2 <- gsub(".*_vs_", "", contrast)
  p <- ggplot(df, aes(x, y)) + geom_vline(xintercept = 0) +
    geom_point(aes(col = significant)) + geom_text(data = data.frame(),
                                                   aes(x = c(Inf, -Inf), y = c(-Inf, -Inf), hjust = c(1, 0), vjust = c(-1, -1), label = c(name1, name2), size = 5,
                                                       fontface = "bold")) + labs(title = contrast, x = expression(log[2] ~ "Fold change")) + theme_DEP1() + theme(legend.position = "none") +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "grey"))
  if (add_names) {
    p <- p + ggrepel::geom_text_repel(data = filter(df, significant),
                                      aes(label = name), size = label_size, box.padding = unit(0.1, "lines"), point.padding = unit(0.1, "lines"),
                                      segment.size = 0.5) }
  if (adjusted) { p <- p + labs(y = expression(-log[10] ~ "Adjusted p-value"))
  } else if(BHadjusted){ p <- p + labs(y = expression(-log[10] ~ "BH Adjusted p-value"))
  } else { p <- p + labs(y = expression(-log[10] ~ "P-value")) }
  if (plot) { return(p) } else {
    df <- df %>% select(name, x, y, significant) %>% arrange(desc(x))
    colnames(df)[c(1, 2, 3)] <- c("protein", "log2_fold_change", "p_value_-log10")
    if (adjusted) { colnames(df)[3] <- "adjusted_p_value_-log10" }
    if (BHadjusted) { colnames(df)[3] <- "BH_adjusted_p_value_-log10" }
    return(df) } }









