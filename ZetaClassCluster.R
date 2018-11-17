# /usr/usc/R/3.4.4/bin/Rscript to run on the USC cluster
require(zetadiv)
require(plyr)
require(dplyr)

#This script analyzes zeta diversity patterns for community groups defined by presence/absence data clustered to class level.
setwd("~/Desktop/CALeDNA/ASVbyClass") #Change on the cluster to ~/panfs/CALeDNA/ASVbyClass

#Get all of the metagenomic data tables, split by class, within a primer set.
Files16S <- list.files(pattern="asv_16S")
Files18S <- list.files(pattern="asv_18S")
FilesCO1 <- list.files(pattern="asv_CO1")
FilesPITS <- list.files(pattern="asv_PITS")
FilesFITS <- list.files(pattern="asv_FITS")

#Create a presence/absence OTU table, clustered by class, within a primer set.
OTUClass <- data.frame()
for(File16S in Files16S){
  OTURaw <- read.table(File16S, header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE) #Individual OTU table by class.
  if(nrow(OTURaw)>1){
    OTURaw <- OTURaw[,-c(1)]
    siteNames <- colnames(OTURaw)
    OTURaw <- sapply(OTURaw,as.numeric)
    OTURow <- as.data.frame(t(as.data.frame(colSums(OTURaw))))
    className <- strsplit(strsplit(File16S,"_")[[1]][3],".csv")[[1]][1]
    rownames(OTURow) <- className
    colnames(OTURow) <- siteNames
    OTUClass <- rbind(OTUClass,OTURow)
  }
}
OTUClass[OTUClass > 0] <- 1
