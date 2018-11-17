# /usr/usc/R/3.4.4/bin/Rscript to run on the USC cluster
require(zetadiv)
require(plyr)
require(dplyr)

#This script analyzes zeta diversity patterns for community groups defined by presence/absence data clustered to class level.
#CALeDNA site metadata.
setwd("~/Desktop/CALeDNA") #Change on the cluster to ~/panfs/CALeDNA
metadata <- read.table("Final_metadata.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Get the cluster ID of all of the sample sites where locations have sufficient numbers of samples.
clusteredSites <- metadata[,c("MatchName","Zeta_4ID")]
clusteredSites <- subset(clusteredSites,is.na(clusteredSites$Zeta_4ID)==FALSE)

setwd("~/Desktop/CALeDNA/ASVbyClass") #Change on the cluster to ~/panfs/CALeDNA/ASVbyClass

#Get all of the metagenomic data tables, split by class, within a primer set.
Files16S <- list.files(pattern="asv_16S")
Files18S <- list.files(pattern="asv_18S")
FilesCO1 <- list.files(pattern="asv_CO1")
FilesPITS <- list.files(pattern="asv_PITS")
FilesFITS <- list.files(pattern="asv_FITS")
FilesAll <- list(Files16S,Files18S,FilesCO1,FilesPITS,FilesFITS)

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

#Create subset presence/absence tables for each geographic cluster.
zetaNum=4 #Highest order of zeta diversity.
for(clusterID in unique(clusteredSites$Zeta_4ID)){
  localCluster <- subset(clusteredSites,clusteredSites$Zeta_4ID==clusterID)
  localClusterOTUs <- OTUClass[,colnames(OTUClass) %in% localCluster$MatchName]
  zetaDecay <- Zeta.decline.ex(localClusterOTUs,order=1:zetaNum,rescale=TRUE,plot=TRUE)
  dat <- data.frame()
  dat[1,1] <- clusterID
  dat[1,2] <- zetaDecay$zeta.val[zetaNum]
  dat[1,3] <- zetaDecay$zeta.val.sd[zetaNum]
  dat[1,4] <- zetaDecay$zeta.exp$coefficients[1]
  dat[1,5] <- zetaDecay$zeta.exp$coefficients[2]
  dat[1,6] <- zetaDecay$aic$AIC[1]
  dat[1,7] <- zetaDecay$zeta.pl$coefficients[1]
  dat[1,8] <- zetaDecay$zeta.pl$coefficients[2]
  dat[1,9] <- zetaDecay$aic$AIC[2]
  print(dat)
}
