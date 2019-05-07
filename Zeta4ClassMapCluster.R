# /usr/usc/R/3.4.4/bin/Rscript to run on the USC cluster
require(zetadiv)
require(plyr)

#This script creates the zeta diversity decay parameters for geographically clustered sites
#with OTUs defined at the family level

##To generate the input data for mapping.
wd <- "~/cmb/CALeDNA" #Change on the cluster to ~/panfs/CALeDNA
setwd(wd)
metadata <- read.table("Final_metadata.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
metadata$clust <- as.factor(as.character(metadata$clust))
#Get the groups for the type of environmental factors.
factorGroups <- read.table("Metadata_explanation.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Unique environmental factors are as follows:
uniqueFactors <- c("Location","Topology","Habitat","BioClim","Soil Properties","Vegetation")
#Get the cluster ID of all of the sample sites where locations have sufficient numbers of samples.
clusteredSites <- metadata[,c("loc","MatchName","Zeta_4ID")]
clusteredSites <- subset(clusteredSites,is.na(clusteredSites$Zeta_4ID)==FALSE)

#Get all of the metagenomic data tables, clustered by family, within a primer set.
primerList <- c("16S","18S","PITS","FITS","CO1")
zetaNum <- 4 #The maximum order to calculate zeta diversity
zetaLocation <- data.frame()
for(primer in primerList){
  primerFiles <- list.files(path=paste(wd,"/ASVbyClass",sep=""),full.names=TRUE,pattern=paste("asv_",primer,sep=""))
  OTUClass <- data.frame()
  for(primerFile in primerFiles){
    OTURaw <- read.table(primerFile, header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE) #Individual OTU table by class.
    if(nrow(OTURaw)>1){
      OTURaw <- OTURaw[,-c(1)]
      siteNames <- colnames(OTURaw)
      OTURaw <- sapply(OTURaw,as.numeric)
      OTURow <- as.data.frame(t(as.data.frame(colSums(OTURaw))))
      className <- strsplit(strsplit(primerFile,"_")[[1]][3],".csv")[[1]][1]
      rownames(OTURow) <- className
      colnames(OTURow) <- siteNames
      OTUClass <- rbind(OTUClass,OTURow)
    }
  }
  #Create OTU tables by primer group, clustered to class.
  OTUClass[OTUClass > 0] <- 1
  OTUClass <- as.data.frame(t(OTUClass))
  #Calculate zeta diversity decay, up to a specified order, for the communities within each geographic cluster
  # with community presence/absences defined at the family level.
  for(clusterID in unique(clusteredSites$Zeta_4ID)){
    uniqueCluster <- clusteredSites[clusteredSites$Zeta_4ID==clusterID,]
    DataLocationSubset <- OTUClass[rownames(OTUClass) %in% uniqueCluster$MatchName,]
    if(nrow(DataLocationSubset) >= zetaNum){
      zetaDecay <- Zeta.decline.ex(DataLocationSubset,order=1:zetaNum,rescale=FALSE,plot=FALSE)
      zetaDecayScaled <- Zeta.decline.ex(DataLocationSubset,order=1:zetaNum,rescale=TRUE,plot=FALSE)
      loc <- unique(uniqueCluster$loc)
      dat <- data.frame()
      dat[1,1] <- primerFile
      dat[1,2] <- clusterID
      dat[1,3] <- primer 
      dat[1,4] <- zetaDecay$zeta.val[zetaNum]
      dat[1,5] <- zetaDecay$zeta.val.sd[zetaNum]
      dat[1,6] <- zetaDecayScaled$zeta.val[zetaNum]
      dat[1,7] <- zetaDecayScaled$zeta.val.sd[zetaNum]
      dat[1,8] <- zetaDecay$zeta.exp$coefficients[1]
      dat[1,9] <- zetaDecay$zeta.exp$coefficients[2]
      dat[1,10] <- zetaDecay$aic$AIC[1]
      dat[1,11] <- zetaDecay$zeta.pl$coefficients[1]
      dat[1,12] <- zetaDecay$zeta.pl$coefficients[2]
      dat[1,13] <- zetaDecay$aic$AIC[2]
      zetaLocation <- rbind(zetaLocation,dat)
      print(dat)
    }
  }
}
colnames(zetaLocation) <- c("file","clusterID","primer",paste("zeta",zetaNum,sep=""),paste("zeta",zetaNum,"sd",sep=""),paste("zeta",zetaNum,"scaled",sep=""),paste("zeta",zetaNum,"sdscaled",sep=""),"ExpIntercept","ExpExp","ExpAIC","PLIntercept","PLExp","PLAIC")
zetaLocation[is.na(zetaLocation)] <- NA
write.table(zetaLocation,"Zeta4ClassMapData.txt",quote=FALSE,sep="\t",row.names = FALSE)
