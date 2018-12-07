# /usr/usc/R/3.4.4/bin/Rscript to run on the USC cluster
require(zetadiv)
require(plyr)
require(dplyr)

#This script analyzes zeta diversity patterns for community groups defined by presence/absence data clustered to family level.
#CALeDNA site metadata.
wd <- "~/Desktop/CALeDNA" #Change on the cluster to ~/panfs/CALeDNA
setwd(wd)
metadata <- read.table("Final_metadata.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
metadata$clust <- as.factor(as.character(metadata$clust))
#Get the groups for the type of environmental factors.
factorGroups <- read.table("Metadata_explanation.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Unique environmental factors are as follows:
uniqueFactors <- c("Location","Topology","Habitat","BioClim","Soil Properties","Vegetation")
#Get the cluster ID of all of the sample sites where locations have sufficient numbers of samples.
clusteredSites <- metadata[,c("MatchName","Zeta_4ID")]
clusteredSites <- subset(clusteredSites,is.na(clusteredSites$Zeta_4ID)==FALSE)

#Get all of the metagenomic data tables, clustered by family, within a primer set.
primerList <- c("16S","18S","PITS","FITS","CO1")
zetaAnalysis <- data.frame()
for(primer in primerList){
  primerFile <- paste("asv_fam_",primer,"_all.csv",sep="")
  OTUFamily <- data.frame()
  #Individual OTU table by family.
  OTURaw <- read.table(primerFile, header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
  if(nrow(OTURaw)>1){
    familyNames <- OTURaw[,1]
    OTURaw <- OTURaw[,-c(1)]
    siteNames <- colnames(OTURaw)
    OTURaw <- sapply(OTURaw,as.numeric)
    OTUFamily <- as.data.frame(t(as.data.frame(OTURaw)))
    colnames(OTUFamily) <- familyNames
    #Create OTU tables by primer group, clustered to family.
    OTUFamily[OTUFamily > 0] <- 1
    #Determine the contributions to the variation in zeta diversity due to
    #geographic separation, environmental factor group, or some unknown factor.
    for(uniqueFactor in uniqueFactors){
      factorSubset <- subset(factorGroups,Category==uniqueFactor | Category=="Coordinates")
      metadataSubset <- metadata
      rownames(metadataSubset) <- metadataSubset$MatchName
      metadataSubset <- subset(metadataSubset,is.na(metadataSubset$Zeta_4ID)==FALSE) #Keep sites with a geographic cluster ID.
      metadataSubset <- metadataSubset[,colnames(metadataSubset) %in% factorSubset$`Column name` | colnames(metadataSubset) %in% c("Latitude","Longitude")] #Keep factors by group, along with GPS coordinates.
      metadataSubset <- metadataSubset[complete.cases(metadataSubset),] #Remove rows with missing data.
      data.xy <- metadataSubset[,c("Latitude","Longitude")] #Create location data frame.
      metadataSubset <- metadataSubset[,-which(names(metadataSubset) %in% c("Latitude","Longitude")),drop=FALSE] #Remove GPS coordinates from factor group for zeta diversity analysis.
      #Force metadata factors into numeric or factor types.
      metadataSubset[sapply(metadataSubset, is.character)] <- lapply(metadataSubset[sapply(metadataSubset, is.character)], as.factor)
      metadataSubset[sapply(metadataSubset, is.integer)] <- lapply(metadataSubset[sapply(metadataSubset, is.integer)], as.numeric)
      #Subset the OTU table by the remaining sample sites
      data.OTU <- OTUFamily[which(row.names(OTUFamily) %in% row.names(metadataSubset)),]
      #Determine contibutions of variance to zeta diversity from environmental factors
      zetaFactor <- Zeta.msgdm(data.spec=data.OTU,data.env=metadataSubset,xy=data.xy,order=4,method.glm="glm.fit.cons",distance.type="ortho")
      zetaVars <-  Zeta.varpart(zetaFactor,method.glm="glm.fit.cons")
      dat <- data.frame()
      dat[1,1] <- primer
      dat[1,2] <- uniqueFactor
      dat[1,3] <- nrow(metadataSubset)
      dat[1,4] <- zetaVars$`Adjusted Rsq`[6]
      dat[1,5] <- zetaVars$`Adjusted Rsq`[4]
      dat[1,6] <- zetaVars$`Adjusted Rsq`[7]
      zetaAnalysis <- rbind(zetaAnalysis,dat)
      print(paste(dat[1,1],dat[1,2],dat[1,3]))
      print(paste("Proportion variation by factor group",zetaVars$`Adjusted Rsq`[6]))
      print(paste("Proportion variation by separation distance",zetaVars$`Adjusted Rsq`[4]))
      print(paste("Proportion variation by unknown factors",zetaVars$`Adjusted Rsq`[7]))
    }
  }
}
colnames(zetaAnalysis) <- c("Primer","FactorGroup","NumSamples","VarFactor","VarDistance","VarUnknown")
write.table(zetaAnalysis,"CALeDNAZeta4Family.txt",quote=FALSE,sep="\t",row.names = FALSE)

#Create subset presence/absence tables for each geographic cluster.
#This is useful for downstream analysis such as mapping.
zetaNum=4 #Highest order of zeta diversity.
zetaCluster <- data.frame()
for(clusterID in unique(clusteredSites$Zeta_4ID)){
  localCluster <- subset(clusteredSites,clusteredSites$Zeta_4ID==clusterID)
  localClusterOTUs <- OTUFamily[,colnames(OTUFamily) %in% localCluster$MatchName]
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
  zetaCluster <- rbind(zetaCluster,dat)
}
