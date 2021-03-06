# /usr/usc/R/3.4.4/bin/Rscript to run on the USC cluster
require(zetadiv)
require(plyr)
require(ggmap)
require(maps)
require(mapview)
require(mapdata)
require(munsell)
require(leaflet)
require(devtools)
require(webshot)
require(viridis)

#This script creates the zeta diversity decay parameters for geographically clustered sites
#with OTUs defined at the family level

##To generate the input data for mapping.
wd <- "~/Desktop/CALeDNA" #Change on the cluster to ~/panfs/CALeDNA
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
    #Calculate zeta diversity decay, up to a specified order, for the communities within each geographic cluster
    # with community presence/absences defined at the family level.
    for(clusterID in unique(clusteredSites$Zeta_4ID)){
      uniqueCluster <- clusteredSites[clusteredSites$Zeta_4ID==clusterID,]
      DataLocationSubset <- OTUFamily[rownames(OTUFamily) %in% uniqueCluster$MatchName,]
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
}
colnames(zetaLocation) <- c("file","clusterID","primer",paste("zeta",zetaNum,sep=""),paste("zeta",zetaNum,"sd",sep=""),paste("zeta",zetaNum,"scaled",sep=""),paste("zeta",zetaNum,"sdscaled",sep=""),"ExpIntercept","ExpExp","ExpAIC","PLIntercept","PLExp","PLAIC")
zetaLocation[is.na(zetaLocation)] <- NA
write.table(zetaLocation,"Zeta4FamilyMapData.txt",quote=FALSE,sep="\t",row.names = FALSE)

##To create the maps of zeta diversity.
#Choose a primer to subset zeta diversity results.
zetaLocation <- read.table("Zeta4FamilyMapData.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
primer <- "16S" #Choose a primer
zetaPrimerSubset <- zetaLocation[which(zetaLocation$primer==primer),]
#Merge zeta diversity data with sample site metadata.
zetaAnalysis <- left_join(zetaPrimerSubset,metadata,by=c("clusterID"="Zeta_4ID"))

#To generate map of data for a given zeta diversity parameter in California.
dev.off()
MapCoordinates <- zetaAnalysis
#Map data.
CalMap = leaflet(MapCoordinates) %>% 
  addTiles()
ColorScale <- colorNumeric(palette=plasma(10),domain=MapCoordinates$zeta4scaled)
CalMap %>% addCircleMarkers(color = ~ColorScale(zeta4scaled), fill = TRUE,radius=0.1,fillOpacity = 0.1) %>% 
  setView(median(MapCoordinates$Longitude),median(MapCoordinates$Latitude),zoom=5) %>%
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addLegend("topright", pal=ColorScale,values=~zeta4scaled,title=paste(primer,"Scaled &#950;<sub>4</sub></b><br>family diversity"))
