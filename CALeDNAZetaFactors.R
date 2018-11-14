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

setwd("~/Desktop/CALeDNA") #Change on the cluster to ~/panfs/CALeDNA

metadata <- read.table("Final_metadata.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE) #CALeDNA site metadata.
#Add in the number of samples per location in the metadata.
locationFrequency <- as.data.frame(table(metadata$loc))
colnames(locationFrequency) <- c("loc","Freq")
metadata <- join(metadata,locationFrequency,by=c("loc"))
names(metadata)[names(metadata) == "loc"] <- "location"

#If the file CALeDNAZetaDecay.txt was already generated on the cluster using Zeta4eDNA.R read it in here.
zetaLocation <- read.table("CALeDNAZetaDecay.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE)

#Choose a primer to subset zeta diversity results.
primerList <- unique(zetaLocation$primer)
primer <- "18S"
zetaPrimerSubset <- zetaLocation[which(zetaLocation$primer==primer),]
#Merge zeta diversity data with sample site metadata.
zetaAnalysis <- join(zetaPrimerSubset,metadata,by=c("location"))

#To generate map of data for a given zeta diversity parameter in California.
dev.off()
MapCoordinates <- zetaAnalysis
#Map data.
CalMap = leaflet(MapCoordinates) %>% 
  addTiles()
ColorScale <- colorNumeric(palette=rainbow(10),domain=MapCoordinates$PLExp)
CalMap %>% addCircleMarkers(color = ~ColorScale(PLExp), fill = TRUE,radius=0.1,fillOpacity = 0.1) %>% 
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addLegend("topright", pal=ColorScale,values=~PLExp,title=paste(primer,"Power-law exponent"))

#Comparison of zeta diversity parameters between communities derived using different primer sets.
zetaPrimerComparison <- data.frame()
i=0
for(primer in primerList){
  zetaPrimerSubset <- zetaLocation[which(zetaLocation$primer==primer),]
  zetaPrimerSubset <- zetaPrimerSubset[,c("location","zeta4","PLExp")]
  zetaPrimerSubset <- zetaPrimerSubset[order(zetaPrimerSubset$location),]
  colnames(zetaPrimerSubset) <- c("location",paste(primer,"zeta4",sep=""),paste(primer,"PLExp",sep=""))
  i=i+1
  if(i==1){
    zetaPrimerComparison <- zetaPrimerSubset
  } else{
    zetaPrimerComparison <- join(zetaPrimerComparison,zetaPrimerSubset,by=c("location"))
  }
}

#Regression between network parameters.
require(Hmisc)
require(corrplot)
require("PerformanceAnalytics")
#Each significance level is associated to a symbol : p-values(0, 0.001, 0.01, 0.05, 0.1, 1) <=> symbols(“***”, “**”, “*”, “.”, " “)
chart.Correlation(zetaPrimerComparison[,c("16Szeta4","18Szeta4","CO1zeta4","FITSzeta4","PITSzeta4")], histogram=FALSE, method="spearman")
chart.Correlation(zetaPrimerComparison[,c("16SPLExp","18SPLExp","CO1PLExp","FITSPLExp","PITSPLExp")], histogram=FALSE, method="spearman")
