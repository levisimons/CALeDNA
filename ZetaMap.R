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

#If the file CALeDNAZetaDecay.txt was already generated on the cluster read it in here.
zetaLocation <- read.table("CALeDNAZetaDecay.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE)

#Choose a primer to subset zeta diversity results.
primerList <- unique(zetaLocation$primer)
primer <- "16S"
zetaPrimerSubset <- zetaLocation[which(zetaLocation$primer==primer),]
#Merge zeta diversity data with sample site metadata.
zetaAnalysis <- join(zetaPrimerSubset,metadata,by=c("location"))

#To generate map of data for a given zeta diversity parameter in California.
dev.off()
MapCoordinates <- zetaAnalysis
#Map data.
CalMap = leaflet(MapCoordinates) %>% 
  addTiles()
ColorScale <- colorNumeric(palette=rainbow(10),domain=MapCoordinates$zeta4scaled)
CalMap %>% addCircleMarkers(color = ~ColorScale(zeta4scaled), fill = TRUE,radius=0.1,fillOpacity = 0.1) %>% 
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addLegend("topright", pal=ColorScale,values=~zeta4scaled,title=paste(primer,"Scaled Zeta_4"))
