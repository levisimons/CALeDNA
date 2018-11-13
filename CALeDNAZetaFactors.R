require(zetadiv)
require(plyr)

setwd("~/Desktop/CALeDNA") #Change on the cluster to ~/panfs/CALeDNA

fileList <- list.files(pattern="asv_deco_dedup") #Get all of the metagenomic data tables to calculate zeta_4 diversity for given areas.
file <- fileList[5]

metadata <- read.table("Final_metadata.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE) #CALeDNA site metadata.
#Add in the number of samples per location in the metadata.
locationFrequency <- as.data.frame(table(metadata$loc))
colnames(locationFrequency) <- c("loc","Freq")
metadata <- join(metadata,locationFrequency,by=c("loc"))

#Read in OTU table.
DataRaw <- read.table(file, header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
#Get primer name.
primer <- strsplit(file,"_")
primer <- primer[[1]][length(primer[[1]])]
primer <- strsplit(primer,".csv")
primer <- primer[[1]][1]
#Get OTU names
taxa <- DataRaw[,1]
#Transpose the OTU table and merge it with the sample metadata
DataT <- as.data.frame(t(DataRaw[,-c(1)]))
colnames(DataT) <- taxa
DataT$MatchName <- rownames(DataT)
MergedData <-join(DataT,metadata,by=c("MatchName"))
#Filter out areas with less than a certain number of a samples.
zetaNum=4
data.spec <- subset(MergedData,Freq>=zetaNum)
data.spec <- data.spec[,colnames(data.spec) %in% taxa]
#Create sample location data frame.
data.xy <- subset(MergedData,Freq>=zetaNum)
data.xy <- data.xy[,c("Latitude","Longitude")]
#Create sample factor data frame.
data.env <- subset(MergedData,Freq>=zetaNum)
data.env <- data.env[,c("elev","Slope","greenness")]

zetaFactors <- Zeta.msgdm(data.spec=data.spec,data.env=data.env,xy=data.xy,order=2,sam=nrow(data.spec),reg.type="glm",distance.type="ortho",normalize=FALSE)
