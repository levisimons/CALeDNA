# /usr/usc/R/3.4.4/bin/Rscript to run on the USC cluster
require(zetadiv)
require(dplyr)

setwd("~/Desktop/CALeDNA") #Change on the cluster to ~/panfs/CALeDNA

fileList <- list.files(pattern="Oct25_2018_zeta4.csv") #Get all of the metagenomic data tables to calculate zeta_4 diversity for given areas.

metadata <- read.table("Transect_all_283_samples_metadata_Zeta4.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE) #CALeDNA site metadata.
#Add in the number of samples per location in the metadata.
locationFrequency <- as.data.frame(table(metadata$loc))
colnames(locationFrequency) <- c("loc","Freq")
metadata <- join(metadata,locationFrequency,by=c("loc"))

for(file in fileList){
  #print(file)
  #Read in OTU table.
  DataRaw <- read.table(file, header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
  #Get OTU names
  taxa <- DataRaw[,1]
  #Transpose the OTU table and merge it with the sample metadata
  DataT <- as.data.frame(t(DataRaw[,-c(1)]))
  colnames(DataT) <- taxa
  DataT$sum.taxonomy <- rownames(DataT)
  DataPlus <-join(DataT,metadata,by=c("sum.taxonomy"))
  #Filter out locations with less than 4 samples.
  DataPlus <- subset(DataPlus,Freq>=4)
  for(location in unique(DataPlus$loc)){
    #print(paste(file,location))
    DataLocationSubset <- subset(DataPlus, loc==location)
    #Now remove metadata elements.
    clipColumn <- which(colnames(DataLocationSubset)=="sum.taxonomy")
    DataLocationSubset <- DataLocationSubset[,-c(clipColumn:ncol(DataLocationSubset))]
    #Make a presence/absence matrix
    DataLocationSubset[DataLocationSubset > 0] <- 1
    zetaN <- Zeta.order.ex(DataLocationSubset,order=3,rescale=TRUE)
    if(zetaN$zeta.val==0 | is.na(zetaN$zeta.val)==TRUE){
      zetaN_scaled=0
    } else{
      zeta1 <- Zeta.order.ex(DataLocationSubset,order=1,rescale=TRUE)
      zetaN_scaled = zetaN$zeta.val/zeta1$zeta.val
    }
    print(paste(file,location,zetaN$zeta.val,zetaN$zeta.val.sd,zetaN_scaled))
   }
}




