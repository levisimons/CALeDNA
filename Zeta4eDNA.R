# /usr/usc/R/3.4.4/bin/Rscript to run on the USC cluster
require(zetadiv)
require(plyr)

setwd("~/Desktop/CALeDNA") #Change on the cluster to ~/panfs/CALeDNA

fileList <- list.files(pattern="asv_deco_dedup") #Get all of the metagenomic data tables to calculate zeta_4 diversity for given areas.

metadata <- read.table("Final_metadata.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE) #CALeDNA site metadata.
#Add in the number of samples per location in the metadata.
locationFrequency <- as.data.frame(table(metadata$loc))
colnames(locationFrequency) <- c("loc","Freq")
metadata <- join(metadata,locationFrequency,by=c("loc"))

zetaLocation <- data.frame()
zetaNum=4

for(file in fileList){
  #print(file)
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
  DataPlus <-join(DataT,metadata,by=c("MatchName"))
  #Filter out locations with less than 4 samples.
  DataPlus <- subset(DataPlus,Freq>=zetaNum)
  for(location in unique(DataPlus$loc)){
    #print(paste(file,location))
    DataLocationSubset <- subset(DataPlus, loc==location)
    #Now remove metadata elements.
    clipColumn <- which(colnames(DataLocationSubset)=="MatchName")
    DataLocationSubset <- DataLocationSubset[,-c(clipColumn:ncol(DataLocationSubset))]
    #Make a presence/absence matrix
    DataLocationSubset[DataLocationSubset > 0] <- 1
    #Calculate zeta decay parameters and a particular order of zeta.
    zetaDecayScaled <- Zeta.decline.ex(DataLocationSubset,order=1:zetaNum,rescale=TRUE,plot=FALSE)
    zetaDecay <- Zeta.decline.ex(DataLocationSubset,order=1:zetaNum,rescale=FALSE,plot=FALSE)
    dat <- data.frame()
    dat[1,1] <- file
    dat[1,2] <- location
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
colnames(zetaLocation) <- c("file","location","primer",paste("zeta",zetaNum,sep=""),paste("zeta",zetaNum,"sd",sep=""),paste("zeta",zetaNum,"scaled",sep=""),paste("zeta",zetaNum,"sdscaled",sep=""),"ExpIntercept","ExpExp","ExpAIC","PLIntercept","PLExp","PLAIC")
zetaLocation[is.na(zetaLocation)] <- NA

#write.table(zetaLocation,"CALeDNAZetaDecay.txt",quote=FALSE,sep="\t",row.names = FALSE)

#If the file CALeDNAZetaDecay.txt was already generated on the cluster read it in here.
zetaLocation <- read.table("CALeDNAZetaDecay.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
