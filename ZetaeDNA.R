# /usr/usc/R/3.4.4/bin/Rscript to run on the USC cluster
require(zetadiv)

setwd("~/Desktop/CALeDNA") #Change on the cluster to ~/panfs/CALeDNA

fileList <- list.files(pattern=".csv")

for(file in fileList){
  print(file)
  filename <- gsub(".csv","",file)
  #Read in OTU table.
  DataRaw <- read.table(paste(filename,".csv",sep=""), header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
  #Get OTU names
  taxa <- DataRaw[,1]
  #Make a presence/absence matrix
  DataPA <- as.data.frame(t(DataRaw[,-c(1)]))
  colnames(DataPA) <- taxa
  DataPA[DataPA > 0] <- 1
  #Calculate zeta diversity decay parameters and plots.
  dat <- data.frame()
  zetaDecay <- Zeta.decline.ex(DataPA[,1:100],orders=1:10,plot=TRUE)
  dat[1,1] <- zetaDecay$zeta.exp$coefficients[1] #Zeta diversity exponential decay intercept.
  dat[1,2] <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
  dat[1,3] <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
  dat[1,4] <- zetaDecay$zeta.pl$coefficients[1] #Zeta diversity power law decay intercept.
  dat[1,5] <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
  dat[1,6] <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
  colnames(dat) <- c("ExpIntercept","ExpExP","ExpAIC","PLIntercept","PLExp","PLAIC")
  write.table(dat,paste(filename,"ZetaModels.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
  print(dat)
  
  dat <- as.data.frame(zetaDecay$zeta.val)
  colnames(dat) <- c("Zeta_N")
  write.table(dat,paste(filename,"ZetaValues.txt",sep=""),quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
  
  file.rename("Rplots.pdf",paste(filename,"ZetaModels.pdf",sep=""))
  dev.off()
}

