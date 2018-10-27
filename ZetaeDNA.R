# /usr/usc/R/3.4.4/bin/Rscript to run on the USC cluster
require(zetadiv)

setwd("~/Desktop/CALeDNA")

Data16SRaw <- read.table("16S_table_nomin_Oct25_2018.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)

taxa <- Data16SRaw[,1]

Data16S <- as.data.frame(t(Data16SRaw[,-c(1)]))
colnames(Data16S) <- taxa

Data16S[Data16S > 0] <- 1

dat <- data.frame()
zetaDecay <- Zeta.decline.ex(Data16S[,1:300],orders=1:10,plot=FALSE)
dat[1,1] <- zetaDecay$zeta.exp$coefficients[1] #Zeta diversity exponential decay intercept.
dat[1,2] <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
dat[1,3] <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
dat[1,4] <- zetaDecay$zeta.pl$coefficients[1] #Zeta diversity power law decay intercept.
dat[1,5] <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
dat[1,6] <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
print(dat)
