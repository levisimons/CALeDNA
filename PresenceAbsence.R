#Generating a presence/absence matrix, aggregated by family, for CALeDNA data.
require(dplyr)

setwd("~/Desktop/CALeDNA")
#Read in OTU counts, aggregated by family.
OTUInput <- read.table("ALL_fam_taxon_richness.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
OTUSummed <- aggregate(. ~ sum.taxonomy, OTUInput, sum)
#Sum counts across duplicate taxa entries.
taxaNames <- as.data.frame(OTUSummed[,1])
#Rename column for taxaNames.
colnames(taxaNames)[colnames(taxaNames)=="OTUSummed[, 1]"] <- "Family"
#Convert to presence/absence
PA <- OTUSummed[,-c(1)]
PA[PA>0] <- 1
#Add names back in.
PA <- cbind(taxaNames,PA)

write.table(PA,"CALeDNAFamilyPA.txt",quote=FALSE,sep="\t",row.names = TRUE)
