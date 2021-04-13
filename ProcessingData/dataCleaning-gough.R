library(qtl)
# load functions to extract marker information
source("getMarkerInfo.R")

#read F2 intercross data after downloading from the url
gouf<-read.cross("csv",file = "../goughF2/goughF2.csv",genotypes = c("GG","GW","WW"),alleles = c("G","W"))

############### quick analysis to check missing values ####################
par(mfrow=c(2,1), las=1)
hist(nmissing(gouf), main="No. missing genotypes by mouse",
           xlab="No. missing genotypes", breaks=140)
rug(nmissing(gouf), col="blue")
hist(nmissing(gouf, what="mar"), main="No. missing genotypes by marker",
             xlab="No. missing genotypes", breaks=50)
rug(nmissing(gouf, what="mar"), col="blue")

plotMissing(gouf)
par(mfrow=c(2,1), las=1)
plot(ntyped(gouf), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(gouf,"mar"), ylab="No. typed individuals", main="No. genotypes by marker")

#
#check if any column have missing values
table(colSums(is.na(gouf$pheno)))
##########################################


### produce genotype probability and phenotype data including missing values
##impute genotypes, compute corresponding genotype probabilities
go.sim<-sim.geno(gouf,step=1,n.draw=16)
goF<-calc.genoprob(go.sim,step = 1)

#extract genotype & prob data excluding `X` chromosome
impg_auto<-pull.draws(go.sim,chr = "-X")
gp_auto<-pull.genoprob(goF,chr="-X")

## get marker info
ginfo<-lapply(goF$geno,getmap)
markerinfo<-getChr(ginfo)
write.csv(markerinfo,file = "../processedData/goughF2_markers.csv")

#write genotype & prob files
write.csv(impg_auto[,,7],file = "../processedData/goughF2_imp_genotypes_auto.csv",row.names=FALSE)
write.csv(gp_auto,file = "../processedData/goughF2_imp_genoprob.csv",row.names = FALSE)


##for 2d genescan 
go2d.sim<-sim.geno(gouf,step=0,stepwidth = "fixed",n.draws = 16)
go2d<-calc.genoprob(go2d.sim,stepwidth = "fixed",step = 0)
gp_2d<-pull.genoprob(go2d,chr="-X")
ginfo<-lapply(go2d.sim$geno,getmap)
markerinfo<-getChr(ginfo)
write.csv(gp_2d,file = "../processedData/goughF2_genoprob_2d_step0_auto.csv",row.names = FALSE)
write.csv(markerinfo,file = "../processedData/goughF2_markerinfo_2d_step0.csv")

###



# extract phenotype data (missing values included)
phen<-pull.pheno(goF)

### impute phenotype data
library(mice)
library(lattice)
wght<-phen[,1:16]
wght_imp<-mice(wght,m=15,defaultMethod = "pmm",seed = 10)
pheno<-complete(wght_imp)
data<-cbind(phen$sex,pheno)
write.csv(data,file="../processedData/goughF2_sex_imp_16weight.csv",row.names = FALSE)


