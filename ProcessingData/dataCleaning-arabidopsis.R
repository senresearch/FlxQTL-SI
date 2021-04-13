# data imputation for arabidopsis data


library(qtl)
# load functions to extract marker information
source("getMarkerInfo.R")

## for 1-d genome scan
# impute genotypes and extract the imputed. 
# genotypes: a =1 (Italian parent). b=2(Swedish parent)
data1<-read.cross("csv",dir="../rawData",file="agren2013_fullGenoMatrix4qtl_withParents.csv",genotypes = c("a","b"),crosstype = "riself")


## adding pseudo markers (at every 1cM between two markers) to genotype data  

#imputation 
data.sim<-sim.geno(data1,step=1,n.draw=16)
data_imp<-pull.draws(data.sim)

marlabels<-lapply(data.sim$geno,getmap)
#markers<-unlist(markerlabels)
Markers<-getChr(marlabels)

write.csv(Markers,file="../processedData/marlabels_agren.csv")
#pick any imputed genotype dataset
write.csv(data_imp[,,5],file = "../processedData/fullrank_imput.csv",row.names = F)




## for 2-d genome scan
#generate pseudo markers at every 5cM (less refined)
mul.sim<-sim.geno(data1,step=5,n.draws = 16)
mul.imp<-pull.draws(mul.sim)
write.csv(mul.imp[,,7],file="../processedData/agren_gene_imp.csv",row.names = F)

#get marker information :marker names, chromosomes, positions
marker.mul<-lapply(mul.sim$geno, getmap)
mul.marlabel<-rbind(cbind(rep(1,length(marker.mul$`1`)),marker.mul$`1`),cbind(rep(2,length(marker.mul$`2`)),marker.mul$`2`),
                    cbind(rep(3,length(marker.mul$`3`)),marker.mul$`3`),cbind(rep(4,length(marker.mul$`4`)),marker.mul$`4`),
                    cbind(rep(5,length(marker.mul$`5`)),marker.mul$`5`))

write.csv(mul.marlabel,file="../processedData/agren_markers_multipleqtl.csv")


# extract phenotype data (missing values included)
phen<-pull.pheno(data1)

### impute phenotype data
library(mice)
library(lattice)
fitness<-phen[,1:6]
fitness_imp<-mice(fitness,defaultMethod = "pmm",seed = 10)
pheno<-complete(fitness_imp)
write.csv(pheno,file="../processedData/pheno2013_imp.csv",row.names = FALSE)


