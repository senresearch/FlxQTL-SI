
## imputation and cleaning based on genotypes data 
library(qtl)

# load functions to extract marker information
source("getMarkerInfo.R")

#load files 
crs1<-get(load("../36_site_yr_fl50/sexAveraged.cross.rda",verbose = T))


#fill compute genotype probabilities 
cross_imputed <- calc.genoprob(crs1, step=1, error.prob=0.001)

#get marker info
marker.sw<-lapply(cross_imputed$geno,getmap.prob)

#final marker info
marker.info<-getChrom(marker.sw) # 6118 markers
write.csv(marker.info,"../processedData/switchgrass36site_yr_markers.csv",quote = F,row.names = T)

#pull out genotype probabilities
geno_imp<-pull.genoprob(cross_imputed) # 750 x (6118 x 4 genotypes)
write.csv(geno_imp,"../processedData/switchgrass36site_yr_geno_prob_imputed.csv",quote = F,row.names = F)

