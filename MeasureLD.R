setwd("~/MADP")
library(rJava)
library(rTASSEL)
options(java.parameters = c("-Xmx15g","Xms15g"))


#load hapmap and phenotype data
genoPathHMP1 <- "C:/Users/natek/Documents/MADP/Genotypes/MADP_IMAGE_geno.hmp.txt"
genoPathHMP1
phenoPath<- "C:/Users/natek/Documents/MADP/Genotypes/MADP_image_pheno_RGB_EV.txt"

pheno<-read.table("MADP_image_pheno_RGB_EV.txt",header=TRUE)


tasGenoHMP<-rTASSEL::readGenotypeTableFromPath(path=genoPathHMP1)
tasGenoHMP
tasPheno<-rTASSEL::readPhenotypeFromPath(path=phenoPath)
tasPheno

tasGenoPheno <- rTASSEL::readGenotypePhenotype(
  genoPathOrObj = tasGenoHMP,
  phenoPathDFOrObj = tasPheno
)
tasGenoPheno

#calculate kinship matrix
tasKin <- rTASSEL::kinshipMatrix(tasObj = tasGenoPheno)
tasKinRMat <- rTASSEL::kinshipToRMatrix(tasKin)
tasKinRMat[1:5, 1:5]

#calculate distance matrix
tasDist <- rTASSEL::distanceMatrix(tasObj = tasGenoPheno)
tasDistRMat <- rTASSEL::distanceToRMatrix(tasDist)
tasDistRMat[1:5, 1:5]

# Calculate MLM
tasMLM1 <- rTASSEL::assocModelFitter(
  tasObj = tasGenoPheno,             # <- our prior TASSEL object
  formula = RGB_eigenvector ~ .,          # <- Phenotype Data
  fitMarkers = TRUE,                 # <- set this to TRUE for GLM
  kinship = NULL,                  # <- our prior kinship object
  fastAssociation = FALSE
)

##find critical r2 value
r2vals<-tasMLM$FastAssociation$r2

rvals <- r2vals^(0.5)
rvals[1:25]

#######################################################################################
#USed TASSEL GUI to build LD table


LDtable<-read.csv("LD.csv")
summary(LDtable)
LDtable1<-as.data.frame(LDtable,use="complete.obs")

use="complete.obs"
#critical R2 value defined by transforming the values by taking the square root, and using the 95th percentile as threshold


LDtable[1:400,]

plot1<-ggplot(LDtable, aes(x=Dist_bp,y=R.2)) + geom_point() 


R2vals <- LDtable$R.2
Rvals <- R2vals^(0.5)
R2vals[225:460]

Rvals1 <- list()

for (x in Rvals){
  if (x == "NaN"){
    print("")
  } else Rvals1<-c(Rvals1,x)
}



#remove nas


Rvals_N<-list()

for (x in Rvals){
  if (x == 1){
    Rvals_N<-c(Rvals_N,x)
  } else print("")
}


Rthreshold <- quantile(Rvals,.95,na.rm=TRUE)
Rthreshold


