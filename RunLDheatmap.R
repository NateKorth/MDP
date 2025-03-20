##Library packages (install.packages() if you need to)
library(LDheatmap)
library(VariantAnnotation)
library(viridis)
library(gaston)
#!!!!!Make sure Beanline column has identical identifiers to SNP file
#Set working directory to location of OTU table and SNP file
setwd("~/MADP/LD")

#Get genomic information from beans:
#VCF<-readVcf("../Genotypes/MADP_geno_imputed.maf_05.vcf")
#VCF$

#in bins
SNPofINT<-read.csv("./SNPsOInterest1MB.1.csv")
#all chromosomes
#SNPofINT<-read.csv("./ALLSigSNPsInChrm.csv")

#Import a vcf file that contains just the SNPs you want to plot:
#all snps in a region:
snps1 <- genotypeToSnpMatrix(readVcf("./SNPsetsOUT1MB/MADP_BIN1_SNPs.vcf"))
SNP1<-SNPofINT$X1


#Only Significant SNPs:
snps1 <- genotypeToSnpMatrix(readVcf("./SNPsOfInterestOUT/MADP_BIN1_SNPs.vcf"))
SNP1<-SNPofINT$X1


#All Significant SNPs
snps1 <- genotypeToSnpMatrix(readVcf("./ALLSNPsinCHRMsOUT/MADP_BIN1_SNPs.vcf"))
SNP1<-SNPofINT$X1

#Smaller Window Around Bins:
snps1 <- genotypeToSnpMatrix(readVcf("./SNPsSetV2OUT/MADP_BIN1_SNPs.vcf"))
SNP1<-SNPofINT$X1

#define a color palette:
rgb.palette <- colorRampPalette(rev(c("white", "firebrick3")), space = "rgb")

#read in a map file
map<-data.table::fread("MADP_V6.geno.map")

pos<-as.data.frame(snps1$map$snp.names)
colnames(pos)<-"SNP"
pos<-merge(pos,map,by="SNP",all.x=T)
pos<-pos[order(pos[3]),]

#Generate heatmap:
LDheatmap(snps1$genotypes,flip = T, genetic.distances = pos$POS,distances="physical",SNP.name = SNP1,title="CHR1",color = rgb.palette(4))

#how to pull all snps from a specific region: 
#######
#get SNPs in each location of interest
map1<-subset(map,CHROM=="1")
map1<-subset(map1,POS<=42036457)
map1<-subset(map1,POS>=41025998)
map1<-map1[,1]
write.table(map1,"MADP_BIN1_SNPs.csv", row.names = F,quote = F,col.names = F)

#After you have this table you can use it to subset a vcf using vcf tools (there's probably a way to do it in R but I don't know it off the top of my head)

