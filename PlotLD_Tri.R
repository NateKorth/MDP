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

snps1 <- genotypeToSnpMatrix(readVcf("./SNPsetsOUT1MB/MADP_BIN1_SNPs.vcf"))
SNP1<-SNPofINT$X1
snps2 <- genotypeToSnpMatrix(readVcf("./SNPsetsOUT1MB/MADP_BIN2_SNPs.vcf"))
SNP2<-SNPofINT$X2
snps3 <- genotypeToSnpMatrix(readVcf("./SNPsetsOUT1MB/MADP_BIN3_SNPs.vcf"))
SNP3<-SNPofINT$X3
snps4 <- genotypeToSnpMatrix(readVcf("./SNPsetsOUT1MB/MADP_BIN4_SNPs.vcf"))
SNP4<-SNPofINT$X4
snps5 <- genotypeToSnpMatrix(readVcf("./SNPsetsOUT1MB/MADP_BIN5_SNPs.vcf"))
SNP5<-SNPofINT$X5
snps6 <- genotypeToSnpMatrix(readVcf("./SNPsetsOUT1MB/MADP_BIN6_SNPs.vcf"))
SNP6<-SNPofINT$X6
snps7 <- genotypeToSnpMatrix(readVcf("./SNPsetsOUT1MB/MADP_BIN7_SNPs.vcf"))
SNP7<-SNPofINT$X7



#Only Significant SNPs
snps1 <- genotypeToSnpMatrix(readVcf("./SNPsOfInterestOUT/MADP_BIN1_SNPs.vcf"))
SNP1<-SNPofINT$X1
snps2 <- genotypeToSnpMatrix(readVcf("./SNPsOfInterestOUT/MADP_BIN2_SNPs.vcf"))
SNP2<-SNPofINT$X2
snps3 <- genotypeToSnpMatrix(readVcf("./SNPsOfInterestOUT/MADP_BIN3_SNPs.vcf"))
SNP3<-SNPofINT$X3
snps4 <- genotypeToSnpMatrix(readVcf("./SNPsOfInterestOUT/MADP_BIN4_SNPs.vcf"))
SNP4<-SNPofINT$X4
snps5 <- genotypeToSnpMatrix(readVcf("./SNPsOfInterestOUT/MADP_BIN5_SNPs.vcf"))
SNP5<-SNPofINT$X5
snps6 <- genotypeToSnpMatrix(readVcf("./SNPsOfInterestOUT/MADP_BIN6_SNPs.vcf"))
SNP6<-SNPofINT$X6
snps7 <- genotypeToSnpMatrix(readVcf("./SNPsOfInterestOUT/MADP_BIN7_SNPs.vcf"))
SNP7<-SNPofINT$X7

#All Significant SNPs
snps1 <- genotypeToSnpMatrix(readVcf("./ALLSNPsinCHRMsOUT/MADP_BIN1_SNPs.vcf"))
SNP1<-SNPofINT$X1
snps2 <- genotypeToSnpMatrix(readVcf("./ALLSNPsinCHRMsOUT/MADP_BIN2_SNPs.vcf"))
SNP2<-SNPofINT$X3
snps3 <- genotypeToSnpMatrix(readVcf("./ALLSNPsinCHRMsOUT/MADP_BIN3_SNPs.vcf"))
SNP3<-SNPofINT$X5
snps4 <- genotypeToSnpMatrix(readVcf("./ALLSNPsinCHRMsOUT/MADP_BIN4_SNPs.vcf"))
SNP4<-SNPofINT$X4
snps5 <- genotypeToSnpMatrix(readVcf("./ALLSNPsinCHRMsOUT/MADP_BIN5_SNPs.vcf"))
SNP5<-SNPofINT$X7
snps6 <- genotypeToSnpMatrix(readVcf("./ALLSNPsinCHRMsOUT/MADP_BIN6_SNPs.vcf"))
SNP6<-SNPofINT$X6

#Smaller Window Around Bins:
snps1 <- genotypeToSnpMatrix(readVcf("./SNPsSetV2OUT/MADP_BIN1_SNPs.vcf"))
SNP1<-SNPofINT$X1
snps2 <- genotypeToSnpMatrix(readVcf("./SNPsSetV2OUT/MADP_BIN2_SNPs.vcf"))
SNP2<-SNPofINT$X2
snps3 <- genotypeToSnpMatrix(readVcf("./SNPsSetV2OUT/MADP_BIN3_SNPs.vcf"))
SNP3<-SNPofINT$X3
snps4 <- genotypeToSnpMatrix(readVcf("./SNPsSetV2OUT/MADP_BIN4.2_SNPs_trim.vcf"))
SNP4<-SNPofINT$X4
snps5 <- genotypeToSnpMatrix(readVcf("./SNPsSetV2OUT/MADP_BIN5_SNPs.vcf"))
SNP5<-SNPofINT$X5
snps6 <- genotypeToSnpMatrix(readVcf("./SNPsSetV2OUT/MADP_BIN6.2_SNPs_trim.vcf"))
SNP6<-SNPofINT$X6
snps7 <- genotypeToSnpMatrix(readVcf("./SNPsSetV2OUT/MADP_BIN7_SNPs.vcf"))
SNP7<-SNPofINT$X7


  
#Bed<-read.bed.matrix("MADP_maf_0.05", bed = paste(basename, ".bed", sep=""),  fam = paste(basename, ".fam", sep=""), bim = paste(basename, ".bim", sep=""), verbose = getOption("gaston.verbose",TRUE))

#Bed<-read.bed.matrix("MADP_maf_0.05", verbose = getOption("gaston.verbose",TRUE))

rgb.palette <- colorRampPalette(rev(c("white", "firebrick3")), space = "rgb")

map<-data.table::fread("MADP_V6.geno.map")

pos<-as.data.frame(snps1$map$snp.names)
colnames(pos)<-"SNP"
pos<-merge(pos,map,by="SNP",all.x=T)
pos<-pos[order(pos[3]),]
LDheatmap(snps1$genotypes,flip = T, genetic.distances = pos$POS,distances="physical",SNP.name = SNP1,title="CHR1",color = rgb.palette(4))

pos<-as.data.frame(snps2$map$snp.names)
colnames(pos)<-"SNP"
pos<-merge(pos,map,by="SNP",all.x=T)
pos<-pos[order(pos[3]),]
LDheatmap(snps2$genotypes,flip = T, genetic.distances = pos$POS,distances="physical",SNP.name = SNP2,title="CHR3",color = rgb.palette(4))

pos<-as.data.frame(snps3$map$snp.names)
colnames(pos)<-"SNP"
pos<-merge(pos,map,by="SNP",all.x=T)
pos<-pos[order(pos[3]),]
LDheatmap(snps3$genotypes,flip = T, genetic.distances = pos$POS,distances="physical",SNP.name = SNP3,title="CHR5",color = rgb.palette(4))

pos<-as.data.frame(snps4$map$snp.names)
colnames(pos)<-"SNP"
pos<-merge(pos,map,by="SNP",all.x=T)
pos<-pos[order(pos[3]),]
LDheatmap(snps4$genotypes,flip = T, genetic.distances = pos$POS,distances="physical",SNP.name = SNP4,title="CHR6",color = rgb.palette(4))

pos<-as.data.frame(snps5$map$snp.names)
colnames(pos)<-"SNP"
pos<-merge(pos,map,by="SNP",all.x=T)
pos<-pos[order(pos[3]),]
LDheatmap(snps5$genotypes,flip = T, genetic.distances = pos$POS,distances="physical",SNP.name = SNP5,title="CHR7a",color = rgb.palette(4))

pos<-as.data.frame(snps6$map$snp.names)
colnames(pos)<-"SNP"
pos<-merge(pos,map,by="SNP",all.x=T)
pos<-pos[order(pos[3]),]
LDheatmap(snps6$genotypes,flip = T, genetic.distances = pos$POS,distances="physical",SNP.name = SNP6,title="CHR7b",LDmeasure="r",color = rgb.palette(4))

pos<-as.data.frame(snps7$map$snp.names)
colnames(pos)<-"SNP"
pos<-merge(pos,map,by="SNP",all.x=T)
pos<-pos[order(pos[3]),]
LDheatmap(snps7$genotypes,flip = T, genetic.distances = pos$POS,distances="physical",SNP.name = SNP7,title="CHR8",color = rgb.palette(4))

#pos<-as.data.frame(snps8$map$snp.names)
#colnames(pos)<-"SNP"
#pos<-merge(pos,map,by="SNP",all.x=T)
#LDheatmap(snps8$genotypes,flip = T, genetic.distances = pos$POS,distances="physical",SNP.name = pos$POS,title="CHR8: 20.8-22.3 MB",color = rgb.palette(200))

#pos<-as.data.frame(snps9$map$snp.names)
#colnames(pos)<-"SNP"
#pos<-merge(pos,map,by="SNP",all.x=T)
#LDheatmap(snps9$genotypes,flip = T, genetic.distances = pos$POS,distances="physical",SNP.name = pos$POS,title="CHR8: 52.5-54.5 MB",color = rgb.palette(200))

#pos<-as.data.frame(snps10$map$snp.names)
#colnames(pos)<-"SNP"
#pos<-merge(pos,map,by="SNP",all.x=T)
#LDheatmap(snps10$genotypes,flip = T, genetic.distances = pos$POS,distances="physical",SNP.name = pos$POS,title="CHR10: 52.5-54.5 MB",color = rgb.palette(200))

#pos<-as.data.frame(snps11$map$snp.names)
#colnames(pos)<-"SNP"
#pos<-merge(pos,map,by="SNP",all.x=T)
#LDheatmap(snps11$genotypes,flip = T, genetic.distances = pos$POS,distances="physical",SNP.name = pos$POS,title="CHR11: 50.6-52.2 MB",color = rgb.palette(200))


	
	
	
	
	
	


#######
#get SNPs in each location of interest
map1<-subset(map,CHROM=="1")
map1<-subset(map1,POS<=42036457)
map1<-subset(map1,POS>=41025998)
map1<-map1[,1]
write.table(map1,"MADP_BIN1_SNPs.csv", row.names = F,quote = F,col.names = F)

map2<-subset(map,CHROM=="3")
map2<-subset(map2,POS<=52306019)
map2<-subset(map2,POS>=51484051)
map2<-map2[,1]
write.table(map2,"MADP_BIN2_SNPs.csv", row.names = F,quote = F,col.names = F)

map3<-subset(map,CHROM=="5")
map3<-subset(map3,POS<=9615026)
map3<-subset(map3,POS>=9206068)
map3<-map3[,1]
write.table(map3,"MADP_BIN3_SNPs.csv", row.names = F,quote = F,col.names = F)

map4<-subset(map,CHROM=="6")
map4<-subset(map4,POS<=3700000)
map4<-subset(map4,POS>=0)
map4<-map4[,1]
write.table(map4,"MADP_BIN4.2_SNPs.csv", row.names = F,quote = F,col.names = F)

map5<-subset(map,CHROM=="7")
map5<-subset(map5,POS<=9661392)
map5<-subset(map5,POS>=9484401)
map5<-map5[,1]
write.table(map5,"MADP_BIN5_SNPs.csv", row.names = F,quote = F,col.names = F)

map6<-subset(map,CHROM=="7")
map6<-subset(map6,POS<=24500000)
map6<-subset(map6,POS>=12800000)
map6<-map6[,1]
write.table(map6,"MADP_BIN6.2_SNPs.csv", row.names = F,quote = F,col.names = F)

map7<-subset(map,CHROM=="8")
map7<-subset(map7,POS<=3838416)
map7<-subset(map7,POS>=3777523)
map7<-map7[,1]
write.table(map7,"MADP_BIN7_SNPs.csv", row.names = F,quote = F,col.names = F)


#Pull just SNPs of intrest
#write.table(SNP1,"MADP_BIN1_SNPs.csv", row.names = F,quote = F,col.names = F)
#write.table(SNP2,"MADP_BIN2_SNPs.csv", row.names = F,quote = F,col.names = F)
#write.table(SNP3,"MADP_BIN3_SNPs.csv", row.names = F,quote = F,col.names = F)
#write.table(SNP4,"MADP_BIN4_SNPs.csv", row.names = F,quote = F,col.names = F)
#write.table(SNP5,"MADP_BIN5_SNPs.csv", row.names = F,quote = F,col.names = F)
#write.table(SNP6,"MADP_BIN6_SNPs.csv", row.names = F,quote = F,col.names = F)
#write.table(SNP7,"MADP_BIN7_SNPs.csv", row.names = F,quote = F,col.names = F)



