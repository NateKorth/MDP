##Library packages (install.packages() if you need to)
library(sommer)
library(ggplot2)
library(readxl)
library(data.table)
library(tidyverse)
library(biomformat)
library(metagenomeSeq)
library(phyloseq)
library(ape)
library(reshape2)
library(vegan)
library(superheat)


#!!!!!Make sure Beanline column has identical identifiers to SNP file
#Set working directory to location of OTU table and SNP file
setwd("~/MADP")

##insert OTU table here:
df <- read_excel("Analysis/770/MADPS770Samples_RarefiedGenusAbsoluteFilterTable.xlsx",sheet="Abs_Counts_filter")



#Subset data by subject
df770 <- as.data.frame(subset(df, Subject=="S770"))
df770 <- subset(df770, ControlType=="Sample")
df770 <- subset(df770, SampleLost=="No")



#remove Baseline and Controls:377:382
table.alpha<-df770[-c(2:23)]
table.alpha0<-as.data.frame(t(table.alpha))

#remove ID column and add it as column names
A.OTU<-table.alpha0[-c(1)]
OTU_IDs<-as.data.frame(table.alpha0[1])
rownames(A.OTU)<-OTU_IDs[,1]

A.OTU[1:3,358:370]

###################################################################
#read in Silva taxonomy | before, in excel seperate taxonomy into columns
taxa<-read.delim("./exported/UFMu_AOX_1-8_silva_taxonomy_FULL.tsv.txt",stringsAsFactors = FALSE,sep="\t")

OTUdata<-AnnotatedDataFrame(taxa)
rownames(OTUdata)<-taxa$Feature.ID

#normalize data
#load metadata
Meta<-load_phenoData("UFMu_AOX_1-8_metadata2.tsv.txt",sep="\t")
#name the columns in Meta data file
colnames(Meta)<-c("Line","Plate#","Rep","Type","Subject","Mutant","DialysisBatch","LineM","LineA")

#subset metadata by sample type and subject 
Ameta<-subset(Meta, Type=="WholeSeed"|Type=="Pericarp" | Type=="FBB16" | Type=="FBB0")


#Match names in meta data to OTU table
ordA<- match(colnames(A.OTU), rownames(Ameta))
Ameta <- Ameta[ordA,]


#insert subset into phenodataframe:
phenotypeDataA<-AnnotatedDataFrame(Ameta)

#Make MetagenomeSeq object
A.data<-newMRexperiment(A.OTU, phenoData = phenotypeDataA, featureData = OTUdata)

#trim that shit yo: taxa only present in 20% of taxa removed
A.datatrim<-filterData(A.data,present = 143)

#try skipping this step?
#normalize by CSS:
pA<-cumNormStatFast(A.datatrim)
A.datatrim<-cumNorm(A.datatrim,p=pA)

#log=TRUE -log2 transform
UFMU_A_normal_log<-as.data.frame(MRcounts(A.datatrim,norm=TRUE,log=TRUE))

#trimmed relative abundance OTU table:
A.relabun.tr<-as.data.frame(MRcounts(A.datatrim,norm=TRUE,log=FALSE))
A.relabun.trim<-sweep(A.relabun.tr,2,colSums(A.relabun.tr),"/")

#get features present in subsets
AIDs<-as.data.frame(rownames(UFMU_A_normal_log))
names(AIDs)<-"Feature.ID"

#subset taxonomy table to only include taxa in OTU table:
A_taxa<-taxa[,c("Feature.ID","Family","Genus")]
A_taxa<-A_taxa[which(A_taxa$Feature.ID %in% AIDs$Feature.ID),]

#define varibles of model: Line, Type, Subject
ALine<-pData(A.datatrim)$Line
ALineM<-pData(A.datatrim)$LineM
ASubject<-pData(A.datatrim)$Subject
AType<-pData(A.datatrim)$Type

#define normalization factor
Anorm.factor <- normFactors(A.datatrim)
Anorm.factor <- log2(Anorm.factor/median(Anorm.factor)+1)

#define model
Amod<- model.matrix(~ALine+ASubject+Anorm.factor)
AmodM<- model.matrix(~ALineM+ASubject+Anorm.factor)
AmodT<- model.matrix(~AType+ASubject+Anorm.factor)

#settings can be same for all samples:
settings<-zigControl(maxit=10, verbose=TRUE)

#fit model:
Afit<-fitZig(obj=A.datatrim,mod=Amod,useCSSoffset =TRUE,control=settings)
AfitM<-fitZig(obj=A.datatrim,mod=AmodM,useCSSoffset =TRUE,control=settings)
AfitT<-fitZig(obj=A.datatrim,mod=AmodT,useCSSoffset =TRUE,control=settings)

Acoefs<-MRcoefs(Afit,group=2,number=200)
AcoefsM<-MRcoefs(AfitM,group=2,number=200)
AcoefsT<-MRcoefs(AfitT,group=2,number=200)

#filter out non-significant OTUs
A.sig.otus<-Acoefs[-which(Acoefs$pvalues>=.1),]
A.sig.otus.M<-AcoefsM[-which(AcoefsM$pvalues>=.1),]
A.sig.otus.T<-AcoefsT[-which(AcoefsM$pvalues>=.2),]

#make table of significant OTUs
A.sig.otus1<-as.data.frame(rownames(A.sig.otus.M))

names(A.sig.otus1)<-"OTU"

#create new row with OTU names to compare and filter based on significance
UFMU_A_normal_log$OTU<-rownames(UFMU_A_normal_log)

#subset all taxa tables based on which OTUs were significant:
A.sig.table <- UFMU_A_normal_log[which(UFMU_A_normal_log$OTU %in% A.sig.otus1$OTU),]

# make list of OTUs 1-3005
OTUs<-as.data.frame(paste("ASV","",seq(1,3005, by=1)))
taxa$taxa<-OTUs
names(taxa$taxa)<-"taxa"

#make list of significant taxa genera:
A.sig.taxa<-taxa[,c("taxa","Genus","Feature.ID")]
A.sig.taxa<-A.sig.taxa[which(A.sig.taxa$Feature.ID %in% A.sig.otus1$OTU),]

#make list of significant IDs to make labels for plots
A.sig.IDS<-as.data.frame(A.sig.taxa)

#Significant taxa relative abundance table:
A.relabun.trim$OTU<-rownames(A.relabun.trim)

A.sig.table<-UFMU_A_normal_log[which(UFMU_A_normal_log$OTU %in% A.sig.IDS$Feature.ID),]
A.sig.table.RA<-A.relabun.trim[which(A.relabun.trim$OTU %in% A.sig.IDS$Feature.ID),]
#
names(A.sig.IDS)<-c("taxa","Genus")
AtaxaIDs<-A.sig.taxa$taxa
A.sig.IDS$IDS<-paste(AtaxaIDs$taxa,"-",A.sig.taxa$Genus)
rownames(A.sig.table)<-paste(A.sig.IDS$IDS)
rownames(A.sig.table.RA)<-paste(A.sig.IDS$IDS)


##Plot specific taxa by ID
A.1.RA<-A.sig.table.RA
A.1.RA$OTU<-NULL
A.1.RA2<-as.data.frame(t(A.1.RA))
ALine1<-as.data.frame(ALine)
A.1.RA2$Line<-ALine1$ALine
ALineM1<-as.data.frame(ALineM)
A.1.RA2$MutantType<-ALineM1$ALineM
AType1<-as.data.frame(AType)
A.1.RA2$Type<-AType1$AType
ASubject1<-as.data.frame(ASubject)
A.1.RA2$Subject<-ASubject1$ASubject
A.1.RA2$Line<-factor(A.1.RA2$Line, levels= c("Md0","Md16","D3","M1","W22","B2","C2","M2","N1","A5","D4","E4","F2","H2","U5","W3","X1"))
A.1.RA2$MutantType<-factor(A.1.RA2$MutantType, levels= c("Md0","Md16","WTC","PAL","Methylglyoxal_Detox","Xylan_Backbone","BAHD","PRX"))

A.1.<-A.sig.table
A.1.$OTU<-NULL
A.1.2<-as.data.frame(t(A.1.))
ALine1<-as.data.frame(ALine)
A.1.2$Line<-ALine1$ALine
ALineM1<-as.data.frame(ALineM)
A.1.2$MutantType<-ALineM1$ALineM
AType1<-as.data.frame(AType)
A.1.2$Type<-AType1$AType
ASubject1<-as.data.frame(ASubject)
A.1.2$Subject<-ASubject1$ASubject
A.1.2$Line<-factor(A.1.2$Line, levels= c("Md0","Md16","D3","M1","W22","B2","C2","M2","N1","A5","D4","E4","F2","H2","U5","W3","X1"))
A.1.2$MutantType<-factor(A.1.2$MutantType, levels= c("Md0","Md16","WTC","PAL","Methylglyoxal_Detox","Xylan_Backbone","BAHD","PRX"))

#
#Subset based on Sample Type:
A.1.2W<-subset(A.1.2, Type=="WholeSeed" | Type=="FBB16" | Type=="FBB0")
A.1.2P<-subset(A.1.2, Type=="Pericarp" | Type=="FBB16" | Type=="FBB0")
A.1.RA2W<-subset(A.1.RA2, Type=="WholeSeed" | Type=="FBB16" | Type=="FBB0")
A.1.RA2P<-subset(A.1.RA2, Type=="Pericarp" | Type=="FBB16" | Type=="FBB0")
#All samples Together###########################################################################################################
#
ggplot(data=A.1.RA2, aes(x=Line,y=`ASV  1 - Bacteroides`,colour=Subject,shape=Type))+
  geom_jitter(width=0.07)+scale_colour_manual(values = c("Tan","burlywood4","Brown4","Orange","darkgreen","cornflowerblue","darkorchid","Red"))+
  ylab("Relative Abundance")+theme_classic()+theme(text=element_text(size=30,family = "serif"))+ggtitle("ASV 1 - Bacteroides")##+theme(legend.position ="none")
#




























##############Read in SNP file in .hmp format and convert to numeric_row format for h2 
#First transpose datafile
GD <- read.delim("Genotypes/469-geno-200k-SNPs-imputed.hmp.txt")

GD[1:5,1:5]

GD1 <- GD[-c(2:11)]
GD2 <- as.data.frame(as.matrix(t(GD1)))
colnames(GD2) <- rownames(GD1)
rownames(GD2) <- colnames(GD1)
##Check data

GD2.1<-as.data.frame(rownames(GD2))
GD2$BeanLine<-GD2.1$`rownames(GD2)`

GD2.2[1:5,1:5]
GD2.2<-GD2[which(GD2$BeanLine %in% df770$BeanLine ),]
GD2.2$BeanLine<-NULL

#convert to numeric
GD3 <- atcg1234(GD2.2, ploidy=2)
GD3$M[1:5,1:5]
#remove extra line at top, make element of list into dataframe
GD5 <- GD3$M
GD5[1:200,1:2]



GD6<-as.data.frame(rownames(GD5))
###################Make covariate matrices A=additive matrix. D=Dominance matrix. E=Epistatic matrix 
A <- A.mat(GD5)
D <- D.mat(GD5)
E <- E.mat(GD5)

#Function to generate Broad sense heritability + plot it 
####Matrix will only work in ids are the same in the matrix file and OTU table!!!
##Order does not have to be the same
#accepts replicates
#Broad sense heritability

################Make LMM Model with additive and dominance matrix add identical column to bean line
###MAKE SURE BeanLine identifiers are IDENTICAL in matrix and df 
df770.log$BeanLine<-df770$BeanLine
df770.log$idd<-df770$BeanLine

#import PCs
PCs<-read.delim("Analysis/770/PCs/MADP_geno_imputed_10PCA1.txt.gapit")

#subset lines in otu table that match line names in GD
df770.1<-df770[which(df770$BeanLine %in% GD6$`rownames(GD5)`),]
#identify line names that don't match:
df770.2<-df770[which(!df770$BeanLine %in% GD6$`rownames(GD5)`),]
#make sure IDs match in PCs
df770.2<-df770[which(!df770$BeanLine %in% PCs$BeanLine),]
df770.2<-PCs[which(!PCs$BeanLine %in% df770$BeanLine),]
#Get IDs:
df770vals<-df770[-c(1:4,6:23)]
df770meds<-df770vals %>% group_by(BeanLine) %>% summarise_each(funs(median))

#remove outliers and make df for single taxa:
#log 2 transform
#range value determines how many standard deviations from mean to remove, we're starting with 3

#Prevotella9
df_Prevotella9<-as.data.frame(df770$Prevotella9)
df_Prevotella9<-log(df_Prevotella9,2)
df_Prevotella9$BeanLine<-df770$BeanLine
df_Prevotella9$idd<-df770$BeanLine
df_Prevotella9$Batch<-as.character(df770$DialysisBatch)
df_Prevotella9$Plate<-as.character(df770$Plate)
df_Prevotella9$Race<-(df770$BeanRace)
names(df_Prevotella9)<-c("Prevotella9","BeanLine","idd","Batch","Plate","Race")
Prevotella9_out<-boxplot(df_Prevotella9$Prevotella9,range=3)$out
df_Prevotella9<-df_Prevotella9[-which(df_Prevotella9$Prevotella9 %in% Prevotella9_out),]
hist(df_Prevotella9$Prevotella9)
#calc_heritability
fit <- mmer(Prevotella9~1, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D), rcov=~units, data=df_Prevotella9)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[3,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_list <- data.frame("Prevotella9", h2 )
names(h2_list)<-c("Genus","h2")
#outputBLUPs for GWAS:
fit1 <- mmer(Prevotella9~1, random=~BeanLine+vs(BeanLine,Gu=A) + vs(idd,Gu=D)+Batch+Plate+Race, rcov=~units, data=df_Prevotella9,tolparinv = 7)
#summary(fit1)
Prevotella9BLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(Prevotella9BLUPS)<-df770meds$BeanLine
#export csv for mapping:
write.csv(Prevotella9BLUPS,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Prevotella9_BLUPS",quote = FALSE,row.names=TRUE)


###Enterobacteriaceae
df_Enterobacteriaceae<-as.data.frame(df770$Enterobacteriaceae)
df_Enterobacteriaceae<-log(df_Enterobacteriaceae,2)
df_Enterobacteriaceae$BeanLine<-df770$BeanLine
df_Enterobacteriaceae$idd<-df770$BeanLine
df_Enterobacteriaceae$Batch<-as.character(df770$DialysisBatch)
df_Enterobacteriaceae$Plate<-as.character(df770$Plate)
df_Enterobacteriaceae$Race<-(df770$BeanRace)
names(df_Enterobacteriaceae)<-c("Enterobacteriaceae","BeanLine","idd","Batch","Plate","Race")
Enterobacteriaceae_out<-boxplot(df_Enterobacteriaceae$Enterobacteriaceae,range=3)$out
df_Enterobacteriaceae<-df_Enterobacteriaceae[-which(df_Enterobacteriaceae$Enterobacteriaceae %in% Enterobacteriaceae_out),]
hist(df_Enterobacteriaceae$Enterobacteriaceae)
#calc_heritability
fit <- mmer(Enterobacteriaceae~1, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D), rcov=~units, data=df_Enterobacteriaceae)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[3,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_list <- data.frame("Enterobacteriaceae", h2 )
names(h2_list)<-c("Genus","h2")
#outputBLUPs for GWAS:
fit1 <- mmer(Enterobacteriaceae~1, random=~BeanLine+vs(BeanLine,Gu=A) + vs(idd,Gu=D)+Batch+Plate+Race, rcov=~units, data=df_Enterobacteriaceae,tolparinv = 7)
#summary(fit1)
EnterobacteriaceaeBLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(EnterobacteriaceaeBLUPS)<-df770meds$BeanLine
#export csv for mapping:
write.csv(EnterobacteriaceaeBLUPS,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Enterobacteriaceae_BLUPS",quote = FALSE,row.names=TRUE)

###Succinivibrio
df_Succinivibrio<-as.data.frame(df770$Succinivibrio)
df_Succinivibrio<-log(df_Succinivibrio,2)
df_Succinivibrio$BeanLine<-df770$BeanLine
df_Succinivibrio$idd<-df770$BeanLine
df_Succinivibrio$Batch<-as.character(df770$DialysisBatch)
df_Succinivibrio$Plate<-as.character(df770$Plate)
df_Succinivibrio$Race<-(df770$BeanRace)
names(df_Succinivibrio)<-c("Succinivibrio","BeanLine","idd","Batch","Plate","Race")
Succinivibrio_out<-boxplot(df_Succinivibrio$Succinivibrio,range=3)$out
df_Succinivibrio<-df_Succinivibrio[-which(df_Succinivibrio$Succinivibrio %in% Succinivibrio_out),]
hist(df_Succinivibrio$Succinivibrio)
#calc_heritability
fit <- mmer(Succinivibrio~1, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D), rcov=~units, data=df_Succinivibrio)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[3,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_list <- data.frame("Succinivibrio", h2 )
names(h2_list)<-c("Genus","h2")
#outputBLUPs for GWAS:
fit1 <- mmer(Succinivibrio~1, random=~BeanLine+vs(BeanLine,Gu=A) + vs(idd,Gu=D)+Batch+Plate+Race, rcov=~units, data=df_Succinivibrio,tolparinv = 7)
#summary(fit1)
SuccinivibrioBLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(SuccinivibrioBLUPS)<-df770meds$BeanLine
#export csv for mapping:
write.csv(SuccinivibrioBLUPS,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Succinivibrio_BLUPS",quote = FALSE,row.names=TRUE)

###Dialister
df_Dialister<-as.data.frame(df770$Dialister)
df_Dialister<-log(df_Dialister,2)
df_Dialister$BeanLine<-df770$BeanLine
df_Dialister$idd<-df770$BeanLine
df_Dialister$Batch<-as.character(df770$DialysisBatch)
df_Dialister$Plate<-as.character(df770$Plate)
df_Dialister$Race<-(df770$BeanRace)
names(df_Dialister)<-c("Dialister","BeanLine","idd","Batch","Plate","Race")
Dialister_out<-boxplot(df_Dialister$Dialister,range=3)$out
df_Dialister<-df_Dialister[-which(df_Dialister$Dialister %in% Dialister_out),]
hist(df_Dialister$Dialister)
#calc_heritability
fit <- mmer(Dialister~1, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D), rcov=~units, data=df_Dialister)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[3,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_list <- data.frame("Dialister", h2 )
names(h2_list)<-c("Genus","h2")
#outputBLUPs for GWAS:
fit1 <- mmer(Dialister~1, random=~BeanLine+vs(BeanLine,Gu=A) + vs(idd,Gu=D)+Batch+Plate+Race, rcov=~units, data=df_Dialister,tolparinv = 7)
#summary(fit1)
DialisterBLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(DialisterBLUPS)<-df770meds$BeanLine
#export csv for mapping:
write.csv(DialisterBLUPS,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Dialister_BLUPS",quote = FALSE,row.names=TRUE)

###Bacteroides
df_Bacteroides<-as.data.frame(df770$Bacteroides)
df_Bacteroides<-log(df_Bacteroides,2)
df_Bacteroides$BeanLine<-df770$BeanLine
df_Bacteroides$idd<-df770$BeanLine
df_Bacteroides$Batch<-as.character(df770$DialysisBatch)
df_Bacteroides$Plate<-as.character(df770$Plate)
df_Bacteroides$Race<-(df770$BeanRace)
names(df_Bacteroides)<-c("Bacteroides","BeanLine","idd","Batch","Plate","Race")
Bacteroides_out<-boxplot(df_Bacteroides$Bacteroides,range=3)$out
df_Bacteroides<-df_Bacteroides[-which(df_Bacteroides$Bacteroides %in% Bacteroides_out),]
hist(df_Bacteroides$Bacteroides)
#calc_heritability
fit <- mmer(Bacteroides~1, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D), rcov=~units, data=df_Bacteroides)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[3,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_list <- data.frame("Bacteroides", h2 )
names(h2_list)<-c("Genus","h2")
#outputBLUPs for GWAS:
fit1 <- mmer(Bacteroides~1, random=~BeanLine+vs(BeanLine,Gu=A) + vs(idd,Gu=D)+Batch+Plate+Race, rcov=~units, data=df_Bacteroides,tolparinv = 7)
#summary(fit1)
BacteroidesBLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(BacteroidesBLUPS)<-df770meds$BeanLine
#export csv for mapping:
write.csv(BacteroidesBLUPS,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Bacteroides_BLUPS",quote = FALSE,row.names=TRUE)

###Coprococcus3
df_Coprococcus3<-as.data.frame(df770$Coprococcus3)
df_Coprococcus3<-log(df_Coprococcus3,2)
df_Coprococcus3$BeanLine<-df770$BeanLine
df_Coprococcus3$idd<-df770$BeanLine
df_Coprococcus3$Batch<-as.character(df770$DialysisBatch)
df_Coprococcus3$Plate<-as.character(df770$Plate)
df_Coprococcus3$Race<-(df770$BeanRace)
names(df_Coprococcus3)<-c("Coprococcus3","BeanLine","idd","Batch","Plate","Race")
Coprococcus3_out<-boxplot(df_Coprococcus3$Coprococcus3,range=3)$out
df_Coprococcus3<-df_Coprococcus3[-which(df_Coprococcus3$Coprococcus3 %in% Coprococcus3_out),]
hist(df_Coprococcus3$Coprococcus3)
#calc_heritability
fit <- mmer(Coprococcus3~1, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D), rcov=~units, data=df_Coprococcus3)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[3,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_list <- data.frame("Coprococcus3", h2 )
names(h2_list)<-c("Genus","h2")
#outputBLUPs for GWAS:
fit1 <- mmer(Coprococcus3~1, random=~BeanLine+vs(BeanLine,Gu=A) + vs(idd,Gu=D)+Batch+Plate+Race, rcov=~units, data=df_Coprococcus3,tolparinv = 7)
#summary(fit1)
Coprococcus3BLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(Coprococcus3BLUPS)<-df770meds$BeanLine
#export csv for mapping:
write.csv(Coprococcus3BLUPS,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Coprococcus3_BLUPS",quote = FALSE,row.names=TRUE)

###Anaerostipes
df_Anaerostipes<-as.data.frame(df770$Anaerostipes)
df_Anaerostipes<-log(df_Anaerostipes,2)
df_Anaerostipes$BeanLine<-df770$BeanLine
df_Anaerostipes$idd<-df770$BeanLine
df_Anaerostipes$Batch<-as.character(df770$DialysisBatch)
df_Anaerostipes$Plate<-as.character(df770$Plate)
df_Anaerostipes$Race<-(df770$BeanRace)
names(df_Anaerostipes)<-c("Anaerostipes","BeanLine","idd","Batch","Plate","Race")
Anaerostipes_out<-boxplot(df_Anaerostipes$Anaerostipes,range=3)$out
df_Anaerostipes<-df_Anaerostipes[-which(df_Anaerostipes$Anaerostipes %in% Anaerostipes_out),]
hist(df_Anaerostipes$Anaerostipes)
#calc_heritability
fit <- mmer(Anaerostipes~1, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D), rcov=~units, data=df_Anaerostipes)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[3,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_list <- data.frame("Anaerostipes", h2 )
names(h2_list)<-c("Genus","h2")
#outputBLUPs for GWAS:
fit1 <- mmer(Anaerostipes~1, random=~BeanLine+vs(BeanLine,Gu=A) + vs(idd,Gu=D)+Batch+Plate+Race, rcov=~units, data=df_Anaerostipes,tolparinv = 7)
#summary(fit1)
AnaerostipesBLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(AnaerostipesBLUPS)<-df770meds$BeanLine
#export csv for mapping:
write.csv(AnaerostipesBLUPS,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Anaerostipes_BLUPS",quote = FALSE,row.names=TRUE)

###Sutterella
df_Sutterella<-as.data.frame(df770$Sutterella)
df_Sutterella<-log(df_Sutterella,2)
df_Sutterella$BeanLine<-df770$BeanLine
df_Sutterella$idd<-df770$BeanLine
df_Sutterella$Batch<-as.character(df770$DialysisBatch)
df_Sutterella$Plate<-as.character(df770$Plate)
df_Sutterella$Race<-(df770$BeanRace)
names(df_Sutterella)<-c("Sutterella","BeanLine","idd","Batch","Plate","Race")
Sutterella_out<-boxplot(df_Sutterella$Sutterella,range=3)$out
df_Sutterella<-df_Sutterella[-which(df_Sutterella$Sutterella %in% Sutterella_out),]
hist(df_Sutterella$Sutterella)
#calc_heritability
fit <- mmer(Sutterella~1, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D), rcov=~units, data=df_Sutterella)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[3,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_list <- data.frame("Sutterella", h2 )
names(h2_list)<-c("Genus","h2")
#outputBLUPs for GWAS:
fit1 <- mmer(Sutterella~1, random=~BeanLine+vs(BeanLine,Gu=A) + vs(idd,Gu=D)+Batch+Plate+Race, rcov=~units, data=df_Sutterella,tolparinv = 7)
#summary(fit1)
SutterellaBLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(SutterellaBLUPS)<-df770meds$BeanLine
#export csv for mapping:
write.csv(SutterellaBLUPS,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Sutterella_BLUPS",quote = FALSE,row.names=TRUE)

###Lachnospiraceae
df_Lachnospiraceae<-as.data.frame(df770$Lachnospiraceae)
df_Lachnospiraceae<-log(df_Lachnospiraceae,2)
df_Lachnospiraceae$BeanLine<-df770$BeanLine
df_Lachnospiraceae$idd<-df770$BeanLine
names(df_Lachnospiraceae)<-c("Lachnospiraceae","BeanLine","idd")
Lachnospiraceae_out<-boxplot(df_Lachnospiraceae$Lachnospiraceae,range=2.5)$out
df_Lachnospiraceae<-df_Lachnospiraceae[-which(df_Lachnospiraceae$Lachnospiraceae %in% Lachnospiraceae_out),]
hist(df_Lachnospiraceae$Lachnospiraceae)
#calc_heritability
fit <- mmer(Lachnospiraceae~1, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D), rcov=~units, data=df_Lachnospiraceae)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[3,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_df <- as.data.frame(cbind("Lachnospiraceae",h2))
names(h2_df)<-c("Genus","h2")
h2_list<-rbind(h2_list,data.frame(h2_df))
#outputBLUPs for GWAS:
fit1 <- mmer(Lachnospiraceae~1, random=~BeanLine, rcov=~units, data=df_Lachnospiraceae)
LachnospiraceaeBLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(LachnospiraceaeBLUPS)<-df770meds$BeanLine
#export csv for mapping:
write.csv(LachnospiraceaeBLUPS,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Lachnospiraceae_BLUPS",quote = FALSE,row.names=TRUE)

###Prevotellaceae_uncultured
df_Prevotellaceae_uncultured<-as.data.frame(df770$Prevotellaceae_uncultured)
df_Prevotellaceae_uncultured<-log(df_Prevotellaceae_uncultured,2)
df_Prevotellaceae_uncultured$BeanLine<-df770$BeanLine
df_Prevotellaceae_uncultured$idd<-df770$BeanLine
names(df_Prevotellaceae_uncultured)<-c("Prevotellaceae_uncultured","BeanLine","idd")
Prevotellaceae_uncultured_out<-boxplot(df_Prevotellaceae_uncultured$Prevotellaceae_uncultured,range=2.5)$out
df_Prevotellaceae_uncultured<-df_Prevotellaceae_uncultured[-which(df_Prevotellaceae_uncultured$Prevotellaceae_uncultured %in% Prevotellaceae_uncultured_out),]
hist(df_Prevotellaceae_uncultured$Prevotellaceae_uncultured)
#calc_heritability
fit <- mmer(Prevotellaceae_uncultured~1, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D), rcov=~units, data=df_Prevotellaceae_uncultured)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[3,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_df <- as.data.frame(cbind("Prevotellaceae_uncultured",h2))
names(h2_df)<-c("Genus","h2")
h2_list<-rbind(h2_list,data.frame(h2_df))
#outputBLUPs for GWAS:
fit1 <- mmer(Prevotellaceae_uncultured~1, random=~BeanLine, rcov=~units, data=df_Prevotellaceae_uncultured)
Prevotellaceae_unculturedBLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(Prevotellaceae_unculturedBLUPS)<-df770meds$BeanLine
#export csv for mapping:
write.csv(Prevotellaceae_unculturedBLUPS,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Prevotellaceae_uncultured_BLUPS",quote = FALSE,row.names=TRUE)


###Coprococcus1
df_Coprococcus1<-as.data.frame(df770$Coprococcus1)
df_Coprococcus1<-log(df_Coprococcus1,2)
df_Coprococcus1$BeanLine<-df770$BeanLine
df_Coprococcus1$idd<-df770$BeanLine
names(df_Coprococcus1)<-c("Coprococcus1","BeanLine","idd")
Coprococcus1_out<-boxplot(df_Coprococcus1$Coprococcus1,range=2.5)$out
df_Coprococcus1<-df_Coprococcus1[-which(df_Coprococcus1$Coprococcus1 %in% Coprococcus1_out),]
hist(df_Coprococcus1$Coprococcus1)
#calc_heritability
fit <- mmer(Coprococcus1~1, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D), rcov=~units, data=df_Coprococcus1)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[3,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_df <- as.data.frame(cbind("Coprococcus1",h2))
names(h2_df)<-c("Genus","h2")
h2_list<-rbind(h2_list,data.frame(h2_df))
#outputBLUPs for GWAS:
fit1 <- mmer(Coprococcus1~1, random=~BeanLine, rcov=~units, data=df_Coprococcus1)
Coprococcus1BLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(Coprococcus1BLUPS)<-df770meds$BeanLine
#export csv for mapping:
write.csv(Coprococcus1BLUPS,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Coprococcus1_BLUPS",quote = FALSE,row.names=TRUE)

###Dorea
df_Dorea<-as.data.frame(df770$Dorea)
df_Dorea<-log(df_Dorea,2)
df_Dorea$BeanLine<-df770$BeanLine
df_Dorea$idd<-df770$BeanLine
names(df_Dorea)<-c("Dorea","BeanLine","idd")
Dorea_out<-boxplot(df_Dorea$Dorea,range=2.5)$out
df_Dorea<-df_Dorea[-which(df_Dorea$Dorea %in% Dorea_out),]
hist(df_Dorea$Dorea)
#calc_heritability
fit <- mmer(Dorea~1, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D), rcov=~units, data=df_Dorea)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[3,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_df <- as.data.frame(cbind("Dorea",h2))
names(h2_df)<-c("Genus","h2")
h2_list<-rbind(h2_list,data.frame(h2_df))
#outputBLUPs for GWAS:
fit1 <- mmer(Dorea~1, random=~BeanLine, rcov=~units, data=df_Dorea)
DoreaBLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(DoreaBLUPS)<-df770meds$BeanLine
#export csv for mapping:
write.csv(DoreaBLUPS,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Dorea_BLUPS",quote = FALSE,row.names=TRUE)

###Blautia
df_Blautia<-as.data.frame(df770$Blautia)
df_Blautia<-log(df_Blautia,2)
df_Blautia$BeanLine<-df770$BeanLine
df_Blautia$idd<-df770$BeanLine
names(df_Blautia)<-c("Blautia","BeanLine","idd")
Blautia_out<-boxplot(df_Blautia$Blautia,range=2.5)$out
df_Blautia<-df_Blautia[-which(df_Blautia$Blautia %in% Blautia_out),]
hist(df_Blautia$Blautia)
#calc_heritability
fit <- mmer(Blautia~1, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D), rcov=~units, data=df_Blautia)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[3,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_df <- as.data.frame(cbind("Blautia",h2))
names(h2_df)<-c("Genus","h2")
h2_list<-rbind(h2_list,data.frame(h2_df))
#outputBLUPs for GWAS:
fit1 <- mmer(Blautia~1, random=~BeanLine, rcov=~units, data=df_Blautia)
BlautiaBLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(BlautiaBLUPS)<-df770meds$BeanLine
#export csv for mapping:
write.csv(BlautiaBLUPS,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Blautia_BLUPS",quote = FALSE,row.names=TRUE)

###Roseburia
df_Roseburia<-as.data.frame(df770$Roseburia)
df_Roseburia<-log(df_Roseburia,2)
df_Roseburia$BeanLine<-df770$BeanLine
df_Roseburia$idd<-df770$BeanLine
names(df_Roseburia)<-c("Roseburia","BeanLine","idd")
Roseburia_out<-boxplot(df_Roseburia$Roseburia,range=2.5)$out
df_Roseburia<-df_Roseburia[-which(df_Roseburia$Roseburia %in% Roseburia_out),]
hist(df_Roseburia$Roseburia)
#calc_heritability
fit <- mmer(Roseburia~1, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D), rcov=~units, data=df_Roseburia)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[3,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_df <- as.data.frame(cbind("Roseburia",h2))
names(h2_df)<-c("Genus","h2")
h2_list<-rbind(h2_list,data.frame(h2_df))
#outputBLUPs for GWAS:
fit1 <- mmer(Roseburia~1, random=~BeanLine, rcov=~units, data=df_Roseburia)
RoseburiaBLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(RoseburiaBLUPS)<-df770meds$BeanLine
#export csv for mapping:
write.csv(RoseburiaBLUPS,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Roseburia_BLUPS",quote = FALSE,row.names=TRUE)

###Fournierella
df_Fournierella<-as.data.frame(df770$Fournierella)
df_Fournierella<-log(df_Fournierella,2)
df_Fournierella$BeanLine<-df770$BeanLine
df_Fournierella$idd<-df770$BeanLine
names(df_Fournierella)<-c("Fournierella","BeanLine","idd")
Fournierella_out<-boxplot(df_Fournierella$Fournierella,range=2.5)$out
df_Fournierella<-df_Fournierella[-which(df_Fournierella$Fournierella %in% Fournierella_out),]
hist(df_Fournierella$Fournierella)
#calc_heritability
fit <- mmer(Fournierella~1, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D), rcov=~units, data=df_Fournierella)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[3,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_df <- as.data.frame(cbind("Fournierella",h2))
names(h2_df)<-c("Genus","h2")
h2_list<-rbind(h2_list,data.frame(h2_df))
#outputBLUPs for GWAS:
fit1 <- mmer(Fournierella~1, random=~BeanLine, rcov=~units, data=df_Fournierella)
FournierellaBLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(FournierellaBLUPS)<-df770meds$BeanLine
#export csv for mapping:
write.csv(FournierellaBLUPS,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Fournierella_BLUPS",quote = FALSE,row.names=TRUE)

###Bifidobacterium
df_Bifidobacterium<-as.data.frame(df770$Bifidobacterium)
df_Bifidobacterium<-log(df_Bifidobacterium,2)
df_Bifidobacterium$BeanLine<-df770$BeanLine
df_Bifidobacterium$idd<-df770$BeanLine
names(df_Bifidobacterium)<-c("Bifidobacterium","BeanLine","idd")
Bifidobacterium_out<-boxplot(df_Bifidobacterium$Bifidobacterium,range=2.5)$out
df_Bifidobacterium<-df_Bifidobacterium[-which(df_Bifidobacterium$Bifidobacterium %in% Bifidobacterium_out),]
hist(df_Bifidobacterium$Bifidobacterium)
#calc_heritability
fit <- mmer(Bifidobacterium~1, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D), rcov=~units, data=df_Bifidobacterium)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[3,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_df <- as.data.frame(cbind("Bifidobacterium",h2))
names(h2_df)<-c("Genus","h2")
h2_list<-rbind(h2_list,data.frame(h2_df))
#outputBLUPs for GWAS:
fit1 <- mmer(Bifidobacterium~1, random=~BeanLine, rcov=~units, data=df_Bifidobacterium)
BifidobacteriumBLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(BifidobacteriumBLUPS)<-df770meds$BeanLine
#export csv for mapping:
write.csv(BifidobacteriumBLUPS,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Bifidobacterium_BLUPS",quote = FALSE,row.names=TRUE)

###Clostridium
df_Clostridium<-as.data.frame(df770$Clostridium)
df_Clostridium<-log(df_Clostridium,2)
df_Clostridium$BeanLine<-df770$BeanLine
df_Clostridium$idd<-df770$BeanLine
names(df_Clostridium)<-c("Clostridium","BeanLine","idd")
Clostridium_out<-boxplot(df_Clostridium$Clostridium,range=2.5)$out
#df_Clostridium<-df_Clostridium[-which(df_Clostridium$Clostridium %in% Clostridium_out),]
hist(df_Clostridium$Clostridium)
#calc_heritability
fit <- mmer(Clostridium~1, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D), rcov=~units, data=df_Clostridium)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[3,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_df <- as.data.frame(cbind("Clostridium",h2))
names(h2_df)<-c("Genus","h2")
h2_list<-rbind(h2_list,data.frame(h2_df))
#outputBLUPs for GWAS:
fit1 <- mmer(Clostridium~1, random=~BeanLine, rcov=~units, data=df_Clostridium)
ClostridiumBLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(ClostridiumBLUPS)<-df770meds$BeanLine
#export csv for mapping:
write.csv(ClostridiumBLUPS,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Clostridium_BLUPS",quote = FALSE,row.names=TRUE)

###Lachnospira
df_Lachnospira<-as.data.frame(df770$Lachnospira)
df_Lachnospira<-log(df_Lachnospira,2)
df_Lachnospira$BeanLine<-df770$BeanLine
df_Lachnospira$idd<-df770$BeanLine
names(df_Lachnospira)<-c("Lachnospira","BeanLine","idd")
Lachnospira_out<-boxplot(df_Lachnospira$Lachnospira,range=2.5)$out
df_Lachnospira<-df_Lachnospira[-which(df_Lachnospira$Lachnospira %in% Lachnospira_out),]
hist(df_Lachnospira$Lachnospira)
#calc_heritability
fit <- mmer(Lachnospira~1, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D), rcov=~units, data=df_Lachnospira)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[3,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_df <- as.data.frame(cbind("Lachnospira",h2))
names(h2_df)<-c("Genus","h2")
h2_list<-rbind(h2_list,data.frame(h2_df))
#outputBLUPs for GWAS:
fit1 <- mmer(Lachnospira~1, random=~BeanLine, rcov=~units, data=df_Lachnospira)
LachnospiraBLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(LachnospiraBLUPS)<-df770meds$BeanLine
#export csv for mapping:
write.csv(LachnospiraBLUPS,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Lachnospira_BLUPS",quote = FALSE,row.names=TRUE)

###Faecalibacterium
df_Faecalibacterium<-as.data.frame(df770$Faecalibacterium)
df_Faecalibacterium<-log(df_Faecalibacterium,2)
df_Faecalibacterium$BeanLine<-df770$BeanLine
df_Faecalibacterium$idd<-df770$BeanLine
names(df_Faecalibacterium)<-c("Faecalibacterium","BeanLine","idd")
Faecalibacterium_out<-boxplot(df_Faecalibacterium$Faecalibacterium,range=2.5)$out
df_Faecalibacterium<-df_Faecalibacterium[-which(df_Faecalibacterium$Faecalibacterium %in% Faecalibacterium_out),]
hist(df_Faecalibacterium$Faecalibacterium)
#calc_heritability
fit <- mmer(Faecalibacterium~1, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D), rcov=~units, data=df_Faecalibacterium)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[3,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_df <- as.data.frame(cbind("Faecalibacterium",h2))
names(h2_df)<-c("Genus","h2")
h2_list<-rbind(h2_list,data.frame(h2_df))
#outputBLUPs for GWAS:
fit1 <- mmer(Faecalibacterium~1, random=~BeanLine, rcov=~units, data=df_Faecalibacterium)
FaecalibacteriumBLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(FaecalibacteriumBLUPS)<-df770meds$BeanLine
#export csv for mapping:
write.csv(FaecalibacteriumBLUPS,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Faecalibacterium_BLUPS",quote = FALSE,row.names=TRUE)

#h2_list_r3<-h2_list
#h2_list_r2.5<-h2_list
h2_list_r2.5.1<-h2_list

h2_round <- round(as.numeric(h2_list_r2.5$h2),3) 
h2_list_r2.5$h2<-as.numeric(h2_list_r2.5$h2)

map <- ggplot(h2_list_r2.5,aes(x=Genus,y=h2))+geom_point(size=13,color="White")+theme(axis.text.x = element_text(angle=45, hjust=1))+theme(text=element_text(size=20,family = "serif"))+
  ggtitle("Common Bean Heritability of Subject770 bacterial response")+scale_y_continuous(name="h2",limits=c(0,0.5))

plot(map)

map+geom_text(aes(label=h2_round))











