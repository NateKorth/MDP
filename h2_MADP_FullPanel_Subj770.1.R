##Library packages (install.packages() if you need to)
library(sommer)
library(ggplot2)
library(readxl)
library(data.table)
library(tidyverse)

#!!!!!Make sure Beanline column has identical identifiers to SNP file
#Set working directory to location of OTU table and SNP file
setwd("~/MADP")

##insert OTU table here:
df <- read_excel("Analysis/770/MADPS770Samples_RarefiedGenusAbsoluteFilterTable.xlsx",sheet="RelAbundGenotypes")



#Subset data by subject
df770 <- as.data.frame(subset(df, Subject=="S770"))
df770 <- subset(df770, ControlType=="Sample")
df770 <- subset(df770, SampleLost=="No")





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
fit <- mmer(Prevotella9~Race, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D)+Plate, rcov=~units, data=df_Prevotella9)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[4,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_list <- data.frame("Prevotella9", h2 )
names(h2_list)<-c("Genus","h2")
#outputBLUPs for GWAS:
fit1 <- mmer(Prevotella9~Race, random=~vs(BeanLine,Gu=A)+Plate, rcov=~units, data=df_Prevotella9,tolparinv = 7)
summary(fit1)
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
#df_Enterobacteriaceae<-df_Enterobacteriaceae[-which(df_Enterobacteriaceae$Enterobacteriaceae %in% Enterobacteriaceae_out),]
hist(df_Enterobacteriaceae$Enterobacteriaceae)
#calc_heritability
fit <- mmer(Enterobacteriaceae~Race, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D)+Plate, rcov=~units, data=df_Enterobacteriaceae)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[4,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_df <- as.data.frame(cbind("Enterobacteriaceae",h2))
names(h2_df)<-c("Genus","h2")
h2_list<-rbind(h2_list,data.frame(h2_df))
names(h2_list)<-c("Genus","h2")
#outputBLUPs for GWAS:
fit1 <- mmer(Enterobacteriaceae~1, random=~BeanLine+vs(BeanLine,Gu=A) + vs(idd,Gu=D)+Race+Plate+Batch, rcov=~units, data=df_Enterobacteriaceae,tolparinv = 7)
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
fit <- mmer(Succinivibrio~Race, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D)+Plate, rcov=~units, data=df_Succinivibrio)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[4,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_df <- as.data.frame(cbind("Succinivibrio",h2))
names(h2_df)<-c("Genus","h2")
h2_list<-rbind(h2_list,data.frame(h2_df))
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
fit <- mmer(Dialister~Race, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D)+Plate, rcov=~units, data=df_Dialister)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[4,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_df <- as.data.frame(cbind("Dialister",h2))
names(h2_df)<-c("Genus","h2")
h2_list<-rbind(h2_list,data.frame(h2_df))
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
fit <- mmer(Bacteroides~Race, random=~ vs(BeanLine,Gu=A) + vs(idd,Gu=D)+Plate, rcov=~units, data=df_Bacteroides)
VC <- summary(fit)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VE <- VC[4,1]
h2 <- (VA+VD)/(VA+VD+VE)
h2_df <- as.data.frame(cbind("Bacteroides",h2))
names(h2_df)<-c("Genus","h2")
h2_list<-rbind(h2_list,data.frame(h2_df))
names(h2_list)<-c("Genus","h2")
#outputBLUPs for GWAS:
fit1 <- mmer(Bacteroides~1, random=~BeanLine+vs(BeanLine,Gu=A) + vs(idd,Gu=D)+Batch+Plate+Race, rcov=~units, data=df_Bacteroides,tolparinv = 7)
#summary(fit1)
BacteroidesBLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(BacteroidesBLUPS)<-df770meds$BeanLine
#export csv for mapping:
write.csv(BacteroidesBLUPS,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Bacteroides_BLUPS",quote = FALSE,row.names=TRUE)
 





################### Make masterheritability list

h765 <- subset(h2_765,Subject==765)
h766 <- subset(h2_766,Subject==766)
h767 <- subset(h2_767,Subject==767)
h768 <- subset(h2_768,Subject==768)

###bind all subjects
h2_Master <- rbind(h765,h766,h767,h768)

###Round to the hundreth to put values on heatmap
h2_round <- round(h2_Master$h2,3) 

h2<-h2_round$

########################### 
#Make Heatmap

map <- ggplot(h2_Master, aes(Subject, Genus))+geom_tile(aes(fill=h2), colour="white")+scale_fill_gradient(low="white", high="red")
map+geom_text(aes(label=h2_round))

plot(map)


#######################################################################################################
##Working the mixed linear model 
#Usage mmer(taxa~1, random=~ factor of interest1+factor of interest2, 
#where f1:f2 tests interactions and vs(spl2D(Row,Column) tests 2 dimensional efx, rcov=~units(cov),data)

fit7 <- mmer(Faecalibacterium~1, random=~ vs(BeanLine,Gu=A)+vs(idd,Gu=D)+Subject+BeanEater+Gender+SujectOrigin+Subject:BeanLine+ BeanRace + MarketClass+vs(spl2D(Row,Column)),rcov=~units, data=df)

#####calculating variance. IF you change model variance estimates must be altered
summary(fit7)$varcomp
VC <- summary(fit7)$varcomp
VA <- VC[1,1]
VD <- VC[2,1]
VLine <-VA+VD
VLine<-VC[1,1]
VSubject <- VC[3,1]
VBeaneater<- VC[4,1]
VGender<-VC[5,1] 
VSorigin<-VC[6,1] 
VSxL<-VC[7,1] 
VRace<-VC[8,1]
VMC<-VC[9,1]
VPosition<-VC[10,1]
Vunexplained <-VC[11,1]

##Make df of all effects of interest and their corresponding variance components
vars1<-data.frame(group=c("Line","Subject","BeanEater","Gender","Sorigin","SxL","BRace","BMC","Position","Unexplained"),
                  value=c(VLine,VSubject,VBeaneater,VGender,VSorigin,VSxL,VRace,VMC,VPosition,Vunexplained))

bp<- ggplot(vars1, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")
pie<-bp+coord_polar("y",start=0)
pie





