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
library(dplyr)
#library("qgraph")
library(corrplot)
library(factoextra)


#!!!!!Make sure Beanline column has identical identifiers to SNP file
#Set working directory to location of OTU table and SNP file
setwd("~/MADP")

##insert OTU table here:
df <- read_excel("Analysis/770/MADPS770Samples_RarefiedGenusAbsoluteFilterTable.xlsx",sheet="MADPS770_subset_abscounts_T")
df1 <- read_excel("Analysis/770/MADPS770Samples_RarefiedGenusAbsoluteFilterTable.xlsx",sheet="Abs_Counts_filter")
df2 <- read_excel("Analysis/770/MADPS770Samples_RarefiedGenusAbsoluteFilterTable.xlsx",sheet="RelAbundGenotypes")



df2.1<-df2[-c(1:23)]

df3<-as.data.frame(t(df2.1))
#df2<-df2[-c(1)]
names(df3)<-df2$BeanLine



df3[21:25,1:5]

df3.2<-df3[which(rownames(df3) %in% df$'#OTU ID'),]

df770<-df[-c(1)]
#rownames(df770)<-df$`#OTU ID`

df770_T<-as.data.frame(t(df3.2))
df770_T.1<-df770_T
df770_T.1$BeanLine<-df2$BeanLine

df770meds<-df770_T.1 %>% group_by(BeanLine) %>% summarise_each(funs(median))

corrmatrix1<-cor(df770,method="pearson")
corrmatrix2<-cor(df770_T,method="pearson")


row.names(corrmatrix1)<-corrmatrix1$X.OTU.ID
cormatrix1<-cormatrix1.W[-c(1)]
#colnames(cormatrix1)<-row.names(cormatrix1)

NamNum1<-as.data.frame(colnames(cormatrix1))
colnames(NamNum1)<-"taxa"

##
rownames(corrmatrix1)<-rownames(NamNum1)
colnames(corrmatrix1)<-rownames(NamNum1)


Graph_pcor <- qgraph(corrmatrix2, graph = "cor",layout="spring",graph=cor, threshold=0.3,title="MADP_770_trim",details=TRUE, theme="colorblind",label.scale=FALSE,label.scale.equal=TRUE)


corrplot(corrmatrix2,type="upper",order="hclust",method="color",tl.col="black",diag=FALSE)#,title="MAPD Subject770 RelativeAbundance Correlation",mar=c(0,0,3,0))


heatmap(as.matrix(corrmatrix2),scale="column")#, col = cm.colors(256))

###############################################################
 #################
#Variation within lines per taxa:

#var770<-with(df770_T.1, tapply(df770_T.1[,1],df770_T.1$BeanLine,boxplot,range=0.5))

#trimoutliers<-function(Line){
#  for (i in rep(1:22, each =1))
#  var770<-with(df770_T.1, tapply(df770_T.1[,i],df770_T.1$BeanLine,boxplot,range=2))
#  varout<-var770$paste0(Line)
#  df770_T.3<-as.data.frame(df770_T.1[-which((df770_T.1[,i] %in% varout)& df770_T.1$BeanLine==paste0(",Line,")),])
#}

#trimoutliers(A_55,1)

#varout<-var770$A_55$out
#df770_T.3<-as.data.frame(df770_T.1[-which((df770_T.1[,1] %in% varout)& df770_T.1$BeanLine=="A_55"),])

#Prevotella9_out<-boxplot(df_Prevotella9$Prevotella9,range=3)$out
#df_Prevotella9<-df_Prevotella9[-which(df_Prevotella9$Prevotella9 %in% Prevotella9_out),]

#for (i in rep(2:22, each = 1))
#{
#  VP770<-as.data.frame(with(df770_T.1, tapply(df770_T.1[,i],df770_T.1$BeanLine,var)))
#  var770<-as.data.frame(cbind(var770,VP770))
#}

#trimoutliers<-function(Line){
#  for (i in rep(1:22, each =1))
#    var770<-with(df770_T.1, tapply(df770_T.1[,i],df770_T.1$BeanLine,boxplot,range=2))
#  varout<-var770$paste0(Line)
#  df770_T.3<-as.data.frame(df770_T.1[-which((df770_T.1[,i] %in% varout)& df770_T.1$BeanLine==paste0(",Line,")),])
#}
GD <- read.delim("Genotypes/469-geno-200k-SNPs-imputed.hmp.txt")
GD1 <- GD[-c(2:11)]
GD2 <- as.data.frame(as.matrix(t(GD1))) 
colnames(GD2) <- rownames(GD1)
rownames(GD2) <- colnames(GD1)

GD2.1<-as.data.frame(rownames(GD2))
GD2$BeanLine<-GD2.1$`rownames(GD2)`

GD2.2<-GD2[which(GD2$BeanLine %in% df2$BeanLine ),]
GD2.2$BeanLine<-NULL

#convert to numeric
GD3 <- atcg1234(GD2.2, ploidy=2)
GD3$M[1:5,1:5]
#remove extra line at top, make element of list into dataframe
GD5 <- GD3$M
GD5[1:20,1:2]

GD6<-as.data.frame(rownames(GD2.2))
names(GD6)<-"BeanLine"

A <- A.mat(GD5)
D <- D.mat(GD5)
#E <- E.mat(GD5)
###################################



df770_T.2<-df770_T.1
df770_T.2$BeanLine<-NULL
microPCs<-prcomp(df770_T.2,scale. = TRUE, center = TRUE)
microPCs.1<-as.data.frame(microPCs$x)
#microPCs.1T<-as.data.frame(t(microPCs.1))

#Faecalibacterium
df_Faecalibacterium<-as.data.frame(df770_T.1$Faecalibacterium)
df_Faecalibacterium<-log(df_Faecalibacterium,2)
df_Faecalibacterium$BeanLine<-df770_T.1$BeanLine
df_Faecalibacterium$idd<-df2$BeanLine
df_Faecalibacterium$Plate<-as.character(df2$Plate)
df_Faecalibacterium$Race<-(df2$BeanRace)
df_Faecalibacterium$PC1<-microPCs.1$PC1
df_Faecalibacterium$PC2<-microPCs.1$PC2
df_Faecalibacterium$PC3<-microPCs.1$PC3
df_Faecalibacterium$PC4<-microPCs.1$PC4
df_Faecalibacterium$PC5<-microPCs.1$PC5
names(df_Faecalibacterium)<-c("Faecalibacterium","BeanLine","idd","Plate","Race","PC1","PC2","PC3","PC4","PC5")
varFaecalibacterium<-with(df_Faecalibacterium, tapply(df_Faecalibacterium$Faecalibacterium,df_Faecalibacterium$BeanLine,boxplot,range=1.8))
df_Faecalibacterium.1<-df_Faecalibacterium
for (name in GD6$BeanLine){
  name1<-paste0("varFaecalibacterium$",name,"$out",sep="")
  name2<-noquote(paste0(name1))
  varout<-eval(parse(text=name2))
  if (is_empty(varout)==TRUE){
    print("nothing")
  } else {
    df_Faecalibacterium.1<-as.data.frame(df_Faecalibacterium.1[-which((df_Faecalibacterium.1$Faecalibacterium %in% varout)& df_Faecalibacterium.1$BeanLine %in% name),])
  }
}
Faecalibacterium_out<-boxplot(df_Faecalibacterium.1$Faecalibacterium,range=6)$out
df_Faecalibacterium.1<-df_Faecalibacterium.1[-which(df_Faecalibacterium.1$Faecalibacterium %in% Faecalibacterium_out),]
hist(df_Faecalibacterium.1$Faecalibacterium)
fit1 <- mmer(Faecalibacterium~cbind(PC1,PC2,PC3,PC4), random=~ vs(BeanLine,Gu=A)+vs(idd,Gu=D)+Plate, rcov=~units, data=df_Faecalibacterium.1,tolparinv = 7)
summary(fit1)
FaecalibacteriumBLUPS<-as.data.frame(randef(fit1)$'u:BeanLine')
rownames(FaecalibacteriumBLUPS)<-df770meds$BeanLine
#order phenotype to match genotype
ordA<- match(GD6$BeanLine,rownames(FaecalibacteriumBLUPS))
FaecalibacteriumBLUPS.1 <- as.data.frame(FaecalibacteriumBLUPS[ordA,])
rownames(FaecalibacteriumBLUPS.1)<-GD6$BeanLine
names(FaecalibacteriumBLUPS.1)<-"Faecalibacterium"
hist(FaecalibacteriumBLUPS.1$Faecalibacterium)
#Scales values to between 0 and 1:
#FaecalibacteriumBLUPS.1N<-(FaecalibacteriumBLUPS.1[1]-min(FaecalibacteriumBLUPS.1[1]))/(max(FaecalibacteriumBLUPS.1[1])-min(FaecalibacteriumBLUPS.1[1]))
#hist(FaecalibacteriumBLUPS.1N$Faecalibacterium)
#export csv for mapping:
write.csv(FaecalibacteriumBLUPS.1,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Faecalibacterium_BLUPS8",quote = FALSE,row.names=TRUE)



#Prevotella9
df_Prevotella9<-as.data.frame(df770_T.1$Prevotella9)
df_Prevotella9<-log(df_Prevotella9,2)
df_Prevotella9$BeanLine<-df770_T.1$BeanLine
df_Prevotella9$idd<-df2$BeanLine
df_Prevotella9$Plate<-as.character(df2$Plate)
df_Prevotella9$Race<-(df2$BeanRace)
df_Prevotella9$PC1<-microPCs.1$PC1
df_Prevotella9$PC2<-microPCs.1$PC2
df_Prevotella9$PC3<-microPCs.1$PC3
df_Prevotella9$PC4<-microPCs.1$PC4
df_Prevotella9$PC5<-microPCs.1$PC5
names(df_Prevotella9)<-c("Prevotella9","BeanLine","idd","Plate","Race","PC1","PC2","PC3","PC4","PC5")
varPrevotella9<-with(df_Prevotella9, tapply(df_Prevotella9$Prevotella9,df_Prevotella9$BeanLine,boxplot,range=0.9))
df_Prevotella9.1<-df_Prevotella9
for (name in GD6$BeanLine){
  name1<-paste0("varPrevotella9$",name,"$out",sep="")
  name2<-noquote(paste0(name1))
  varout<-eval(parse(text=name2))
  if (is_empty(varout)==TRUE){
   print("nothing")
  } else {
   df_Prevotella9.1<-as.data.frame(df_Prevotella9.1[-which((df_Prevotella9.1$Prevotella9 %in% varout)& df_Prevotella9.1$BeanLine %in% name),])
  }
}
Prevotella9_out<-boxplot(df_Prevotella9.1$Prevotella9,range=4)$out
df_Prevotella9.1<-df_Prevotella9.1[-which(df_Prevotella9.1$Prevotella9 %in% Prevotella9_out),]
hist(df_Prevotella9.1$Prevotella9)
fit1 <- mmer(Prevotella9~Race, random=~ vs(BeanLine,Gu=A)+vs(Plate)+PC1+PC2+PC3, rcov=~units, data=df_Prevotella9.1,tolparinv = 7)
summary(fit1)
Prevotella9BLUPS<-as.data.frame(randef(fit1)$'u:BeanLine')
rownames(Prevotella9BLUPS)<-df770meds$BeanLine
#order phenotype to match genotype
ordA<- match(GD6$BeanLine,rownames(Prevotella9BLUPS))
Prevotella9BLUPS.1 <- as.data.frame(Prevotella9BLUPS[ordA,])
rownames(Prevotella9BLUPS.1)<-GD6$BeanLine
names(Prevotella9BLUPS.1)<-"Prevotella9"
hist(Prevotella9BLUPS.1$Prevotella9)
#export csv for mapping:
write.csv(Prevotella9BLUPS.1,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Prevotella9_BLUPS",quote = FALSE,row.names=TRUE)


#Bifidobacterium
#Bifidobacterium
df_Bifidobacterium<-as.data.frame(df770_T.1$Bifidobacterium)
df_Bifidobacterium<-log(df_Bifidobacterium,2)
df_Bifidobacterium$BeanLine<-df770_T.1$BeanLine
df_Bifidobacterium$idd<-df2$BeanLine
df_Bifidobacterium$Plate<-as.character(df2$Plate)
df_Bifidobacterium$Race<-(df2$BeanRace)
names(df_Bifidobacterium)<-c("Bifidobacterium","BeanLine","idd","Plate","Race")
varBifidobacterium<-with(df_Bifidobacterium, tapply(df_Bifidobacterium$Bifidobacterium,df_Bifidobacterium$BeanLine,boxplot,range=0.9))
df_Bifidobacterium.1<-df_Bifidobacterium
for (name in GD6$BeanLine){
  name1<-paste0("varBifidobacterium$",name,"$out",sep="")
  name2<-noquote(paste0(name1))
  varout<-eval(parse(text=name2))
  if (is_empty(varout)==TRUE){
    print("nothing")
  } else {
    df_Bifidobacterium.1<-as.data.frame(df_Bifidobacterium.1[-which((df_Bifidobacterium.1$Bifidobacterium %in% varout)& df_Bifidobacterium.1$BeanLine %in% name),])
  }
}
Bifidobacterium_out<-boxplot(df_Bifidobacterium.1$Bifidobacterium,range=5)$out
df_Bifidobacterium.1<-df_Bifidobacterium.1[-which(df_Bifidobacterium.1$Bifidobacterium %in% Bifidobacterium_out),]
hist(df_Bifidobacterium.1$Bifidobacterium)
fit1 <- mmer(Bifidobacterium~Race, random=~vs(BeanLine,Gu=A)+Plate, rcov=~units, data=df_Bifidobacterium.1,tolparinv = 7)
summary(fit1)
BifidobacteriumBLUPS<-as.data.frame(randef(fit1)$'u:BeanLine')
rownames(BifidobacteriumBLUPS)<-df770meds$BeanLine
#order phenotype to match genotype
ordA<- match(GD6$BeanLine,rownames(BifidobacteriumBLUPS))
BifidobacteriumBLUPS.1 <- as.data.frame(BifidobacteriumBLUPS[ordA,])
rownames(BifidobacteriumBLUPS.1)<-GD6$BeanLine
names(BifidobacteriumBLUPS.1)<-"Bifidobacterium"
hist(BifidobacteriumBLUPS.1)
#export csv for mapping:
write.csv(BifidobacteriumBLUPS.1,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Bifidobacterium_BLUPS",quote = FALSE,row.names=TRUE)


#Coprococcus3
df_Coprococcus3<-as.data.frame(df770_T.1$Coprococcus3)
df_Coprococcus3<-log(df_Coprococcus3,2)
df_Coprococcus3$BeanLine<-df770_T.1$BeanLine
df_Coprococcus3$idd<-df2$BeanLine
df_Coprococcus3$Plate<-as.character(df2$Plate)
df_Coprococcus3$Race<-(df2$BeanRace)
names(df_Coprococcus3)<-c("Coprococcus3","BeanLine","idd","Plate","Race")
varCoprococcus3<-with(df_Coprococcus3, tapply(df_Coprococcus3$Coprococcus3,df_Coprococcus3$BeanLine,boxplot,range=0.9))
df_Coprococcus3.1<-df_Coprococcus3
for (name in GD6$BeanLine){
  name1<-paste0("varCoprococcus3$",name,"$out",sep="")
  name2<-noquote(paste0(name1))
  varout<-eval(parse(text=name2))
  if (is_empty(varout)==TRUE){
    print("nothing")
  } else {
    df_Coprococcus3.1<-as.data.frame(df_Coprococcus3.1[-which((df_Coprococcus3.1$Coprococcus3 %in% varout)& df_Coprococcus3.1$BeanLine %in% name),])
  }
}
Coprococcus3_out<-boxplot(df_Coprococcus3.1$Coprococcus3,range=4)$out
df_Coprococcus3.1<-df_Coprococcus3.1[-which(df_Coprococcus3.1$Coprococcus3 %in% Coprococcus3_out),]
hist(df_Coprococcus3.1$Coprococcus3)
fit1 <- mmer(Coprococcus3~Race, random=~vs(BeanLine,Gu=A)+Plate, rcov=~units, data=df_Coprococcus3.1,tolparinv = 7)
summary(fit1)
Coprococcus3BLUPS<-as.data.frame(randef(fit1)$'u:BeanLine')
rownames(Coprococcus3BLUPS)<-df770meds$BeanLine
#order phenotype to match genotype
ordA<- match(GD6$BeanLine,rownames(Coprococcus3BLUPS))
Coprococcus3BLUPS.1 <- as.data.frame(Coprococcus3BLUPS[ordA,])
rownames(Coprococcus3BLUPS.1)<-GD6$BeanLine
names(Coprococcus3BLUPS.1)<-"Coprococcus3"
hist(Coprococcus3BLUPS.1)
#export csv for mapping:
write.csv(Coprococcus3BLUPS.1,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Coprococcus3_BLUPS",quote = FALSE,row.names=TRUE)
##########################
#PCAs
df770_T.2<-df770_T.1
df770_T.2$BeanLine<-NULL
PCA<- prcomp(df770_T.2, scale=TRUE)
PCA
ev<-as.data.frame(PCA$x)

summary(PCA)

#PC1
df_PC1<-as.data.frame(ev$PC1)
#df_PC1<-log(df_PC1,2)
df_PC1$PC2<-ev$PC2
df_PC1$PC3<-ev$PC3
df_PC1$BeanLine<-df770_T.1$BeanLine
df_PC1$idd<-df2$BeanLine
df_PC1$Plate<-as.character(df2$Plate)
df_PC1$Race<-(df2$BeanRace)
names(df_PC1)<-c("PC1","PC2","PC3","BeanLine","idd","Plate","Race")
varPC1<-with(df_PC1, tapply(df_PC1$PC1,df_PC1$BeanLine,boxplot,range=0.9))
df_PC1.1<-df_PC1
for (name in GD6$BeanLine){
  name1<-paste0("varPC1$",name,"$out",sep="")
  name2<-noquote(paste0(name1))
  varout<-eval(parse(text=name2))
  if (is_empty(varout)==TRUE){
    print("nothing")
  } else {
    df_PC1.1<-as.data.frame(df_PC1.1[-which((df_PC1.1$PC1 %in% varout)& df_PC1.1$BeanLine %in% name),])
  }
}
PC1_out<-boxplot(df_PC1.1$PC1,range=3)$out
df_PC1.1<-df_PC1.1[-which(df_PC1.1$PC1 %in% PC1_out),]
hist(df_PC1.1$PC1)
fit1 <- mmer(PC1~Race, random=~vs(BeanLine,Gu=A)+Plate, rcov=~units, data=df_PC1.1,tolparinv = 7)
summary(fit1)
PC1BLUPS<-as.data.frame(randef(fit1)$'u:BeanLine')
rownames(PC1BLUPS)<-df770meds$BeanLine
#order phenotype to match genotype
ordA<- match(GD6$BeanLine,rownames(PC1BLUPS))
PC1BLUPS.1 <- as.data.frame(PC1BLUPS[ordA,])
rownames(PC1BLUPS.1)<-GD6$BeanLine
names(PC1BLUPS.1)<-"PC1"
hist(PC1BLUPS.1)
#export csv for mapping:
write.csv(PC1BLUPS.1,"./Analysis/770/Phenotypes/BLUPS/MADP_770_PC1_BLUPS",quote = FALSE,row.names=TRUE)

#PC2
df_PC2<-as.data.frame(ev$PC2)
#df_PC2<-log(df_PC2,2)
df_PC2$PC2<-ev$PC2
df_PC2$PC3<-ev$PC3
df_PC2$BeanLine<-df770_T.1$BeanLine
df_PC2$idd<-df2$BeanLine
df_PC2$Plate<-as.character(df2$Plate)
df_PC2$Race<-(df2$BeanRace)
names(df_PC2)<-c("PC2","PC2","PC3","BeanLine","idd","Plate","Race")
varPC2<-with(df_PC2, tapply(df_PC2$PC2,df_PC2$BeanLine,boxplot,range=0.9))
df_PC2.1<-df_PC2
for (name in GD6$BeanLine){
  name1<-paste0("varPC2$",name,"$out",sep="")
  name2<-noquote(paste0(name1))
  varout<-eval(parse(text=name2))
  if (is_empty(varout)==TRUE){
    print("nothing")
  } else {
    df_PC2.1<-as.data.frame(df_PC2.1[-which((df_PC2.1$PC2 %in% varout)& df_PC2.1$BeanLine %in% name),])
  }
}
PC2_out<-boxplot(df_PC2.1$PC2,range=3)$out
df_PC2.1<-df_PC2.1[-which(df_PC2.1$PC2 %in% PC2_out),]
hist(df_PC2.1$PC2)
fit1 <- mmer(PC2~Race, random=~vs(BeanLine,Gu=A)+Plate, rcov=~units, data=df_PC2.1,tolparinv = 7)
summary(fit1)
PC2BLUPS<-as.data.frame(randef(fit1)$'u:BeanLine')
rownames(PC2BLUPS)<-df770meds$BeanLine
#order phenotype to match genotype
ordA<- match(GD6$BeanLine,rownames(PC2BLUPS))
PC2BLUPS.1 <- as.data.frame(PC2BLUPS[ordA,])
rownames(PC2BLUPS.1)<-GD6$BeanLine
names(PC2BLUPS.1)<-"PC2"
hist(PC2BLUPS.1)
#export csv for mapping:
write.csv(PC2BLUPS.1,"./Analysis/770/Phenotypes/BLUPS/MADP_770_PC2_BLUPS",quote = FALSE,row.names=TRUE)