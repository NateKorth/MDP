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
library("qgraph")
library(corrplot)


#!!!!!Make sure Beanline column has identical identifiers to SNP file
#Set working directory to location of OTU table and SNP file
setwd("~/MADP")

##insert OTU table here:
df <- read_excel("Analysis/770/MADPS770Samples_RarefiedGenusAbsoluteFilterTable.xlsx",sheet="MADPS770_subset_abscounts_T")
df1 <- read_excel("Analysis/770/MADPS770Samples_RarefiedGenusAbsoluteFilterTable.xlsx",sheet="Abs_Counts_filter")


df770<-df[-c(1)]
rownames(df770)<-df$`#OTU ID`

df770_T<-as.data.frame(t(df770))

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


Graph_pcor <- qgraph(corrmatrix2, graph = "pcor",layout="spring",graph=cor, threshold=0.3,title="MADP_770_trim",details=TRUE, theme="colorblind",label.scale=FALSE,label.scale.equal=FALSE)


corrplot(corrmatrix2,type="upper",order="hclust",method="color",tl.col="black",diag=FALSE)


heatmap(as.matrix(corrmatrix2),scale="column")#, col = cm.colors(256))

###############################################################
df770_T$BeanLine<-df1$Genotype
GD <- read.delim("Genotypes/469-geno-200k-SNPs-imputed.hmp.txt")

GD1[1:5,1:45]

GD1 <- GD[-c(2:11)]
GD2 <- as.data.frame(as.matrix(t(GD1))) 
colnames(GD2) <- rownames(GD1)
rownames(GD2) <- colnames(GD1)
##Check data


GD2.1<-as.data.frame(rownames(GD2))
GD2$BeanLine<-GD2.1$`rownames(GD2)`

GD2.2[1:5,1:5]
GD2.2<-GD2[which(GD2$BeanLine %in% df770_T$BeanLine ),]
GD2.2$BeanLine<-NULL

#convert to numeric
GD3 <- atcg1234(GD2.2, ploidy=2)
GD3$M[1:5,1:5]
#remove extra line at top, make element of list into dataframe
GD5 <- GD3$M
GD5[1:200,1:2]

GD6<-as.data.frame(rownames(GD2.2))

###################Make covariate matrices A=additive matrix. D=Dominance matrix. E=Epistatic matrix 
A <- A.mat(GD5)
D <- D.mat(GD5)
E <- E.mat(GD5)
###

#subset lines in otu table that match line names in GD
df770_T.1<-df770_T[which(df770_T$BeanLine %in% GD6$`rownames(GD2.2)`),]
#identify line names that don't match:
df770.2<-df770_T[which(!df770_T$BeanLine %in% GD6$`rownames(GD2.2)`),]
df1.1<-df1[which(df1$Genotype %in% GD6$`rownames(GD2.2)`),]

#Get IDs:
df770vals<-df770_T.1#[-c(1:4,6:23)]
df770meds<-df770vals %>% group_by(BeanLine) %>% summarise_each(funs(median))

#########################################################################################
#Variation within lines per taxa:

var770<-with(df770_T.1, tapply(df770_T.1[,1],df770_T.1$BeanLine,boxplot,range=0.5))

trimoutliers<-function(Line){
  for (i in rep(1:22, each =1))
  var770<-with(df770_T.1, tapply(df770_T.1[,i],df770_T.1$BeanLine,boxplot,range=2))
  varout<-var770$paste0(Line)
  df770_T.3<-as.data.frame(df770_T.1[-which((df770_T.1[,i] %in% varout)& df770_T.1$BeanLine==paste0(",Line,")),])
}

trimoutliers(A_55,1)


varout<-var770$A_55$out
df770_T.3<-as.data.frame(df770_T.1[-which((df770_T.1[,1] %in% varout)& df770_T.1$BeanLine=="A_55"),])


for (name in df770_T.1$BeanLine)
{
var1<-var770$name$out

}



#Prevotella9_out<-boxplot(df_Prevotella9$Prevotella9,range=3)$out
#df_Prevotella9<-df_Prevotella9[-which(df_Prevotella9$Prevotella9 %in% Prevotella9_out),]



for (i in rep(2:22, each = 1))
{
  VP770<-as.data.frame(with(df770_T.1, tapply(df770_T.1[,i],df770_T.1$BeanLine,var)))
  var770<-as.data.frame(cbind(var770,VP770))
}

trimoutliers<-function(Line){
  for (i in rep(1:22, each =1))
    var770<-with(df770_T.1, tapply(df770_T.1[,i],df770_T.1$BeanLine,boxplot,range=2))
  varout<-var770$paste0(Line)
  df770_T.3<-as.data.frame(df770_T.1[-which((df770_T.1[,i] %in% varout)& df770_T.1$BeanLine==paste0(",Line,")),])
}

#Prevotella9
df_Prevotella9<-as.data.frame(df770_T.1$Prevotella9)
#df_Prevotella9<-log(df_Prevotella9,2)
df_Prevotella9$BeanLine<-df770_T.1$BeanLine
df_Prevotella9$idd<-df1.1$Genotype
#df_Prevotella9$Batch<-as.character(df1$DialysisBatch)
df_Prevotella9$Plate<-as.character(df1.1$Plate)
df_Prevotella9$Race<-(df1.1$BeanRace)
names(df_Prevotella9)<-c("Prevotella9","BeanLine","idd","Plate","Race")

varPrevotella9<-with(df_Prevotella9, tapply(df_Prevotella9$Prevotella9,df_Prevotella9$BeanLine,boxplot,range=.01))
df_Prevotella9.1<-df_Prevotella9

for (name in df_Prevotella9$BeanLine)
  {
  line1<-paste()
  varout<-varPrevotella9$name$out
  if(is.null(varout)==TRUE){
    next
  }
  df_Prevotella9.1<-as.data.frame(df_Prevotella9.1[-which((df_Prevotella9.1$Prevotella9 %in% varout)& df_Prevotella9.1$BeanLine==line1),])
}

line1<-paste("","Zorro","",sep="")
  
for(i in 1:10){
  str<-paste("modelCheck(var",i,"_d.bug)",sep="")
  print(str)
}


varout<-varPrevotella9$Zorro$out
df_Prevotella9.1<-as.data.frame(df_Prevotella9[-which((df_Prevotella9$Prevotella9 %in% varout) & df_Prevotella9$BeanLine==line1),])





#Prevotella9_out<-boxplot(df_Prevotella9$Prevotella9,range=3)$out
#df_Prevotella9<-df_Prevotella9[-which(df_Prevotella9$Prevotella9 %in% Prevotella9_out),]

hist(df_Prevotella9$Prevotella9)
fit1 <- mmer(Prevotella9~Race, random=~BeanLine+Plate, rcov=~units, data=df_Prevotella9,tolparinv = 7)
summary(fit1)
Prevotella9BLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(Prevotella9BLUPS)<-df770meds$BeanLine
#order phenotype to match genotype
ordA<- match(GD6$`rownames(GD2.2)`,rownames(Prevotella9BLUPS))
Prevotella9BLUPS.1 <- as.data.frame(Prevotella9BLUPS[ordA,])
rownames(Prevotella9BLUPS.1)<-GD6$`rownames(GD2.2)`
names(Prevotella9BLUPS.1)<-"Prevotella9"
#export csv for mapping:
write.csv(Prevotella9BLUPS.1,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Prevotella9_BLUPS",quote = FALSE,row.names=TRUE)

#Bifidobacterium
df_Bifidobacterium<-as.data.frame(df770_T.1$Bifidobacterium)
#df_Bifidobacterium<-log(df_Bifidobacterium,2)
df_Bifidobacterium$BeanLine<-df1.1$Genotype
df_Bifidobacterium$idd<-df1.1$Genotype
#df_Bifidobacterium$Batch<-as.character(df1$DialysisBatch)
df_Bifidobacterium$Plate<-as.character(df1.1$Plate)
df_Bifidobacterium$Race<-(df1.1$BeanRace)
names(df_Bifidobacterium)<-c("Bifidobacterium","BeanLine","idd","Plate","Race")
Bifidobacterium_out<-boxplot(df_Bifidobacterium$Bifidobacterium,range=3)$out
df_Bifidobacterium<-df_Bifidobacterium[-which(df_Bifidobacterium$Bifidobacterium %in% Bifidobacterium_out),]
hist(df_Bifidobacterium$Bifidobacterium)
fit1 <- mmer(Bifidobacterium~Race, random=~BeanLine+Plate, rcov=~units, data=df_Bifidobacterium,tolparinv = 7)
summary(fit1)
BifidobacteriumBLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(BifidobacteriumBLUPS)<-df770meds$BeanLine
#order phenotype to match genotype
ordA<- match(GD6$`rownames(GD2.2)`,rownames(BifidobacteriumBLUPS))
BifidobacteriumBLUPS.1 <- as.data.frame(BifidobacteriumBLUPS[ordA,])
rownames(BifidobacteriumBLUPS.1)<-GD6$`rownames(GD2.2)`
names(BifidobacteriumBLUPS.1)<-"Bifidobacterium"
#export csv for mapping:
write.csv(BifidobacteriumBLUPS.1,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Bifidobacterium_BLUPS",quote = FALSE,row.names=TRUE)

#Coprococcus3
df_Coprococcus3<-as.data.frame(df770_T.1$Coprococcus3)
#df_Coprococcus3<-log(df_Coprococcus3,2)
df_Coprococcus3$BeanLine<-df1.1$Genotype
df_Coprococcus3$idd<-df1.1$Genotype
#df_Coprococcus3$Batch<-as.character(df1$DialysisBatch)
df_Coprococcus3$Plate<-as.character(df1.1$Plate)
df_Coprococcus3$Race<-(df1.1$BeanRace)
names(df_Coprococcus3)<-c("Coprococcus3","BeanLine","idd","Plate","Race")
Coprococcus3_out<-boxplot(df_Coprococcus3$Coprococcus3,range=3)$out
#df_Coprococcus3<-df_Coprococcus3[-which(df_Coprococcus3$Coprococcus3 %in% Coprococcus3_out),]
hist(df_Coprococcus3$Coprococcus3)
fit1 <- mmer(Coprococcus3~Race, random=~BeanLine+Plate, rcov=~units, data=df_Coprococcus3,tolparinv = 7)
summary(fit1)
Coprococcus3BLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(Coprococcus3BLUPS)<-df770meds$BeanLine
#order phenotype to match genotype
ordA<- match(GD6$`rownames(GD2.2)`,rownames(Coprococcus3BLUPS))
Coprococcus3BLUPS.1 <- as.data.frame(Coprococcus3BLUPS[ordA,])
rownames(Coprococcus3BLUPS.1)<-GD6$`rownames(GD2.2)`
names(Coprococcus3BLUPS.1)<-"Coprococcus3"
#export csv for mapping:
write.csv(Coprococcus3BLUPS.1,"./Analysis/770/Phenotypes/BLUPS/MADP_770_Coprococcus3_BLUPS",quote = FALSE,row.names=TRUE)

##########################
#PCAs
df770_T.2<-df770_T.1
df770_T.2$BeanLine<-NULL
PCA<- prcomp(df770_T.2)
ev<-as.data.frame(PCA$x)


#PC1
df_PC1<-as.data.frame(ev$PC1)
#df_PC1<-log(df_PC1,2)
df_PC1$BeanLine<-df1.1$Genotype
df_PC1$idd<-df1.1$Genotype
#df_PC1$Batch<-as.character(df1$DialysisBatch)
df_PC1$Plate<-as.character(df1.1$Plate)
df_PC1$Race<-(df1.1$BeanRace)
names(df_PC1)<-c("PC1","BeanLine","idd","Plate","Race")
PC1_out<-boxplot(df_PC1$PC1,range=4)$out
df_PC1<-df_PC1[-which(df_PC1$PC1 %in% PC1_out),]
hist(df_PC1$PC1)
fit1 <- mmer(PC1~Race, random=~BeanLine+Plate, rcov=~units, data=df_PC1,tolparinv = 7)
summary(fit1)
PC1BLUPS<-as.data.frame(randef(fit1)$BeanLine)
rownames(PC1BLUPS)<-df770meds$BeanLine
#order phenotype to match genotype
ordA<- match(GD6$`rownames(GD2.2)`,rownames(PC1BLUPS))
PC1BLUPS.1 <- as.data.frame(PC1BLUPS[ordA,])
rownames(PC1BLUPS.1)<-GD6$`rownames(GD2.2)`
names(PC1BLUPS.1)<-"PC1"
#export csv for mapping:
write.csv(PC1BLUPS.1,"./Analysis/770/Phenotypes/BLUPS/MADP_770_PC1_BLUPS",quote = FALSE,row.names=TRUE)

