library(ggplot2)
library(dplyr)

setwd("~/MADP")

data1<-read.csv("CulmulativeGWAS.FarmCPU.csv")

Color<-cbind(data1[2:3],"Color",data1$RGB_eigenvector.FarmCPU)
names(Color)<-c("Chromosome","BP","Phenotype","Pvalue")
Color<-Color[-which(Color$Pvalue>0.0001),]

PC1<-cbind(data1[2:3],"PC1",data1$PC1.FarmCPU)
names(PC1)<-c("Chromosome","BP","Phenotype","Pvalue")
PC1<-PC1[-which(PC1$Pvalue>0.0001),]

PC2<-cbind(data1[2:3],"PC2",data1$PC2.FarmCPU)
names(PC2)<-c("Chromosome","BP","Phenotype","Pvalue")
PC2<-PC2[-which(PC2$Pvalue>0.0001),]

Bifidobacterium<-cbind(data1[2:3],"Bifidobacterium",data1$Bifidobacterium.FarmCPU)
names(Bifidobacterium)<-c("Chromosome","BP","Phenotype","Pvalue")
Bifidobacterium<-Bifidobacterium[-which(Bifidobacterium$Pvalue>0.0001),]

data2<-rbind(Color,PC1,PC2,Bifidobacterium)

don<-data2 %>%
  group_by(Chromosome) %>%
  summarise(chr_len=max(BP)) #%>%
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(data2, ., by=c("Chromosome"="Chromosome")) %>%
  arrange(Chromosome, BP) %>%
  mutate( BPcum=BP+tot)

axisdf=don %>% group_by(Chromosome) %>% summarise(center=(max(BPcum)+min(BPcum))/2)

ggplot(don,aes(x=BPcum,y=-log10(Pvalue)))+
  geom_point( aes(color=as.factor(Phenotype)), alpha=0.9, size=2.5) +
  scale_x_continuous( label = axisdf$Chromosome, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, .2) ) +
  scale_color_manual(values = rep(c("Darkgreen","Black","Purple","Blue"), 22 ))+theme_classic()+
  theme(panel.grid.mor.x = element_line())

  theme_bw()

  theme(legend.position="none",panel.border = element_blank(),
        panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
