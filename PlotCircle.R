setwd("~/MADP/Analysis/ManuscriptFigures/")

Faecali<-read.csv("FaecaliSNPs.csv")
Faecali1<-cbind(Faecali[1:3],Faecali[5])
mels<-Faecali$SNP[c(7742:7867,32597:32763,32764:32861,44952:45081,45083:45181,55540:56487,56489:57109,62484:64313,64315:64412,66534:68114,72230:72274)]
mels2<-Faecali$SNP[c(7742:7867,32597:32861,44952:45181,55540:57109,62484:64412,66534:68114,72230:72274)]

pcent<-Faecali$SNP[c(1513:7462,9858:16778,23260:29576,34706:40719,43436:54372,55541:59730,64382:70975,73743:86662,91446:92195,98227:109214,112779:128358)]
library(rMVP)

faecali<-read.csv("FaecaliSNPs.csv")
#mels<-faecali$SNP[c(7743,7866,32596,32860,44951,45180,55540,57109,62484,64412,66534,68114,72231,72273)]

MVP.Report(Faecali1,plot.type="c",r=5,col=c("grey36","grey72"),pch=3,multracks=TRUE,chr.labels=paste("C",c(1:11),sep=""),
           threshold=c(3.18e-7),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("blue1"),signal.line=0,signal.col=c("blue1"),chr.den.col=c("blue","yellow","chocolate2"),cir.legend=TRUE,
           bin.size=2.5e6,outward=TRUE,file.type="jpg",memo="",band=0.4,dpi=300,highlight=mels2,highlight.col = ("orangered2"),highlight.cex = 1,cir.band=0.3)

#MVP.Report(Faecali,plot.type="m",r=4.5,col=c("blue","green","red"),pch=3,multracks=TRUE,chr.labels=paste("C",c(1:11),sep=""),
          # threshold=c(3.18e-7),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("blue1"),signal.line=0,signal.col=c("blue1"),chr.den.col=c("blue","yellow","chocolate2"),cir.legend=TRUE,
          # bin.size=2.5e6,outward=TRUE,file.type="jpg",memo="",band=0.5,dpi=300,highlight=mels,highlight.col = ("orangered2"),highlight.cex = 1,cir.band=0.2)


#MVP.Report(Faecali1,plot.type="c",r=0.5,col=c("grey30","grey60"),pch=3,multracks=TRUE,chr.labels=paste("C",c(1:11),sep=""),
#           threshold=c(3.18e-7),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red",
#          "blue"),signal.line=1,signal.col=c("red"),chr.den.col=c("darkgreen","yellow","red"),cir.legend=TRUE,
#           bin.size=1.5e6,outward=TRUE,file.type="jpg",memo="",band=0.5,dpi=300,highlight=NULL,cir.band=0.2)


#Plot shapes:

MelD<-read.csv("Mel_input.csv")

ggplot(MelD, aes(x=x,y=y))+geom_point(aes(fill=Family),colour="black",pch=21,size=18)+theme_classic()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        text=element_text(size=14,family = "sans"),plot.margin=unit(c(1,0,0,0),"cm"))+
  xlim(-3,3)+ylim(0,4)
  


####Make a plot for the book
Faecali<-read.csv("FaecaliSNPs.csv")
Faecali1<-cbind(Faecali[1:3],Faecali[7])
SNPs<-as.numeric(nrow(Faecali1))
thresh<-0.05/SNPs
MVP.Report(Faecali1,plot.type=c("m"),col=c("grey36","grey72"),chr.labels=paste("C",c(1:11),sep=""),
           threshold=thresh,threshold.lty=2,threshold.col=c("black"),signal.col=c("blue1"),chr.den.col=NULL,
           file.type="jpg",memo="",band=0.4,dpi=600)

 