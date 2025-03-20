library(sommer)
library(bigmemory)
library(biganalytics)


df1 <- read.csv("MADP/SB_5_8_subset_v2.csv")

df3 <- subset(df2, Subject=="S765")
df4 <- subset(df2, Subject=="S766")

sommer_geth2 <- function(genus){
  f <- formula(paste0(genus, ' ~ 1'))
  fit <- mmer(f, random= ~BeanLine+Subject+Subject:BeanLine+Subject:Gender+Subject:SujectOrigin+BeanLine:MarketClass+BeanLine:BeanRace+vs(spl2D(Row,Column)), data=df1)
  ###vc <- summary(fit)$varcomp
  n.env <- length(levels(df1$Subject))
  out <- pin(fit, formula(paste0(genus, ' ~ V1/ (V1 + V3/n.env + V8/(3*n.env))'))) 
  return(out)
}



fit1 <- mmer(Bacteroides ~1, random = ~BeanLine+Subject+Subject:BeanLine+Subject:Gender+Subject:SujectOrigin+BeanLine:MarketClass+BeanLine:BeanRace+vs(spl2D(Row,Column)), data=df2)
summary(fit1)$varcomp
vc<-summary(fit1)$varcomp
n.env <- length(levels(df1$Subject))
pin(fit1, h2 ~ V1/ (V1 + V3/n.env + V8/(3*n.env)))

#*Rare Taxa (only 1 or 2 subjects)
sommer_geth2("Bifidobacterium")
sommer_geth2("Biff_rel")
sommer_geth2("Akkermansia")
sommer_geth2("Collinsella")
sommer_geth2("Prevotella")
hist(df1$Bacteroides)
sommer_geth2("Bacteroides")

sommer_geth2("Enterococcus")
sommer_geth2("Streptococcus")
sommer_geth2("Clostridium")
sommer_geth2("Blautia")

##############################
#h2 calc based on relative abundance
fit1 <- mmer(Bifidobacterium~1, random= ~BeanLine+Subject+Subject:Gender+Subject:SujectOrigin+vs(spl2D(Row,Column)), data=df2) 
summary(fit1)$varcomp
###vc <- summary(fit)$varcomp


df2 <- read.csv("MADP/SB_5_8_Relabd.csv")

sommer_geth2_ra <- function(genus){
  f <- formula(paste0(genus, ' ~ 1'))
  fit <- mmer(f, random= ~BeanLine+BeanRace+Subject+Subject:Gender+Subject:SujectOrigin+Subject:BeanEater+vs(spl2D(Row,Column)), data=df2)
  vc <- summary(fit)
  n.env <- length(levels(df2$Subject))
  out <- pin(fit, formula(paste0(genus, ' ~ (V1+V2)/(V1+V2+ V3/n.env + V8/(3*n.env))'))) 
  return(vc)
}

sommer_geth2_ra7 <- function(genus){
  f <- formula(paste0(genus, ' ~ 1'))
  fit <- mmer(f, random= ~BeanLine, data=df4)
  vc <- summary(fit)
  n.env <- length(levels(df4$Subject))
  out <- pin(fit, formula(paste0(genus, ' ~ V1/(V1+V2/3)'))) 
  return(out)
}

fit2 <- mmer(Odoribacter~1, ~BeanLine, data=df3)
hist(df3$Clostridium)

hist(df2$Bifidobacterium)
sommer_geth2_ra6("Bifidobacterium")
sommer_geth2_ra6("Odoribacter")
sommer_geth2_ra6("Clostridium")
sommer_geth2_ra6("Blautia")
sommer_geth2_ra6("Dorea")
sommer_geth2_ra6("Roseburia")
sommer_geth2_ra6("Ruminococcus")

hist(df2$Faecalibacterium)
sommer_geth2_ra("Faecalibacterium")

hist(df4$Akkermansia)
 sommer_geth2_ra7("Akkermansia")

####Matrix "singular" when using BeanRace in mmer
sommer_geth2_ra3("Sutterella")
sommer_geth2_ra3("Escherichia")
sommer_geth2_ra3("Coprococcus")
sommer_geth2_ra3("Parabacteroides")

hist(df2$Sutterella)

####incompatible matrix dimension?
### Bacteroides won't accept spacial dimensions?
sommer_geth2_ra4("Bacteroides")

hist(df2$Phascolarctobacterium)
sommer_geth2_ra5("Phascolarctobacterium")

sommer_geth2_ra3 <- function(genus){
  f <- formula(paste0(genus, ' ~ 1'))
  fit <- mmer(f, random= ~BeanLine+Subject+Subject:Gender+Subject:SujectOrigin+Subject:BeanEater+vs(spl2D(Row,Column)), data=df2)
  vc <- summary(fit)$varcomp
  n.env <- length(levels(df2$Subject))
  out <- pin(fit, formula(paste0(genus, ' ~ V1/(V1+ V2/n.env + V7/(3*n.env))'))) 
  return(out)
}

sommer_geth2_ra4 <- function(genus){
  f <- formula(paste0(genus, ' ~ 1'))
  fit <- mmer(f, random= ~BeanLine+Subject+Subject:Gender+Subject:SujectOrigin+Subject:BeanEater, data=df2)
  vc <- summary(fit)$varcomp
  n.env <- length(levels(df2$Subject))
  out <- pin(fit, formula(paste0(genus, ' ~ V1/(V1+ V2/n.env + V6/(3*n.env))'))) 
  return(out)
}


sommer_geth2_ra5 <- function(genus){
  f <- formula(paste0(genus, ' ~ 1'))
  fit <- mmer(f, random= ~BeanLine+Subject+Gender+SujectOrigin+BeanEater+vs(spl2D(Row,Column)), data=df2)
  vc <- summary(fit)$varcomp
  n.env <- length(levels(df2$Subject))
  out <- pin(fit, formula(paste0(genus, ' ~ V1/(V1+ V2/n.env + V7/(3*n.env))'))) 
  return(vc)
}