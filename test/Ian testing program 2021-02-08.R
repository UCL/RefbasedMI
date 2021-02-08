# Ian's testing program for Rmimix
# Main methods - J2R, CIR, CR, MAR - after discontinuation

# Install mimix if required
if(!require(mimix)) {
  if(!require(devtools)) install.packages('devtools') 
  library(devtools) 
  install_github("UCL/mimix")
}

# Load mimix
library(mimix)
packageVersion("mimix")
library(mice)
packageVersion("mice")

setwd("C:/ado/ian/test_Rmimix")

# Open data
load("C:/ado/ian/Rmimix/data/asthma.RData")
asthma$treat[1:200]<-3 # creates a 3rd arm
asthma$base2 <- asthma$base^2 # creates a 2nd covariate
asthma$fev<-asthma$fev*1000
head(asthma)

# explore data (using skills learned from Michelle)
library(tidyverse)
asthma %>% count(treat)
asthma %>% filter(!is.na(fev)) %>% 
  group_by(treat, time) %>% 
  summarise(n=n(), fevmean=mean(fev), fevsd=sd(fev))

# J2R
impJ2R1 <- mimix(data="asthma",
  covar=c("base","base2"),
  depvar="fev",
  treatvar="treat",
  idvar="id",
  timevar="time",
  M=2,
  reference=1,
  method="J2R",
  seed=101,
  prior="jeffreys",
  burnin=1000,
  bbetween=NULL,
  methodvar=NULL
)
J2R1<-impJ2R1 %>% filter(treat==1) %>% filter(.imp>0) %>% filter(patt>0)

# CR
impCR1 <- mimix(data="asthma",
               covar=c("base","base2"),
               depvar="fev",
               treatvar="treat",
               idvar="id",
               timevar="time",
               M=2,
               reference=1,
               method="CR",
               seed=101,
               prior="jeffreys",
               burnin=1000,
               bbetween=NULL,
               methodvar=NULL
)
CR1<-impCR1 %>% filter(treat==1) %>% filter(.imp>0) %>% filter(patt>0)
max(abs(J2R1-CR1))

# CIR
impCIR1 <- mimix(data="asthma",
               covar=c("base","base2"),
               depvar="fev",
               treatvar="treat",
               idvar="id",
               timevar="time",
               M=2,
               reference=1,
               method="CIR",
               seed=101,
               prior="jeffreys",
               burnin=1000,
               bbetween=NULL,
               methodvar=NULL
)
CIR1<-impCIR1 %>% filter(treat==1) %>% filter(.imp>0) %>% filter(patt>0)
max(abs(J2R1-CIR1))

# same again with reference = active

# J2R
impJ2R2 <- mimix(data="asthma",
                covar=c("base","base2"),
                depvar="fev",
                treatvar="treat",
                idvar="id",
                timevar="time",
                M=2,
                reference=2,
                method="J2R",
                seed=101,
                prior="jeffreys",
                burnin=1000,
                bbetween=NULL,
                methodvar=NULL
)

# CR
impCR2 <- mimix(data="asthma",
               covar=c("base","base2"),
               depvar="fev",
               treatvar="treat",
               idvar="id",
               timevar="time",
               M=2,
               reference=2,
               method="CR",
               seed=101,
               prior="jeffreys",
               burnin=1000,
               bbetween=NULL,
               methodvar=NULL
)

# CIR
impCIR2 <- mimix(data="asthma",
                covar=c("base","base2"),
                depvar="fev",
                treatvar="treat",
                idvar="id",
                timevar="time",
                M=2,
                reference=2,
                method="CIR",
                seed=101,
                prior="jeffreys",
                burnin=1000,
                bbetween=NULL,
                methodvar=NULL
)

J2R1<-impJ2R1 %>% filter(treat==1) %>% filter(.imp>0) %>% filter(patt>0)
CR1<-impCR1 %>% filter(treat==1) %>% filter(.imp>0) %>% filter(patt>0)
CIR1<-impCIR1 %>% filter(treat==1) %>% filter(.imp>0) %>% filter(patt>0)
max(abs(J2R1-CR1))
max(abs(J2R1-CIR1))

## show similar results (imputed values and treatment effect) when used with a different seed 

# CIR
impCIR2A <- mimix(data="asthma",
                 covar=c("base","base2"),
                 depvar="fev",
                 treatvar="treat",
                 idvar="id",
                 timevar="time",
                 M=100,
                 reference=2,
                 method="CIR",
                 seed=1037,
                 prior="jeffreys",
                 burnin=1000,
                 bbetween=NULL,
                 methodvar=NULL
)

# CIR
impCIR2B <- mimix(data="asthma",
                 covar=c("base","base2"),
                 depvar="fev",
                 treatvar="treat",
                 idvar="id",
                 timevar="time",
                 M=100,
                 reference=2,
                 method="CIR",
                 seed=4501,
                 prior="jeffreys",
                 burnin=1000,
                 bbetween=NULL,
                 methodvar=NULL
)

fitA <- pool(with(as.mids(impCIR2A), lm(fev.12~as.factor(treat)+base+base2)))
summary(fitA)
sqrt(fitA$pooled$b/fitA$pooled$m)

fitB <- pool(with(as.mids(impCIR2B), lm(fev.12~as.factor(treat)+base+base2)))
summary(fitB)
sqrt(fitB$pooled$b/fitB$pooled$m)



# interim missings
intJ2R1 <- impJ2R1 %>% filter(.id==5051|.id==5115|.id==5333) %>% filter(.imp==1)
intCIR2 <- impCIR2 %>% filter(.id==5051|.id==5115|.id==5333) %>% filter(.imp==1)
intJ2R1
intCIR2
# 5051 at time 8, 12 are post-discontinuation -> should differ
# all others are obs/interim -> should agree



# DELTA
# CIR
impCIRDELTA <- mimix(data="asthma",
                 covar=c("base","base2"),
                 depvar="fev",
                 treatvar="treat",
                 idvar="id",
                 timevar="time",
                 M=2,
                 reference=2,
                 method="CIR",
                 seed=101,
                 prior="jeffreys",
                 burnin=1000,
                 bbetween=NULL,
                 methodvar=NULL,
                 delta=c(1,2,3,4),dlag=c(1,2,3,4)
)
# dropout after first visit should have delta's = cum(1*2,2*3,3*4)=2,8,20
# dropout after 2nd visit should have delta's = cum(1*3,2*4)=3,11
# dropout after 3rd visit should have delta = 1*4=4
# CORRECT!

frame=cbind(impCIR2,impCIRDELTA[,4:7]-impCIR2[,4:7])