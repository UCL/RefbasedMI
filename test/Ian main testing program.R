#####################################################################
# Ian's main testing program for Rmimix
# Revised 25may2021
# Revised 8feb2021
#####################################################################


# Install mimix if required
if(!require(mimix)) {
  if(!require(devtools)) install.packages('devtools') 
  library(devtools) 
  install_github("UCL/mimix")
}

# Load mimix
library(RefBasedMI)
packageVersion("RefBasedMI")
library(mice)
packageVersion("mice")

setwd("C:/ado/ian/Rmimix/test")

# Open and modify data
load("C:/ado/ian/Rmimix/data/asthma.RData")
asthma$treat[1:200]<-3 # creates a 3rd arm
asthma$base2 <- asthma$base^2 # creates a 2nd covariate
asthma$fev<-asthma$fev*1000
head(asthma)

# Explore data 
library(tidyverse)
asthma %>% count(treat)
asthma %>% filter(!is.na(fev)) %>% 
  group_by(treat, time) %>% 
  summarise(n=n(), fevmean=mean(fev), fevsd=sd(fev))

#####################################################################
# Main methods - J2R, CIR, CR, MAR - after discontinuation
#####################################################################

# J2R
impJ2R1 <- RefBasedMI(data=asthma,
  covar=c("base","base2"),
  depvar=fev,
  treatvar=treat,
  idvar=id,
  timevar=time,
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
# TEST
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



#####################################################################
# interim missings and deltas
#####################################################################
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




#####################################################################
# TEST CAUSAL ROUTINES (was Ian_test_causal.R)
#####################################################################
# SETTINGS FOR THIS TEST
K0=0.1
K1=0.9
tweak=3

# MODIFY ASTHMA DATA
asthma2=asthma
asthma2$fev=(asthma2$fev+tweak*(asthma2$time==4)*(asthma2$treat==2))

library(tidyverse)
# ggplot(data=asthma, aes(x=time,y=fev,group=id,colour=treat))+
#   geom_line()

# J2R - just to get data in right shape and to set up useful variables
orig <- mimix(data="asthma",
  covar=c("base"),
  depvar="fev",
  treatvar="treat",
  idvar="id",
  timevar="time",
  M=1,
  reference=1,
  method="J2R",
  seed=101,
  prior="jeffreys",
  burnin=1000,
  bbetween=NULL,
  methodvar=NULL
) %>% filter(.imp==0) 
nonmon = (is.na(orig$fev.2)&!is.na(orig$fev.4)| 
                  is.na(orig$fev.4)&!is.na(orig$fev.8) | 
                  is.na(orig$fev.8)&!is.na(orig$fev.12))
pattern = (1000*!is.na(orig$fev.2)) + (100*!is.na(orig$fev.4))+ (10*!is.na(orig$fev.8)) + (!is.na(orig$fev.12))
orig[nonmon,]

# causal, original data
causal1 <- mimix(data="asthma",
                 covar=c("base"),
                 depvar="fev",
                 treatvar="treat",
                 idvar="id",
                 timevar="time",
                 M=1,
                 reference=1,
                 method="causal",
                 seed=101,
                 prior="jeffreys",
                 burnin=1000,
                 bbetween=NULL,
                 methodvar=NULL,
                 K0=K0,K1=K1
) %>% filter(.imp==1) %>% select("fev.2","fev.4","fev.8","fev.12")

# causal, same K, modified data
causal1m <- mimix(data="asthma2",
                 covar=c("base"),
                 depvar="fev",
                 treatvar="treat",
                 idvar="id",
                 timevar="time",
                 M=1,
                 reference=1,
                 method="causal",
                 seed=101,
                 prior="jeffreys",
                 burnin=1000,
                 bbetween=NULL,
                 methodvar=NULL,
                 K0=K0,K1=K1
) %>% filter(.imp==1) %>% select("fev.2","fev.4","fev.8","fev.12")

# key comparison
diff=causal1m-causal1
compare = cbind(orig,pattern,diff)

orig[pattern==1100 & orig$treat==2,]
causal1[pattern==1100 & orig$treat==2,]
causal1m[pattern==1100 & orig$treat==2,]
diff[pattern==1100 & orig$treat==2,]
# each row should read 0,1,0.0625,0.00390625
diff[pattern==1000 & orig$treat==2,]
# each row should read 0,0,0,0
diff[pattern==1110 & orig$treat==2,]
# each row should read 0,1,0,0

# express as tests: check changes happen correctly
abs(diff[pattern==1100 & orig$treat==2,"fev.8"]   - tweak*K0*K1^4) < 1E-10
abs(diff[pattern==1100 & orig$treat==2,"fev.12"]  - tweak*K0*K1^8) < 1E-10
abs(diff[pattern==1110 & orig$treat==2,"fev.4"]   - tweak)     < 1E-10

# check changes don't happen where they shouldn't
abs(diff[pattern==1110 & orig$treat==2,"fev.12"]  - 0) < 1E-10
abs(diff[pattern==1000 & orig$treat==2,"fev.8"]   - 0) < 1E-10

# further tests needed for non-monotone cases

#####################################################################
# END OF TESTS
#####################################################################
