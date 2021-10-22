#####################################################################
# Ian's main testing program for RefBasedMI
# 22/10/2021: corrected sortorder test
# Revised 16jul2021
# Revised 25may2021
# Revised 8feb2021
#####################################################################


# Install mimix if required
if(!require(RefBasedMI)) {
  if(!require(devtools)) install.packages('devtools') 
  library(devtools) 
  install_github("UCL/RefBasedMI")
}

if(!require(tidyverse)) install.packages('tidyverse') 
  
# Load mimix
library(RefBasedMI)
packageVersion("RefBasedMI")
library(mice)
packageVersion("mice")

setwd("C:/ado/ian/RefBasedMI/test")

# Open and modify data
load("C:/ado/ian/RefBasedMI/data/asthma.RData")
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

# simple regression
lm(fev~factor(treat)+base+base2, data=asthma)

# Set up test results holder
start=TRUE
test=as.data.frame(start)

#####################################################################
# Main methods - J2R, CIR, CR, MAR - after discontinuation
#####################################################################

# MAR
impMAR <- RefBasedMI(data=asthma,
                      covar=c(base,base2),
                      depvar=fev,
                      treatvar=treat,
                      idvar=id,
                      timevar=time,
                      M=2,
                      method="MAR",
                      seed=101,
                      prior="jeffreys",
                      burnin=1000,
                      bbetween=NULL,
                      methodvar=NULL
)
MAR1<-impMAR %>% filter(treat==1) %>% filter(.imp>0) %>% select(-treat)
MAR2<-impMAR %>% filter(treat==2) %>% filter(.imp>0) %>% select(-treat)

# J2R
impJ2R1 <- RefBasedMI(data=asthma,
                      covar=c(base,base2),
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
J2R1<-impJ2R1 %>% filter(treat==1) %>% filter(.imp>0) %>% select(-treat)

# TEST: J2R1 should = MAR in reference group
test$J2R1eqMAR <- max(abs(J2R1-MAR1))==0

# CR
impCR1 <- RefBasedMI(data=asthma,
               covar=c(base,base2),
               depvar=fev,
               treatvar=treat,
               idvar=id,
               timevar=time,
               M=2,
               reference=1,
               method="CR",
               seed=101,
               prior="jeffreys",
               burnin=1000,
               bbetween=NULL,
               methodvar=NULL
)
CR1<-impCR1 %>% filter(treat==1) %>% filter(.imp>0) %>% select(-treat)
# TEST: CR1 should = MAR in reference group
test$CR1eqMAR <- max(abs(CR1-MAR1))==0

# CIR
impCIR1 <- RefBasedMI(data=asthma,
                     covar=c(base,base2),
                     depvar=fev,
                     treatvar=treat,
                     idvar=id,
                     timevar=time,
                     M=2,
                     reference=1,
                     method="CIR",
                     seed=101,
                     prior="jeffreys",
                     burnin=1000,
                     bbetween=NULL,
                     methodvar=NULL
)
CIR1<-impCIR1 %>% filter(treat==1) %>% filter(.imp>0) %>% select(-treat)
# TEST: CIR1 should = MAR in reference group
test$CIR1eqMAR <- max(abs(CIR1-MAR1))<1E-12


#####################################################################
# same again with reference = active
#####################################################################

# J2R
impJ2R2 <- RefBasedMI(data=asthma,
                      covar=c(base,base2),
                      depvar=fev,
                      treatvar=treat,
                      idvar=id,
                      timevar=time,
                      M=2,
                      reference=2,
                      method="J2R",
                      seed=101,
                      prior="jeffreys",
                      burnin=1000,
                      bbetween=NULL,
                      methodvar=NULL
)
J2R2<-impJ2R2 %>% filter(treat==2) %>% filter(.imp>0) %>% select(-treat)
# TEST: J2R1 should = MAR in reference group
test$J2R2eqMAR <- max(abs(J2R2-MAR2))==0

# CR
impCR2 <- RefBasedMI(data=asthma,
                     covar=c(base,base2),
                     depvar=fev,
                     treatvar=treat,
                     idvar=id,
                     timevar=time,
                     M=2,
                     reference=2,
                     method="CR",
                     seed=101,
                     prior="jeffreys",
                     burnin=1000,
                     bbetween=NULL,
                     methodvar=NULL
)
CR2<-impCR2 %>% filter(treat==2) %>% filter(.imp>0) %>% select(-treat)
# TEST: CR2 should = MAR in reference group
test$CR2eqMAR <- max(abs(CR2-MAR2))==0

# CIR
impCIR2 <- RefBasedMI(data=asthma,
                      covar=c(base,base2),
                      depvar=fev,
                      treatvar=treat,
                      idvar=id,
                      timevar=time,
                      M=2,
                      reference=2,
                      method="CIR",
                      seed=101,
                      prior="jeffreys",
                      burnin=1000,
                      bbetween=NULL,
                      methodvar=NULL
)
CIR2<-impCIR2 %>% filter(treat==2) %>% filter(.imp>0) %>% select(-treat)
# TEST: CIR2 should = MAR in reference group
test$CIR2eqMAR <- max(abs(CIR2-MAR2))<1E-12


# TEST: no missing values after imputation
test$nomiss <- sum(is.na(c(MAR1,MAR2, J2R1,J2R2,CR1,CR2,CIR1,CIR2)))==0


# TEST: same sort order
head(asthma)
head(impCIR2 %>% filter(.imp==1))
test$sortorder=sum((asthma %>% select(id,time))!=(impCIR2 %>% filter(.imp==1) %>% select(id,time)))==0


#####################################################################
## show similar results (imputed values and treatment effect) 
## when used with a different seed 
#####################################################################

# CIR
impCIR2A <- RefBasedMI(data=asthma,
                       covar=c(base,base2),
                       depvar=fev,
                       treatvar=treat,
                       idvar=id,
                       timevar=time,
                       M=100,
                       reference=1,
                       method="CIR",
                       seed=1037,
                       prior="jeffreys",
                       burnin=1000,
                       bbetween=NULL,
                       methodvar=NULL
)

# CIR
impCIR2B <- RefBasedMI(data=asthma,
                       covar=c(base,base2),
                       depvar=fev,
                       treatvar=treat,
                       idvar=id,
                       timevar=time,
                       M=100,
                       reference=1,
                       method="CIR",
                       seed=4501,
                       prior="jeffreys",
                       burnin=1000,
                       bbetween=NULL,
                       methodvar=NULL
)

fitA <- pool(with(as.mids(impCIR2A %>% filter(time==12)), 
                  lm(fev~as.factor(treat)+base+base2,subset=(time==12))
                  ))
summary(fitA)
sqrt(fitA$pooled$b/fitA$pooled$m) # Monte Carlo errors

fitB <- pool(with(as.mids(impCIR2B %>% filter(time==12)), 
                  lm(fev~as.factor(treat)+base+base2,subset=(time==12))
                  ))
summary(fitB)
sqrt(fitB$pooled$b/fitB$pooled$m) # Monte Carlo errors

# Monte Carlo Z statistics for the differences between runs
MCZ <- (fitB$pooled$estimate-fitA$pooled$estimate) / 
  sqrt(fitA$pooled$b/fitA$pooled$m+fitB$pooled$b/fitB$pooled$m)
MCZ
test$diffseed <- max(abs(MCZ))<3 & max(abs(MCZ))>0.1

#####################################################################
# interim missings should be imputed the same way by different methods
# but final missings shouldn't
#####################################################################
intJ2R1 <- impJ2R1 %>% filter(.id==5051|.id==5115|.id==5333) %>% 
  filter(.imp==1) %>% arrange(id,time)
intCIR2 <- impCIR2 %>% filter(.id==5051|.id==5115|.id==5333) %>% 
  filter(.imp==1) %>% arrange(id,time)
intorig <- asthma  %>% filter(id==5051|id==5115|id==5333) %>% arrange(id,time)
summ <- intJ2R1 %>% select(id,time)
summ$orig <- intorig$fev
summ$J2R <-intJ2R1$fev 
summ$CIR <-intCIR2$fev 
summ$diff <-intJ2R1$fev-intCIR2$fev 
summ
# 5051 at time 8, 12 are post-discontinuation -> should differ
# all others are obs/interim -> should agree
test$interim <- 
summ %>% filter((id==5051 & time>4) == (diff==0)) %>% select(diff) %>% summarise(interim=n()) ==0
test

#####################################################################
# Deltas
#####################################################################


# DELTA
# CIR
impCIRDELTA <- RefBasedMI(data=asthma,
                 covar=c(base,base2),
                 depvar=fev,
                 treatvar=treat,
                 idvar=id,
                 timevar=time,
                 M=2,
                 reference=2,
                 method="CIR",
                 seed=101,
                 prior="jeffreys",
                 burnin=1000,
                 bbetween=NULL,
                 methodvar=NULL,
                 delta=c(1,2,3,4),
                 dlag=c(1,0,0,0)
)

compare <- impCIR2 %>% 
  inner_join(impCIRDELTA,by=c("id","time",".id",".imp","base","base2","treat")) %>% 
  filter(.imp==1)
compare$diff <- compare$fev.y-compare$fev.x
compare %>% group_by(time) %>% summarise(meandiff=mean(diff))
test$deltaworks <- compare %>% summarise(deltaworks=mean(diff)) !=0



#####################################################################
# TEST CAUSAL ROUTINES (was Ian_test_causal.R)
#####################################################################
# SETTINGS FOR THIS TEST
K0=1
K1=0.9
tweak=3

# MODIFY ASTHMA DATA
asthmamod=asthma
asthmamod$fev = asthmamod$fev+tweak*(asthmamod$time==4)*(asthmamod$treat==3)

library(tidyverse)
# ggplot(data=asthma, aes(x=time,y=fev,group=id,colour=treat))+
#   geom_line()

# causal, original data
causal1 <- RefBasedMI(data=asthma,
                 covar=c("base"),
                 depvar=fev,
                 treatvar=treat,
                 idvar=id,
                 timevar=time,
                 M=1,
                 reference=1,
                 method="causal",
                 seed=101,
                 prior="jeffreys",
                 burnin=1000,
                 bbetween=NULL,
                 methodvar=NULL,
                 K0=K0,K1=K1
) %>% filter(.imp==1) 

# causal, same K, modified data
causal1mod <- RefBasedMI(data=asthmamod,
                      covar=c("base"),
                      depvar=fev,
                      treatvar=treat,
                      idvar=id,
                      timevar=time,
                      M=1,
                      reference=1,
                      method="causal",
                      seed=101,
                      prior="jeffreys",
                      burnin=1000,
                      bbetween=NULL,
                      methodvar=NULL,
                      K0=K0,K1=K1
) %>% filter(.imp==1) 

# monotone pattern OXXX: id=5017, treat=3
asthma %>% filter(id==5017)
# monotone pattern OOXX: id=5030, treat=3
asthma %>% filter(id==5030)
# monotone pattern OOOX: id=5074, treat=3
asthma %>% filter(id==5074)
# non-monotone pattern XOXX: id=5051, treat=3
asthma %>% filter(id==5051)

asthma %>% inner_join(asthmamod,by=c("id","time")) %>% filter(id==5030)

compare <- causal1 %>% 
  inner_join(causal1mod,by=c("id","time")) %>% 
  inner_join(asthma,by=c("id","time")) %>% 
  select(id,time,treat,fev,fev.x,fev.y)
compare$diff <- compare$fev.y-compare$fev.x
err=c()

# monotone pattern OXXX: id=5017, treat=3: unaffected by tweak
err=rbind(err, sum(compare %>% filter(id==5017) %>% select(diff) - 0))

# monotone pattern OOXX: id=5030, treat=3: fully affected by tweak
err=rbind(err, sum(compare %>% filter(id==5030) %>% filter(time==2) %>% select(diff) - 0))
err=rbind(err, sum(compare %>% filter(id==5030) %>% filter(time==4) %>% select(diff) - tweak))
err=rbind(err, sum(compare %>% filter(id==5030) %>% filter(time==8) %>% select(diff) - tweak*K0*K1^4))
err=rbind(err, sum(compare %>% filter(id==5030) %>% filter(time==12) %>% select(diff) - tweak*K0*K1^8))

# monotone pattern OOOX: id=5074, treat=3: imputed data unaffected by tweak since causal model relates to last obs time
err=rbind(err, sum(compare %>% filter(id==5074) %>% filter(time==4) %>% select(diff) - tweak))
err=rbind(err, sum(compare %>% filter(id==5074) %>% filter(time!=4) %>% select(diff) - 0))

# non-monotone pattern XOXX: id=5051, treat=3: fully affected by tweak
err=rbind(err, sum(compare %>% filter(id==5051) %>% filter(time==2) %>% select(diff) - 0))
err=rbind(err, sum(compare %>% filter(id==5051) %>% filter(time==4) %>% select(diff) - tweak))
err=rbind(err, sum(compare %>% filter(id==5051) %>% filter(time==8) %>% select(diff) - tweak*K0*K1^4))
err=rbind(err, sum(compare %>% filter(id==5051) %>% filter(time==12) %>% select(diff) - tweak*K0*K1^8))

test$causal <- max(abs(err))>1E-12

#####################################################################
# END OF TESTS: NOW PRINT SUMMARY
#####################################################################

t(test)