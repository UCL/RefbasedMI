#####################################################################
# Main testing program for RefBasedMI
# 31/3/2023: assumed to run in "test" directory
# 6/3/2023: increased permitted error from 1E-12 to 1E-8; corrected direction of causal test.
# 7/7/2022: corrected the 11/12/13 and 2/4/6 tests
# 14/6/2022: improved check of comparable results; improved tests of 11/12/13 and 2/4/6; added overall check of any tests being missed
# 29/4/2022: added test of treat=11/12/13 and 2/4/6
# 20/4/2022: updated linear model analysis
# 4/1/2022: changed .id to id
# 25/10/2021: added check that data imputed with no outcomes or baselines
# 22/10/2021: corrected sortorder test
# Revised 16jul2021
# Revised 25may2021
# Revised 8feb2021
#####################################################################

# Note date of testing
print(date())

# Install RefBasedMI - do this each time to run latest version
if(!require(devtools)) install.packages('devtools') 
library(devtools) 
if("package:RefBasedMI" %in% search()) detach("package:RefBasedMI", unload=TRUE) 
install_github("UCL/RefBasedMI",ref="dev",force=TRUE)

if(!require(tidyverse)) install.packages('tidyverse') 
  
# Load RefBasedMI
library(RefBasedMI)
packageVersion("RefBasedMI")
library(mice)
packageVersion("mice")

# Open and modify data
load("../data/asthma.RData")
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
asthma %>% filter(is.na(fev)) %>% 
  group_by(treat, time) %>% 
  summarise(n=n())
asthma %>% 
  group_by(treat, time) %>% 
  summarise(n=n())

# simple regression
lm(fev~factor(treat)+base+base2, data=asthma)

# Set up test results holder
start=TRUE
test=as.data.frame(start)

#####################################################################
# Main methods - J2R, CIR, CR, MAR - after discontinuation
#####################################################################

# MAR
print("Main methods - MAR")
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

# Check data are useable by MI
micefit <- with(data = as.mids(impMAR), lm(fev ~ factor(treat) + base, subset=(time==12)))
summary(pool(micefit))

# Check results are comparable - REQUIRES VISUAL CHECK
lmfit <- lm(fev ~ factor(treat) + base, data=asthma, subset=(time==12))
summary(lmfit)
# and that raw and imputed data have comparable structure
str(asthma)
str(impMAR)

# J2R
print("Main methods - J2R")
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
print("Main methods - CR")
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
print("Main methods - CIR")
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
test$CIR1eqMAR <- max(abs(CIR1-MAR1))<1E-8


#####################################################################
# same again with reference = active
#####################################################################

# J2R
print("J2R with reference = active")
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
print("CR with reference = active")
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
print("CIR with reference = active")
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
test$CIR2eqMAR <- max(abs(CIR2-MAR2))<1E-8


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
print("CIR with M=100, run 1")
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
print("CIR with M=100, run 2")
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
intJ2R1 <- impJ2R1 %>% filter(id==5051|id==5115|id==5333) %>% 
  filter(.imp==1) %>% arrange(id,time)
intCIR2 <- impCIR2 %>% filter(id==5051|id==5115|id==5333) %>% 
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
print("Delta-adjustment")
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
  inner_join(impCIRDELTA,by=c("id","time",".imp","base","base2","treat")) %>% 
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
print("Causal, original data")
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
print("Causal, modified data")
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

test$causal <- max(abs(err))<1E-8

#####################################################################
# DOES THE PROGRAM IMPUTE IF WE HAVE NO OBSERVED DATA?
# Added 25oct2021
#####################################################################

print("DOES THE PROGRAM IMPUTE IF WE HAVE NO OBSERVED DATA?")

asthma2 <- asthma
asthma2[asthma$id==5001,"fev"] <- NA

impnooutcomes <- RefBasedMI(data=asthma2,
                   depvar=fev,
                   treatvar=treat,
                   idvar=id,
                   timevar=time,
                   M=2,
                   method="J2R",
                   reference=1,
                   seed=101,
                   prior="jeffreys",
                   burnin=1000,
                   bbetween=NULL,
                   methodvar=NULL,
                   delta=c(1,2,3,4), dlag=c(.1,.1,.1,.05)
)

test$nooutcomes <- 
  (impnooutcomes %>% filter(id==5001) %>% filter(.imp==0) %>% summarise(nmiss=sum(is.na(fev)))) > 
  (impnooutcomes %>% filter(id==5001) %>% filter(.imp>0) %>% summarise(nmiss=sum(is.na(fev)))) 



#####################################################################
# DOES PROGRAM WORK AND REPORT SAME RESULTS IF TREATMENT IS 11/12/13?
#####################################################################

# J2R
asthmaten<-asthma
asthmaten$treat=asthma$treat+10
print("DOES PROGRAM WORK AND REPORT SAME RESULTS IF TREATMENT IS 11/12/13?")
impJ2R1ten <- RefBasedMI(data=asthmaten,
                         covar=c(base,base2),
                         depvar=fev,
                         treatvar=treat,
                         idvar=id,
                         timevar=time,
                         M=2,
                         reference=11,
                         method="J2R",
                         seed=101,
                         prior="jeffreys",
                         burnin=1000,
                         bbetween=NULL,
                         methodvar=NULL
)
# test correct #rows 
test$treat111213 <- nrow(impJ2R1ten) == nrow(impJ2R1)
# ... and if so test correct results
if(test$treat111213) 
  test$treat111213 <- 
  sum(impJ2R1$treat!=impJ2R1ten$treat-10)==0 & 
  sum(impJ2R1$fev[impJ2R1$.imp>0] != impJ2R1ten$fev[impJ2R1ten$.imp>0])==0

#####################################################################
# DOES PROGRAM WORK AND REPORT SAME RESULTS IF TREATMENT IS 2/4/6?
#####################################################################

# J2R
asthmatwo<-asthma
asthmatwo$treat=asthma$treat*2
print("DOES PROGRAM WORK AND REPORT SAME RESULTS IF TREATMENT IS 2/4/6?")
impJ2R1two <- RefBasedMI(data=asthmatwo,
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
# test correct #rows 
test$treat246 <- nrow(impJ2R1two) == nrow(impJ2R1)
# ... and if so test correct results
if(test$treat246) 
  test$treat246 <- 
  sum(impJ2R1$treat!=impJ2R1two$treat/2)==0 & 
  sum(impJ2R1$fev[impJ2R1$.imp>0] != impJ2R1two$fev[impJ2R1two$.imp>0])==0


#####################################################################
# END OF TESTS: NOW PRINT SUMMARY
#####################################################################

test$COMPLETE <- ncol(test)==16
print(t(test))
