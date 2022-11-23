# Ian's testing program for RefBasedMI: to compare with Stata
# 8jul2020: 3 arms, 2 covariates, added causal vs CIR
# 23nov2022: update to RefBasedMI syntax

# Install RefBasedMI if required
if(!require(RefBasedMI)) {
  if(!require(devtools)) install.packages('devtools') 
  library(devtools) 
  install_github("UCL/RefBasedMI")
}

# Load RefBasedMI
library(RefBasedMI)
packageVersion("RefBasedMI")
library(mice)
packageVersion("mice")

setwd("C:/ado/ian/RefBasedMI")

# Open data
load("data/asthma.RData")
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

# MAR
impMAR1 <- RefBasedMI(data=asthma,
                      covar=c("base","base2"),
                      depvar=fev,
                      treatvar=treat,
                      idvar=id,
                      timevar=time,
                      M=100,
                      method="MAR",
                      seed=101,
                      prior="jeffreys",
                      burnin=1000,
                      bbetween=NULL,
                      methodvar=NULL
)

# J2R
impJ2R1 <- RefBasedMI(data=asthma,
                      covar=c("base","base2"),
                      depvar=fev,
                      treatvar=treat,
                      idvar=id,
                      timevar=time,
                      M=100,
                      reference=1,
                      method="J2R",
                      seed=101,
                      prior="jeffreys",
                      burnin=1000,
                      bbetween=NULL,
                      methodvar=NULL
)

# CR
impCR1 <- RefBasedMI(data=asthma,
               covar=c("base","base2"),
               depvar=fev,
               treatvar=treat,
               idvar=id,
               timevar=time,
               M=100,
               reference=1,
               method="CR",
               seed=101,
               prior="jeffreys",
               burnin=1000,
               bbetween=NULL,
               methodvar=NULL
)

# CIR
impCIR1 <- RefBasedMI(data=asthma,
               covar=c("base","base2"),
               depvar=fev,
               treatvar=treat,
               idvar=id,
               timevar=time,
               M=100,
               reference=1,
               method="CIR",
               seed=101,
               prior="jeffreys",
               burnin=1000,
               bbetween=NULL,
               methodvar=NULL
)

# same again with reference = active

# J2R
impJ2R2 <- RefBasedMI(data=asthma,
                covar=c("base","base2"),
                depvar=fev,
                treatvar=treat,
                idvar=id,
                timevar=time,
                M=100,
                reference=2,
                method="J2R",
                seed=101,
                prior="jeffreys",
                burnin=1000,
                bbetween=NULL,
                methodvar=NULL
)

# CR
impCR2 <- RefBasedMI(data=asthma,
               covar=c("base","base2"),
               depvar=fev,
               treatvar=treat,
               idvar=id,
               timevar=time,
               M=100,
               reference=2,
               method="CR",
               seed=101,
               prior="jeffreys",
               burnin=1000,
               bbetween=NULL,
               methodvar=NULL
)

# CIR
impCIR2 <- RefBasedMI(data=asthma,
                covar=c("base","base2"),
                depvar=fev,
                treatvar=treat,
                idvar=id,
                timevar=time,
                M=100,
                reference=2,
                method="CIR",
                seed=101,
                prior="jeffreys",
                burnin=1000,
                bbetween=NULL,
                methodvar=NULL
)

# analyse
CC <- lm(data=asthma, formula=fev~as.factor(treat)+base+base2, subset=(time==12))
CC2 <- summary(CC)$coefficients
colnames(CC2)<-c("Estimate","StdError","t","P")

fit <- with(as.mids(impMAR1), lm(fev~as.factor(treat)+base+base2), subset=(time==12))
impMAR1res <- summary(pool(fit))

fit <- with(as.mids(impJ2R1), lm(fev~as.factor(treat)+base+base2), subset=(time==12))
impJ2R1res <- summary(pool(fit))
fit <- with(as.mids(impCR1), lm(fev~as.factor(treat)+base+base2), subset=(time==12))
impCR1res <- summary(pool(fit))
fit <- with(as.mids(impCIR1), lm(fev~as.factor(treat)+base+base2), subset=(time==12))
impCIR1res <- summary(pool(fit))
fit <- with(as.mids(impJ2R2), lm(fev~as.factor(treat)+base+base2), subset=(time==12))
impJ2R2res <- summary(pool(fit))
fit <- with(as.mids(impCR2), lm(fev~as.factor(treat)+base+base2), subset=(time==12))
impCR2res <- summary(pool(fit))
fit <- with(as.mids(impCIR2), lm(fev~as.factor(treat)+base+base2), subset=(time==12))
impCIR2res <- summary(pool(fit))


# output
write.dta(as.data.frame(CC2),file="test/RvsStata/R_CC.dta")
write.dta(impMAR1res,file="test/RvsStata/R_MAR_ref1.dta")
write.dta(impMAR1res,file="test/RvsStata/R_MAR_ref2.dta")
write.dta(impJ2R1res,file="test/RvsStata/R_J2R_ref1.dta")
write.dta(impCR1res, file="test/RvsStata/R_CR_ref1.dta")
write.dta(impCIR1res,file="test/RvsStata/R_CIR_ref1.dta")
write.dta(impJ2R2res,file="test/RvsStata/R_J2R_ref2.dta")
write.dta(impCR2res, file="test/RvsStata/R_CR_ref2.dta")
write.dta(impCIR2res,file="test/RvsStata/R_CIR_ref2.dta")
