# Ian's testing program for Rmimix: to compare with Stata
# 8jul2020: 3 arms, 2 covariates, added causal vs CIR

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
impCR1 <- mimix(data="asthma",
               covar=c("base","base2"),
               depvar="fev",
               treatvar="treat",
               idvar="id",
               timevar="time",
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

# same again with reference = active

# J2R
impJ2R2 <- mimix(data="asthma",
                covar=c("base","base2"),
                depvar="fev",
                treatvar="treat",
                idvar="id",
                timevar="time",
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
impCR2 <- mimix(data="asthma",
               covar=c("base","base2"),
               depvar="fev",
               treatvar="treat",
               idvar="id",
               timevar="time",
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

# analyse
fit <- with(as.mids(impJ2R1), lm(fev.12~as.factor(treat)+base+base2))
summary(pool(fit))
fit <- with(as.mids(impCR1), lm(fev.12~as.factor(treat)+base+base2))
summary(pool(fit))
fit <- with(as.mids(impCIR1), lm(fev.12~as.factor(treat)+base+base2))
summary(pool(fit))
fit <- with(as.mids(impJ2R2), lm(fev.12~as.factor(treat)+base+base2))
summary(pool(fit))
fit <- with(as.mids(impCR2), lm(fev.12~as.factor(treat)+base+base2))
summary(pool(fit))
fit <- with(as.mids(impCIR2), lm(fev.12~as.factor(treat)+base+base2))
summary(pool(fit))
