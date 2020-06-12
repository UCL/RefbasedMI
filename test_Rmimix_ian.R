# Ian's testing program for Rmimix
# 3jun2020

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
head(asthma)

# J2R
impJ2R <- mimix(data="asthma",
  covar="base",
  depvar="fev",
  treatvar="treat",
  idvar="id",
  timevar="time",
  M=2,
  refer=1,
  method="J2R",
  seed=101,
  prior="jeffreys",
  burnin=1000,
  bbetween=NULL,
  methodvar=NULL
)

# analyse
impdata <- as.mids(impJ2R)
fit <- with(impdata, lm(fev.12~treat+base))
summary(pool(fit))

# J2R + delta
impJ2Rdelta <- mimix(data="asthma",
  covar="base",
  depvar="fev",
  treatvar="treat",
  idvar="id",
  timevar="time",
  M=2,
  refer=1,
  method="J2R",
  seed=101,
  prior="jeffreys",
  burnin=1000,
  bbetween=NULL,
  methodvar=NULL,
  delta=c(2,2,2,2),
  dlag=c(1,-.5,-.25,-.125)
)
sum(is.na(subset(impJ2Rdelta,.imp>0)))

# evaluate delta
compare=cbind(impJ2R,impJ2Rdelta-impJ2R)
subset(compare,.id==5017)
subset(compare,.id==5030)
subset(compare,.id==5051)
subset(compare,.id==5115)

# Causal
impcausal1 <- mimix(data="asthma",
                    covar=c("base"),
                    depvar="fev",
                    treatvar="treat",
                    idvar="id",
                    timevar="time",
                    M=2,
                    refer=1,
                    method="causal",
                    seed=101,
                    prior="jeffreys",
                    burnin=1000,
                    K0=0, 
                    K1=0
)
# verify causal with k0=k1=0 is J2R
sum(impcausal1!=impJ2R,na.rm=TRUE)

# Causal
impcausal2 <- mimix(data="asthma",
   covar=c("base"),
   depvar="fev",
   treatvar="treat",
   idvar="id",
   timevar="time",
   M=2,
   refer=1,
   method="causal",
   seed=101,
   prior="jeffreys",
   burnin=1000,
   K0=1,
   K1=1
)

sum(is.na(subset(impcausal2,.imp>0)))

# evaluate delta
compare=cbind(impcausal2,impcausal2-impJ2R)
subset(impcausal2,.id==5017)
subset(impJ2R,.id==5017)
subset(impcausal2,.id==5115)
subset(impJ2R,.id==5115)

