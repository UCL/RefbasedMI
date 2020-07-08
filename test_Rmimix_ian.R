# Ian's testing program for Rmimix
# 3jun2020
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
asthma$treat[1:200]<-3
asthma$base2 <- asthma$base^2
asthma$fev<-asthma$fev*1000
head(asthma)
table(asthma$treat)


# J2R
impJ2R <- mimix(data="asthma",
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

# analyse
impdata <- as.mids(impJ2R)
fit <- with(impdata, lm(fev.12~as.factor(treat)+base+base2))
summary(pool(fit))

# CR
impCR <- mimix(data="asthma",
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
impCIR <- mimix(data="asthma",
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

# J2R + delta
impJ2Rdelta <- mimix(data="asthma",
  covar=c("base","base2"),
  depvar="fev",
  treatvar="treat",
  idvar="id",
  timevar="time",
  M=2,
  reference = 1,
  method="J2R",
  seed=101,
  prior="jeffreys",
  burnin=1000,
  bbetween=NULL,
  methodvar=NULL,
  delta=c(2,5,8,11),
  dlag=c(4,3,2,1)
)
# these delta's should make imputed values increase as follows after discontinuation:
# discontinuation before first visit: 2*4, 2*4+5*3, ... = 8, 23, 39, 50
# discontinuation before 2nd visit: 5*4, 5*4+8*3, ... = 20, 44, 66
# discontinuation before 3rd visit: 8*4, 8*4+11*3 = 32, 65
# discontinuation before 4th visit: 11*4 = 44
sum(is.na(subset(impJ2Rdelta,.imp>0)))

# evaluate delta
compare=cbind(impJ2R,impJ2Rdelta-impJ2R)
subset(compare,.id==5017)
subset(compare,.id==5030)
subset(compare,.id==5051)
subset(compare,.id==5115)

# Causal
impcausal1 <- mimix(data="asthma",
                    covar=c("base","base2"),
                    depvar="fev",
                    treatvar="treat",
                    idvar="id",
                    timevar="time",
                    M=2,
                    reference=1,
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
   covar=c("base","base2"),
   depvar="fev",
   treatvar="treat",
   idvar="id",
   timevar="time",
   M=2,
   reference=1,
   method="causal",
   seed=101,
   prior="jeffreys",
   burnin=1000,
   K0=1,
   K1=1
)
# verify causal with k0=k1=1 is CIR
sum(impcausal2!=impCIR,na.rm=TRUE)
compare=cbind(impcausal2,impcausal2-impCIR)
subset(compare,.id==5017)
subset(compare,.id==5030)
subset(compare,.id==5051)
subset(compare,.id==5115)

subset(impCIR,.id==5051)
subset(impcausal2,.id==5051)



# evaluate difference
compare=cbind(impcausal2,impcausal2-impJ2R)
subset(compare,.id==5017)
subset(compare,.id==5030)
subset(compare,.id==5051)
subset(compare,.id==5115)
