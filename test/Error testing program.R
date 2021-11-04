#####################################################################
# Error testing program for RefBasedMI
# IW 4/11/2021
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

# non-existent covariate
impJ2R1 <- RefBasedMI(data=asthma,
                      covar=c(base,base3),
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
# Error in names(get("data"))[[grep(paste0("^", scovar[[i]], "$"), names(get("data")))]] : 
#  attempt to select less than one element in get1index
# this is runmimix.R line 111


# non-existent depvar
impJ2R1 <- RefBasedMI(data=asthma,
                      covar=c(base,base2),
                      depvar=fev55,
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
# Error in `[.data.frame`(get("data"), c(idvar, depvar, timevar, treatvar,  : 
# undefined columns selected 
# similar error for treatvar, idvar, timevar
