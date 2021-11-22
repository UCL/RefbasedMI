# Test RefBasedMI with artificial data with 12 times, 4 groups and no covariate
# IW 22nov2021
# twelvetimes20.csv  has  20 indiviuals / group
# twelvetimes200.csv has 200 indiviuals / group

# Problem 1: if you use the data set with 20 / group,
#   you get "MCMC procedure aborted" in each imputation, 
#   yet a final "mcmcNorm Loop finished" and it soldiers on
# I want it to abort!

# Problem 2: with 200 / group the mcmcloop appears to succeed but then get the error:
# Imputing interim missing values under MAR:
#   
#   Error in `colnames<-`(`*tmp*`, value = idvar) : 
#   attempt to set 'colnames' on an object with less than two dimensions 

library(tidyverse)
t12 <- read_csv("twelvetimes200.csv")
summary(t12)


# J2R
impJ2R1 <- RefBasedMI(data=t12,
                      depvar=yvar,
                      treatvar=group,
                      idvar=myid,
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


