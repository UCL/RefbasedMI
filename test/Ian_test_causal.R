# Ian's program for testing causal routine in Rmimix
# Concept: considering outcomes Y2,Y4,Y8,Y12
#   increasing all values of Y4 by 1 should make changes:
#     in patterns 0100 and 1100
#       all imputed values of Y8 by K0*K1^4
#       all imputed values of Y12 by K0*K1^8
#     in other patterns: no change
# IW 3dec2020

### load mimix AVOIDING THE PACKAGE
source("C:\\ado\\ian\\Rmimix\\R\\Runmimix.R")
source("C:\\ado\\ian\\Rmimix\\R\\proprocess.R")
source("C:\\ado\\ian\\Rmimix\\R\\utilities.R")

library(mice)
packageVersion("mice")

setwd("C:/ado/ian/test_Rmimix")

# SETTINGS FOR THIS TEST
K0=0.1
K1=0.9
tweak=3

# Open data
load("C:/ado/ian/Rmimix/data/asthma.RData")
head(asthma)
table(asthma$treat)

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
