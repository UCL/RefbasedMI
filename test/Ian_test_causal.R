# Ian's program for testing causal routine in Rmimix

library(mice)
packageVersion("mice")

setwd("C:/ado/ian/test_Rmimix")

# Open data
load("C:/ado/ian/Rmimix/data/asthma.RData")
head(asthma)
table(asthma$treat)

# J2R
impJ2R <- mimix(data="asthma",
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
)

head(impJ2R[impJ2R$.imp==1])