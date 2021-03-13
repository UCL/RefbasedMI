# Ian's program to test what happens with no observed data
# 12mar2021

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

# make id=5001 wholly missing
asthma[1:4,"fev"]<-NA

methods<-c("MAR","J2R","CR","CIR","LMCF")
imp1 <- list()
imp2 <- list()
for(i in 1:length(methods)){
  print(methods[i])
  imp1[[i]] <- mimix(data=asthma,
                     covar="base",
                     depvar=fev,
                     treatvar=treat,
                     idvar=id,
                     timevar=time,
                     M=1,
                     reference=1,
                     method=methods[i],
                     seed=101,
                     prior=jeffreys,
                     burnin=1000,
                     bbetween=NULL,
                     methodvar=NULL
  )
  imp2[[i]] <- mimix(data=asthma,
                     covar="base",
                     depvar=fev,
                     treatvar=treat,
                     idvar=id,
                     timevar=time,
                     M=1,
                     reference=2,
                     method=methods[i],
                     seed=101,
                     prior=jeffreys,
                     burnin=1000,
                     bbetween=NULL,
                     methodvar=NULL
  )
}
print("Ref=1")
for(i in 1:length(methods)){
  print(methods[i])
  print(imp1[[i]] %>% filter(.id==5001))
}
print("Ref=2")
for(i in 1:length(methods)){
  print(methods[i])
  print(imp2[[i]] %>% filter(.id==5001))
}
# LMCF doesn't impute no-obs reference-arm individuals under MAR
# All methods don't impute no-obs active-arm individuals under MAR
