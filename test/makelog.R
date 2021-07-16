# make log file
setwd("C:/ado/ian/Rmimix/test")
sink("Ian main testing program.log", append=FALSE, split=TRUE)
source("Ian main testing program.R")
sink()
