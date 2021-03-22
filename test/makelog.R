# make log file
setwd("C:/ado/ian/Rmimix/test")
sink("allmissing_R.log", append=FALSE, split=TRUE)
source("allmissing_R.R")
sink()
