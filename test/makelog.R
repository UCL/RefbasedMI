# make log file
setwd("C:/ado/ian/RefBasedMI/test")
sink("main testing program.log", append=FALSE, split=TRUE)
source("main testing program.R")
sink()

# this omits lots of R output
# an alternative is
#    R CMD BATCH "main testing program.R"
# but that fails in a way the main prog doesn't