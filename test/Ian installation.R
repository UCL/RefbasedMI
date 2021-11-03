# Commands needed to install RefBasedMI from github
# IW 19nov2020 revised 30mar2021 and 22oct2021
# need to have previously installed Rtools

# RefBasedMI
devtools::install_github("UCL/RefBasedMI",ref="master",force=TRUE,dependencies=FALSE)
packageVersion("RefBasedMI")
library(RefBasedMI)

### recreate help (locally)
setwd("C:/ado/ian/RefbasedMI")
roxygen2::roxygenise()

help(RefBasedMI)

# norm2
# install_url('https://cran.r-project.org/src/contrib/Archive/norm2/norm2_2.0.3.tar.gz')
library(norm2)


### ALTERNATIVE WAY AVOIDING THE PACKAGE
# source("C:\\ado\\ian\\Rmimix\\R\\Runmimix.R")
# source("C:\\ado\\ian\\Rmimix\\R\\proprocess.R")
# source("C:\\ado\\ian\\Rmimix\\R\\utilities.R")


