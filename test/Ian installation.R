# Commands needed to install mimix from github
# IW 19nov2020 revised 30mar2021
# need to have previously installed Rtools

# mimix
if(!require(devtools)) install.packages('devtools') 
library(devtools) 
install_github("UCL/RefBasedMI")
# 12mar2021: Kevin has put latest package on developer branch
# install_github("UCL/mimix", ref="developer")
# couldn't update package: vectrs, mice

# norm2
# install_url('https://cran.r-project.org/src/contrib/Archive/norm2/norm2_2.0.3.tar.gz')
library(norm2)

# remember that mimix sessions need to include:
library(mimix)


### ALTERNATIVE WAY AVOIDING THE PACKAGE

source("C:\\ado\\ian\\Rmimix\\R\\Runmimix.R")
source("C:\\ado\\ian\\Rmimix\\R\\proprocess.R")
source("C:\\ado\\ian\\Rmimix\\R\\utilities.R")
