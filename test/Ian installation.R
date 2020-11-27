# Commands needed to install mimix from github
# IW 19nov2020
# need to have previously installed Rtools
if(!require(devtools)) install.packages('devtools') 
library(devtools) 
install_url('https://cran.r-project.org/src/contrib/Archive/norm2/norm2_2.0.3.tar.gz')
install_github("UCL/mimix")
# couldn't update package: vectrs, mice

# remember that mimix sessions need to include:
library(mimix)


### ALTERNATIVE WAY AVOIDING THE PACKAGE

source("C:\\ado\\ian\\Rmimix\\R\\Runmimix.R")
source("C:\\ado\\ian\\Rmimix\\R\\proprocess.R")
source("C:\\ado\\ian\\Rmimix\\R\\utilities.R")
