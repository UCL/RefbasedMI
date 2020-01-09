



#for select
library(dplyr)

#pivot_wider
library(tidyr)

#emnorm
library(norm2)

#cholsolve
library(sparseinv) 

#stat.desc
library(pastecs)

#library(plyr)
#read data
# need haven
library(haven)




# need also input refer
#read data
mxdata <-read.csv("./asthma.csv")
#mxdata <-read.csv(pathdata)

testread <-function(data) {
mxdata <-read.csv(paste0("./",data))
                  }
testkm <- testread("asthma.csv")