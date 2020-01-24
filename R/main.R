#####################################################################################  
# Rmimix.R                                                                          #
# R program to mimic stata program mimix                                            #
# ie reference based imputation                                                     #
# Note 1st part is to set up a summary table based on missing data pattern-         #
# mg  mimix_group                                                                   #
# reflects the pattern and treatment group configuration of the raw data            #
# then acts as a looping mechanism                                                  #
# norm2 is used as MCMC multivariate normal                                         #
#                                                                                   #
# calls functions listed in functions.R file                                        #
# function preprodata prepares data and and finds mg (the mimix group)              # 
# function Runmimix  performs major analysis                                        #
# the required packages as listed in utilities file                                 #
# this version 6/1/2020                                                             #
# v0.0.2                                                                              #
# Author : Kevin McGrath                                                            #
#####################################################################################


# to check home directory  
#getwd()
#source("N:/Documents/GitHub/mimix/mimixR/functions.R")

# remove existing files
rm(list = ls())

# file refers to functions called from main program  
source("functions.R")

#read the data file - in csv format
mxdata<- readdata("asthma.csv")

# or directly from github
mxdata<-read.table("http://raw.githubusercontent.com/UCL/mimix/master/data/asthma.csv",header=TRUE,sep=",",fileEncoding = "UTF-8-BOM")

# acupuncture data
mxdata<-read.table("http://raw.githubusercontent.com/UCL/mimix/master/data/accupuncture.csv",header=TRUE,sep=",",fileEncoding = "UTF-8-BOM")
# next2 below just for testing whther no typ eworks instead of char for treat

#save treat col as want o recode as numeric
mxdata$treatcopy<- mxdata$treat
mxdata$treat<-(mxdata$treatcopy=="accupuncture")*1
mxdata$treat<-(mxdata$treatcopy=="control")*1+1
# at moment hard coded to use treat as argument , so rename any teratmnet variable in input data to treat.
kmargs <- list("head","treat","id","time","head_base",1000,2,"J2R",101)
mimix_outputlist <- Runmimix(kmargs)

# do need to introduce a covar argument to cope when more than 1 covariate variable.  

# Assign list of input parameters 
kmargs <- list("fev","treat","id","time","base",10000,2,"J2R",201)


# run main program outputting list containing the M imputed data sets  

#mimix_outputlist=do.call('Runmimix', kmargs)
mimix_outputlist <- Runmimix(kmargs)
#mimix_outputlist<-Runmimix("fev","treat","id","time","base",100,2,"J2R",101)


# for program timings
#system.time(do.call('Runmimix', kmargs))

# save list of imputed data
mata_all_newlist <- mimix_outputlist[1]
# pattern matching
mg <- (mimix_outputlist[2])
# Number of imputations
M <- unlist(mimix_outputlist[3])
#method chosen
meth <- unlist(mimix_outputlist[4])



# produce summary for individual
analyselist(meth,"5456")
analyselist(meth,"100")

#system.time(analyselist(meth,"5456"))

# to handle and combine the outputted muliple data sets see for example
#https://rdrr.io/cran/norm2/man/miInference.html
# but instead of mcmcResult$imp.list just create new list from the saved list output

# convert imputed data list into combined data set
# dimension of data set, nrows in pattern times no imputations, 
dimlist <- (nrow(mg[[1]])*M)
# extract from nested list 
mata_all_newlist[[1]][[dimlist]]
# combine into data set containing M imputed datasets 
mata_all_newData1x <- do.call(rbind,mata_all_newlist[[1]])


# then sort into M data sets and maybe split into M lists 
km1x<-mata_all_newData1x[order(mata_all_newData1x$II,mata_all_newData1x$SNO),]
# to get the list
kmlist1x <- split(km1x,km1x$II)
# so has M elements in list
# can obtain a list of coefficients and their se's from a regression
# declare list for estimates 
est.list <- as.list(NULL)
# declare lists for se's 
std.err.list <- as.list( NULL )
for( m in 1:M ){
  mod<-lm(fev12~as.factor(treat)+base,data=kmlist1x[[m]] )
  
  est.list[[m]] <- coef(summary(mod))[,1]
  std.err.list[[m]] <- coef(summary(mod))[,2] }
## combine the results by rules of Barnard and Rubin (1999)
## with df.complete = 27, because a paired t-test in a dataset with
## N=28 cases has 27 degrees of freedom
miResult <- miInference(est.list, std.err.list, df.complete=182)
print(miResult)



# Alternatively 
# this from CRAN vignette/amelia.pdf
b.out<-NULL
se.out<-NULL
for( m in 1:M ){
  ols.out<-lm(fev12~treat+base,data=testkmlist1k[[m]] )
  b.out<-rbind(b.out,ols.out$coef)
  se.out<-rbind(se.out,coef(summary(ols.out))[,2])
}
library(Amelia)
combined.results<-mi.meld(q=b.out,se=se.out)
# to write as Stat file
#write.amelia(obj=mata_all_newData1k,file.stem = "outdata",format = "dta")
print(combined.results)

# sign tests compard with Stata
pttestf(10000,10000,1.614,0.378,1.626,0.369)
pttestf(1000,1000,11.912,11.116,12.34,10.4)

#for program timings
end_time <- proc.time()
time_taken <- end_time - start_time
print(paste("Time taken:", time_taken[1]))

system.time(analyselist("5137")) 

 
  






