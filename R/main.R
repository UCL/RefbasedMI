#####################################################################################
# Rmimix.R
# R program to mimic stata program mimix
# ie reference based imputation 
# Note 1st part is to set up a summary table based on missing data pattern- mimix_group
# reflects the pattern and treatment group configuration of the raw data
# then acts as a lookup table to provide the looping mechanism 
# norm2 is used as MCMC multivariate normal  
# this version 6/1/2020
# v0.2
# Author : Kevin McGrath
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

# either save Stata data file dirertory into csv format or use
# haven  to read Stata file directly 
#library(haven)
#mxdata <- read_dta("asthma.dta") 


# Assign list of input parameters 
kmargs <- list("fev","treat","id","time","base",10000,2,"J2R",301)

# run main program outputting list containing the M imputed data sets  

mimix_outputlist=do.call('Runmimix', kmargs)

# for program timings
#system.time(do.call('Runmimix', kmargs))

# list of imputed data
mata_all_newlist <- mimix_outputlist[1]
# pattern matching
mg <- (mimix_outputlist[2])
# Number of imputations
M <- unlist(mimix_outputlist[3])
#method chosen
meth <- unlist(mimix_outputlist[4])


#if run below then outputs the data set to screen so better to run as above 
#Runmimix(kmargs)
# equiv to
#Runmimix("fev","treat","id","time","base","J2R",10000,2)

# produce summary for individual
analselist(meth,"5456")

system.time(analselist(meth,"5456"))

# to handle and combine the outputted muliple data sets see for example
#https://rdrr.io/cran/norm2/man/miInference.html
# but instead of mcmcResult$imp.list just create new list from the saved list output


dimlist <- (nrow(mg[[1]])*M)
mata_all_newlist[[1]][[dimlist]]
mata_all_newData1k <- do.call(rbind,mata_all_newlist[[1]])
# then sort into M data sets and maybe split into M lists !?
testkm1k<-mata_all_newData1k[order(mata_all_newData1k$II,mata_all_newData1k$SNO),]
# to get the list
testkmlist1k <- split(testkm1k,testkm1k$II)
# so has M elemets in list
est.list <- as.list(NULL)
std.err.list <- as.list( NULL )
for( m in 1:M ){
  yimp <- testkmlist1k[[m]]  # one imputed dataset
  diff <- yimp[,"fev2"] 
  est.list[[m]] <- mean(diff)
  std.err.list[[m]] <- sqrt( var(diff) / length(diff) ) }
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
pttestf(10000,10000,1.040,0.385,1.051,0.341)


#for program timings
end_time <- proc.time()
time_taken <- end_time - start_time
print(paste("Time taken:", time_taken[1]))

system.time(analselist("5137")) 

 
  






