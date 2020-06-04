#### program showing example run for mimix program ######   

# remove existing files
rm(list = ls())

# if not install from github use source ("/functions.R") 
# eg source("dir-path/functions.R")

# read csv data file directly from github
# asthma data
mxdata<-read.table("http://raw.githubusercontent.com/UCL/mimix/master/inst/extdata/asthma.csv",header=TRUE,sep=",",fileEncoding = "UTF-8-BOM")

#recode treatment groups to 1,2..
mxdata$treat<-ifelse(mxdata$treat==2,1,ifelse(mxdata$treat==3,2,9))


mimix_output<-mimix("mxdata",c("base"),"fev","treat","id","time",100,1,"J2R",101,"jeffreys",1000,NULL,NULL)   # replace by ridge c("ridge","0.5") rids warnings


# convert imputed dataset to mids format (from mice package)
impdata= as.mids(mimix_output, .id="SNO",.imp="II")

# fit specified model to each imputed data set and pool results together (Rubin's rules), functions from mice package
fit<-with(impdata, lm(fev.12~treat+base))
summary(pool(fit))


# alternaive method based on miInfernce from norm2 
regressimp(mimix_output, "fev.12~treat+base")

# to examine estimates for individual ids
varlist <- c("fev.2","fev.4","fev.8","fev.12","base")
# produce summary for individual
analyselist(5099,mimix_output,varlist)





###########################################################################################################################################
#anti-depressant data 
#mxdata<-read.table("http://raw.githubusercontent.com/UCL/mimix/master/data/SASantidep.csv",header=TRUE,sep=",",fileEncoding = "UTF-8-BOM")
#anti-depressant data with methdovar & refernecevar
mxdata<-read.table("http://raw.githubusercontent.com/UCL/mimix/master/data/inst/extdata/SASantidepmethodvar.csv",header=TRUE,sep=",",fileEncoding = "UTF-8-BOM")


#recode treatment variable to numeric

mxdata$TREATMENT.NAME<-ifelse(mxdata$TREATMENT.NAME=="PLACEBO",1,ifelse(mxdata$TREATMENT.NAME=="DRUG",2,9))


# mxdata <- subset(mxdata, POOLED.INVESTIGATOR !="999") 

# try running methodvar, note cant hanle "NULL" 
mimix_outputlist <- Runmimix(c("basval"),"change","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",2,NULL,NULL,101,c("jeffreys"),1000,NULL,c("methodvar","referencevar")) 



mimix_outputlist <- Runmimix(c("basval"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",1000,"DRUG","J2R",101,c("jeffreys")) #c("ridge","0.5"))
impdatasets<-getimpdatasets(mimix_outputlist)
analyselist("2721",mimix_outputlist,c("change4","change5","change6","change7","basval"))
regressimp(impdatasets,"change4~TREATMENT.NAME+basval")

#try with
with(impdatasets,lm("change4~TREATMENT.NAME+basval"))
#4623


# dont need to introduce a covar argument to cope when more than 1 covariate variable.  
# Assign list of input parameters 
#kmargs <- list("fev","treat","id","time",covar("treat","base"),1,2,"J2R",201)


# run main program outputting list containing the M imputed data sets  

#mimix_outputlist=do.call('Runmimix', kmargs)

#mimix_outputlist<-Runmimix("fev","treat","id","time","base",100,2,"J2R",101)

# for program timings
#system.time(do.call('Runmimix', kmargs))



impdatasets<-getimpdatasets(mimix_outputlist)
regressimp(impdatasets,"fev12~treat+base")

# produce summary for individual
analyselist(meth,"5456",c("fev2","fev4","fev8","fev12","base"))
analyselist("100",mimix_outputlist,c("head3","head12","sex","head_base"))
analyselist(meth,"4623",c("HAMD17.TOTAL4","HAMD17.TOTAL5","HAMD17.TOTAL6","HAMD17.TOTAL7"))

#system.time(analyselist(meth,"5456"))

# to handle and combine the outputted muliple data sets see for example
#https://rdrr.io/cran/norm2/man/miInference.html
# but instead of mcmcResult$imp.list just create new list from the saved list output

# convert imputed data list into combined data set
# dimension of data set, nrows in pattern times no imputations, 
dimlist <- (nrow(mg[[1]])*M)
# extract from nested list 
#head(mata_all_newlist[[1]][[dimlist]])
# combine into data set containing M imputed datasets 
mata_all_newData1x <- do.call(rbind,mata_all_newlist[[1]][[dimlist]])


# then sort into M data sets and split into M lists 
# ort by imputation no and patient id (SNO)
impdatasets<-mata_all_newData1x[order(mata_all_newData1x$II,mata_all_newData1x$SNO),]
# to get the list
implist1x <- split(impdatasets,impdatasets$II)
# so has M elements in list
# can obtain a list of coefficients and their se's from a regression
# declare list for estimates 
est.list <- as.list(NULL)
# declare lists for se's 
std.err.list <- as.list( NULL )
for( m in 1:M ){
  #mod<-lm(fev12~as.factor(treat)+base,data=kmlist1x[[m]] )
  #mod<-lm(head12~head_base+sex,data=implist1x[[m]] )
  mod<-lm(HAMD17.TOTAL7~basval+HAMD17.TOTAL6,data=implist1x[[m]] )
  est.list[[m]] <- coef(summary(mod))[,1]
  std.err.list[[m]] <- coef(summary(mod))[,2] }
## combine the results by rules of Barnard and Rubin (1999)
## with df.complete = 27, because a paired t-test in a dataset with
## N=28 cases has 27 degrees of freedom
miResult <- miInference(est.list, std.err.list, df.complete=801)
print(miResult)

#mixed model long data repeated measures

# HAMD17 is repeatd measure, treatment between subjects effect
library(lme4)
lmeModel = lmer(HAMD17~treat*time+ (1|))
# need to reshape to long data format to run lmer.
# so reshape each list element of the imputed data sets then use apply?

#test<-reshape(mimix_outputlist[[1]][1],

# sort by SNO and time
#tail(test[order(test$id,test$SNO,test$II),],10)
impdatalong <- test[order(test$id,test$SNO,test$II),]
implonglist1x <- split(impdatalong,impdatalong$II)
# declare lists for model estimates and their se's 
est.list <- as.list(NULL)
std.err.list <- as.list( NULL )

for( m in 1:M ){
  lm_mod <-  lmer(HAMD17~TREATMENT+basval+POOLED.INVESTIGATOR+(1|SNO),data=implonglist1x[[m]])
  est.list[[m]] <- coef(summary(lm_mod))[,1]
  std.err.list[[m]] <- coef(summary(lm_mod))[,2]
}
## combine the results by rules of Barnard and Rubin (1999)
## with df.complete = 27, because a paired t-test in a dataset with
## N=28 cases has 27 degrees of freedom
miResult <- miInference(est.list, std.err.list, df.complete=171)
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


pttestf(200,200,1.302,0.623,1.367,0.573)
