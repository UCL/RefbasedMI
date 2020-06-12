####################################### anti-depressant data ##################################################################
antidepressant<-read.table("http://raw.githubusercontent.com/UCL/mimix/master/inst/extdata/SASantidepmethodvar.csv",header=TRUE,sep=",",fileEncoding = "UTF-8-BOM")
#recode covariates to numeric
antidepressant$PATIENT.SEX <- as.numeric(antidepressant$PATIENT.SEX)
antidepressant$TREATMENT.NAME <- as.numeric(antidepressant$TREATMENT.NAME)
# 2= drug, 1= placebo
#recode treatment variable to numeric

antidepressant$PATIENT.SEX <- as.numeric(antidepressant$PATIENT.SEX)
# NOTE in input data pooled.investigator has 999 outlier

############################################# Test Causal method, when Kd = 0 same as J2R, Kd=1 same as CIR
impantiCausal_2NoDelta <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"Causal",101,c("jeffreys"),1000,NULL,NULL,,,1,0.7)
impantiCausal_2 <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"Causal",101,c("jeffreys"),1000,NULL,NULL,c(0.5,0.5,1,1 ),c(-3,-4,-5,-6),1,2)
impantiCausal_1 <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",5,1,"Causal",101,c("jeffreys"),1000,NULL,NULL,,1)
impantiJ2R <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"J2R",101,c("jeffreys"),1000,NULL,NULL,,1)
impantiCIR <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"CIR",101,c("jeffreys"),1000,NULL,NULL,,1)
impantiCR <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"CR",101,c("jeffreys"),1000,NULL,NULL,,1)
impantiLMCF <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"LMCF",101,c("jeffreys"),1000,NULL,NULL,,1)
impantiMAR <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"MAR",101,c("jeffreys"),1000,NULL,NULL,,1)
impantiCausal <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"Causal",101,c("jeffreys"),1000,NULL,NULL,,,1,0)



varlist <-c("HAMD17.TOTAL.4","HAMD17.TOTAL.5","HAMD17.TOTAL.6","HAMD17.TOTAL.7")
analyselist("3618",impantiCausal_2,varlist)
analyselist("2104",impantiCausal_2,varlist)
analyselist("2230",impantiCausal_2,varlist)
analyselist("1513",impantiCausal_2,varlist)
analyselist("2104",impantiCausal_2NoDelta,varlist)
analyselist("2230",impantiCausal_2NoDelta,varlist)
analyselist("1513",impantiCausal_2NoDelta,varlist)
analyselist("2721",impantiCausal_2NoDelta,varlist)
analyselist("2721",impantiCausal_1,varlist)
analyselist("2721",impantiJ2R,varlist)
analyselist("2104",impantiCIR,varlist)
analyselist("2230",impantiCIR,varlist)
analyselist("1513",impantiCIR,varlist)
analyselist("2104",impantiCR,varlist)
analyselist("2230",impantiCR,varlist)
analyselist("1513",impantiCR,varlist)
analyselist("2104",impantiLMCF,varlist)
analyselist("2230",impantiLMCF,varlist)
analyselist("1513",impantiLMCF,varlist)
analyselist("2104",impantiMAR,varlist)
analyselist("2230",impantiMAR,varlist)
analyselist("1513",impantiMAR,varlist)
analyselist("3618",impantiCausal,varlist)
analyselist("2104",impantiCausal,varlist)
analyselist("2230",impantiCausal,varlist)
analyselist("1513",impantiCausal,varlist)

library(mice)
impdata= as.mids(impantiCausal)
fit<-with(data= as.mids(impantiCausal_2),expr = lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+PATIENT.SEX))
summary(pool(fit))          


fit<-with(data= as.mids(impantiCausal_2),expr = lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+PATIENT.SEX))
summary(pool(fit))          


fit<-with(data= as.mids(impantiCausal_2NoDelta), expr = lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+PATIENT.SEX))
summary(pool(fit))

fit<-with(data= as.mids(impantiJ2R), expr = lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+PATIENT.SEX))
summary(pool(fit))
fit<-with(data= as.mids(impantiCIR), expr = lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+PATIENT.SEX))
summary(pool(fit))
fit<-with(data= as.mids(impantiCR), expr = lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+PATIENT.SEX))
summary(pool(fit))
fit<-with(data= as.mids(impantiLMCF), expr = lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+PATIENT.SEX))
summary(pool(fit))
fit<-with(data= as.mids(impantiMAR), expr = lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+PATIENT.SEX))
summary(pool(fit))

impantiJ2R_3covar <- mimix("antidepressant",c("basval","POOLED.INVESTIGATOR","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",5,1,"J2R",101,c("jeffreys"),1000,NULL,NULL,,0.5,c(0.5,0.5,1,1),c(0.1,0.1,0.2,0.2))
varlist <-c("HAMD17.TOTAL.4","HAMD17.TOTAL.5","HAMD17.TOTAL.6","HAMD17.TOTAL.7")
analyselist("2721",impantiJ2R_3covar,varlist)
# NOTE pooled.investigator has a 999 outlier in it
fit<-with(data= as.mids(impantiJ2R_3covar), expr = lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+PATIENT.SEX))
summary(pool(fit))



###################### INDIVIDUAL-SPECIFIC method (plus delta)   ##############################################################
impantiIndivDt <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",5,1,NULL,101,c("jeffreys"),1000,NULL,c("methodvar","referencevar"),c(0.5,0.5,1,1),c(-3,-4,-5,-6))

impantiIndiv <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",5,1,,101,c("jeffreys"),1000,NULL,c("methodvar","referencevar"),)
varlist <-c("HAMD17.TOTAL.4","HAMD17.TOTAL.5","HAMD17.TOTAL.6","HAMD17.TOTAL.7")
analyselist("2721",impantiIndivDt,varlist)
analyselist("2721",impantiIndiv,varlist)
fit<-with(data= as.mids(impantiIndivDt, .id="SNO",.imp="II"),expr=lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+PATIENT.SEX))
summary(pool(fit))




## method with Causal option so temporary add Causal into methodvar
##antidepressant$methodvar <- ifelse(antidepressant$PATIENT.NUMBER=="1507", "Causal", antidepressant$methodvar)
impantiIndivCausal <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",5,1,,101,c("jeffreys"),1000,NULL,c("methodvar","referencevar"),NULL,NULL,1,0.5)
# to find all .id==3411
impantiIndivCausal[which(impantiIndivCausal$.id=="3411"),]


#################################   with/wout covars
impantiIndivNOcov <- mimix("antidepressant",,"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",5,1,,101,c("jeffreys"),1000,NULL,c("methodvar","referencevar"))
analyselist("2721",impantiIndivNOcov,varlist)
fit<-with(data= as.mids(impantiIndivNOcov, .id="SNO",.imp="II"),expr=lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+PATIENT.SEX))
summary(pool(fit))

pttestf<- function(n1,n2,mn1,sd1,mn2,sd2) {
  pttest = stats::pt((((mn1 - mn2) -0) /sqrt(sd1^2/n1+sd2^2/n2)),(n1+n2-2))
  return(pttest)
  
  # save into package     
  save(antidepressant, file="data/antidepressant.RData")  
  
  # check non missins are not imputed !
  
}