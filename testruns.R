# testruns and examples for mimix program

######################################### asthma data-set testing J2R method #############################################
# testing all options, m=5 imputations , method = J2R, referEnce treatment gp =1,seed = 101,  Delta specified
impdatasetJ2R<-mimix("asthma",c("base"),"fev","treat","id","time",200,1,"J2R",101,"jeffreys",1000,NULL,NULL,NULL,c(0.5,0.5,1,1 ),NULL )

# convert imputed dataset to mids format (from mice package)
impdatamids= as.mids(impdatasetJ2R)

#  converting mids type not  pmm to mimix and  and predictormatrix not relevant 
impdatamids$method[impdatamids$method %in% "pmm"] <- "mimix"
impdatamids$predictorMatrix<-"N/A"

# fit specified model to each imputed data set and pool results together (Rubin's rules), functions from mice package
fit<-with(impdatamids, lm(fev.12~treat+base))
summary(pool(fit))

# previously regressimp(impdatasetJ2R[impdatasetJ2R$II!="0",],"fev.12~treat+base")

# to examine estimates for individual ids
varlist <- c("fev.2","fev.4","fev.8","fev.12")
# print estimates for individuals  , note 5051 interim example
analyselist(5017,impdatasetJ2R,varlist)
analyselist(5127,impdatasetJ2R,varlist)
analyselist(5099,impdatasetJ2R,varlist)
analyselist(5051,impdatasetJ2R,varlist)
analyselist(5094,impdatasetJ2R,varlist)
analyselist(5333,impdatasetJ2R,varlist)
analyselist(5456,impdatasetJ2R,varlist)


# select 1 id from each pattern group 
analyselist(5017,impdatasetJ2jeffreys,varlist)
analyselist(5127,impdatasetJ2jeffreys,varlist)
analyselist(5099,impdatasetJ2jeffreys,varlist)
analyselist(5051,impdatasetJ2jeffreys,varlist)
analyselist(5094,impdatasetJ2jeffreys,varlist)
analyselist(5333,impdatasetJ2jeffreys,varlist)

#id=5017 (pattern OXXX) gets exactly the same imputations as J2R but should be different, 
#id= 5115 (pattern OOXO) i.e. only interim missing
#patt 0  treat 1  5001  0000
#        treat 2  5003  0000
#patt 14 treat 1  5017  0xxx
#patt 14 treat 2  5127  0xxx
#patt 12 treat 1  5044 00xx
#patt 12 treat 2  5099 00xx
#patt 13 treat 1  5051 x0xx
#patt 13       2
#patt 8        1  5074 000x
#patt 8        2  5094 000x
#patt 4 treat  2  5115 00x0
#patt 7 treat  1  5333 xxx0



##################################################### no delta option #########################################
impdatasetJ2R_Nd<-mimix("asthma",c("base"),"fev","treat","id","time",200,1,"J2R",101,"jeffreys",1000,NULL,NULL,, )
analyselist(5099,impdatasetJ2R_Nd,varlist)
fit<-with(data= as.mids(impdatasetJ2R_Nd, .id="SNO",.imp="II"), expr = lm(fev.12~treat))
summary(pool(fit))

################################################ differnt priors ################################################
impdatasetJ2jeffreys<-(mimix("asthma",c("base"),"fev","treat","id","time",200,1,"J2R",101,"jeffreys",1000,NULL,NULL,NULL) )
analyselist(5017,impdatasetJ2jeffreys,varlist)
analyselist(5127,impdatasetJ2jeffreys,varlist)
analyselist(5099,impdatasetJ2jeffreys,varlist)
analyselist(5051,impdatasetJ2jeffreys,varlist)
analyselist(5094,impdatasetJ2jeffreys,varlist)
analyselist(5333,impdatasetJ2jeffreys,varlist)


impdatasetJ2Rridge<-(mimix("asthma",c("base"),"fev","treat","id","time",5,1,"J2r",101,"ridge",1000,NULL,NULL,NULL,c(0.5,0.5,1,1 ) ) )
analyselist(5017,impdatasetJ2Rridge,varlist)



###################################### NO COVAR #######################################################################
# test for no covar option, create new data set asthmafev0,  converts base variable into initial depvar (here fev.0)
tmp<-(unique(asthma[,c("id","treat","base")]))
tmp$time<-0
tmp$fev <- tmp$base
tmp2<-(rbind(tmp,asthma))
asthmafev0<-tmp2[order(tmp2[,1]),]

# NOw use asthmafev0 data set with  no covars
imptestdatasetJ2R5<-(mimix("asthmafev0",NULL,"fev","treat","id","time",5,1,"J2R",101,"jeffreys",1000,NULL,NULL,NULL,0.5 ))
varlistnull <- c("fev.2","fev.4","fev.8","fev.12","fev.0")
analyselist(5099,imptestdatasetJ2R5,varlistnull)
###########################################################################################################################





#####################################Testing different methods ############################################################


#test all the  methods at least run ok

######################################### MAR method jeffreys default
impdatasetMAR<-(mimix("asthma",c("base"),"fev","treat","id","time",5,1,"MAR",101,,1000,NULL,NULL,NULL ) )
# ridge prior
impdatasetMARridge<-(mimix("asthma",c("base"),"fev","treat","id","time",2,1,"MAR",101,"ridge",1000,NULL,NULL,NULL ) )
analyselist(5099,impdatasetMAR,varlist)

######################################### CIR
impdatasetCIR<-(mimix("asthma",c("base"),"fev","treat","id","time",5,1,"CIR",101,,1000,NULL,NULL, ) )
impdatasetCIRridge<-(mimix("asthma",c("base"),"fev","treat","id","time",5,1,"CIR",101,"ridge",1000,NULL,NULL,NULL,c(0.5,0.5,1,1 ) ) )
analyselist(5099,impdatasetCIR,varlist)
analyselist(5456,impdatasetCIR,varlist)
analyselist(5456,impdatasetCIRridge,varlist)
analyselist(5051,impdatasetCRridge,varlist)

######################################## CR
  impdatasetCR<-(mimix("asthma",c("base"),"fev","treat","id","time",5,1,"CR",101,,1000,NULL,NULL,NULL ) )
impdatasetCRridge<-(mimix("asthma",c("base"),"fev","treat","id","time",5,1,"CR",101,"ridge",1000,NULL,NULL,NULL, ) )
impdatasetCRridgeDelt<-(mimix("asthma",c("base"),"fev","treat","id","time",100,1,"CR",101,"ridge",1000,NULL,NULL,NULL,c(0.5,0.5,1,1),c(0.5,0.5,2,2)  ) )
analyselist(5051,impdatasetCRridge,varlist)
analyselist(5051,impdatasetCRridgeDelt,varlist)

analyselist(5099,impdatasetCR,varlist)
analyselist(5017,impdatasetCR,varlist)
analyselist(5099,impdatasetCRridge,varlist)
analyselist(5017,impdatasetCRridge,varlist)

########################################## LMCF
impdatasetLMCF<-(mimix("asthma",c("base"),"fev","treat","id","time",5,1,"LMCF",101,,1000,NULL,NULL,NULL ) )
impdatasetLMCFridgedelt<-(mimix("asthma",c("base"),"fev","treat","id","time",5,1,"LMCF",101,"ridge",1000,NULL,NULL,NULL,c(0.5,0.5,1,1),c(0.5,0.5,2,2) ) )
analyselist(5456,impdatasetLMCF,varlist)
analyselist(5099,impdatasetLMCF,varlist)
analyselist(5456,impdatasetLMCFridge,varlist)
analyselist(5099,impdatasetLMCFridge,varlist)
analyselist(5051,impdatasetLMCFridge,varlist)
analyselist(5051,impdatasetLMCFridgedelt,varlist)
analyselist(5115,impdatasetLMCFridgedelt,varlist)


##########################################  Causal method
# ERROR!!! induced because kd not specified so give errr msg
impdatasetCausalridge<-(mimix("asthma",c("base"),"fev","treat","id","time",5,1,"Causal",101,"ridge",1000,NULL,NULL,NULL,K0=1,K1=1 ) )
impdatasetCausal_11<-(mimix("asthma",c("base"),"fev","treat","id","time",5,2,"Causal",101,,1000,NULL,NULL,NULL,K0=1,K1=1 ) )
impdatasetJ2R<-(mimix("asthma",c("base"),"fev","treat","id","time",5,1,"J2R",101,,1000,NULL,NULL,NULL,K0=1,K1=1 ) )
impdatasetJ2Rdelta<-(mimix("asthma",c("base"),"fev","treat","id","time",100,1,"J2R",101,,1000,NULL,NULL,NULL,delta=c(1,5,2,3),K0=1,K1=1 ) )


analyselist(5456,impdatasetCausal,varlist)
analyselist(5611,impdatasetCausal,varlist)
analyselist(5099,impdatasetCausal,varlist)
analyselist(5017,impdatasetCausal,varlist)
analyselist(5115,impdatasetCausal,varlist)

analyselist(5611,impdatasetCausal_11,varlist)
analyselist(5017,impdatasetCausal_11,varlist)
analyselist(5115,impdatasetCausal_11,varlist)
analyselist(5611,impdatasetCausal_00,varlist)
analyselist(5017,impdatasetCausal_00,varlist)
analyselist(5115,impdatasetCausal_00,varlist)

analyselist(5456,impdatasetJ2R,varlist)
analyselist(5611,impdatasetJ2R,varlist)
analyselist(5099,impdatasetJ2R,varlist)
analyselist(5017,impdatasetJ2R,varlist)
analyselist(5115,impdatasetJ2R,varlist)
analyselist(5051,impdatasetJ2R,varlist)
analyselist(5017,impdatasetJ2Rdelta,varlist)
analyselist(5115,impdatasetJ2Rdelta,varlist)
analyselist(5051,impdatasetJ2Rdelta,varlist)
analyselist(5017,impdatasetJ2Rnodelta,varlist)
analyselist(5115,impdatasetJ2Rnodelta,varlist)
analyselist(5051,impdatasetJ2Rnodelta,varlist)

#convert to mids object (in mice package)
#impdatamids<-as.mids(impdatasets)
#impdatamids$method[IndivDtmids$method %in% "pmm"] <- "mimix"
#impdatamids$predictorMatrix<-"N/A"

library(mice)
fit<-with(data= as.mids(impdatasetJ2Rdelta), expr = lm(fev.12~treat+base))
summary(pool(fit))


# multiply outcomes by 10 to check same treatment effect 19/06
impdatasetJ2Rdelta$fev.4 <- impdatasetJ2Rdelta$fev.4*10
impdatasetJ2Rdelta$fev.8 <- impdatasetJ2Rdelta$fev.8*10
impdatasetJ2Rdelta$fev.12 <- impdatasetJ2Rdelta$fev.12*10
impdatasetJ2Rdelta$base <- impdatasetJ2Rdelta$base*10
fit<-with(data= as.mids(impdatasetJ2Rdelta), expr = lm(fev.12~treat+base))
summary(pool(fit))


# testimp  function (defined below)
################################## testing non-missing data not affected by imputation
# and all missing values are imputed 
# original data set with missing values in long format so select out time point
asthmafev2<-(subset(asthma[,c("fev","time")],time==2))
asthmafev4<-(subset(asthma[,c("fev","time")],time==4))
asthmafev8<-(subset(asthma[,c("fev","time")],time==8))
asthmafev12<-(subset(asthma[,c("fev","time")],time==12))

testimp(asthmafev2,impdatasetJ2R,5,"fev.2","fev")
testimp(asthmafev4,impdatasetJ2R,5,"fev.4","fev")
testimp(asthmafev8,impdatasetJ2R,5,"fev.8","fev")
testimp(asthmafev12,impdatasetJ2R,5,"fev.12","fev")


# original data set with missing values in long format so select out time point
antidepHamd17vis5<-(subset(antidepressant[,c("HAMD17.TOTAL","VISIT.NUMBER")],VISIT.NUMBER==5))
antidepHamd17vis6<-(subset(antidepressant[,c("HAMD17.TOTAL","VISIT.NUMBER")],VISIT.NUMBER==6))
antidepHamd17vis7<-(subset(antidepressant[,c("HAMD17.TOTAL","VISIT.NUMBER")],VISIT.NUMBER==7))
antidepHamd17vis4<-(subset(antidepressant[,c("HAMD17.TOTAL","VISIT.NUMBER")],VISIT.NUMBER==4))

testimp(antidepHamd17vis4,impantiJ2R,100,"HAMD17.TOTAL.4","HAMD17.TOTAL")
testimp(antidepHamd17vis5,impantiJ2R,100,"HAMD17.TOTAL.5","HAMD17.TOTAL")
testimp(antidepHamd17vis6,impantiJ2R,100,"HAMD17.TOTAL.6","HAMD17.TOTAL")
testimp(antidepHamd17vis7,impantiJ2R,100,"HAMD17.TOTAL.7","HAMD17.TOTAL")

#  define test function 
testimp <- function(dataf,imputdata,M,var,TOTALvar)
{
  #expand orig dataset same length as imputed 
  dummydataf<-dataf
  #concatenate as much as M times
  for ( i in 1:M) {
    dummydataf<-rbind( dummydataf,dataf)
  }
  # browser()  
  datacomp<-cbind(dummydataf,subset(imputdata,select=c(var,".imp")))
  datacompna<-(subset(datacomp,!is.na(datacomp[,TOTALvar])))
  # this compares original values with imputed dat set  should be identical!! 
  cat("\n",paste(var,"non-missing = ",(nrow(datacompna)/(M+1))))
  if (all(datacompna[,1]==datacompna[,3]) ) {cat("\nnon-missing data unchanged")} 
  
  # need to check the nas have imputed values so se;ect out rows with nas
  datacompjustna<-subset(datacomp,(is.na(datacomp[,TOTALvar]) & .imp!=0))
  # check this should be 0 if all imputed.
  cat("\n", "number nas = ",paste(nrow(datacompjustna)/M))
  if (sum(is.na(datacompjustna[,3]))==0) {cat("\n",paste(var,"all imputed"))} 
  #if (all.equal(datacompna[,1],datacompna[,3]) ) {cat("original data unchanged")}   
  # return(dummydataf)  
}


##################################### test whether causal(0,0) gives same as J2R.
impantiCausal_2NoDelta00 <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"Causal",101,c("jeffreys"),1000,NULL,NULL,NULL,,,0,0)
impantiJ2R <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"J2R",101,c("jeffreys"),1000,NULL,NULL,NULL,,1)
#check whether same results 
all.equal(impantiCausal_2NoDelta00,impantiJ2R)


varlist <-c("HAMD17.TOTAL.4","HAMD17.TOTAL.5","HAMD17.TOTAL.6","HAMD17.TOTAL.7")
analyselist("3618",impantiCausal_2,varlist)
analyselist("2104",impantiCausal_2,varlist)
analyselist("2230",impantiCausal_2,varlist)
analyselist("3618",impantiJ2R ,varlist)
analyselist("2104",impantiJ2R ,varlist)
analyselist("2230",impantiJ2R ,varlist)

################################################ test whether causal(1,1) gives same as CIR?.
impantiCausal_2NoDelta11 <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"Causal",101,c("jeffreys"),1000,NULL,NULL,NULL,,,1,1)
impantiCIR <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"CIR",101,c("jeffreys"),1000,NULL,NULL,NULL,,1)
all.equal(impantiCausal_2NoDelta11,impantiCIR)



analyselist(5099,impdatasetCausalridge,varlist)


fit<-with(data= as.mids(impdatasetCIRridge, .id="SNO",.imp="II"), expr = lm(fev.12~treat+base))
summary(pool(fit))
fit<-with(data= as.mids(impdatasetCIRridge), expr = lm(fev.12~treat+base))
summary(pool(fit))

####################################### acupuncture data ##################################################################

#recode treatment groups to 1,2..
#acupuncture$treat<-ifelse(acupuncture$treat=="control",1,ifelse(acupuncture$treat=="acupuncture",2,9))

impacupJ2R_3covar <- mimix(c("head_base","sex","age"),"head","treat","id","time",5,1,data="acupuncture","J2R",201,c("jeffreys"),1000,NULL,NULL)
varlist <- c("head.3","head.12","head_base","sex","age")
analyselist("100",impacupJ2R_3covar,varlist)

fit<-with(data= as.mids(impacupJ2R_3covar, .id="SNO",.imp="II"), expr = lm(head.12~treat+head_base+age+sex))
summary(pool(fit))

#regressimp(impacupJ2R_3covar[impacupJ2R_3covar$II!="0",],"head.12~treat+head_base+age+sex")

###################################### J2R no covars
impacupJ2R_nocovars <- mimix(NULL,"head","treat","id","time",5,1,data="acupuncture","J2R",201,c("jeffreys"),1000,NULL,NULL)
varlist <- c("head.3","head.12")
analyselist("100",impacupJ2R_nocovars,varlist)
fit<-with(data= as.mids(impacupJ2R_nocovars, .id="SNO",.imp="II"), expr = lm(head.12~treat))
summary(pool(fit))


###################################### LMCF
impacupLMCF_2covar <- mimix(c("head_base","sex"),"head","treat","id","time",5,1,data="acupuncture","LMCF",201,,1000,NULL,NULL)
varlist <- c("head.3","head.12","head_base","sex")
analyselist("100",impacupLMCF_2covar,varlist)


###################################### Causal
impacupCausal <- mimix(c("head_base","sex"),"head","treat","id","time",5,1,data="acupuncture","Cuasal",201,c("jeffreys"),1000,NULL,NULL,NULL,2)
varlist <- c("head.3","head.12","head_base","sex")
analyselist("100",impacupCausal,varlist)



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
impantiCausal_2NoDelta <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"Causal",101,c("jeffreys"),1000,NULL,NULL,NULL,,,1,0.7)
impantiCausal_2 <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"Causal",101,c("jeffreys"),1000,NULL,NULL,NULL,c(0.5,0.5,1,1 ),c(-3,-4,-5,-6),1,2)
impantiCausal_1 <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",5,1,"Causal",101,c("jeffreys"),1000,NULL,NULL,,1)
impantiJ2R <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"J2R",101,c("jeffreys"),1000,NULL,NULL,,1)
impantiCIR <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"CIR",101,c("jeffreys"),1000,NULL,NULL,NULL)
impantiCR <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"CR",101,c("jeffreys"),1000,NULL,NULL,NULL,NULL)
impantiLMCF <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"LMCF",101,c("jeffreys"),1000,NULL,NULL,,1)
impantiMAR <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"MAR",101,c("jeffreys"),1000,NULL,NULL,,1)
impantiCausal <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"Causal",101,c("jeffreys"),1000,NULL,NULL,NULL,,,1,0)



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


impdata= as.mids(impantiJ2R)
fit<-with(data= as.mids(impantiJ2R),expr = lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+PATIENT.SEX))
summary(pool(fit))          

impantiJ2Rlist<-subset(impantiJ2R,impantiJ2R$.imp!=0)
b.out<-NULL
se.out<-NULL
for( m in 1:M ){
  ols.out<-lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+PATIENT.SEX,data=impantiJ2Rlist )
  b.out<-rbind(b.out,ols.out$coef)
  se.out<-rbind(se.out,coef(summary(ols.out))[,2])
}
combined.results<-mi.meld(q=b.out,se=se.out)
print(combined.results)
# same ests and ses as mice ! just not p values
print(combined.results)
# not works on non amelia class but can output to stata anyway! 
write.amelia(impantiJ2R, separate = FALSE, file.stem = "outdata", format = "dta")
missmap(impantiJ2R)
missmap(subset(impantiJ2R,impantiJ2R$.imp==0))

library(mitools)
antidep.implist<-imputationList(impantiJ2Rlist)
# not works
fitmi<-with(antidep.implist,lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+PATIENT.SEX))

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

impantiJ2R_3covar <- mimix("antidepressant",c("basval","POOLED.INVESTIGATOR","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",5,1,"J2R",101,c("jeffreys"),1000,NULL,NULL,NULL,,0.5,c(0.5,0.5,1,1),c(0.1,0.1,0.2,0.2))
varlist <-c("HAMD17.TOTAL.4","HAMD17.TOTAL.5","HAMD17.TOTAL.6","HAMD17.TOTAL.7")
analyselist("2721",impantiJ2R_3covar,varlist)
# NOTE pooled.investigator has a 999 outlier in it
fit<-with(data= as.mids(impantiJ2R_3covar), expr = lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+PATIENT.SEX))
summary(pool(fit))

mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"Causal",101,c("jeffreys"),1000,NULL,NULL,NULL,,1,0.7)

###################### INDIVIDUAL-SPECIFIC method (plus delta)   ##############################################################
impantiIndivDt <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",5,NULL,NULL,101,c("jeffreys"),1000,NULL,"methodcol","referencecol",c(0.5,0.5,1,1),c(-3,-4,-5,-6))

impantiIndiv <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",5,1,,101,c("jeffreys"),1000,NULL,"methodvar","referencevar",)
varlist <-c("HAMD17.TOTAL.4","HAMD17.TOTAL.5","HAMD17.TOTAL.6","HAMD17.TOTAL.7")
analyselist("2721",impantiIndivDt,varlist)
analyselist("2721",impantiIndiv,varlist)
fit<-with(data= as.mids(impantiIndivDt, .id="SNO",.imp="II"),expr=lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+PATIENT.SEX))
summary(pool(fit))

# check original values unchanged



## method with Causal option so temporary add Causal into methodvar
##antidepressant$methodvar <- ifelse(antidepressant$PATIENT.NUMBER=="1507", "Causal", antidepressant$methodvar)
impantiIndivCausal <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",5,1,,101,c("jeffreys"),1000,NULL,"methodcol","referencecol",NULL,NULL,1,0.5)
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
save(antidepressant, file="data/antidepressant.RData",version=2)  
save(asthma, file="data/asthma.RData",version=2)   
# check non missins are not imputed !
  
}



#' mimix: A package for Reference-based imputation for longitudinal clinical trials with protocol deviation
#' 
#' similar to the Stata mimix function
#' 
#' The mimix package contains the functions preprodata and preproIndivdata to 
#'  process long longitudinal data into wide data format
#'  
#'  Also the function Addelta to add delta adjustment to the imputed estimates 
#' @docType package
#' @name mimix
NULL     

# so no almost reday for mice!
#fit<-with(data= as.mids(impanticausalun, .id="SNO",.imp="II"), expr = lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+POOLED.INVESTIGATOR+PATIENT.SEX))
#summary(pool(fit))

#convert to mids object (in mice package)
#impdatamids<-as.mids(impdatasets)
#impdatamids$method[IndivDtmids$method %in% "pmm"] <- "mimix"
#impdatamids$predictorMatrix<-"N/A"


# test whether causal(0,0) gives same as J2R.
impantiCausal_2NoDelta00 <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"Causal",101,c("jeffreys"),1000,NULL,NULL,NULL,,,0,0)
varlist <-c("HAMD17.TOTAL.4","HAMD17.TOTAL.5","HAMD17.TOTAL.6","HAMD17.TOTAL.7")
analyselist("3618",impantiCausal_2,varlist)
analyselist("2104",impantiCausal_2,varlist)
analyselist("2230",impantiCausal_2,varlist)
impantiJ2R <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"J2R",101,c("jeffreys"),1000,NULL,NULL,NULL,,1)
analyselist("3618",impantiJ2R ,varlist)
analyselist("2104",impantiJ2R ,varlist)
analyselist("2230",impantiJ2R ,varlist)
#check whether same results 
all.equal(impantiCausal_2NoDelta00,impantiJ2R)

# original data set with missing values in long format so select out time point
antidepHamd17vis5<-(subset(antidepressant[,c("HAMD17.TOTAL","VISIT.NUMBER")],VISIT.NUMBER==5))
antidepHamd17vis6<-(subset(antidepressant[,c("HAMD17.TOTAL","VISIT.NUMBER")],VISIT.NUMBER==6))
antidepHamd17vis7<-(subset(antidepressant[,c("HAMD17.TOTAL","VISIT.NUMBER")],VISIT.NUMBER==7))
antidepHamd17vis4<-(subset(antidepressant[,c("HAMD17.TOTAL","VISIT.NUMBER")],VISIT.NUMBER==4))


# then concatenate as much as M times
testimp <- function(dataf,imputdata,M,var,TOTALvar)
    {
     dummydataf<-dataf
     for ( i in 1:M) {
      dummydataf<-rbind( dummydataf,dataf)
     }
     # browser()  
   datacomp<-cbind(dummydataf,subset(imputdata,select=c(var,".imp")))
   datacompna<-(subset(datacomp,!is.na(datacomp[,TOTALvar])))
   # this compares original values with imputed dat set  should be identical!! 
   cat("\n",paste(var,"non-missing = ",(nrow(datacompna)/(M+1))))
   if (all(datacompna[,1]==datacompna[,3]) ) {cat("\nnon-missing data unchanged")} 
   
   # need to check the nas have imputed values so se;ect out rows with nas
   datacompjustna<-subset(datacomp,(is.na(datacomp[,TOTALvar]) & .imp!=0))
   # check this should be 0 if all imputed.
    cat("\n", "number nas = ",paste(nrow(datacompjustna)/M))
   if (sum(is.na(datacompjustna[,3]))==0) {cat("\n",paste(var,"all imputed"))} 
  #if (all.equal(datacompna[,1],datacompna[,3]) ) {cat("original data unchanged")}   
     # return(dummydataf)  
}
testimp(antidepHamd17vis4,impantiJ2R,100,"HAMD17.TOTAL.4","HAMD17.TOTAL")
testimp(antidepHamd17vis5,impantiJ2R,100,"HAMD17.TOTAL.5","HAMD17.TOTAL")
testimp(antidepHamd17vis6,impantiJ2R,100,"HAMD17.TOTAL.6","HAMD17.TOTAL")
testimp(antidepHamd17vis7,impantiJ2R,100,"HAMD17.TOTAL.7","HAMD17.TOTAL")

# now take this to compare with imputed data set same column  
datacomp<-cbind(concatfun(antidepHamd17vis5,100),subset(impantiJ2R,select=c("HAMD17.TOTAL.5",".imp")))
# now exclude rows with na's 
# and compare cols after taking out NAs
datacompna<-(subset(datacomp,!is.na(HAMD17.TOTAL)))
# this compares original values with imputed dat set  should be identical!! 
if (all.equal(datacompna[,1],datacompna[,3]) ) {cat("original data unchanged")}   
# need to check the nas have imputed values so se;ect out rows with nas
datacompjustna<-(subset(datacomp,is.na(HAMD17.TOTAL) & .imp!=0))
# check this should be 0 if all imputed.
if (sum(is.na(datacompjustna[,3]))==0) {cat("all imputed")} 




all.equal(datacompjustna[,1],datacompjustna[,3])      

comparefun <- function(datac) {
      ifelse(datac[1] != datac[3]), {
        cat("\nNot equal")
      } 
}
  

all.equal(antidepHamd17vis5,subset(impantiJ2R,select=c("HAMD17.TOTAL.5"))
datavisit5<-cbind(antidepHamd17vis5,subset(impantiJ2R,select=c("HAMD17.TOTAL.5")))

# now test differnce between ttal and total.5
testdiff<-function(data,testvar)  {
   all.equal(data[,tstdiff],
             
# check orig values same after imputation
#origdata<-as.data.frame(subset(impdatasetJ2Rdelta,.imp==0))
#identical(all.equal(origdata,subset(impdatasetJ2Rdelta,.imp==1)),TRUE)
# for( m in 1:M ) {
#impdata_m<-as.data.frame(subset(impdatasetJ2Rdelta,.imp==m))
# cat("\n",identical(all.equal(origdata$fev.2,impdata_m$fev.2),TRUE))
# cat("\n",identical(all.equal(origdata$fev.4,impdata_m$fev.4),TRUE))
             }             
             
