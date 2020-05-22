# testruns and examples for mimix program

######################################### asthma data-set testing J2R method #############################################
# testing all options, m=5 imputations , method = J2R, referEnce treatment gp =1,seed = 101,  Delta specified
impdatasetJ2R<-mimix("asthma",c("base"),"fev","treat","id","time",5,1,"J2R",101,"jeffreys",1000,NULL,NULL,c(0.5,0.5,1,1 ),NULL )

# convert imputed dataset to mids format (from mice package)
impdata= as.mids(impdatasetJ2R, .id="SNO",.imp="II")

# fit specified model to each imputed data set and pool results together (Rubin's rules), functions from mice package
fit<-with(impdata, lm(fev.12~treat+base))
summary(pool(fit))

# previously regressimp(impdatasetJ2R[impdatasetJ2R$II!="0",],"fev.12~treat+base")

# to examine estimates for individual ids
varlist <- c("fev.2","fev.4","fev.8","fev.12","base")
# print estimates for individuals
analyselist(5099,impdatasetJ2R,varlist)
analyselist(5456,impdatasetJ2R,varlist)

##################################################### no delta option #########################################
impdatasetJ2R_Nd<-mimix("asthma",c("base"),"fev","treat","id","time",1000,1,"J2R",301,"jeffreys",1000,NULL,NULL,, )
analyselist(5099,impdatasetJ2R_Nd,varlist)
fit<-with(data= as.mids(impdatasetJ2R_Nd, .id="SNO",.imp="II"), expr = lm(fev.12~treat))
summary(pool(fit))

################################################ differnt priors ################################################
impdatasetJ2jeffreys<-(mimix("asthma",c("base"),"fev","treat","id","time",5,1,"J2R",101,"jeffreys",1000,NULL,NULL,NULL,0.5 ) )
analyselist(5099,impdatasetJ2jeffreys,varlist)


impdatasetJ2Rridge<-(mimix("asthma",c("base"),"fev","treat","id","time",5,1,"J2r",101,"ridge",1000,NULL,NULL,c(0.5,0.5,1,1 ) ) )
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
impdatasetCIRridge<-(mimix("asthma",c("base"),"fev","treat","id","time",5,1,"CIR",101,"ridge",1000,NULL,NULL, ) )
analyselist(5099,impdatasetCIR,varlist)
analyselist(5456,impdatasetCIR,varlist)
analyselist(5456,impdatasetCIRridge,varlist)

######################################## CR
impdatasetCR<-(mimix("asthma",c("base"),"fev","treat","id","time",5,1,"CR",101,,1000,NULL,NULL,NULL ) )
impdatasetCRridge<-(mimix("asthma",c("base"),"fev","treat","id","time",5,1,"CR",101,"ridge",1000,NULL,NULL,NULL ) )
analyselist(5099,impdatasetCR,varlist)
analyselist(5017,impdatasetCR,varlist)
analyselist(5099,impdatasetCRridge,varlist)
analyselist(5017,impdatasetCRridge,varlist)

########################################## LMCF
impdatasetLMCF<-(mimix("asthma",c("base"),"fev","treat","id","time",5,1,"LMCF",101,,1000,NULL,NULL,NULL ) )
impdatasetLMCFridge<-(mimix("asthma",c("base"),"fev","treat","id","time",5,1,"LMCF",101,"ridge",1000,NULL,NULL,NULL ) )
analyselist(5456,impdatasetLMCF,varlist)
analyselist(5099,impdatasetLMCF,varlist)
analyselist(5456,impdatasetLMCFridge,varlist)
analyselist(5099,impdatasetLMCFridge,varlist)


##########################################  Causal method
# ERROR!!! induced because kd not specified so give errr msg
impdatasetCausalridge<-(mimix("asthma",c("base"),"fev","treat","id","time",5,1,"Causal",101,"ridge",1000,NULL,NULL,NULL, ) )
impdatasetCausal<-(mimix("asthma",c("base"),"fev","treat","id","time",5,1,"Causal",101,,1000,NULL,NULL,NULL,0.5 ) )

analyselist(5099,impdatasetCausal,varlist)



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

#recode treatment variable to numeric
#antidepressant$TREATMENT.NAME<-ifelse(antidepressant$TREATMENT.NAME=="PLACEBO",1,ifelse(antidepressant$TREATMENT.NAME=="DRUG",2,99))
#recode covarates to numeric
antidepressant$PATIENT.SEX <- as.numeric(antidepressant$PATIENT.SEX)
# NOTE in input data pooled.investigator has 999 outlier

############################################# Test Causal method, when Kd = 0 same as J2R, Kd=1 same as CIR
impantiCausal_2NoDelta <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",5,1,"Causal",101,c("jeffreys"),1000,NULL,NULL,,2)
impantiCausal_2 <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",5,1,"Causal",101,c("jeffreys"),1000,NULL,NULL,c(0.5,0.5,1,1 ),2)
impantiCausal_1 <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",5,1,"Causal",101,c("jeffreys"),1000,NULL,NULL,,1)
impantiJ2R <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",5,1,"J2R",101,c("jeffreys"),1000,NULL,NULL,,1)
impantiCIR <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",5,1,"CIR",101,c("jeffreys"),1000,NULL,NULL,,1)


varlist <-c("HAMD17.TOTAL.4","HAMD17.TOTAL.5","HAMD17.TOTAL.6","HAMD17.TOTAL.7")
analyselist("2721",impantiCausal_2,varlist)
analyselist("2721",impantiCausal_2NoDelta,varlist)
analyselist("2721",impantiCausal_1,varlist)
analyselist("2721",impantiJ2R,varlist)
analyselist("2721",impantiCIR,varlist)
fit<-with(data= as.mids(impantiCausal_2, .id="SNO",.imp="II"), expr = lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+PATIENT.SEX))
summary(pool(fit))


impantiJ2R_3covar <- mimix("antidepressant",c("basval","POOLED.INVESTIGATOR","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",5,1,"J2R",101,c("jeffreys"),1000,NULL,NULL,,0.5)
varlist <-c("HAMD17.TOTAL.4","HAMD17.TOTAL.5","HAMD17.TOTAL.6","HAMD17.TOTAL.7")
analyselist("2721",impantiJ2R_3covar,varlist)
# NOTE pooled.investigator has a 999 outlier in it



###################### INDIVIDUAL-SPECIFIC method (plus delta)   ##############################################################
impantiIndivDt <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",5,1,NULL,101,c("jeffreys"),1000,NULL,c("methodvar","referencevar"),c(0.5,0.5,1,1 ))

impantiIndiv <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",5,1,,101,c("jeffreys"),1000,NULL,c("methodvar","referencevar"),)
varlist <-c("HAMD17.TOTAL.4","HAMD17.TOTAL.5","HAMD17.TOTAL.6","HAMD17.TOTAL.7")
analyselist("2721",impantiIndivDt,varlist)
analyselist("2721",impantiIndiv,varlist)
fit<-with(data= as.mids(impantiIndivDt, .id="SNO",.imp="II"),expr=lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+PATIENT.SEX))
summary(pool(fit))



#################################   with/wout covars
impantiIndivNOcov <- mimix("antidepressant",,"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",5,1,,101,c("jeffreys"),1000,NULL,c("methodvar","referencevar"))
analyselist("2721",impantiIndivNOcov,varlist)
fit<-with(data= as.mids(impantiIndivNOcov, .id="SNO",.imp="II"),expr=lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+PATIENT.SEX))
summary(pool(fit))







