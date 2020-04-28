# testruns and examples for mimix program

######################################### asthma data-set testing all methods #############################################

impdatasetJ2R<-(mimix("asthma",c("base"),"fev","treat","id","time",2,1,"J2R",101,"jeffreys",1000,NULL,NULL,c(0.5,0.5,1,1 ) ))
varlist <- c("fev.2","fev.4","fev.8","fev.12","base")
analyselist(5017,impdatasetJ2R,varlist)
regressimp(impdatasetJ2R,"fev.12~treat+base")


impdatasetJ2Rridge<-(mimix("asthma",c("base"),"fev","treat","id","time",2,1,"J2R",101,"ridge",1000,NULL,NULL,NULL ) )
analyselist(5017,impdatasetJ2Rridge,varlist)
regressimp(impdatasetJ2Rridge,"fev.12~treat+base")


impdatasetJ2Rridge<-(mimix("asthma",c("base"),"fev","treat","id","time",2,1,"J2R",101,"ridge",1000,NULL,NULL,c(0.5,0.5,1,1 ) ) )
analyselist(5017,impdatasetJ2Rridge,varlist)
regressimp(impdatasetJ2Rridge,"fev.12~treat+base")

#test all the  methods at least run ok

impdatasetMARridge<-(mimix("asthma",c("base"),"fev","treat","id","time",2,1,"MAR",101,"ridge",1000,NULL,NULL,NULL ) )
analyselist(5017,impdatasetMARridge,varlist)
regressimp(impdatasetMARridge,"fev.12~treat+base")

impdatasetCIRridge<-(mimix("asthma",c("base"),"fev","treat","id","time",2,1,"CIR",101,"ridge",1000,NULL,NULL,NULL ) )
analyselist(5017,impdatasetCIRridge,varlist)
regressimp(impdatasetCIRridge,"fev.12~treat+base")

impdatasetCRridge<-(mimix("asthma",c("base"),"fev","treat","id","time",2,1,"CR",101,"ridge",1000,NULL,NULL,NULL ) )
analyselist(5017,impdatasetCRridge,varlist)
regressimp(impdatasetCRridge,"fev.12~treat+base")

impdatasetLMCFridge<-(mimix("asthma",c("base"),"fev","treat","id","time",2,1,"LMCF",101,"ridge",1000,NULL,NULL,NULL ) )
analyselist(5017,impdatasetLMCFridge,varlist)
regressimp(impdatasetLMCFridge,"fev.12~treat+base")


impdatasetCausalridge<-(mimix("asthma",c("base"),"fev","treat","id","time",2,1,"Causal",101,"ridge",1000,NULL,NULL,NULL ) )
# Causal asks for interactive  Kd value
analyselist(5017,impdatasetCausalridge,varlist)
regressimp(impdatasetCausalridge,"fev.12~treat+base")


####################################### acupuncture data ##################################################################


#recode treatment groups to 1,2..
#acupuncture$treat<-ifelse(acupuncture$treat=="control",1,ifelse(acupuncture$treat=="acupuncture",2,9))

impacupJ2R <- mimix("acupuncture",c("head_base","sex","age"),"head","treat","id","time",100,1,"J2R",201,c("jeffreys"),1000,NULL,NULL)
varlist <- c("head.3","head.12","head_base","sex","age")
analyselist("151",impacupJ2R,varlist)
regressimp(impacupJ2R,"head.12~treat+head_base+age+sex")





####################################### antidepressant data ##################################################################

#recode treatment variable to numeric
#antidepressant$TREATMENT.NAME<-ifelse(antidepressant$TREATMENT.NAME=="PLACEBO",1,ifelse(antidepressant$TREATMENT.NAME=="DRUG",2,99))

#recode covarates to numeric
antidepressant$PATIENT.SEX <- as.numeric(antidepressant$PATIENT.SEX)

impanticausal <- mimix("antidepressant",c("basval","POOLED.INVESTIGATOR","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"Causal",101,c("jeffreys"),1000,NULL,NULL)

varlist <-c("HAMD17.TOTAL.4","HAMD17.TOTAL.5","HAMD17.TOTAL.6","HAMD17.TOTAL.7")
analyselist("2721",impanticausal,varlist)
regressimp(impanticausal,"HAMD17.TOTAL.7~TREATMENT.NAME+basval+POOLED.INVESTIGATOR+PATIENT.SEX")



# individual specific method (plus delta)
impantiIndiv <- mimix("antidepressant",c("basval","POOLED.INVESTIGATOR","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,NULL,101,c("jeffreys"),1000,NULL,c("methodvar","referencevar"),c(0.5,0.5,1,1 ))
analyselist("2721",impantiIndiv,varlist)
regressimp(impantiIndiv,"HAMD17.TOTAL.7~TREATMENT.NAME+basval+POOLED.INVESTIGATOR+PATIENT.SEX")








