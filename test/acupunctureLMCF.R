###########################################################################
# test program to  compare R on acupunture data focusing on LmFCF option  ~
# uses acupunture data and when covariate acts as base response           #
###########################################################################
 # uses analyselist fun
analyselist <-function(id,datlist,varlist) {
  #browser()
  datano <- subset(datlist,id==datlist$.id)
  cat(paste0("\ncase = ",id))
  cat(paste0("\n treatarm = ",subset(datano$treat,datano$.imp==0),"\n"))
  # numbers denote the descriptive stats to display
  t(round(pastecs::stat.desc(datano)[,varlist],8)[c(1,9,13,4,8,5),])
}

#  
# create time=0 base response  value
# ie convert basval to response  value at visit 0

testacup<-acupuncture 
testacup3<-(subset(testacup,testacup$time==3))
testacup3$head<- testacup3$head_base
testacup3$time <- 0
testacupunt3<-(rbind(testacup3,testacup))
head(testacupunt3[order(testacupunt3[,"id"],testacupunt3[,"time"]),])
# need to order  sts4Dpatt[order(sts4Dpatt[,treatvar],sts4Dpatt[,"patt"]),]
acupunturetime0<- testacupunt3[order(testacupunt3[,"id"],testacupunt3[,"time"]),]

# control is 1 acupuncture is 2
# 0702 hav established lmcf ok when no covar ie same as stata nd sas so try covar as base value ie iuse acupunturetime0
# compare stata 3/2/20
# ref 1 same as 0 in stata
# testted no diff in reg group for lmcf! 4/02/20
impdatasetLMCFAcupBase<-mimix("acupuncture",c("head_base"),"head","treat","id","time",10,1,"LMCF",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,K0=1,K1=0.5,mle=0  )
sink(file("analyseLMCFacupuntBase.csv"),append=TRUE)
fit<-with(data= as.mids(impdatasetLMCFAcupBase), expr = lm(head.12~treat+head_base))
#sink(file("mimixheadBaseLMCFAcup.csv"),append=TRUE)
summary(pool(fit))
varlist <- c("head.3","head.12");
analyselist(104,impdatasetLMCFAcupBase,varlist)
analyselist(333,impdatasetLMCFAcupBase,varlist)
analyselist(386,impdatasetLMCFAcupBase,varlist)
analyselist(787,impdatasetLMCFAcupBase,varlist)
analyselist(435,impdatasetLMCFAcupBase,varlist)
analyselist(697,impdatasetLMCFAcupBase,varlist)
analyselist(100,impdatasetLMCFAcupBase,varlist)
analyselist(101,impdatasetLMCFAcupBase,varlist)
sink()

impdatasetLMCFAcupnoBase<-mimix("acupuncture",,"head","treat","id","time",1000,2,"LMCF",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,K0=1,K1=1,mle=0  )
sink(file("LMCFacupuntnoBase.csv"),append=TRUE)
fit<-with(data= as.mids(impdatasetLMCFAcupnoBase), expr = lm(head.12~treat))
summary(pool(fit))
varlist <- c("head.3","head.12");
analyselist(104,impdatasetLMCFAcupnoBase,varlist)
analyselist(333,impdatasetLMCFAcupnoBase,varlist)
analyselist(386,impdatasetLMCFAcupnoBase,varlist)
analyselist(787,impdatasetLMCFAcupnoBase,varlist)
analyselist(435,impdatasetLMCFAcupnoBase,varlist)
analyselist(697,impdatasetLMCFAcupnoBase,varlist)
analyselist(100,impdatasetLMCFAcupnoBase,varlist)
analyselist(101,impdatasetLMCFAcupnoBase,varlist)
sink()
impdatasetMARAcupnoBase<-mimix("acupuncture",,"head","treat","id","time",1000,,"MAR",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,K0=1,K1=1,mle=0  )
#sink(file("analyseLMCFacupuntBase.csv"),append=TRUE)
fit<-with(data= as.mids(impdatasetMARAcupnoBase), expr = lm(head.12~treat))
summary(pool(fit))
varlist <- c("head.3","head.12");
analyselist(104,impdatasetMARAcupnoBase,varlist)
analyselist(333,impdatasetMARAcupnoBase,varlist)
analyselist(386,impdatasetMARAcupnoBase,varlist)
analyselist(787,impdatasetMARAcupnoBase,varlist)
analyselist(435,impdatasetMARAcupnoBase,varlist)
analyselist(697,impdatasetMARAcupnoBase,varlist)
analyselist(100,impdatasetMARAcupnoBase,varlist)
analyselist(101,impdatasetMARAcupnoBase,varlist)



impdatasetCausalACcup09<-mimix("acupunturetime0",NULL,"head","treat","id","time",1000,2,"Causal",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,K0=1,K1=0.9,mle=0  )
library(mice)
sink(file("analyseCausalJ2Rref2.csv"),append=TRUE)
fit<-with(as.mids(impdatasetCausalACcup09), lm(head.12~treat))
summary(pool(fit))
varlist <- c("head.3","head.12")
analyselist(100,impdatasetCausalACcup09,varlist)
analyselist(101,impdatasetCausalACcup09,varlist)
analyselist(104,impdatasetCausalACcup09,varlist)
analyselist(105,impdatasetCausalACcup09,varlist)
analyselist(138,impdatasetCausalACcup09,varlist)
analyselist(139,impdatasetCausalACcup09,varlist)
sink()

impdatasetJ2RCcupBase<-mimix("acupuncture",c("head_base","sex","age"),"head","treat","id","time",1000,1,"Causal",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,K0=1,K1=0.5,mle=0  )

fit<-with(as.mids(impdatasetJ2RCcupBase), lm(head.12~treat))
summary(pool(fit))

sink()
# compare stata 3/2/20
# ref arm 1 same as 0 in stata
sink(file("analyseJ2AcupRef.csv"),append=TRUE)
impdatasetCIRAcupBaseref1<-mimix("acupuncture",c("head_base"),"head","treat","id","time",10,1,"CIR",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,K0=1,K1=0.5,mle=0  )
library(mice)
sink(file("mimixJ2AcupBaseref1.csv"),append=TRUE)
fit<-with(data= as.mids(impdatasetCIRAcupBaseref1), expr = lm(head.12~treat+head_base))
summary(pool(fit))
varlist <- c("head.3","head.12");
analyselist(104,impdatasetJ2RAcupBaseref1,varlist)
analyselist(333,impdatasetJ2RAcupBaseref1,varlist)
analyselist(386,impdatasetJ2RAcupBaseref1,varlist)
analyselist(787,impdatasetJ2RAcupBaseref1,varlist)
analyselist(435,impdatasetJ2RAcupBaseref1,varlist)
analyselist(697,impdatasetJ2RAcupBaseref1,varlist)
analyselist(100,impdatasetJ2RAcupBaseref1,varlist)
analyselist(101,impdatasetJ2RAcupBaseref1,varlist)
sink()


# save so to red in stata
write.csv(acupunturetime0,'acupunturetime0.csv')

# compare with  time0  ie head_base is resp at time 0
sink(file("time0LMCFAcupRef1.csv"),append=TRUE)
#bbetween=100 doesnt make not much differnce  
# so try deleting all those cases with all missing values
#lapply(acupuncture, if  

#impdatasetLMCFAcup0<-mimix("acupunturetime0",,"head","treat","id","time",1000,1,"LMCF",54321,"jeffreys",1000,bbetween=100,NULL,NULL,NULL,NULL,K0=1,K1=1,mle=0  )
impdatasetLMCFAcup0<-mimix("acupunturetime0",,"head","treat","id","time",1000,2,"LMCF",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,K0=1,K1=1,mle=0  )
library(mice)
sink(file("mimixtime0LMCFAcupt0.csv"),append=TRUE)
fit<-with(data= as.mids(impdatasetLMCFAcup0), expr = lm(head.12~treat))
summary(pool(fit))
fit<-with(data= as.mids(impdatasetLMCFAcup0), expr = lm(head.12~treat+head_base))
summary(pool(fit))
varlist <- c("head.3","head.12");
analyselist(104,impdatasetLMCFAcup0,varlist)
analyselist(333,impdatasetLMCFAcup0,varlist)
analyselist(386,impdatasetLMCFAcup0,varlist)
analyselist(787,impdatasetLMCFAcup0,varlist)
analyselist(435,impdatasetLMCFAcup0,varlist)
analyselist(697,impdatasetLMCFAcup0,varlist)
analyselist(100,impdatasetLMCFAcup0,varlist)
analyselist(101,impdatasetLMCFAcup0,varlist)
sink()

# so now trea effects closer but case 100 still large in r !  try deleting all missing va;ues ! next    



# in stata ref grp 1=acupuncture, 0=control, so in r equiv to 1  is 2
impdatasetCIRAcupBaseref2<-mimix("acupuncture",c("head_base"),"head","treat","id","time",1000,2,"CIR",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,K0=1,K1=1,mle=0  )
sink(file("mimixCIRAcupuntBaseRef2.csv"),append=TRUE)
fit<-with(data= as.mids(impdatasetCIRAcupBaseref2), expr = lm(head.12~treat+head_base))
summary(pool(fit))
varlist <- c("head.3","head.12");
analyselist(104,impdatasetCIRAcupBaseref2,varlist)
analyselist(333,impdatasetCIRAcupBaseref2,varlist)
analyselist(386,impdatasetCIRAcupBaseref2,varlist)
analyselist(787,impdatasetCIRAcupBaseref2,varlist)
analyselist(435,impdatasetCIRAcupBaseref2,varlist)
analyselist(697,impdatasetCIRAcupBaseref2,varlist)
analyselist(100,impdatasetCIRAcupBaseref2,varlist)
analyselist(101,impdatasetCIRAcupBaseref2,varlist)
sink()

impdatasetCIRAcupBaseref1<-mimix("acupuncture",c("head_base"),"head","treat","id","time",1000,1,"CIR",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,K0=1,K1=1,mle=0  )
sink(file("mimixCIRAcupuntBaseRef1.csv"),append=TRUE)
fit<-with(data= as.mids(impdatasetCIRAcupBaseref1), expr = lm(head.12~treat+head_base))
summary(pool(fit))
varlist <- c("head.3","head.12");
analyselist(104,impdatasetCIRAcupBaseref1,varlist)
analyselist(333,impdatasetCIRAcupBaseref1,varlist)
analyselist(386,impdatasetCIRAcupBaseref1,varlist)
analyselist(787,impdatasetCIRAcupBaseref1,varlist)
analyselist(435,impdatasetCIRAcupBaseref1,varlist)
analyselist(697,impdatasetCIRAcupBaseref1,varlist)
analyselist(100,impdatasetCIRAcupBaseref1,varlist)
analyselist(101,impdatasetCIRAcupBaseref1,varlist)
sink()