########################################################## Test program for mimix v009 #################################################  
# for outputs find in rlogs folder  
# to create log files, initiate log header which can then be used to append by using sink
#library(logr)
#log_locn<- "testCRref2seed.log"
#log_open(log_locn)
#sink(paste0("C:/GitHub/mimix/mimixR/log/",log_locn,""),append=TRUE,split=TRUE)
# operations
# to close log  
#sink()
#close_log


#analyselist <-function(id,datlist,varlist) {
  #browser()
#  datano <- subset(datlist,id==datlist$.id)
#  cat(paste0("\ncase = ",id))
#  cat(paste0("\n treatarm = ",subset(datano$treat,datano$.imp==0),"\n"))
  # numbers denote the descriptive stats to display
#  t(round(pastecs::stat.desc(datano)[,varlist],3)[c(1,9,13,4,8,5),])
#}

######################################## compare  methods same imputed values in  ref arm for comparing main methods  J2r CIR CR MAR #####

# set up log header 
log_locn<- "testCRref1test2.log"
log_open(log_locn)
#log_print('impdatasetJ2Rref1 and ,CIR,CR to test reference arm estimates the same')
# can use after close sink()  log_close()
sink(paste0("C:/GitHub/mimix/mimixR/log/",log_locn,""),append=TRUE,split=TRUE)

#cat("
impdatasetJ2Rref2<-mimix("asthma",c("base"),"fev","treat","id","time",10,2,"J2R",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
impdatasetCIRref2<-mimix("asthma",c("base"),"fev","treat","id","time",10,2,"CIR",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
impdatasetCRref2<-mimix("asthma",c("base"),"fev","treat","id","time",10,2,"CR",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
impdatasetMARref2<-mimix("asthma",c("base"),"fev","treat","id","time",10,2,"MAR",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
#\n")
varlist <- c("fev.2","fev.4","fev.8","fev.12")

#cat("
# 5333 an interim
analyselist(5333,impdatasetJ2Rref2,varlist)
analyselist(5333,impdatasetCIRref2,varlist)
analyselist(5333,impdatasetCRref2,varlist)
analyselist(5333,impdatasetMARref2,varlist)
#\n")
analyselist(5051,impdatasetJ2Rref2,varlist)
analyselist(5051,impdatasetCIRref2,varlist)
analyselist(5051,impdatasetCRref2,varlist)
analyselist(5051,impdatasetMARref2,varlist)

#5115 interim
analyselist(5115,impdatasetJ2Rref2,varlist)
analyselist(5115,impdatasetCIRref2,varlist)
analyselist(5115,impdatasetCRref2,varlist)
analyselist(5115,impdatasetMARref2,varlist)


analyselist(5094,impdatasetMARref2,varlist)
analyselist(5094,impdatasetJ2Rref2,varlist)
analyselist(5094,impdatasetCIRref2,varlist)
analyselist(5094,impdatasetCRref2,varlist)

analyselist(5017,impdatasetMARref2,varlist)
analyselist(5017,impdatasetJ2Rref2,varlist)
analyselist(5017,impdatasetCIRref2,varlist)
analyselist(5017,impdatasetCRref2,varlist)


impdatasetJ2Rref1<-mimix("asthma",c("base"),"fev","treat","id","time",10,1,"J2R",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
impdatasetCIRref1<-mimix("asthma",c("base"),"fev","treat","id","time",10,1,"CIR",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
impdatasetCRref1<-mimix("asthma",c("base"),"fev","treat","id","time",10,1,"CR",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
impdatasetMARref1<-mimix("asthma",c("base"),"fev","treat","id","time",10,1,"MAR",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
#cat("
analyselist(5333,impdatasetJ2Rref1,varlist)
analyselist(5333,impdatasetCIRref1,varlist)
analyselist(5333,impdatasetCRref1,varlist)
analyselist(5333,impdatasetMARref1,varlist)

analyselist(5051,impdatasetJ2Rref1,varlist)
analyselist(5051,impdatasetCIRref1,varlist)
analyselist(5051,impdatasetCRref1,varlist)
analyselist(5051,impdatasetMARref1,varlist)

analyselist(5115,impdatasetJ2Rref1,varlist)
analyselist(5115,impdatasetCIRref1,varlist)
analyselist(5115,impdatasetCRref1,varlist)
analyselist(5115,impdatasetMARref1,varlist)

analyselist(5094,impdatasetMARref1,varlist)
analyselist(5094,impdatasetJ2Rref1,varlist)
analyselist(5094,impdatasetCIRref1,varlist)
analyselist(5094,impdatasetCRref1,varlist)

analyselist(5017,impdatasetMARref1,varlist)
analyselist(5017,impdatasetJ2Rref1,varlist)
analyselist(5017,impdatasetCIRref1,varlist)
analyselist(5017,impdatasetCRref1,varlist)

analyselist(5044,impdatasetMARref1,varlist)
analyselist(5044,impdatasetJ2Rref1,varlist)
analyselist(5044,impdatasetCIRref1,varlist)
analyselist(5044,impdatasetCRref1,varlist)

analyselist(5059,impdatasetMARref1,varlist)
analyselist(5059,impdatasetJ2Rref1,varlist)
analyselist(5059,impdatasetCIRref1,varlist)
analyselist(5059,impdatasetCRref1,varlist)
#\n")
sink()
log_close()


#impdatasetJ2RnoD<-mimix("asthma",c("base"),"fev","treat","id","time",3,2,"J2R",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
#impdatasetJ2RnoC<-mimix("asthma",,"fev","treat","id","time",3,2,"J2R",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )

sink("C:/GitHub/mimix/mimixR/log/testCRref1test2.log",append=TRUE,split=TRUE)
varlist <- c("fev.2","fev.4","fev.8","fev.12");

#cat("Some more stuff here...\n")) use cat for sink then uncomment and run cmds


####################################################  test differnt seeds ###################################################################
log_locn<- "testref1seed.log"
log_open(log_locn)

sink(paste0("C:/GitHub/mimix/mimixR/log/",log_locn,""),append=TRUE,split=TRUE)
sink()

#cat('
impdatasetJ2Rref1s101<-mimix("asthma",c("base"),"fev","treat","id","time",1000,1,"J2R",101,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
impdatasetJ2Rref1s54321<-mimix("asthma",c("base"),"fev","treat","id","time",1000,1,"J2R",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
impdatasetCIRref1s101<-mimix("asthma",c("base"),"fev","treat","id","time",1000,1,"CIR",101,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
impdatasetCIRref1s54321<-mimix("asthma",c("base"),"fev","treat","id","time",1000,1,"CIR",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
impdatasetCRref1s101<-mimix("asthma",c("base"),"fev","treat","id","time",1000,1,"CR",101,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
impdatasetCRref1s54321<-mimix("asthma",c("base"),"fev","treat","id","time",1000,1,"CR",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
#\n')


# seed on 5017
impdatasetCIRref1s101<-mimix("asthma",c("base"),"fev","treat","id","time",2000,1,"CIR",101,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
impdatasetCIRref1s54321<-mimix("asthma",c("base"),"fev","treat","id","time",2000,1,"CIR",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
impdatasetCIRref1s999<-mimix("asthma",c("base"),"fev","treat","id","time",2000,1,"CIR",999,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
impdatasetCIRref1s1394<-mimix("asthma",c("base"),"fev","treat","id","time",2000,1,"CIR",1394,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )

analyselist(5017,impdatasetCIRref1s101,varlist)
fit<-with(data= as.mids(impdatasetCIRref1s101), expr = lm(fev.12~treat+base))
summary(pool(fit))
analyselist(5017,impdatasetCIRref1s54321,varlist)
fit<-with(data= as.mids(impdatasetCIRref1s54321), expr = lm(fev.12~treat+base))
summary(pool(fit))
analyselist(5017,impdatasetCIRref1s999,varlist)
fit<-with(data= as.mids(impdatasetCIRref1s999), expr = lm(fev.12~treat+base))
summary(pool(fit))
analyselist(5017,impdatasetCIRref1s1394,varlist)
fit<-with(data= as.mids(impdatasetCIRref1s1394), expr = lm(fev.12~treat+base))
summary(pool(fit))


log_locn<- "testCIRref2test2.log"

#If <- 
log_open(log_locn)

sink(paste0("C:/GitHub/mimix/mimixR/log/",log_locn,""),append=TRUE,split=TRUE)

#cat('
impdatasetCIRref2s101<-mimix("asthma",c("base"),"fev","treat","id","time",2000,2,"CIR",101,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
impdatasetCIRref2s54321<-mimix("asthma",c("base"),"fev","treat","id","time",2000,2,"CIR",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
impdatasetCIRref2s999<-mimix("asthma",c("base"),"fev","treat","id","time",2000,2,"CIR",999,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
impdatasetCIRref2s1394<-mimix("asthma",c("base"),"fev","treat","id","time",2000,2,"CIR",1394,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
#\n')
#cat('
analyselist(5059,impdatasetCIRref1s101,varlist)
analyselist(5059,impdatasetCIRref1s54321,varlist)
analyselist(5059,impdatasetCIRref1s999,varlist)
analyselist(5059,impdatasetCIRref1s1394,varlist)
#\n')
#cat('
analyselist(5017,impdatasetCIRref2s101,varlist)
analyselist(5017,impdatasetCIRref2s54321,varlist)
analyselist(5017,impdatasetCIRref2s999,varlist)
analyselist(5017,impdatasetCIRref2s1394,varlist)
#\n')
#cat('
impdatasetCIRref2s101_1000<-mimix("asthma",c("base"),"fev","treat","id","time",1000,2,"CIR",101,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
impdatasetCIRref2s54321_1000<-mimix("asthma",c("base"),"fev","treat","id","time",1000,2,"CIR",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
impdatasetCIRref2s999_1000<-mimix("asthma",c("base"),"fev","treat","id","time",1000,2,"CIR",999,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
impdatasetCIRref2s1394_1000<-mimix("asthma",c("base"),"fev","treat","id","time",1000,2,"CIR",1394,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
#\n')
#cat('
analyselist(5017,impdatasetCIRref2s101_1000,varlist)
analyselist(5017,impdatasetCIRref2s54321_1000,varlist)
analyselist(5017,impdatasetCIRref2s999_1000,varlist)
analyselist(5017,impdatasetCIRref2s1394_1000,varlist)
#\n')
#cat('
fit<-with(data= as.mids(impdatasetCIRref2s101), expr = lm(fev.12~treat+base))
summary(pool(fit))
fit<-with(data= as.mids(impdatasetCIRref2s101_1000), expr = lm(fev.12~treat+base))
summary(pool(fit))
#\n')
sink()
log_close()

#cat("
library(mice)
fit<-with(data= as.mids(impdatasetJ2Rref1s101), expr = lm(fev.12~treat+base))
summary(pool(fit))
fit<-with(data= as.mids(impdatasetJ2Rref1s54321), expr = lm(fev.12~treat+base))
summary(pool(fit))
#\n")
#cat("
fit<-with(data= as.mids(impdatasetCIRref1s101), expr = lm(fev.12~treat+base))
summary(pool(fit))
fit<-with(data= as.mids(impdatasetCIRref1s54321), expr = lm(fev.12~treat+base))
summary(pool(fit))
#\n")
#cat("
fit<-with(data= as.mids(impdatasetCRref1s101), expr = lm(fev.12~treat+base))
summary(pool(fit))
fit<-with(data= as.mids(impdatasetCRref1s54321), expr = lm(fev.12~treat+base))
summary(pool(fit))
#\n")

#1001 tests show prgram behaves well one discrepancy o seed is CIR 5051 case
# treat effects differ slightly and diff df's though due to mice?
#cat("
analyselist(5051,impdatasetJ2Rref1s101,varlist)
analyselist(5051,impdatasetJ2Rref1s54321,varlist)
analyselist(5333,impdatasetJ2Rref1s101,varlist)
analyselist(5333,impdatasetJ2Rref1s54321,varlist)
analyselist(5094,impdatasetJ2Rref1s101,varlist)
analyselist(5094,impdatasetJ2Rref1s54321,varlist)
analyselist(5115,impdatasetJ2Rref1s101,varlist)
analyselist(5115,impdatasetJ2Rref1s54321,varlist)
analyselist(5017,impdatasetJ2Rref1s101,varlist)
analyselist(5017,impdatasetJ2Rref1s54321,varlist)
analyselist(5044,impdatasetJ2Rref1s101,varlist)
analyselist(5044,impdatasetJ2Rref1s54321,varlist)
analyselist(5059,impdatasetJ2Rref1s101,varlist)
analyselist(5059,impdatasetJ2Rref1s54321,varlist)
#\n")

#cat("
analyselist(5051,impdatasetCIRref1s101,varlist)
analyselist(5051,impdatasetCIRref1s54321,varlist)
analyselist(5333,impdatasetCIRref1s101,varlist)
analyselist(5333,impdatasetCIRref1s54321,varlist)
analyselist(5094,impdatasetCIRref1s101,varlist)
analyselist(5094,impdatasetCIRref1s54321,varlist)
analyselist(5115,impdatasetCIRref1s101,varlist)
analyselist(5115,impdatasetCIRref1s54321,varlist)
analyselist(5017,impdatasetCIRref1s101,varlist)
analyselist(5017,impdatasetCIRref1s54321,varlist)
analyselist(5044,impdatasetCIRref1s101,varlist)
analyselist(5044,impdatasetCIRref1s54321,varlist)
analyselist(5059,impdatasetCIRref1s101,varlist)
analyselist(5059,impdatasetCIRref1s54321,varlist)
#\n")

#cat("
analyselist(5051,impdatasetCRref1s101,varlist)
analyselist(5051,impdatasetCRref1s54321,varlist)
analyselist(5333,impdatasetCRref1s101,varlist)
analyselist(5333,impdatasetCRref1s54321,varlist)
analyselist(5094,impdatasetCRref1s101,varlist)
analyselist(5094,impdatasetCRref1s54321,varlist)
analyselist(5115,impdatasetCRref1s101,varlist)
analyselist(5115,impdatasetCRref1s54321,varlist)
analyselist(5017,impdatasetCRref1s101,varlist)
analyselist(5017,impdatasetCRref1s54321,varlist)
analyselist(5044,impdatasetCRref1s101,varlist)
analyselist(5044,impdatasetCRref1s54321,varlist)
analyselist(5059,impdatasetCRref1s101,varlist)
analyselist(5059,impdatasetCRref1s54321,varlist)
#\n")
sink()
log_close()

##################################################################### test causal #################################################################
log_locn<- "testref2CausalK0K154321.log"
log_open(log_locn)

sink(paste0("C:/GitHub/mimix/mimixR/log/",log_locn,""),append=TRUE,split=TRUE)

#cat('
# compare with J2R K0=0
impdatasetCausalK0<-mimix("asthma",c("base"),"fev","treat","id","time",1000,1,"Causal",101,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,K0=0,K1=0.9,mle=0 )
#cat('
impdatasetref2CausalK0K154321<-mimix("asthma",c("base"),"fev","treat","id","time",1000,2,"Causal",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,K0=1,K1=1,mle=0 )
#\n')
impdatasetJ2Rref1<-mimix("asthma",c("base"),"fev","treat","id","time",1000,1,"J2R",101,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,K0=0,K1=0.9,mle=0 )

#\n')
#cat('
impdatasetJ2Rref1<-mimix("asthma",c("base"),"fev","treat","id","time",1000,1,"J2R",101,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,K0=0,K1=0.9,mle=0 )
#\n')

impdatasetJ2Rref1Causal54321<-mimix("asthma",c("base"),"fev","treat","id","time",1000,1,"Causal",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,K0=1,K1=0.9,mle=0 )

#cat('
fit<-with(data= as.mids(impdatasetCausalK0), expr = lm(fev.12~treat+base))
summary(pool(fit))

#\n')
#cat(' 
fit<-with(data= as.mids(impdatasetref2CausalK0K154321), expr = lm(fev.12~treat+base))
summary(pool(fit))

analyselist(5051,impdatasetref2CausalK0K154321,varlist)
analyselist(5333,impdatasetref2CausalK0K154321,varlist)
analyselist(5094,impdatasetref2CausalK0K154321,varlist)
analyselist(5115,impdatasetref2CausalK0K154321,varlist)
analyselist(5017,impdatasetref2CausalK0K154321,varlist)
analyselist(5044,impdatasetref2CausalK0K154321,varlist)
analyselist(5059,impdatasetref2CausalK0K154321,varlist)
#\n')
sink()
log_close()



analyselist(5051,impdatasetJ2Rref1CausalK0,varlist)
analyselist(5333,impdatasetJ2Rref1CausalK0,varlist)

analyselist(5094,impdatasetJ2Rref1CausalK0,varlist)

analyselist(5115,impdatasetJ2Rref1CausalK0,varlist)

analyselist(5017,impdatasetJ2Rref1CausalK0,varlist)

analyselist(5044,impdatasetJ2Rref1CausalK0,varlist)

analyselist(5059,impdatasetJ2Rref1CausalK0,varlist)

#\n")
sink()

##################################################################### delta  #############################################################

log_locn<- "testref2delta101.log"
log_open(log_locn)

sink(paste0("C:/GitHub/mimix/mimixR/log/",log_locn,""),append=TRUE,split=TRUE)

#cat(' 
impdatasetJ2Rdel<-mimix("asthma",c("base"),"fev","treat","id","time",10,2,"J2r",101,"jeffreys",1000,NULL,NULL,NULL,c(2,2,1.5,1.5),NULL,1,1,mle=0 )

fit<-with(data= as.mids(impdatasetJ2Rdel), expr = lm(fev.12~treat+base))
summary(pool(fit))
varlist <- c("fev.2","fev.4","fev.8","fev.12");
analyselist(5333,impdatasetJ2Rdel,varlist)
analyselist(5051,impdatasetJ2Rdel,varlist)
analyselist(5115,impdatasetJ2Rdel,varlist)
analyselist(5017,impdatasetJ2Rdel,varlist)
analyselist(5044,impdatasetJ2Rdel,varlist)
analyselist(5059,impdatasetJ2Rdel,varlist)
#\n')

#cat(' 
impdatasetCIRdel<-mimix("asthma",c("base"),"fev","treat","id","time",10,2,"CIR",101,"jeffreys",1000,NULL,NULL,NULL,c(2,2,1.5,1.5),NULL,1,1,mle=0 )

fit<-with(data= as.mids(impdatasetCIRdel), expr = lm(fev.12~treat+base))
summary(pool(fit))
varlist <- c("fev.2","fev.4","fev.8","fev.12");
analyselist(5333,impdatasetCIRdel,varlist)
analyselist(5051,impdatasetCIRdel,varlist)
analyselist(5115,impdatasetCIRdel,varlist)
analyselist(5017,impdatasetCIRdel,varlist)
analyselist(5044,impdatasetCIRdel,varlist)
analyselist(5059,impdatasetCIRdel,varlist)
#\n')

#cat(' 
impdatasetCRdel<-mimix("asthma",c("base"),"fev","treat","id","time",10,2,"CR",101,"jeffreys",1000,NULL,NULL,NULL,c(2,2,1.5,1.5),NULL,1,1,mle=0 )

fit<-with(data= as.mids(impdatasetCRdel), expr = lm(fev.12~treat+base))
summary(pool(fit))
varlist <- c("fev.2","fev.4","fev.8","fev.12");
analyselist(5333,impdatasetCRdel,varlist)
analyselist(5051,impdatasetCRdel,varlist)
analyselist(5115,impdatasetCRdel,varlist)
analyselist(5017,impdatasetCRdel,varlist)
analyselist(5044,impdatasetCRdel,varlist)
analyselist(5059,impdatasetCRdel,varlist)
#\n')

#cat(' 
impdatasetLMCFdel<-mimix("asthma",c("base"),"fev","treat","id","time",10,2,"LMCF",101,"jeffreys",1000,NULL,NULL,NULL,c(2,2,1.5,1.5),NULL,1,1,mle=0 )

fit<-with(data= as.mids(impdatasetLMCFdel), expr = lm(fev.12~treat+base))
summary(pool(fit))
varlist <- c("fev.2","fev.4","fev.8","fev.12");
analyselist(5333,impdatasetLMCFdel,varlist)
analyselist(5051,impdatasetLMCFdel,varlist)
analyselist(5115,impdatasetLMCFdel,varlist)
analyselist(5017,impdatasetLMCFdel,varlist)
analyselist(5044,impdatasetLMCFdel,varlist)
analyselist(5059,impdatasetLMCFdel,varlist)
#\n')

#cat(' 
impdatasetMARdel<-mimix("asthma",c("base"),"fev","treat","id","time",10,2,"MAR",101,"jeffreys",1000,NULL,NULL,NULL,c(2,2,1.5,1.5),NULL,1,1,mle=0 )

fit<-with(data= as.mids(impdatasetMARdel), expr = lm(fev.12~treat+base))
summary(pool(fit))
varlist <- c("fev.2","fev.4","fev.8","fev.12");
analyselist(5333,impdatasetMARdel,varlist)
analyselist(5051,impdatasetMARdel,varlist)
analyselist(5115,impdatasetMARdel,varlist)
analyselist(5017,impdatasetMARdel,varlist)
analyselist(5044,impdatasetMARdel,varlist)
analyselist(5059,impdatasetMARdel,varlist)
#\n')
#cat(' 
impdatasetCausaldel<-mimix("asthma",c("base"),"fev","treat","id","time",10,2,"Causal",101,"jeffreys",1000,NULL,NULL,NULL,c(2,2,1.5,1.5),NULL,1,0.9,mle=0 )

fit<-with(data= as.mids(impdatasetCausaldel), expr = lm(fev.12~treat+base))
summary(pool(fit))
varlist <- c("fev.2","fev.4","fev.8","fev.12");
analyselist(5333,impdatasetCausaldel,varlist)
analyselist(5051,impdatasetCausaldel,varlist)
analyselist(5115,impdatasetCausaldel,varlist)
analyselist(5017,impdatasetCausaldel,varlist)
analyselist(5044,impdatasetCausaldel,varlist)
analyselist(5059,impdatasetCausaldel,varlist)
#\n')

sink()
log_close()

cat(' 
impdatasetCRnodel<-mimix("asthma",c("base"),"fev","treat","id","time",10,2,"CR",101,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,1,1,mle=0 )

fit<-with(data= as.mids(impdatasetCRdel), expr = lm(fev.12~treat+base))
summary(pool(fit))
varlist <- c("fev.2","fev.4","fev.8","fev.12");
analyselist(5333,impdatasetCRnodel,varlist)
analyselist(5051,impdatasetCRnodel,varlist)
analyselist(5115,impdatasetCRnodel,varlist)
analyselist(5017,impdatasetCRnodel,varlist)
analyselist(5044,impdatasetCRnodel,varlist)
analyselist(5059,impdatasetCRnodel,varlist)
\n')
#cat(' 
impdatasetLMCFnodel<-mimix("asthma",c("base"),"fev","treat","id","time",10,2,"LMCF",101,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,1,1,mle=0 )

fit<-with(data= as.mids(impdatasetLMCFnodel), expr = lm(fev.12~treat+base))
summary(pool(fit))
varlist <- c("fev.2","fev.4","fev.8","fev.12");
analyselist(5333,impdatasetLMCFnodel,varlist)
analyselist(5051,impdatasetLMCFnodel,varlist)
analyselist(5115,impdatasetLMCFnodel,varlist)
analyselist(5017,impdatasetLMCFnodel,varlist)
analyselist(5044,impdatasetLMCFnodel,varlist)
analyselist(5059,impdatasetLMCFnodel,varlist)
#\n')
#,1000,NULL,NULL,NULL,c(0.5,0.5,1,1),c(1,1,1,1),1,0.6,mle=0)

#cat(' 
impdatasetLMCFdel<-mimix("asthma",c("base"),"fev","treat","id","time",10,2,"LMCF",101,"jeffreys",1000,NULL,NULL,NULL,c(2,2,1.5,1.5),NULL,1,1,mle=0 )

fit<-with(data= as.mids(impdatasetLMCFdel), expr = lm(fev.12~treat+base))
summary(pool(fit))
varlist <- c("fev.2","fev.4","fev.8","fev.12");
analyselist(5333,impdatasetLMCFdel,varlist)
analyselist(5051,impdatasetLMCFdel,varlist)
analyselist(5115,impdatasetLMCFdel,varlist)
analyselist(5017,impdatasetLMCFdel,varlist)
analyselist(5044,impdatasetLMCFdel,varlist)
analyselist(5059,impdatasetLMCFdel,varlist)
#\n')

log_print("analyselist(5051,impdatasetCIRref2s101,varlist),console=TRUE")
log_print(analyselist(5051,impdatasetCIRref2s101,varlist))
log_print("analyselist(5051,impdatasetCIRref2s54321,varlist),console=TRUE")
log_print(analyselist(5051,impdatasetCIRref2s54321,varlist))
log_print("analyselist(5333,impdatasetCRref2s101,varlist),console=TRUE")
log_print(analyselist(5333,impdatasetCIRref2s101,varlist))
log_print("analyselist(5333,impdatasetCRref2s54321,varlist),console=TRUE")
log_print(analyselist(5333,impdatasetCIRref2s54321,varlist))
log_print("analyselist(5094,impdatasetCIRref2s101,varlist),console=TRUE")
log_print(analyselist(5094,impdatasetCIRref2s101,varlist))
log_print("analyselist(5094,impdatasetCIRref2s54321,varlist),console=TRUE")
log_print(analyselist(5094,impdatasetCIRref2s54321,varlist))
log_print("analyselist(5115,impdatasetCIRref2s101,varlist),console=TRUE")
log_print(analyselist(5115,impdatasetCIRref2s101,varlist))
log_print("analyselist(5115,impdatasetCIRref2s54321,varlist),console=TRUE")
log_print(analyselist(5115,impdatasetCIRref2s54321,varlist))
log_print("analyselist(5017,impdatasetCIRref2s101,varlist),console=TRUE")
log_print(analyselist(5017,impdatasetCIRref2s101,varlist))
log_print("analyselist(5017,impdatasetCIRref2s54321,varlist),console=TRUE")
log_print(analyselist(5017,impdatasetCIRref2s54321,varlist))
log_print("analyselist(5044,impdatasetCIRref2s101,varlist),console=TRUE")
log_print(analyselist(5044,impdatasetCIRref2s101,varlist))
log_print("analyselist(5044,impdatasetCIRref2s54321,varlist),console=TRUE")
log_print(analyselist(5044,impdatasetCIRref2s54321,varlist))
log_print("analyselist(5059,impdatasetCIRref2s101,varlist),console=TRUE")
log_print(analyselist(5059,impdatasetCIRref2s101,varlist))
log_print("analyselist(5059,impdatasetCIRref2s54321,varlist),console=TRUE")
log_print(analyselist(5059,impdatasetCIRref2s54321,varlist))
log_close()


log_print("analyselist(5051,impdatasetCRref2s101,varlist),console=TRUE")
log_print(analyselist(5051,impdatasetCRref2s101,varlist))
log_print("analyselist(5051,impdatasetCRref2s54321,varlist),console=TRUE")
log_print(analyselist(5051,impdatasetCRref2s54321,varlist))
log_print("analyselist(5333,impdatasetCRref2s101,varlist),console=TRUE")
log_print(analyselist(5333,impdatasetCRref2s101,varlist))
log_print("analyselist(5333,impdatasetCRref2s54321,varlist),console=TRUE")
log_print(analyselist(5333,impdatasetCRref2s54321,varlist))
log_print("analyselist(5094,impdatasetCRref2s101,varlist),console=TRUE")
log_print(analyselist(5094,impdatasetCRref2s101,varlist))
log_print("analyselist(5094,impdatasetCRref2s54321,varlist),console=TRUE")
log_print(analyselist(5094,impdatasetCRref2s54321,varlist))
log_print("analyselist(5115,impdatasetCRref2s101,varlist),console=TRUE")
log_print(analyselist(5115,impdatasetCRref2s101,varlist))
log_print("analyselist(5115,impdatasetCRref2s54321,varlist),console=TRUE")
log_print(analyselist(5115,impdatasetCRref2s54321,varlist))
log_print("analyselist(5017,impdatasetCRref2s101,varlist),console=TRUE")
log_print(analyselist(5017,impdatasetCRref2s101,varlist))
log_print("analyselist(5017,impdatasetCRref2s54321,varlist),console=TRUE")
log_print(analyselist(5017,impdatasetCRref2s54321,varlist))
log_print("analyselist(5044,impdatasetCRref2s101,varlist),console=TRUE")
log_print(analyselist(5044,impdatasetCRref2s101,varlist))
log_print("analyselist(5044,impdatasetCRref2s54321,varlist),console=TRUE")
log_print(analyselist(5044,impdatasetCRref2s54321,varlist))
log_print("analyselist(5059,impdatasetCRref2s101,varlist),console=TRUE")
log_print(analyselist(5059,impdatasetCRref2s101,varlist))
log_print("analyselist(5059,impdatasetCRref2s54321,varlist),console=TRUE")
log_print(analyselist(5059,impdatasetCRref2s54321,varlist))
log_close()


impdatasetCIRref2<-mimix("asthma",c("base"),"fev","treat","id","time",10,2,"CIR",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
impdatasetCIRref2<-mimix("asthma",c("base"),"fev","treat","id","time",10,2,"CIR",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )



library(logger)

#initiate log header
log_locn<- "testCRref2seed.log"

library(logr)
If <- 
log_open(log_locn)
sink()

log_txt<-file("C:/GitHub/mimix/test.txt",open="a")
log_info("executing {impdatasetJ2Rref2}")

log_con<- file("C:/GitHub/mimix/mimixR/test.log")
cat("write to file",file=log_con)
log_con<-file("test.log",open="a")
close(log_con)
writeLines(Main_Methods,"C:/GitHub/mimix/test.txt")

log_txt<- file("C:/GitHub/mimix/test.txt")


sink("C:/GitHub/mimix/test.txt",append=T)
cat("\nHers an appended lineimpdatasetJ2Rref2<-mimix(\"asthma",c("base"),"fev","treat","id","time\",10,2,\"J2R\",54321,\"jeffreys\",1000,NULL,NULL,NULL,NULL,NULL )
")
sink()

cat('impdatasetJ2RnoD<-mimix("asthma",c("base"),"fev","treat","id","time",1000,2,"J2R",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0)',file=log_txt)
impdatasetJ2RnoD<-mimix("asthma",c("base"),"fev","treat","id","time",1000,2,"J2R",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )

library(mice)

cat('fit<-with(data= as.mids(impdatasetJ2RnoD), expr = lm(fev.12~treat+base))',file=log_txt)
fit<-with(data= as.mids(impdatasetJ2RnoD), expr = lm(fev.12~treat+base))
summary(pool(fit))
varlist <- c("fev.2","fev.4","fev.8","fev.12")

Main_Methods<- c("J2R","CIR","CR","MAR")
for (val in 1:length(Main_Methods)) {
  analyselist(5051,assign(paste0("impdataset",Main_Methods[val],"noD")),varlist)
}

impdatasetJ2Rref2

MMethod <- c("J2R","CIR","CR","MAR")
test_SNOs<-c("5051","5333","5017","5059","5115")

test_fun<-function(X) {
  for (val in 1:length(Main_Methods)) {
    cat(paste0("\nimpdataset",Main_Methods[[val]],"ref2 ",X))
    cat(paste0("\n treatarm = ",subset(impdatasetJ2Rref2$treat,(impdatasetJ2Rref2$.id==X) & (impdatasetJ2Rref2$.imp==0)) ))
    cat("\n")
    print(analyselist(X,get(paste0("impdataset",Main_Methods[[val]],"ref2")),varlist))
  }
}
sapply(test_SNOs,test_fun)
 
for (val in 1:length(Main_Methods)) {
   cat(paste0("\nimpdataset",Main_Methods[[val]],"noD"))
   print(analyselist(5059,get(paste0("impdataset",Main_Methods[[val]],"ref")),varlist))
}

If <- log_open(tmp)
log_txt<- ("C:/GitHub/mimix/test.txt")
If <- log_open(log_txt)
log_txt<-file("C:/GitHub/mimix/test.txt",open="a")

log_print(analyselist(5115,impdatasetJ2RnoD,varlist))
log_print(analyselist(5333,impdatasetJ2RnoD,varlist))
log_print(analyselist(5017,impdatasetJ2RnoD,varlist))
log_print(analyselist(5059,impdatasetJ2RnoD,varlist))


hmc<-(analyselist(5115,impdatasetJ2RnoD,varlist))
log_close
writeLines(readLines(log_txt))


impdatasetCIRnoD<-mimix("asthma",c("base"),"fev","treat","id","time",10,1,"CIR",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL )

fit<-with(data= as.mids(impdatasetCIRnoD), expr = lm(fev.12~treat+base))
summary(pool(fit))
varlist <- c("fev.2","fev.4","fev.8","fev.12");
analyselist(5051,impdatasetCIRnoD,varlist)
analyselist(5115,impdatasetCIRnoD,varlist)
analyselist(5333,impdatasetCIRnoD,varlist)
analyselist(5017,impdatasetCIRnoD,varlist)
analyselist(5059,impdatasetCIRnoD,varlist)


impdatasetCRnoD<-mimix("asthma",c("base"),"fev","treat","id","time",10,1,"CR",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL )

fit<-with(data= as.mids(impdatasetCRnoD), expr = lm(fev.12~treat+base))
summary(pool(fit))
varlist <- c("fev.2","fev.4","fev.8","fev.12");
analyselist(5051,impdatasetCRnoD,varlist)
analyselist(5115,impdatasetCRnoD,varlist)
analyselist(5333,impdatasetCRnoD,varlist)
analyselist(5017,impdatasetCRnoD,varlist)
analyselist(5059,impdatasetCRnoD,varlist)

# try create log header in file first
log_con<-file("testMAR.log",open="a")
sink("C:/GitHub/mimix/testMAR.txt",append=T,split=TRUE)

system.time(
impdatasetJ2RnoD<-mimix("asthma",c("base"),"fev","treat","id","time",1000,1,"J2R",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0  )
)

fit<-with(data= as.mids(impdatasetMARnoD), expr = lm(fev.12~treat+base))
summary(pool(fit))
varlist <- c("fev.2","fev.4","fev.8","fev.12");
analyselist(5051,impdatasetMARnoD,varlist)
analyselist(5115,impdatasetMARnoD,varlist)
analyselist(5333,impdatasetMARnoD,varlist)
analyselist(5017,impdatasetMARnoD,varlist)
analyselist(5059,impdatasetMARnoD,varlist)
sink()


closeAllConnections()

#try other data set
impantiJ2Rtest_NoDelta <- mimix("antidepressant",c("basval"),"change","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,2,"J2R",101,c("jeffreys"),1000,NULL,NULL,NULL,,,K0,K1,mle=0)

#also try the acupunctre data?
load("C:/GitHub/mimix/mimixOld/data/acupuncture.RData")
impdatasetJ2RAcup<-mimix("acupuncture",c("head_base"),"head","treat","id","time",1000,2,"J2R",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0  )

#NOTE ref 1 corresponds ref =2 in stata
impdatasetJ2RCcup<-mimix("acupuncture",,"head","treat","id","time",1000,1,"J2R",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0  )


library(mice)
fit<-with(data= as.mids(impdatasetJ2RAcup), expr = lm(head.12~treat+head_base))
summary(pool(fit))
varlist <- c("head.3","head.12");
analyselist(104,impdatasetJ2RAcup,varlist)
analyselist(333,impdatasetJ2RAcup,varlist)
analyselist(386,impdatasetJ2RAcup,varlist)
analyselist(787,impdatasetJ2RAcup,varlist)
analyselist(435,impdatasetJ2RAcup,varlist)
analyselist(697,impdatasetJ2RAcup,varlist)
analyselist(100,impdatasetJ2RAcup,varlist)
analyselist(101,impdatasetJ2RAcup,varlist)

impdatasetJ2RnoD<-mimix("asthma",c("base"),"fev","treat","id","time",100,2,"J2R",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )
impdatasetJ2RnoC<-mimix("asthma",,"fev","treat","id","time",100,2,"J2R",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,NULL,NULL,mle=0 )



