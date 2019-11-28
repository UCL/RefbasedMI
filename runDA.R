# having preprocessed long data into wide data AND found mimix_group
# next  
# perform DA using norm2 
# by looping thru ntreat
#

# pick up from test<-preprodata("fev","treat","id","time","base")
library(norm2)
library(sparseinv) 
#M<-1

for (val in t(ntreat)) {
  print(paste0("prenormdat",val))
  assign(paste0("prenormdat",val),subset(finaldatS,treat==val))
}  

# need to generate this seq
# fev1,fev2,fev3,fev4,base
# from response &  covar loop
# use the ntime var  

mi_impute <-function(idvar,timevar,depvar,covar) {
  #preprodata(depvar,treatvar,idvar,timevar,covar)
   tst<-pivot_wider(mxdata,id_cols=c(idvar),names_from=timevar,names_prefix=depvar,values_from=depvar)
   # add ont covariates and add on to resp list
   respvars<-c(names(tst[,-1]),covar)
   return(respvars)
}

# so tsts2 is th equiv to mi_impute so try to put in norm2
tst2 <- mi_impute("id","time","fev","base")
# put commas in
tst3<-paste(tst2,collapse = ",")


#try putting a matrix instead of a formula seems ok
#seems to works

##CHECK WITH STATA , MULTIPLE DOESNT KINCK IN YET!!1? 

#*CREATE AN EMPTY MATRIX FOR COLLECTING IMPUTED DATA 	
# just an empty row but take dimensions and col types fron mata_obs
#mata_all_new<-array(data=1,dim=c(1,2+ncol(mata_Obs)))
GI<-c(0)
II<-c(0)
SNO <-mata_Obs[1,1]
dropid<-c("id")
mata_ObsX<-mata_Obs[,!(names(mata_Obs) %in% dropid)] 
mata_all_new <- cbind(GI,II,mata_ObsX,SNO)
mata_all_new[ mata_all_new>=0] <-NA 

#try
#for (val in t(ntreat)) {
#  for(m in 1:M) {
#    #prnormobj <- subset(prenormdat2, select=c(tst2))
#    kmvar=get(paste0("prenormdat",val))
#    testprnormobj<-assign(paste0("prnormobj",val), subset(kmvar, select=c(tst2)))
#    emResultT<-emNorm(testprnormobj,prior = "ridge",prior.df=0.5)
#    mcmcResultT<- mcmcNorm(emResultT,iter=1000,multicycle = 100)
#    assign(paste0("mcmcResultT",val,m),mcmcResultT)
#    assign(paste0("parambeta",val,m),mcmcResultT$param$beta)
#    assign(paste0("paramsigma",val,m),mcmcResultT$param$sigma)
#    print(paste0("parambeta",val,m))
#}                                 
#}
library(tictoc)
#M<-1000
tic("total")
M<-1000
for (val in t(ntreat)) {
    #prnormobj <- subset(prenormdat2, select=c(tst2))
    kmvar=get(paste0("prenormdat",val))
    testprnormobj<-assign(paste0("prnormobj",val), subset(kmvar, select=c(tst2)))
  for(m in 1:M) {  
    #emResultT<-emNorm(testprnormobj,prior = "ridge",prior.df=0.5)
    emResultT<-emNorm(testprnormobj,prior = "jeffreys")
    mcmcResultT<- mcmcNorm(emResultT,iter=1000,multicycle = NULL,prior = "jeffreys")
    assign(paste0("mcmcResultT",val,m),mcmcResultT)
    assign(paste0("parambeta",val,m),mcmcResultT$param$beta)
    assign(paste0("paramsigma",val,m),mcmcResultT$param$sigma)
    print(paste0("parambeta",val,m))
  }                                 
}
toc()


## these are now refrenced so can use in below imputation loop


  
# now go to the START main analysis bit
# need to use look up file to process pattern (by treatment)
# mimix_group<- test[[5]]
# need define missing pattern vector within look-up table
 
#*CREATE AN EMPTY MATRIX FOR COLLECTING IMPUTED DATA 	
# just an empty row but take dimensions and col types fron mata_obs
#mata_all_new<-array(data=1,dim=c(1,2+ncol(mata_Obs)))
GI<-c(0)
II<-c(0)
SNO <-mata_Obs[1,1]
dropid<-c("id")
mata_ObsX<-mata_Obs[,!(names(mata_Obs) %in% dropid)] 
mata_all_new <- cbind(GI,II,mata_ObsX,SNO)
mata_all_new[ mata_all_new>=0] <-NA 

#testing
SigmaRefer11_new = matrix(,nrow=5,ncol=5)
S12_pat11_new =matrix(,nrow=2,ncol=3)

# loop over pattern lookup
#   for ( i in  1:nrow(mg)) 
#testingloop over last pattern  
#for  (i in nrow(mg):nrow(mg))

meth <- 'MAR' 
for (i in 1:nrow(mg))
  {
  #define mata_miss as vector of 1's denoting missing using col names ending i ".1"
  mata_miss <- mg[i,grep("*..1",colnames(mg)),drop=F]
  #mata_miss <- mimix_group[i,c(2,3,4,5)]     #define mata_miss
  #assumes covariate non missing
  mata_miss$covar.1 <-0                       #assuming cov col is non missing
  mata_nonmiss <- (ifelse(mata_miss==0,1,0))  #define mata-nonmiss from miss
  
  # need transform nonmiss,miss to c lists - ie. index the  
  c_mata_miss<-which(mata_miss==1)
  c_mata_nonmiss<-which(mata_nonmiss==1)
  # eg no missing is c(1,2,3,4,5)
  
  # if non missing then have to amend above !!!  
  #S11 will be all of matrix, S12,S22 will be Null
  #if c_mata_miss
  # 1's signify missing values, 0 indicates not missing.
  # use if statement to create suitable mat_obs submat     
  
  
  # count of pattern by treatment
  cnt<- mg$X1[i]

  # treatment grp
  trtgp<- mg$treat[i]
  cat("trtgp = ", trtgp)
  # reference grp
  # refer col is just a constant , might be better than this method 
  refer <- 2 # change this to  input argument!
  # see what chanfe refer makes!!
  #refer <-3 doesnt improve the ests
  # try to mcmcrEsult as 2 matrices  

  # multiple  simulations start here #########  
for ( m in  1:M)  { 
  #*FOR INDIVIDUALS WITH NO MISSING DATA COPY COMPLETE DATA INTO THE NEW DATA MATRIX mata_all_new `m' TIMES
			#if `pat' == 0{ 
  if(length(c_mata_miss)==0 ) {
          st<-mg[i,"X1cum"]-mg[i,"X1"]+1
          en <-mg[i,"X1cum"]
          SNO<-mata_Obs[c(st:en),1]
          mata_new <- mata_Obs[c(st:en),2:ncol(mata_Obs)]
          GI <- array(data=mg[i,"treat"],dim=c(mg[i,"X1"],1))
         #II  no imputations 
          II <- array(data=m,dim=c(mg[i,"X1"],1))
          mata_new=cbind(GI,II,mata_new,SNO)
          #names(mata_all_new)<-names(mata_new)
          mata_all_new=rbind(mata_all_new,mata_new)
  } else {
    #*FOR INDIVIDUALS WITH  MISSING DATA  `m' TIMES  
#} 
#for mimix_group

#assuming for now only 2 treatment groups ! 

# want the treatment grp means and where missing overwrite with refer grp (when differs!)     
#mata_means_trt <-mat_Betas[trt_gp,] 
#mata_means_ref <-mat_Betas[refer,]
print(paste("trtgp= ",trtgp))
print(paste("refer= ",refer))    
  

mata_means_trt <- get(paste0("parambeta",trtgp,m)) 
mata_means_ref <- get(paste0("parambeta",refer,m))

Sigmatrt <- get(paste0("paramsigma",trtgp,m))

if (meth== 'MAR')  {
  
   mata_Means <- mata_means_trt 
   mata_means<-mata_Means[rep(seq(nrow(mata_Means)),each=mg$X1[i]),]
  
   
   
   S11 <-Sigmatrt[c_mata_nonmiss,c_mata_nonmiss]
   S12 <-Sigmatrt[c_mata_nonmiss,c_mata_miss]
   S22 <-Sigmatrt[c_mata_miss,c_mata_miss]
   

}
else if (meth == 'J2R') { 
 

#one way is to element multiply (because 1,0) then add 
mata_means_t <- unlist(mata_means_trt)*mata_nonmiss
mata_means_r <- unlist(mata_means_ref)*mata_miss
# so when all missing  1,1,1, ... then all contributing comes from reference means      
mata_means <- mata_means_r+mata_means_t
# and preserve names   
colnames(mata_means) <- colnames(mata_means_trt)
#replicate to number of rows defined by X1 
mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]
  
############# SIGMA is from paramsigma  the reference group ################

#SigmaRef <- Reduce(rbind,mat_Sigma[1]) 
# this could be the error was list[1] rather than m! , mkaes a bit of difference, still discrepancy with stata 


# do we ever  use SigmaTrt !!?? in j2r? 
SigmaRefer <- get(paste0("paramsigma",refer,m))
# when reading in Stata sigmas
# needs to to ths as tibble will fail in cholesky
#SigmaRefer <- as.matrix(get(paste0("paramsigmaStata",refer,m)))

print(paste0("paramsigma",refer,m))
#SigmaRefer <- Reduce(rbind,res1sigma.list[m])

S11 <-SigmaRefer[c_mata_nonmiss,c_mata_nonmiss]
S12 <-SigmaRefer[c_mata_nonmiss,c_mata_miss]
S22 <-SigmaRefer[c_mata_miss,c_mata_miss]

}
else if (meth=='CR') {
  mata_means_ref <- get(paste0("parambeta",refer,m))
  mata_means_r <- unlist(mata_means_ref)*1
  #mata_means_r <- unlist(mata_means_ref)*mata_miss
  #mata_means <- mata_means_r
  # and preserve names   
  #colnames(mata_means) <- colnames(mata_means_trt)
  #replicate to number of rows defined by X1 
  mata_means <- mata_means_r
  colnames(mata_means) <- colnames(mata_means_ref)
  mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]
  
  SigmaRefer <- get(paste0("paramsigma",refer,m))
  S11 <-SigmaRefer[c_mata_nonmiss,c_mata_nonmiss]
  S12 <-SigmaRefer[c_mata_nonmiss,c_mata_miss]
  S22 <-SigmaRefer[c_mata_miss,c_mata_miss]
}
else if (meth=='CIR') {
  # pre-deviating use mean of trt gp up to last obs time bfore deviating, post-deviating use mean from ref grp 
  
  #mata_means_trt <- get(paste0("parambeta",trtgp,m)) 
  #mata_means_ref <- get(paste0("parambeta",refer,m))
  
  # put equiv to mimix 
  mata_Means <-  mata_means_trt
  MeansC <- mata_means_ref
  
  #might be better to copy mimix algol
  
  mata_means<-CIR_loop(c_mata_miss,mata_Means,MeansC)
  #returns mata_means as single row
  # then duplicate over patt rows
  #replicate to number of rows defined by X1 
  mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]
  
  SigmaRefer <- get(paste0("paramsigma",refer,m))
  # when reading in Stata sigmas
  # needs to to ths as tibble will fail in cholesky
  #SigmaRefer <- as.matrix(get(paste0("paramsigmaStata",refer,m)))
  
  print(paste0("paramsigma",refer,m))
  #SigmaRefer <- Reduce(rbind,res1sigma.list[m])
  
  S11 <-SigmaRefer[c_mata_nonmiss,c_mata_nonmiss]
  S12 <-SigmaRefer[c_mata_nonmiss,c_mata_miss]
  S22 <-SigmaRefer[c_mata_miss,c_mata_miss]
 }  
  else if (meth=='LMCF') { 
    mata_Means <-  mata_means_trt
    # no ref MeansC <- mata_means_ref
    mata_means<-LMCF_loop(c_mata_miss,mata_Means)
    Sigmatrt <- get(paste0("paramsigma",trtgrp,m))
    # when reading in Stata sigmas
    S11 <-Sigmatrt[c_mata_nonmiss,c_mata_nonmiss]
    S12 <-Sigmatrt[c_mata_nonmiss,c_mata_miss]
    S22 <-Sigmatrt[c_mata_miss,c_mata_miss]
 }


  


###################### MNAR IMPUTATION ################
# need insert routine when ALL missing values
# else
###################### MNAR IMPUTATION ################



m1 <- mata_means[c_mata_nonmiss]
m2 <- mata_means[c_mata_miss]


#then mata_obs is the sequential selection of rows according to the mimix_group variable X1 values
# stil goimg thru the rows in mimix_group want to read off the X1 values to create subsets   
#for (l in 1:mimix_group[i,"X1"]) {


# mata_obs <- mimix_d2cols[i:j,]

#need a counter to accumulate j  - easiest way is to ceate a cumulstive col in mimix_group
#sol is to use cumX1 minus X1 ,thats all there is to it
#for ( i in  1:nrow(mimix_group))  {
#print(i)
j <- mg[i,"X1"]
print(paste0(" count in patt = ", j))
k <- mg[i,"X1cum"]
startrow <-(k-j+1)  
stoprow  <-(k)
print("startrow stoprow = ")
print(startrow)
print(stoprow)

#raw1 <- mata_obs[, c_mata_nonmiss]
preraw <- mata_Obs[c(startrow:stoprow),2:ncol(mata_Obs)]
raw1 <- preraw[,c_mata_nonmiss]
print("mata_raw = ")
#print(raw1)

#}
# to here seems ok 8/11/19 
#J2r uses Sigma derived from refernce group

#t=cholsolve(S11,S12)
#cholsolve(A, B) solves AX=B and returns X for symmetric (Hermitian), positive-definite A.  cholsolve() returns a matrix of missing values if A is not positive definite or if A is singular.
# so want solve S11X=S12


  t_mimix =cholsolve(Q=S11,y=S12)   
  conds <-  S22-t(S12)%*%t_mimix
  #  
  #meanval = as.matrix(m2) + as.matrix(raw1 - m1)%*%as.matrix(t_mimix)
  # below for CR but may not work with J2R
  meanval = (m2) + as.matrix(raw1 - m1)%*%as.matrix(t_mimix)
  U <- chol(conds)
  # mg[i,X1] is equiv to Stata counter, miss_count is no. of missing, so 
  miss_count=rowSums(mata_miss)
  Z<-qnorm(matrix(runif( mg[i,"X1"]* miss_count,0,1),mg[i,"X1"],miss_count))
  # check same input parameters for inverse norm gen as in stata 
  print(paste0('mg[i,X1] =',mg[i,"X1"],' miss_count= ',miss_count)) 
 
  mata_y1 = meanval+Z%*%t(U)
  
  #define new matrix from observed,  id column  (the last)
  mata_new <- preraw
 
  # assigning the columns where the missing values  
  if(length(c_mata_miss)==0 ) { mata_new <- mata_Obs[,c(2:nct+1)]
  }else {
  mata_new[,c_mata_miss] <- mata_y1
  }
  # if(length(c_mata_miss)!=0 )
  
  #save U,Z for imp(m) and patt(i)
  # this worked for Z
  
  #assign(paste0("Z11",i,"_imp",m),subset(Z))
  #assign(paste0("U11",i,"_imp",m),subset(U))
  assign(paste0("S12",i,"_imp",m),S12)
  assign(paste0("S22",i,"_imp",m),S22)
  assign(paste0("mata_y1_11",i,"_imp",m),mata_y1)
  assign(paste0("conds11",i,"_imp",m),conds) 
  # no idea why below cuases error!
  assign(paste0("t_mimix11",i,"_imp",m),t_mimix)
  
  #assuming this  from Stata if "`interim'"==""{
  # SNO just id col
  SNO <- mata_Obs[c(startrow:stoprow),1]
  #SNO <- mata_ObsX[,ncol(mata_Obs)]
  # GI treatment grp column 1 (here),II imputation number col, mata_new matrix then SNO is id col.
  GI <- array(data=mg[i,1],dim=c(mg[i,"X1"],1))
  #II  no imputations 
  II <- array(data=m,dim=c(mg[i,"X1"],1))
  
  #doesnt need SNO, as id already in
  #mata_new<-cbind(GI,II,mata_new,SNO)
  mata_new<-cbind(GI,II,mata_new,SNO)
  mata_all_new<-rbind(mata_all_new,mata_new)
  #stata equilv \ is row bind NOt cbind!!?
  }
    
}
} #for M StOP HERE!!

mata_all_new_rmnaMAR <- na.omit(mata_all_new)
mata_all_new_rmnaCIR <- na.omit(mata_all_new)
mata_all_new_rmnaStataSig <- na.omit(mata_all_new) 
# to save the results file
mata_all_new_rmna1k <- na.omit(mata_all_new) 
save (mata_all_new_rmna,file = "kmmata_all_new_rmna181119noiter.RData")
save (mata_all_new_rmnaStataSig,file = "kmmata_all_new_rmna191119Stata.RData")
  
#from ttest.r
#5333 is 1st 3 fields missing, treat=2  
mata_all_new_rmna5333 <- filter(mata_all_new_rmna1k,SNO == "5333")
mata_all_new_rmnaJeff5456 <- filter(mata_all_new_rmna1kJeff,SNO == "5456")
mata_all_new_rmna5333 <- filter(mata_all_new_rmna,SNO == "5333") 
mata_all_new_rmnaref35456 <- filter(mata_all_new_rmnarefer3,SNO == "5456") 

mata_all_new_rmnaStataSig5456 <- filter(mata_all_new_rmnaStataSig,SNO == "5456")

mata_all_new_rmnaCIR5456 <- filter(mata_all_new_rmnaCIR,SNO == "5456") 
mata_all_new_rmnaMAR5456 <- filter(mata_all_new_rmnaMAR,SNO == "5456")
round(stat.desc(mata_all_new_rmnaMAR5456[,c("fev2","fev4","fev8","fev12")]),3)[c(1,9,13,4,8,5),]
print("mata_means = ")
print(mata_means) 

library(pastecs)
for ( val in c("fev2","fev4","fev8","fev12") ) { 
  print(round(stat.desc(mata_all_new_rmnaJeff5456[,val]),3 )[c(1,9,13,4,8,5)])
}
                                       
#prnormobj <- subset(prenormdat2, select=c(base,fev2,fev4,fev8,fev12))
#this not work, but take out the "'s does work, have to use cbind
emResultT<-emNorm(prnormobj,prior = "ridge",prior.df=0.5)
#emResultT<-emNorm(cbind(fev2,fev4,fev8,fev12,base) ~1,data=prenormdat2,prior = "ridge",prior.df=0.5)
#uniqdat<-unique(mxdata[c(idvar,covar)])
#respcov<-cbind(respvars,uniqdat)
# declare matrix/list to hold m result outputs
res.list <-as.list(NULL)

res1beta.list <- as.list(NULL)
res1sigma.list <- as.list(NULL)
res2beta.list <- as.list(NULL)
res2sigma.list <- as.list(NULL)
# set M no. imputations in terms of data sets, 1000 took all morning!
# for now just hard code 2 treats, later generalise
M <- 1000
for(m in 1:M) {
  #emNorm(cbind(fev1,fev2,fev3,fev4,base) ~1,data=preimputemvn1)
  #summary(emNorm(cbind(fev1,fev2,fev3,fev4,base) ~1,data=preimputemvn1))
  #emResult1<-emNorm(cbind(fev1,fev2,fev3,fev4,base) ~1,data=prenormdat2,prior = "ridge",prior.df=0.5)
  emResultT<-emNorm(prnormobj,prior = "ridge",prior.df=0.5)
  #mcmcResult <- mcmcNorm(cbind(fev1,fev2,fev3,fev4,base) ~1,data=preimputemvn1,starting.values = emResult$param)
  mcmcResult1<- mcmcNorm(emResult1,iter=5000,multicycle = 100)
  res1beta.list[[m]] <- mcmcResult1$param$beta
  res1sigma.list[[m]] <-mcmcResult1$param$sigma
  
  # necessary to use ridge because finite fifference problem
  # see  user manual norm2User Guide 7,3 
  emResult2<-emNorm(cbind(fev1,fev2,fev3,fev4,base) ~1,data=preimputemvn2,prior = "ridge",prior.df=0.5)
  
  #mcmcResult <- mcmcNorm(cbind(fev1,fev2,fev3,fev4,base) ~1,data=preimputemvn1,starting.values = emResult$param)
  mcmcResult2<- mcmcNorm(emResult2,iter=5000,multicycle = 100)
  res2beta.list[[m]] <- mcmcResult2$param$beta
  res2sigma.list[[m]] <-mcmcResult2$param$sigma
}


