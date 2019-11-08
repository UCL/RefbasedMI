# having preprocessed long data into wide data AND found mimix_group
# next  
# perform DA using norm2 
# by looping thru ntreat
#

library(norm2)

ntreat<-test[[3]]
finaldatS<-test[[2]]

for (val in t(ntreat)) {
  print(paste0("prenormdat",val))
  assign(paste0("prenormdat",val),subset(finaldatS,treat==val))
}  

# need to generate this seq
# fev1,fev2,fev3,fev4,base
# from response &  covar loop
# use the ntime var  

mi_impute <-function(idvar,timevar,depvar,covar) {
   tst<-pivot_wider(mxdata,id_cols=c(idvar),names_from=timevar,names_prefix=depvar,values_from=depvar)
   # add ont covariates and add on to resp list
   respvars<-c(names(tst[,-1]),covar)
   return(respvars)
}

# so tsts2 is th equiv to mi_impute so try to put in norm2
tst2 <- mi_impute("id","time","fev","base")
tst3<-paste(tst2,collapse = ",")


#try formula
fmla <- as.formula(paste(tst2,collapse=","))

#try putting a matrix instead of a formula?
#seems to work
M<-2
for(m in 1:M) {
for (val in t(ntreat)) {
 #prnormobj <- subset(prenormdat2, select=c(tst2))
  kmvar=get(paste0("prenormdat",val))
  testprnormobj<-assign(paste0("prnormobj",val), subset(kmvar, select=c(tst2)))
 emResultT<-emNorm(testprnormobj,prior = "ridge",prior.df=0.5)
 mcmcResultT<- mcmcNorm(emResultT,iter=5000,multicycle = 100)
 assign(paste0("mcmcResultT",val,m),mcmcResultT)
 assign(paste0("parambeta",val,m),mcmcResultT$param$beta)
 assign(paste0("paramsigma",val,m),mcmcResultT$param$sigma)
 print(paste0("parambeta",val,m))
 }                                 
}
  
# now go to the START main analysis bit
# need to use look up file to process pattern (by treatment)
# mimix_group<- test[[5]]
# need define missing pattern vector within look-up table
mg<-test[[5]]
for ( i in  1:nrow(mg)) 
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
  # try to mcmcrEsult as 2 matrices  
  
#} 
#for mimix_group

#assuming for now only 2 treatment groups ! 

# want the treatment grp means and where missing overwrite with refer grp (when differs!)     
#mata_means_trt <-mat_Betas[trt_gp,] 
#mata_means_ref <-mat_Betas[refer,]

mata_means_trt <- get(paste0("parambeta",trtgp)) 
mata_means_ref <- get(paste0("parambeta",refer))
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

SigmaRefer <- get(paste0("paramsigma",refer))
print(paste0("paramsigma",refer))
#SigmaRefer <- Reduce(rbind,res1sigma.list[m])

S11 <-SigmaRefer[c_mata_nonmiss,c_mata_nonmiss]
S12 <-SigmaRefer[c_mata_nonmiss,c_mata_miss]
S22 <-SigmaRefer[c_mata_miss,c_mata_miss]


#  S11 <- mcmcResult1$param$sigma[c_mata_nonmiss,c_mata_nonmiss]
#  S12 <- mcmcResult1$param$sigma[c_mata_nonmiss,c_mata_miss]
#  S22 <- mcmcResult1$param$sigma[c_mata_miss,c_mata_miss]

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
print(j)
k <- mg[i,"X1cum"]
startrow <-(k-j+1)  
stoprow  <-(k)
print("startrow stoprow = ")
print(startrow)
print(stoprow)

mata_raw <- mata_Obs[c(startrow:stoprow),]
print("mata_raw = ")
print(mata_raw)

#}
# to here seems ok 8/11/19 
#J2r uses Sigma derived from refernce group
library(sparseinv) 
#t=cholsolve(S11,S12)
#cholsolve(A, B) solves AX=B and returns X for symmetric (Hermitian), positive-definite A.  cholsolve() returns a matrix of missing values if A is not positive definite or if A is singular.
# so want solve S11X=S12


  t_mimix =cholsolve(Q=S11,y=S12)   
  conds <-  S22-t(S12)%*%t_mimix
  # perhaps below should be checked? 
  meanval = as.matrix(m2) + as.matrix(raw1 - m1)%*%as.matrix(t_mimix)
  U <- chol(conds)
  # mg[i,X1] is equiv to Stata counter, miss_count is no. of missing, so 
  miss_count=rowSums(mata_miss)
  Z<-qnorm(matrix(runif( mg[i,"X1"]* miss_count,0,1),mg[i,"X1"],miss_count))
  
  print(Z) 
 
  mata_y1 = meanval+Z%*%t(U)
  
  #define new matrix from observed,  id column  (the last)
  mata_new <- mata_obs
 
  # assigning the columns where the missing values  
  mata_new[,c_mata_miss] <- mata_y1
  # if(length(c_mata_miss)!=0 )
  #if(length(c_mata_miss)==0 ) { mata_new <- mata_obs[,c(1:nct)]}
  
  #assuming this ture from Stata if "`interim'"==""{
  # SNO just id col
  SNO <- mata_obs[,ncol(mata_obs)]
  # GI treatment grp column 1 (here),II imputation number col, mata_new matrix then SNO is id col.
  GI <- array(data=mg[i,1],dim=c(mg[i,"X1"],1))
  #II  no imputations 
  II <- array(data=imp,dim=c(mg[i,"X1"],1))
  
  #doesnt need SNO, as id already in
  #mata_new<-cbind(GI,II,mata_new,SNO)
  mata_new<-cbind(GI,II,mata_new)
  mata_all_new<-rbind(mata_all_new,mata_new)
  #stata equilv \ is row bind NOt cbind!!?
}



print("mata_means = ")
print(mata_means) 


                                       
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


