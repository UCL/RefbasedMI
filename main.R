install.packages(c("tidyr","dplyr","plyr","haven","pastecs","norm2","sparseinv"))

#pivot_wider
library(tidyr)

#for select
library(dplyr)

library(plyr)
#read data
# need haven
library(haven)
library(pastecs)
library(norm2)
library(sparseinv) 

#be nice to plot the patt?

# need also input refer

refer <-2
getwd()
source("N:/Documents/GitHub/mimix/mimixR/functions.R")

kmargs <- list("fev","treat","id","time","base",100,2,"J2R")
mata_all_newlist=do.call('Runmimix', kmargs)

mxdata<- readdata("asthma.csv")
Runmimix("fev","treat","id","time","base","J2R",10000,2)
Runmimix(kmargs)

#ry without seed
Runmimix<- function(depvar,treatvar,idvar,timevar,covar,M=1,refer=1,meth=NULL) {
# 
rm(list = ls())
source("functions.R")
#mxdata <-read.csv(paste0("./",data))
#mxdata <-read.csv("./asthma.csv")
# empyty R environmet?
#set.seed(101)
 

set.seed(301)
#set.seed(seedval)
# call preprocesssing with variable names from data set
#testlist<-preprodata("fev","treat","id","time","base",10000,2,"J2R")
#testlist<- preprodata(depvar,treatvar,idvar,timevar,covar,M,refer,meth)

#try repack arguments for call to prepro
#testlist<- preprodata(unlist(kmargs))

 testlist = do.call( preprodata,kmargs)


# returns list from preprodata function
ntreat<-unlist(testlist[[4]])
#stopifnot()
finaldatS<-testlist[[2]]
mg<-testlist[[5]]
# vital to get the mata_obs correctly sorted! so corresponds with mimix_group lookup 
# to be consistent with Stata move the base col after the fevs!
mata_Obs <- testlist[[2]]
M <- testlist[[11]]
refer <- testlist[[12]]
stopifnot(refer %in% ntreat)
meth <- testlist[[13]]

# so tsts2 is th equiv to mi_impute so try to put in norm2
tst2 <- mi_impute("id","time","fev","base")
# put commas in
tst3<-paste(tst2,collapse = ",")

#create input data sets for each tment from which to model 
for (val in t(ntreat)) {
  print(paste0("prenormdat",val))
  assign(paste0("prenormdat",val),subset(finaldatS,treat==val))
}  



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
mata_all_newlist <- vector('list',M*nrow(mg))
#Warning message:
#In Ops.factor(left, right) : ‘>=’ not meaningful for factors



# run the mcmc simulations over the treatment grps 
# create a matrix for param files, a beta  and sigma matrices
#paramMatrix<-matrix(1:(nrow(ntreat)*M),nrow=M,dimnames=list(c(1:M),c(1:nrow(ntreat))))

#paramMatrix<-matrix(,nrow=(nrow(ntreat)*M),ncol=2)

#create  emptylist for each treat and multiple m's
paramBiglist <- vector('list',length(ntreat)*M)
#create 

#paramMatrixT<-matrix(,nrow=2,ncol=M)
start_time <- proc.time()
#instead iof start _time try system_time
#paramBetaMatrixT <- matrix(1,nrow=1,ncol=5)
#paramSigmaMatrixT <- matrix(,nrow=5,ncol=5)
iter<-0
system.time(

  #try ser up as many Result data files as treatments instead of one big file?.
  
  
for (val in t(ntreat)) {
  #prnormobj <- subset(prenormdat2, select=c(tst2))
  #prnormobj <- assign(paste0("prenormdat",val),subset(finaldatS,treat==val))
  kmvar=get(paste0("prenormdat",val))
  prnormobj<-assign(paste0("prnormobj",val), subset(kmvar, select=c(tst2)))
  #create  emptylist for each treat 
 # assign(paste0("paramTRTlist",val),vector('list',M))
  print(paste0("Looping for treatment = ",val," performing mcmcNorm for m = 1 to ",M))
  for(m in 1:M) {  
    #emResultT<-emNorm(testprnormobj,prior = "ridge",prior.df=0.5)
    emResultT<-emNorm(prnormobj,prior = "jeffreys")
    mcmcResultT<- mcmcNorm(emResultT,iter=1000,multicycle = NULL,prior = "jeffreys")
   #try saving parm files  to a matrix instead of indiv parm files
    
    # assign doesnty make much differenc in tghis loop , perhaps efficiency comes from not calling them later?
    
    # keep fo rnow to test againt biglist 
    # teszted ok for j2r
    #assign(paste0("param",val,m),mcmcResultT$param)
     
     iter<-iter+1    
     paramBiglist[[iter]] <- mcmcResultT$param
     
   
     
     #this bit trying matrix
     #if (m==1) { 
     #paramBetaMatrixT = as.matrix(mcmcResultT$param)[,1][[1]]   
     #paramSigmaMatrixT= as.matrix(mcmcResultT$param)[,1][[2]]
     #} else{
    #   paramBetaMatrixT = rbind(paramBetaMatrixT,as.matrix(mcmcResultT$param)[,1][[1]])
    #   paramSigmaMatrixT = rbind(paramSigmaMatrixT,as.matrix(mcmcResultT$param)[,1][[2]]) 
     #} 
     
    # make this more efficient by saving just result 
   #assign(paste0("mcmcResultT",val,m),mcmcResultT)
    # above not wprk later so tr saving param
    
    #assign(paste0("parambeta",val,m),mcmcResultT$param$beta)
    #assign(paste0("paramsigma",val,m),mcmcResultT$param$sigma)
   
    #print(paste0("parambeta",val,m))
    #return(list(paste0("parambeta",val,m),paste0("paramsigma",val,m)))
  }  

}
) # system.time 
print(paste0("mcmcNorm Loop finished, m = ",M))

# try and use lappy instead of loop for M




# can repeat interactively from here
# now loop over the lookup table mg, looping over every pattern
# declare iterate for saving data
m_mg_iter<-0
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
  
  
  # count of pattern by treatment
  cnt<- mg$X1[i]
  
  # treatment grp
  trtgp<- mg$treat[i]
  cat("\ntrtgp = ", trtgp)
  cat("\n Looping within each pattern,no patients = ", cnt)

  #need to convert (relate) treatment group to position in ntreat (create Pindex vector) 
  #something like
  #position(trtgp,"ntreat")
  # in fact this should do it
  trtgpindex<-which(trtgp==ntreat)
  referindex<-which(refer==ntreat) 
  
  # multiple  simulations start here within the pattern loop #########  
  for ( m in  1:M)  { 
    #*FOR INDIVIDUALS WITH NO MISSING DATA COPY COMPLETE DATA INTO THE NEW DATA MATRIX mata_all_new `m' TIMES
    #if `pat' == 0{ 
    m_mg_iter<-m_mg_iter+1
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
      mata_all_newlist[[m_mg_iter]]=mata_new
     # mata_all_new=rbind(mata_all_new,mata_new)
    } else {
      
      #FOR INDIVIDUALS WITH  MISSING DATA  `m' TIMES  
      # dependent on method chosen 
      
      if (meth== 'MAR')  {
        
        #mata_means <- get(paste0("param",trtgp,m))[1]  
        mata_means <- paramBiglist[[M*(trtgpindex-1)+m]][1]
        # mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]
        # convert from list element to matrix
        mata_means <- mata_means[[1]]
        
        
        #Sigmatrt <- get(paste0("param",trtgp,m))[2]
        Sigmatrt <- paramBiglist[[M*(trtgpindex-1)+m]][2]
        S11 <-Sigmatrt[[1]][c_mata_nonmiss,c_mata_nonmiss]
        S12 <-Sigmatrt[[1]][c_mata_nonmiss,c_mata_miss]
        S22 <-Sigmatrt[[1]][c_mata_miss,c_mata_miss]
        
        
      }
      else if (meth == 'J2R') { 
        
        # changed saving the result into  just the param file, list of 2 so can use list index here
       #treatmnets are 1.. M then M+1 ..2M .. etc
        mata_means_trt <- paramBiglist[[M*(trtgpindex-1)+m]][1]
        mata_means_ref <- paramBiglist[[M*(referindex-1)+m]][1]
        
        #mata_means_trt <- get(paste0("param",trtgp,m))[1]
        #mata_means_ref <- get(paste0("param",refer,m))[1]
        #mata_means_trt <- get(paste0("mcmcResultT",trtgp,m,"$param$beta")) 
        #mata_means_ref <- get(paste0("mcmcResultT",refer,m,"$param$beta"))
        
        #mata_means <- get(paste0("parambeta",trtgp,m))
        #one way is to element multiply (because 1,0) then add 
        mata_means_t <- unlist(mata_means_trt)*mata_nonmiss
        mata_means_r <- unlist(mata_means_ref)*mata_miss
        # so when all missing  1,1,1, ... then all contributing comes from reference means      
        mata_means <- mata_means_r+mata_means_t
        # and preserve names   
        colnames(mata_means) <- colnames(mata_means_trt)
        #replicate to number of rows defined by X1 
        #mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]
        
        
        ############# SIGMA is from paramsigma  the reference group ################
        
        
        # do we ever  use SigmaTrt !!?? in j2r? 
        
         SigmaRefer <- paramBiglist[[M*(referindex-1)+m]][2]
        #SigmaRefer <- get(paste0("param",refer,m))[2] 
        
        # get(paste0("mcmcResultT",refer,m,"$param$sigma"))
        # when reading in Stata sigmas
        # needs to to ths as tibble will fail in cholesky
        #SigmaRefer <- as.matrix(get(paste0("paramsigmaStata",refer,m)))
        
        #print(paste0("paramsigma",refer,m))
       # print(paste0("refer,m = ",refer," ",m))   
        #SigmaRefer <- Reduce(rbind,res1sigma.list[m])
        
        # note use of [[1]] as is matrix rathe than list
        S11 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_nonmiss]
        S12 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_miss]
        S22 <-SigmaRefer[[1]][c_mata_miss,c_mata_miss]
        
      }
      else if (meth=='CR') {
        mata_means <- paramBiglist[[M*(referindex-1)+m]][1]
        #mata_means <- get(paste0("param",refer,m))[1]
        # convert from list to matrix
        mata_means <- mata_means[[1]]
        #mata_means_r <- unlist(mata_means_ref)*1
        #mata_means_r <- unlist(mata_means_ref)*mata_miss
        #mata_means <- mata_means_r
        # and preserve names   
        #colnames(mata_means) <- colnames(mata_means_trt)
        #replicate to number of rows defined by X1 
        #mata_means <- mata_means_r
        #colnames(mata_means) <- colnames(mata_means_ref)
        
        #mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]
        
        #SigmaRefer <- get(paste0("param",refer,m))[2]
        SigmaRefer <- paramBiglist[[M*(referindex-1)+m]][2]
        S11 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_nonmiss]
        S12 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_miss]
        S22 <-SigmaRefer[[1]][c_mata_miss,c_mata_miss]
      }
      else if (meth=='CIR') {
        # pre-deviating use mean of trt gp up to last obs time bfore deviating, post-deviating use mean from ref grp 
        
        #mata_means_trt <- get(paste0("parambeta",trtgp,m)) 
        #mata_means_ref <- get(paste0("parambeta",refer,m))
        
        # put equiv to mimix 
        #mata_Means <- get(paste0("param",trtgp,m))[1]
        mata_Means <- paramBiglist[[M*(trtgpindex-1)+m]][1]
        # convert from list to matrix
        mata_Means <- mata_Means[[1]]
        #mata_Means <-  get(paste0("parambeta",trtgp,m))
        #MeansC <-  get(paste0("param",refer,m))[1]
        MeansC <-  paramBiglist[[M*(referindex-1)+m]][1]
        
        #might be better to copy mimix algol
        
        mata_means<-CIR_loop(c_mata_miss,mata_Means,MeansC)
        #returns mata_means as single row
        # then duplicate over patt rows
        #replicate to number of rows defined by X1 
        # mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]
        
        #SigmaRefer <- get(paste0("paramsigma",refer,m))
        
        #SigmaRefer <- get(paste0("param",refer,m))[2] 
        SigmaRefer <- paramBiglist[[M*(referindex-1)+m]][2]
        # when reading in Stata sigmas
        # needs to to ths as tibble will fail in cholesky
        #SigmaRefer <- as.matrix(get(paste0("paramsigmaStata",refer,m)))
        
        #print(paste0("paramsigma",refer,m))
        #SigmaRefer <- Reduce(rbind,res1sigma.list[m])
        
        S11 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_nonmiss]
        S12 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_miss]
        S22 <-SigmaRefer[[1]][c_mata_miss,c_mata_miss]
      }  
      else if (meth=='LMCF') { 
        
        #mata_Means <-  get(paste0("param",trtgp,m))[1]
        mata_Means <- paramBiglist[[M*(trtgpindex-1)+m]][1]
        # convert from list to matrix
        mata_Means <- mata_Means[[1]]
        # no ref MeansC <- mata_means_ref
        mata_means<-LMCF_loop(c_mata_miss,mata_Means)
        #mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]
        
        
        #Sigmatrt <- get(paste0("param",trtgp,m))[2]
        Sigmatrt <- paramBiglist[[M*(trtgpindex-1)+m]][2]
        # when reading in Stata sigmas
        S11 <-Sigmatrt[[1]][c_mata_nonmiss,c_mata_nonmiss]
        S12 <-Sigmatrt[[1]][c_mata_nonmiss,c_mata_miss]
        S22 <-Sigmatrt[[1]][c_mata_miss,c_mata_miss]
      }  #if meth
      
      
      
      # loop still open for row(mg)
      
      ###################### MNAR IMPUTATION ################
      # need insert routine when ALL missing values
      # else
      ###################### MNAR IMPUTATION ################
      
      
      #make sure these are single row vectors! as mistake in LMCF but have to be duplicate rows so add ,s
      #and move after dup fun
      #Error in mata_means[, c_mata_nonmiss] : incorrect number of dimensions
      #m1 <- mata_means[c_mata_nonmiss]
      #m2 <- mata_means[c_mata_miss]
      mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]
      if (!is.null(nrow(mata_means)) )  {
        m1 <- mata_means[,c_mata_nonmiss]
        m2 <- mata_means[,c_mata_miss]
      } else {
        m1 <- mata_means[c_mata_nonmiss]
        m2 <- mata_means[c_mata_miss]
      }
      
      
      #then mata_obs is the sequential selection of rows according to the mimix_group variable X1 values
      
      #need a counter to accumulate j  - easiest way is to ceate a cumulstive col in mimix_group
      
      
      
      j <- mg[i,"X1"]
      
    #for debug  
    #  print(paste0(" count in patt = ", j))
     
      k <- mg[i,"X1cum"]
      startrow <-(k-j+1)  
      stoprow  <-(k)
      
    # for debug
    #  print("startrow stoprow = ")
    #  print(startrow)
    #  print(stoprow)
      
      #raw1 <- mata_obs[, c_mata_nonmiss]
      preraw <- mata_Obs[c(startrow:stoprow),2:ncol(mata_Obs)]
      raw1 <- preraw[,c_mata_nonmiss]
  
    # for debug    
      
    # print("mata_raw = ")
      
      
      t_mimix =cholsolve(Q=S11,y=S12)   
      conds <-  S22-t(S12)%*%t_mimix
      #  
      #meanval = as.matrix(m2) + as.matrix(raw1 - m1)%*%as.matrix(t_mimix)     
      # below for CR but need checks work with J2R
      # need edit m2 as errormsg "data frame with 0 columns and 1 row"
      if (meth=='J2R') {
        meanval = as.matrix(m2) + as.matrix(raw1 - m1)%*%as.matrix(t_mimix)
      } else  {
        meanval = (m2) + as.matrix(raw1 - m1)%*%as.matrix(t_mimix)
      } 
      U <- chol(conds)
      # mg[i,X1] is equiv to Stata counter, miss_count is no. of missing, so 
      miss_count=rowSums(mata_miss)
      # gen erate inverse normal
      Z<-qnorm(matrix(runif( mg[i,"X1"]* miss_count,0,1),mg[i,"X1"],miss_count))
      # check same input parameters for inverse norm gen as in stata 
     
    #for debug   
    #  print(paste0('mg[i,X1] =',mg[i,"X1"],' miss_count= ',miss_count)) 
      
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
      
      #what i this below? 
      #presumablywas  for testing?
      #assign(paste0("Z11",i,"_imp",m),subset(Z))
      #assign(paste0("U11",i,"_imp",m),subset(U))
      # assign(paste0("S12",i,"_imp",m),S12)
      # assign(paste0("S22",i,"_imp",m),S22)
      # assign(paste0("mata_y1_11",i,"_imp",m),mata_y1)
      # assign(paste0("conds11",i,"_imp",m),conds) 
      # no idea why below cuases error!
      # assign(paste0("t_mimix11",i,"_imp",m),t_mimix)
      
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
      
      #this works but bette to pre initialise data structure outsidr loop
      mata_new<-cbind(GI,II,mata_new,SNO)
    
      mata_all_newlist[[m_mg_iter]]=mata_new
      #mata_all_new<-rbind(mata_all_new,mata_new)
      #mata_all_new<- na.omit(mata_all_new)
      #stata equilv \ is row bind NOt cbind!!?
      
    } #if patt ==0r
  } #for row[mg] 
 
} #for M StOP HERE!!

  return(mata_all_newlist)

} # for runmimix test


print(mg)

end_time <- proc.time()
time_taken <- end_time - start_time
print(paste("Time taken:", time_taken[1]))

system.time(analselist("5137")) 

analselist(meth,"5456") 
  
pttestf(10000,10000,1.040,0.385,1.051,0.341)


analse <- function(meth,no)  {
 assign( paste0("mata_all_new_rmna",meth), na.omit(mata_all_new))
 assign( paste0("mata_all_new_rmna",meth,no) , filter(get(paste0("mata_all_new_rmna",meth)),SNO == no))
 print(paste("method= ",meth,"SNO = ",SNO))
 t(round(stat.desc(get(paste0("mata_all_new_rmna",meth,no))[,c("fev2","fev4","fev8","fev12")]),3)[c(1,9,13,4,8,5),])
}

# to handle muliple data sets see
#https://rdrr.io/cran/norm2/man/miInference.html
# but instead of mcmcResult$imp.list just create new list from the saved list output
# mata_all_newlist ?

mata_all_newlist[[110000]]
mata_all_newData10k <- do.call(rbind,mata_all_newlist)
# then sort into M data sets and maybe split into M lists !?
testkm10k<-mata_all_newData10k[order(mata_all_newData10k$II,mata_all_newData10k$SNO),]
# to get the list
testkmlist10k <- split(testkm10k,testkm10k$II)
# so has M elemets in list
est.list <- as.list(NULL)
std.err.list <- as.list( NULL )
for( m in 1:M ){
  yimp <- testkmlist[[m]]  # one imputed dataset
  diff <- yimp[,"fev2"] 
  est.list[[m]] <- mean(diff)
  std.err.list[[m]] <- sqrt( var(diff) / length(diff) ) }
## combine the results by rules of Barnard and Rubin (1999)
## with df.complete = 27, because a paired t-test in a dataset with
## N=28 cases has 27 degrees of freedom
miResult <- miInference(est.list, std.err.list, df.complete=182)
#miResult <- miInference(est.list, std.err.list, confidence=0.95)


# this from CRAN vignette/amelia.pdf
b.out<-NULL
se.out<-NULL
for( m in 1:M ){
   ols.out<-lm(fev12~treat+base,data=testkmlist10k[[m]] )
   b.out<-rbind(b.out,ols.out$coef)
   se.out<-rbind(se.out,coef(summary(ols.out))[,2])
}
library(Amelia)
combined.results<-mi.meld(q=b.out,se=se.out)
write.amelia(obj=mata_all_newData10k,file.stem = "outdata",format = "dta")
print(miResult)



#can observe list of list elements by 
#mata_all_newlist[[101]]
#collapse list into large data set 
#do call obtains a list rather than list of lists, rbind handles elements on list to create a matrix
#then need to sort by imputation to group into M data sets 
mata_all_newData <- do.call(rbind,mata_all_newlist)

df <-mata_all_newData[order(mata_all_newData$II,mata_all_newData$SNO),]
J2r5456<- subset((do.call(rbind,mata_all_newlist)),SNO==5456)
t(round(stat.desc(J2r5456)[,c("fev2","fev4","fev8","fev12")],3)[c(1,9,13,4,8,5),])

#not work
#mata_all_new_unlist <- unlist(mata_all_newlist)
lapply(mata_all_newlist,mean)
#try a loop
subSNOx <- head(mata_all_newlist[[1]],1) 
subSNOx[subSNOx>=0] <-NA
for (i in 1:(M*nrow(mg)))  {
   subSNO<- subset((mata_all_newlist[[i]]),SNO=="5456")
   subSNOx<- rbind(subSNOx,subSNO)
   subSNOx<-na.omit(subSNOx)
}
t(round(stat.desc(subSNOx)[,c("fev2","fev4","fev8","fev12")],3)[c(1,9,13,4,8,5),])
t(round(stat.desc(get(paste0("subSNOx",meth,no))[,c("fev2","fev4","fev8","fev12")]),3)[c(1,9,13,4,8,5),])

# try and combine list elements into data set 
# when m =100, then list 1001..1100 is 3 patients 100 imputation data set



mata_all_new_rmnaJ2R <- na.omit(mata_all_new)
mata_all_new_rmnaJ2R5456 <- filter(mata_all_new_rmnaJ2R,SNO == "5456")
t(round(stat.desc(mata_all_new_rmnaJ2R5456[,c("fev2","fev4","fev8","fev12")]),3)[c(1,9,13,4,8,5),])


analysemeth <- function(meth,no) {
     return(assign(paste0("mata_all_new_rmnaT",meth), (na.omit(mata_all_new) )) )
    # return(datana)
}
 
analyseSNO <- function()  {  
     #datanaSNO<-assign(paste0(datana,no), filter(datana,SNO == no))
     #datana<-assign(paste0(data,"na",meth), na.omit(data))
     return(datana)
    }
analyse(meth,5456)     
     
{
     datanaSNO<-assign(paste0("mata_all_new_rmnaJ2R","5456"), filter(mata_all_new_rmnaJ2R,SNO == 5456))
     descstats<- t(round(stat.desc(datanaSNO[,c("fev2","fev4","fev8","fev12")]),3)[c(1,9,13,4,8,5),])
     return(descstats)
     }
analyse(mata_all_new,meth,5456)  
kmtest<-assign(paste0("mata_all_new_rmna",meth), na.omit(mata_all_new))



