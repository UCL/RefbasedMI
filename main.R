install.packages(c("tidyr","dplyr","plyr","haven","pastecs","norm2","sparseinv"))
library(tidyr)
#pivot_wider

#for select
library(dplyr)
library(plyr)
#read data
# need haven
library(haven)
library(pastecs)
library(norm2)
library(sparseinv) 

# need also input refer

refer <-2
getwd()
source("N:/Documents/GitHub/mimix/mimixR/functions.R")

# 
rm(list = ls())
source("functions.R")
mxdata <-read.csv("./asthma.csv")
# empyty R environmet?

testlist<-preprodata("fev","treat","id","time","base",1000,2,"MAR")

# returns list from preprodata function
ntreat<-testlist[[4]]
finaldatS<-testlist[[2]]
mg<-testlist[[5]]
# vital to get the mata_obs correctly sorted! so corresponds with mimix_group lookup 
# to be consistent with Stata move the base col after the fevs!
mata_Obs <- testlist[[2]]
M <- testlist[[11]]
refer <- testlist[[12]]
meth <- testlist[[13]]

# so tsts2 is th equiv to mi_impute so try to put in norm2
tst2 <- mi_impute("id","time","fev","base")
# put commas in
tst3<-paste(tst2,collapse = ",")

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

#Warning message:
#In Ops.factor(left, right) : ‘>=’ not meaningful for factors



for (val in t(ntreat)) {
  print(paste0("prenormdat",val))
  assign(paste0("prenormdat",val),subset(finaldatS,treat==val))
}  

# run the mcmc simulations over the treatment grps 
for (val in t(ntreat)) {
  #prnormobj <- subset(prenormdat2, select=c(tst2))
  #prnormobj <- assign(paste0("prenormdat",val),subset(finaldatS,treat==val))
  kmvar=get(paste0("prenormdat",val))
  prnormobj<-assign(paste0("prnormobj",val), subset(kmvar, select=c(tst2)))
  for(m in 1:M) {  
    #emResultT<-emNorm(testprnormobj,prior = "ridge",prior.df=0.5)
    emResultT<-emNorm(prnormobj,prior = "jeffreys")
    mcmcResultT<- mcmcNorm(emResultT,iter=1000,multicycle = NULL,prior = "jeffreys")

    
    #try to make this more efficient by saving just result 
   #assign(paste0("mcmcResultT",val,m),mcmcResultT)
    # above not wprk later so tr saving param
    assign(paste0("param",val,m),mcmcResultT$param)
    #assign(paste0("parambeta",val,m),mcmcResultT$param$beta)
    #assign(paste0("paramsigma",val,m),mcmcResultT$param$sigma)
   
    #print(paste0("parambeta",val,m))
    #return(list(paste0("parambeta",val,m),paste0("paramsigma",val,m)))
  }  
  print(paste0("mcmcNorm Loop finished, m = ",M))
}

# can repeat interactively from here
# now loop over the lookup table mg, looping over every pattern
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
  cat("trtgp = ", trtgp)
  cat("multiple  simulations start here within the pattern loop")
  
# multiple  simulations start here within the pattern loop #########  
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
   
    #FOR INDIVIDUALS WITH  MISSING DATA  `m' TIMES  
    # dependent on method chosen 

        if (meth== 'MAR')  {
  
               mata_means <- get(paste0("param",trtgp,m))[1]  
  # mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]
               # convert from list element to matrix
               mata_means <- mata_means[[1]]
            
  
               Sigmatrt <- get(paste0("param",trtgp,m))[2]
               S11 <-Sigmatrt[[1]][c_mata_nonmiss,c_mata_nonmiss]
               S12 <-Sigmatrt[[1]][c_mata_nonmiss,c_mata_miss]
               S22 <-Sigmatrt[[1]][c_mata_miss,c_mata_miss]
  
  
        }
     else if (meth == 'J2R') { 
       
             # changed saving the result into  just the param file ,list of 2 so can use list index here
                                
               mata_means_trt <- get(paste0("param",trtgp,m))[1]
               mata_means_ref <- get(paste0("param",refer,m))[1]
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
               
            SigmaRefer <- get(paste0("param",refer,m))[2] 
          
                # get(paste0("mcmcResultT",refer,m,"$param$sigma"))
  # when reading in Stata sigmas
  # needs to to ths as tibble will fail in cholesky
  #SigmaRefer <- as.matrix(get(paste0("paramsigmaStata",refer,m)))
  
            #print(paste0("paramsigma",refer,m))
            print(paste0("refer,m = ",refer," ",m))   
  #SigmaRefer <- Reduce(rbind,res1sigma.list[m])
  
            # note use of [[1]] as is matrix rathe than list
            S11 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_nonmiss]
            S12 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_miss]
            S22 <-SigmaRefer[[1]][c_mata_miss,c_mata_miss]
  
          }
              else if (meth=='CR') {
              mata_means <- get(paste0("param",refer,m))[1]
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
  
                SigmaRefer <- get(paste0("param",refer,m))[2]
                S11 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_nonmiss]
                S12 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_miss]
                S22 <-SigmaRefer[[1]][c_mata_miss,c_mata_miss]
         }
           else if (meth=='CIR') {
  # pre-deviating use mean of trt gp up to last obs time bfore deviating, post-deviating use mean from ref grp 
  
  #mata_means_trt <- get(paste0("parambeta",trtgp,m)) 
  #mata_means_ref <- get(paste0("parambeta",refer,m))
  
            # put equiv to mimix 
             mata_Means <- get(paste0("param",trtgp,m))[1]
            # convert from list to matrix
             mata_Means <- mata_Means[[1]]
              #mata_Means <-  get(paste0("parambeta",trtgp,m))
              MeansC <-  get(paste0("param",refer,m))[1]
  
  #might be better to copy mimix algol
  
             mata_means<-CIR_loop(c_mata_miss,mata_Means,MeansC)
  #returns mata_means as single row
  # then duplicate over patt rows
  #replicate to number of rows defined by X1 
  # mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]
  
             #SigmaRefer <- get(paste0("paramsigma",refer,m))
             SigmaRefer <- get(paste0("param",refer,m))[2] 
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
                
              mata_Means <-  get(paste0("param",trtgp,m))[1]
              # convert from list to matrix
              mata_Means <- mata_Means[[1]]
              # no ref MeansC <- mata_means_ref
               mata_means<-LMCF_loop(c_mata_miss,mata_Means)
              #mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]
  
              
             Sigmatrt <- get(paste0("param",trtgp,m))[2]
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
          mata_new<-cbind(GI,II,mata_new,SNO)
          mata_all_new<-rbind(mata_all_new,mata_new)
#stata equilv \ is row bind NOt cbind!!?
    
    } #if patt ==0r
  } #for row[mg] 
 
} #for M StOP HERE!!

analse(meth,"5456")
pttestf(1000,1000,1.652,0.430,1.674,0.414)


analse <- function(meth,no)  {
 assign( paste0("mata_all_new_rmna",meth), na.omit(mata_all_new))
 assign( paste0("mata_all_new_rmna",meth,no) , filter(get(paste0("mata_all_new_rmna",meth)),SNO == no))
 print(paste("method= ",meth,"SNO = ",SNO))
 t(round(stat.desc(get(paste0("mata_all_new_rmna",meth,no))[,c("fev2","fev4","fev8","fev12")]),3)[c(1,9,13,4,8,5),])
}


    

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



