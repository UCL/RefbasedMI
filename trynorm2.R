############################################################################
# R program to mimic stata program mimix
# ie reference based imputation 
# Note 1st part is to set up a summary table - mimix_group
# reflects the pattern and treatment group configuration of the raw data
# then acts as a looping mechanism 
# norm2 is used as MCMC multivariate normal  
#########################################

set.seed(101) #so can reproduce results
# initial program reads stata data files 
library(haven)

preimputemvn1 <- read_dta("preimputemvn1.dta")
preimputemvn2 <- read_dta("preimputemvn2.dta")


#note refer 2 is to -placebo but in preimputemvn1

#These data sets already prepared in Stata - need to write code to perform here
# alt read in raw data set 
asthma <- read_dta("asthma.dta")
library(tidyr)

#recode treat, Placebo = 1, active =2
 
asthma$treatrc<-ifelse(asthma$treat==2,1,2)
#convert to wide format 
asthma_wide <-pivot_wider(asthma,id_cols=id,names_from=time,names_prefix="fev",values_from=c(fev))
#need to add cov col, ie  breaks col
#asthma_wide$ pivot_wider(asthma,id_cols=id,names_from=time,values_from=c(base),values_fn = list(breaks=mean))
# can check non missing
#also need to add the treat col 
treatrc<-asthma %>% group_by(id) %>%summarise(treatrc=first(treatrc),baseq=first(base)) 
#bring all the cols together
asthma_wide<-merge(asthma_wide,treatrc,by= "id")
# split into treatment data-sets
#can't have paste on LHS
for (i in 2:3) {
paste0("preimpute","i") <-  subset(asthma_wide,treatq==i)
}

# 1 is placebo
preimpute1 <-  subset(asthma_wide,treatrc==1)
# 2 is active treatment
preimpute2 <-  subset(asthma_wide,treatrc==2)





#Data augmentation da.norm
#find the starting values by em

#prelim.norm2(preimputemvn1)
#em.norm()
library(norm2)
emNorm(cbind(fev1,fev2,fev3,fev4,base) ~1,data=preimputemvn1)
#summary(emNorm(cbind(fev1,fev2,fev3,fev4,base) ~1,data=preimputemvn1))

emResult1<-emNorm(cbind(fev1,fev2,fev3,fev4,base) ~1,data=preimputemvn1)
#set.seed(5321) #so can reproduce results

#mcmcResult <- mcmcNorm(cbind(fev1,fev2,fev3,fev4,base) ~1,data=preimputemvn1,starting.values = emResult$param)

mcmcResult1<- mcmcNorm(emResult1,iter=5000,multicycle = 100)

#desired outputs as below 
#mean matrix
mcmcResult1$param$beta
#var-cov matrixr 
mcmcResult1$param$sigma

#try and make the beta parameters matrices , once working change the above
mat_Betas <-rbind(mcmcResult1$param$beta ,mcmcResult2$param$beta)
# the 1st parameters are  
mat_Betas[1,]
#2nd betas are 
mat_Betas[2,]


#mat_Sigma <-rbind(mcmcResult1$param$sigma ,mcmcResult2$param$sigma)

# store Siigma parametes in matrixlist to enable easy referal
mat_Sigma <-list(mcmcResult1$param$sigma ,mcmcResult2$param$sigma)

#mat_Sigma[1] is placebo, ie the refernce in this case  


########### As a test can take the parameter files from stata

####  1stly just try the means


test1 <- c(2.057, 1.953,2.032,1.936,2.13)
test2 <- c(2.217,2.33,2.312,2.28,2.1)
mat_Betas <-rbind(test1,test2)
# and replicate rows! (later)

#same for sigma , use stata values

xlist <-c(0.468,.271,.328,.299,.25,.377,.452,.27,.34,.607,.28,.214,.233,.204,.29)

# this from stata sigma param
mcmcResult1$param$sigmaTest<-mcmcResult1$param$sigma
mcmcResult1$param$sigmaTest[,1]<-c(0.468,.271,.328,.299,.25)
mcmcResult1$param$sigmaTest[,2]<-c(0.271,.377,.452,.27,.34)
mcmcResult1$param$sigmaTest[,3]<-c(0.328,.452,.607,.28,.214)
mcmcResult1$param$sigmaTest[,4]<-c(0.299,.27,.28,.233,.204)
mcmcResult1$param$sigmaTest[,5]<-c(0.25,.34,.214,.204,.29)

i=1
paste0("mcmcResult",i,"$","param$beta")
names() 
suffix <- seq(1:2)
suffix
mynames <-paste0(mcmcResult,suffix)
mynames







#prelim.norm2(preimputemvn1)
#em.norm()

emNorm(cbind(fev1,fev2,fev3,fev4,base) ~1,data=preimputemvn2)
#summary(emNorm(cbind(fev1,fev2,fev3,fev4,base) ~1,data=preimputemvn1))

emResult2<-emNorm(cbind(fev1,fev2,fev3,fev4,base) ~1,data=preimputemvn2)


#mcmcResult <- mcmcNorm(cbind(fev1,fev2,fev3,fev4,base) ~1,data=preimputemvn1,starting.values = emResult$param)

mcmcResult2<- mcmcNorm(emResult2,iter=5000,multicycle = 100)
mcmcResult2$param
mcmcResult1$param
#  param gives the betas (means) and sigmas (var matrice)

# try and create mimix_group(matrix) data set equilv in R
#read input data
library(haven)
mimix_d2 <- read_dta("mimix_d2.dta")
# mimix_d2 is the data form which msiing pattern found

#need similar to collapse cmd
# cycle over the fev vars to note NAs
data <- mimix_d2
colnames(data)[colnames(data)=="__00000D"] <- "refer"
colnames(data)[colnames(data)=="__000001"] <- "method"
colnames(data)[colnames(data)=="trt_grp"] <- "method"
if (is.na(data$fev1) ) {data$pat <-1}
f <- function(f)
ifelse (is.na(data$fev1) & is.na(data$fev2) & is.na(data$fev3) & is.na(data$fev4)( data$pat <-1,)
data$pat
which(is.na(data$fev1))

if (is.na(data$fev1) ) {data$pat <-1}
f <- function(f)
  ifelse (is.na(f) & is.na(data$fev2) & is.na(data$fev3) & is.na(data$fev4)( data$pat <-1,)

# want to aggregate but after generating a pattern variable
apply(data,1,is.na)

logicna <- (is.na(data$fev1) &  is.na(data$fev2) & is.na(data$fev3) & is.na(data$fev4))  
logic_3 <-


data[!logicna,]
data$patt <- logicna
case when (is.na(data$fev1) &  !is.na(data$fev2) & !is.na(data$fev3) & !is.na(data$fev4))  

#before aggreating would like grouping var, ie define pattern

varlist <- c("fev1","fev2","fev3","fev4")
library(dplyr)
#found better way then this!
data$patt <- ifelse(is.na(data$fev1),7,
                      ifelse(is.na(data$fev2),14,
                      ifelse(is.na(data$fev3),12,
                      ifelse(is.na(data$fev4),4,0 ))) ) 

data$patt <- ifelse(ata$patt <- ifelse(is.na(data$fev1),1,
                                       ifelse(is.na(data$fev2),2,
                                              ifelse(is.na(data$fev3),3,
                                                     ifelse(is.na(data$fev4),4,0 ))) ) ) 
datafev <- data.frame(data$treat,data$fev1,data$fev2,data$fev3,data$fev4,1)
datafev <- data.frame("treat" = data$treat,"fev1" = data$fev1,"fev2" = data$fev2,"fev3" = data$fev3,"fev4" =data$fev4,1,"refer"=data$refer,"method"=data$method)
#convet fevs into 0 and 1s then aggreg should get the pattern
datafev$fev1 <-  ifelse(!is.na(data$fev1),0,1)
datafev$fev2 <-  ifelse(!is.na(data$fev2),0,1)
datafev$fev3 <-  ifelse(!is.na(data$fev3),0,1)
datafev$fev4 <-  ifelse(!is.na(data$fev4),0,1)

sort.datafev <-datafev[order(datafev$treat,datafev$fev1,datafev$fev2,datafev$fev3,datafev$fev4),]
sort.datafev1na <-datafev[order("fev1na",na.last = FALSE),]

#Yes this works!!  tryagg equiv to mimix_group except need to create patt variable same as Stata  
tryagg <- aggregate(sort.datafev,by =list(sort.datafev$treat,sort.datafev$fev1,sort.datafev$fev2,sort.datafev$fev3,sort.datafev$fev4),FUN="sum")
#convert non nas into 0's?
#next stage to obtain   

#this works!  to obtain mimix_group
#but summaris efun migh be nicer as keeps orig var names
tryagg$patttest <- tryagg$Group.2+tryagg$Group.3*2^1+tryagg$Group.4*2^2+tryagg$Group.5*2^3 
mimix_group<-data.frame(tryagg[order(tryagg$Group.1),])
#above works but also need the other variables corresponding to mimix_group
#i.e.
tryagg$refer<-c(2)
tryagg$method<-c(3)
#reset row names 
rownames(mimix_group) <-seq(length=nrow(mimix_group))

tail(data)


#pseudo code from stata program 
#after norm2 runs we need from the parms file, the means and cov matrix  
#mimix_all(i) derived from parms for i treatments  

#mimix_group used after the parms routine and before the method routine

#obtain the means and vars matrices from parms
prefix <- "preimputemvn"
suffix <- seq(1:2)

my.names <-paste0(prefix,suffix)

<- read_dta(paste0("preimputemvn1.dta")

my.names <- read_dta(paste0(prefix,suffix,".dta"))
my.names <- (paste0(prefix,suffix,".dta"))
> my.names
prefix <- "preimputemvn"
suffix <- seq(1:2)

left.names <-paste0(prefix,suffix)
a = c(1,2)
for(i in a) 
{
  #i = 1:2  
  #left <-paste0(prefix,i)
  paste0(prefix,i)   read_dta(my.names[i])
}  

left

}
preimputemvn1

mcmcResult2$param$beta
mcmcResult2$param$sigma
# these are equivalnet to  mata_mean_group`i'_imp`k' and mata_VAR_group`i'_imp`k'
# so use these to obtain  matrices fro method = 3

# loop through the max_indicator var as found form agg cmds 


# so have mimix_group,  mcmcResulti$param$beta  & mcmcResulti$param$sigma

# need to loop thru the treatments (2 in thus case) and 
 
# loop thru the mimix_group matrix  obtaining mata_obs for each loop  
#taking values such as trt_grp , pat is 0,7,8... , counter is no. in grp


#try to obtain S11 and S12 from sigma matrix using pattern from mimix_group

#hard code
mimix_group$refer<-1

mimix_group
#assume covar has no missing values
mata_miss <- mimix_group[1,c(2,3,4,5)]
mata_miss$group.6 <-1
mata_miss
mata_nonmiss <- t(ifelse(mata_miss==0,1,ifelse(mata_miss==1,0,mata_miss))) 
t(mata_nonmiss)

#now the 2nd loop
mata_miss <- mimix_group[2,c(2,3,4,5)]
mata_miss$group.6 <-1
mimix_group$cumX1 <-cumsum(mimix_group$X1)
#add a 2nd col to enable row selection
#not used
#mimix_group$nextcumX1[1] <-1
#mimix_group$nextcumX1 <-(mimix_group$X1) + mimix_group$cumX1 
#declare miss_count as sum across patt
mimix_group$sum_1s <-rowSums(mimix_group[,2:5])


# need to turn into a loop running thru mimix_group
# nrow(mimx_group) corrrsponds to max_indicator

for ( i in  1:nrow(mimix_group)) 
{
  mata_miss <- mimix_group[i,c(2,3,4,5)]
  #assumes covariate non missing
  mata_miss$group.6 <-0
  mata_nonmiss <- (ifelse(mata_miss==0,1,0))
                           
 # immediat-e below for debug 
  print(i)
  print(mata_miss)
  print(mata_nonmiss)
}
#23/9} delte out to increment i for next processing?

# note mcmcResult2 pretty big!

#S11 uses mat_nonmiss, so try here
S11 <- mcmcResult2$param$sigma[mata_nonmiss,mata_nonmiss]

# want to extract from Sigma using mata_miss, nonmiss same as mimix using as indices?
# like  S11 = Sigma[mata_S_nonmiss, mata_S_nonmiss]
# looks like take the elemts from 1st colunm and same elemnts from last row.


mata_miss*mcmcResult2$param$sigma

mcmcResult2$param$sigma[c(1,5),c(2,3,4)]
mcmcResult2$param$sigma[c(1,5),c(1,5)]  # works for S11 nonmiss, nonmiss  
mcmcResult2$param$sigma[c(2,3,4),c(1,5)] # not works !

#so need transform nonmiss,miss to c lists.
c_mata_miss<-which(mata_miss==1)
c_mata_nonmiss<-which(mata_nonmiss==1)
#test
S11 <- mcmcResult1$param$sigma[c_mata_nonmiss,c_mata_nonmiss]
S12 <- mcmcResult1$param$sigma[c_mata_nonmiss,c_mata_miss]
S22 <- mcmcResult1$param$sigma[c_mata_miss,c_mata_miss]


######### START main analysis run HERE ####################################

#CREATE AN EMPTY MATRIX FOR COLLECTING IMPUTED DATA 
#mata_all_new initialises as a single row with ...'s
#then recursively added (later)
mata_all_new <- matrix(nrow=1, ncol=8) #empty matrix
# if any left over from previous
# been using this as starter point  
mata_all_new<-mata_all_new[-c(1:nrow(mata_all_new)),]
#hard code for now
names(mata_all_new) <-names(mata_new2)


#need to loop thru mimix_group ,define  (nrows(mimix_group) equiv to max_indicator)
# NOtE as  cyling thru read off the treatment and reference values from relevant cols.   
for ( i in  1:nrow(mimix_group)) 
{
  mata_miss <- mimix_group[i,c(2,3,4,5)]  #define mata_miss
  #assumes covariate non missing
  mata_miss$group.6 <-0                     #assuming cov col is non missing
  mata_nonmiss <- (ifelse(mata_miss==0,1,0))  #define mata-nonmiss from miss
  
 # need transform nonmiss,miss to c lists.
  c_mata_miss<-which(mata_miss==1)
 # if non missing then have to amend above !!!  
  #S11 will be all of matrix, S12,S22 will be Null
  #if c_mata_miss
 # 1's signify missing values, 0 indicates not missing.
 # use if statement to create suitable mat_obs submat     
  
  c_mata_nonmiss<-which(mata_nonmiss==1)

#j2r requires VAR from reference grp  , 
# as a list elemt need convert to matrix  
# here refernce =1  
  
   SigmaRef <- Reduce(rbind,mat_Sigma[1])
   
   S11 <-SigmaRef[c_mata_nonmiss,c_mata_nonmiss]
   S12 <-SigmaRef[c_mata_nonmiss,c_mata_miss]
   S22 <-SigmaRef[c_mata_miss,c_mata_miss]
                  
   
#  S11 <- mcmcResult1$param$sigma[c_mata_nonmiss,c_mata_nonmiss]
#  S12 <- mcmcResult1$param$sigma[c_mata_nonmiss,c_mata_miss]
#  S22 <- mcmcResult1$param$sigma[c_mata_miss,c_mata_miss]
  
  print("S12 i = " )
  print(S12)
  print(i)
  
# need expand mata_means to same no rows as counter but only after deriving mata_Means
# for  J2r I believe is case of taking relevant treatment betas and overwriting when missing by reference treatment if different     
  
# count of pattern by treatment
  cnt<- mimix_group$X1[i]
# treatment grp
  trt_gp<- mimix_group$Grouptrt[i]
# reference grp
  refer <- mimix_group$refer[i] 
# try to mcmcrEsult as 2 matrices  
   
   #DFbeta<-mcmcResult[refer]1$param$beta
   DFbeta<-mat_Betas[trt_gp,]
   
# want the treatment grp means and where missing overwrite with refer grp (when differs!)     
   mata_means_trt <-mat_Betas[trt_gp,] 
   mata_means_ref <-mat_Betas[refer,]
 
#one way is to element multiply (because 1,0) then add 
   mata_means_t <- mata_means_trt*mata_nonmiss
   mata_means_r <- mata_means_ref*mata_miss
# so when all missing  1,1,1, ... then all contributing comes from reference means      
   mata_means <- mata_means_r+mata_means_t
# and preserve names   
   names(mata_means) <- names(mata_means_trt)

# replicate to number of rows defined by X1 
   mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mimix_group$X1[i]),]
   
   
# these are the replications   
 #  betasX1<- DFbeta[rep(seq(nrow(DFbeta)),each=mimix_group$X1[i]),] 
 #   print("betasX1")
 #   print(betasX1)
    #matrix(rep(mcmcResult1$param$beta,each=10,nrow=10))
# i not incrementing below so comment out !   
#}

#now can go to *MNAR IMPUTATION equlv section:  but within the i loop
#J2R
# uses means treatment group and for missing, replace with refer group means   
   

   


   
   # mat_Betas
   #m1<- betasX1[c_mata_nonmiss]
   #m2<- betasX1[c_mata_miss]

#raw1 is found from mata_obs  fom the sorted data (ie mimix_d2) by subsetting on the pattern 
#also reorder index col from 1st to  last
mimix_d2c <- mimix_d2[,c(2,3:ncol(mimix_d2),1)]
#must use integers to index cols
#select cols from ordered input data
mimix_d2cols <- mimix_d2c[,c(1:4,7,ncol(mimix_d2c))]


#then mata_obs is the sequential selection of rows according to the mimix_group variable X1 values
# stil goimg thru the rows in mimix_group want to read off the X1 values to create subsets   
#for (l in 1:mimix_group[i,"X1"]) {
    j <- mimix_group[i,"X1"]
    print("j=")
    print(j)
   # mata_obs <- mimix_d2cols[i:j,]
  
#need a counter to accumulate j  - easiest way is to ceate a cumulstive col in mimix_group
#sol is to use cumX1 minus X1 ,thats all there is to it
#for ( i in  1:nrow(mimix_group))  {
   #print(i)
   j <- mimix_group[i,"X1"]
   k <- mimix_group[i,"cumX1"]
   startrow <-(k-j+1)  
   stoprow  <-(k)
   print("startrow stoprow = ")
   print(startrow)
   print(stoprow)
   
   
   mata_obs <- mimix_d2cols[c(startrow:stoprow),]
   print("mata_obs = ")
   print(mata_obs)
#} 
   
   #J2r uses Sigma derived from refernce group
   #Mean 
 
#NOTE in stata no t calcs as no missing pattern for 1st 37 obs. so need insert If statement somewhere!
 if(length(c_mata_miss)!=0 ) {
    raw1 <- mata_obs[, c_mata_nonmiss]
    t=cholsolve(S11,S12)
#cholsolve(A, B) solves AX=B and returns X for symmetric (Hermitian), positive-definite A.  cholsolve() returns a matrix of missing values if A is not positive definite or if A is singular.
# so want solve S11X=S12

  library(sparseinv)
  t_mimix =cholsolve(Q=S11,y=S12)   
   conds <-  S22-t(S12)%*%t_mimix
 
   
   U <- chol(conds)
   #stata counter is X1 and misscount is sum of miss pattern (adding 1's )from mimix_group
    mimix_group[i,"X1"]
   #in stata mata uniform generates values of uniform distrib between 0-1 for dimension of matrix given!  
    #inverse of CDF given by qnorm
    #runif generates n uniform distributed random numbers in given interval(in this case 0,1)
    # need a matrix dimension X1 by sum_1s 
    #qnorm takes inverse cdf
  #seed enusre get same result each time
    set.seed(101)
    
    m1 <- mata_means[c_mata_nonmiss]
    m2 <- mata_means[c_mata_miss]
  print("mata_means = ")
  print(mata_means) 
      
    Z<-qnorm(matrix(runif( mimix_group[i,"X1"]* mimix_group[i,"sum_1s"],0,1),mimix_group[i,"X1"],mimix_group[i,"sum_1s"]))
  meanval = as.matrix(m2) + as.matrix(raw1 - m1)%*%as.matrix(t_mimix)
  mata_y1 = meanval+Z%*%t(U)

  # betas correctly programmed as used stata values in test, difference is in t, ie cholsolve, ie S11 S12 matrices
  
    
  print(mata_y1)         
  # define mata_new
  # mata_new 1st consists of mat_obs (same as is used mata_nonmiss, seems unnecessary calculation)
  # then fill in the missing columns with the values obtained from mata_y1   
  # so I assume mata_y1 are the imputed values!!
  # I think is just ok to add the missing columns to  mata_obs thoug Stat does bit extra
  # nct = no cov + ntime 
  # in gthis instance cheat by hard coding 5
  nct  <-5 
 #define new matrix from observed  
  mata_new <- mata_obs[,c(1:nct)]
 #filling in the missing values (columns) 
  mata_new[,c_mata_miss] <- mata_y1
}  # if(length(c_mata_miss)!=0 )
   if(length(c_mata_miss)==0 ) { mata_new <- mata_obs[,c(1:nct)]}
     
   #imp is the imputation number (this instance we only do 1)
  #trt_grp is control or tretament in 1st col here 
  imp <-1
  SNO <- mata_obs[,ncol(mata_obs)]
  # GI treatment grp column 1 (here),II imputation number col, mata_new matrix then SNO is id col.
  GI <- array(data=mimix_group[i,1],dim=c(mimix_group[i,"X1"],1))
  II <- array(data=imp,dim=c(mimix_group[i,"X1"],1))
  mata_new2<-cbind(GI,II,mata_new,SNO)
  # add on in 1st and 2nd colsmata
  #names(mata_new2) <- NULL
  
  
  mata_all_new <-bind_rows(mata_all_new,mata_new2)
  print(mata_all_new)
 }

   
   
