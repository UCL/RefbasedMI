# to obtain mimix_group more elegantly than in trynorm2.R

#instead of agg used pivot_wide to get preimputei dfs 
#convert to wide format 
asthma_wide <-pivot_wider(asthma,id_cols=id,names_from=time,names_prefix="fev",values_from=c(fev))
#need to add cov col, ie  breaks col
#asthma_wide$ pivot_wider(asthma,id_cols=id,names_from=time,values_from=c(base),values_fn = list(breaks=mean))
# can check non missing
#also need to add the treat col 
treatq<-asthma %>% group_by(id) %>%summarise(treatq=first(treat),baseq=first(base)) 
#bring all the cols together
asthma_wide<-merge(asthma_wide,treatq,by= "id")
# split into treatment data-sets
#can't have paste on LHS
for (i in 2:3) {
  paste0("preimpute","i") <-  subset(asthma_wide,treatq==i)
}



#create the patt variable, powers of 2  after creating the 0,1 rows!

#create the 0,0,1.. vectors, 1 signifying missing
#this does create new matrix cols for all cols , some will be redundant but less likely to miss-code 
#but complicated to refer 
#preimpute1$Missdummy <-  ifelse(!is.na(preimpute1),0,1)
#easier to recode indiviually
asthma_wide$dummyfev2 <-  ifelse(!is.na(asthma_wide$fev2),0,1)
asthma_wide$dummyfev4 <-  ifelse(!is.na(asthma_wide$fev4),0,1)
asthma_wide$dummyfev8 <-  ifelse(!is.na(asthma_wide$fev8),0,1)
asthma_wide$dummyfev12 <-  ifelse(!is.na(asthma_wide$fev12),0,1)

#create psttern flag
asthma_wide$pattflag <- (asthma_wide$dummyfev2+asthma_wide$dummyfev4*2^1+asthma_wide$dummyfev8*2^2+asthma_wide$dummyfev12*2^3 )

#sort by id,treat and pattflag
asthma_wide<-asthma_wide[order(asthma_wide$treat,asthma_wide$pattflag),]

# to obtain mimix_group
#mimix_group is summary by the pattflag and treatments,in this case 11 combinations , 
# so add a column of 1's to asthma_wide to act for sum
asthma_wide$X1 <-1
Myagg <-aggregate(asthma_wide[c("pattflag","treatq","dummyfev2","dummyfev4","dummyfev8","dummyfev12","X1")],by=list(asthma_wide$pattflag,asthma_wide$treatq),FUN="sum")
# note Group.1, Group.2 refers to  pattflag, treatq - and can retrieve dummy values by testing if 0 then 0 else 1.
# need to convert dummyfev vars back to dummy's
Myagg$dummyfev2 <- ifelse(Myagg$dummyfev2==0,0,1)
Myagg$dummyfev4 <- ifelse(Myagg$dummyfev4==0,0,1)
Myagg$dummyfev8 <- ifelse(Myagg$dummyfev8==0,0,1)
Myagg$dummyfev12 <- ifelse(Myagg$dummyfev12==0,0,1)

# as want to add sum_1s  and cumX1 cols.
Myagg$cumX1 <-cumsum(Myagg$X1)
Myagg$sum_1s <-rowSums(ifelse(Myagg[,5:8]==0,0,1))


#Myagg is the mimix_group data set. 

#tryagg <- aggregate(sort.datafev,by =list(sort.datafev$treat,sort.datafev$fev1,sort.datafev$fev2,sort.datafev$fev3,sort.datafev$fev4),FUN="sum")
#trysumm <-summarize(asthma_wide,by=list(pattflag),n=n())



preimpute2 <-  subset(asthma_wide,treatq==2)
preimpute3 <-  subset(asthma_wide,treatq==3)



# nrow(mimx_group) corrrsponds to max_indicator
# testing routine below    
for ( i in  1:nrow(Myagg)) 
{
  mata_miss <- Myagg[,c(5,6,7,8)]
  #assumes covariate non missing
  mata_miss$group.6 <-0
  mata_nonmiss <- (ifelse(mata_miss==0,1,0))
  
  # immediat-e below for debug 
  print(i)
  print(mata_miss)
  print(mata_nonmiss)
}


# need transform nonmiss,miss to c lists.
c_mata_miss<-which(mata_miss==1)
c_mata_nonmiss<-which(mata_nonmiss==1)
#test
#S11 <- mcmcResult1$param$sigma[c_mata_nonmiss,c_mata_nonmiss]
#S12 <- mcmcResult1$param$sigma[c_mata_nonmiss,c_mata_miss]
#S22 <- mcmcResult1$param$sigma[c_mata_miss,c_mata_miss]


# Data Augmentation - Using norm2
#find the starting values by em

prelim.norm2(preimputemvn1)
em.norm()

emNorm(cbind(fev1,fev2,fev3,fev4,base) ~1,data=preimputemvn1)
#summary(emNorm(cbind(fev1,fev2,fev3,fev4,base) ~1,data=preimputemvn1))

emResult1<-emNorm(cbind(fev1,fev2,fev3,fev4,base) ~1,data=preimputemvn1)
set.seed(5321) #so can reproduce results

#mcmcResult <- mcmcNorm(cbind(fev1,fev2,fev3,fev4,base) ~1,data=preimputemvn1,starting.values = emResult$param)

mcmcResult1<- mcmcNorm(emResult1,iter=5000,multicycle = 100)

#desired outputs as below 
#mean matrix
mcmcResult1$param$beta
#var-cov matrixr 
mcmcResult1$param$sigma



