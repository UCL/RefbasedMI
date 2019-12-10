

# trying to write a subroutine for this part of data preprocessing
A <- function(x) {  ifelse(!is.na(x),0,1) }


# this section to convert from long to wide data format
preprodata<- function(depvar,treatvar,idvar,timevar,covar,M,refer,meth)  {
  #extract relevant vars
 
  fevdata<-select(mxdata,idvar,depvar,timevar)
  # only want to widen the dose var, ie fev, so take the other variables
  #reshape to wide and assign new names 
  # prefix must be supplied from input argument rather than hard coded, this canbe hard coded
  sts4<-pivot_wider(fevdata,id_cols=c(idvar),names_from=timevar,names_prefix=depvar,values_from=depvar)
  #sts4 is just response data, can join later treat and covar cols from finaldat 
  #assumes no NA for these vars should add in checking routines !!
  uniqdat<-unique(mxdata[c(idvar,covar,treatvar)])
  
  
  finaldat<- merge(sts4,uniqdat,by=idvar)
  # now try and sort on treatvar, doesnt work so instead of sorting , select on treat
 
  
  # need find no. ntreat  to loop over
  #ntreat<-unique(unique(mxdata$treatvar))
  ntreatcol<-(select(mxdata,treatvar))
  ntreat <- unique(ntreatcol)
  ntimecol<-(select(mxdata,timevar))
  ntime<-unique(ntimecol)
  
  
  #in order to aggregate by pattern need create dummy vars
  
  # for sts4 array drop 1st col and create newvars for rest
  STSdummy<-apply(sts4[,2:ncol(sts4)],MARGIN=2,FUN=A)
  # append to names
  colnames(STSdummy) <- paste(colnames(STSdummy),'.1')
  # merge back on to data 
  sts4D<-cbind(finaldat,STSdummy)
  # create powrs of 2
  pows2 <- sapply(1:ncol(STSdummy),function(i) STSdummy[,i]*2^(i-1))
  #need to add up to find patt
  patt <- rowSums(pows2)
  #sts4Dpatt<-cbind(sts4D,patt)
  sts4Dpatt<-cbind(sts4D,patt)
  
  # then sort by patt( last col) and treat and find cumulative sequence  AND split on treatment var!
  # so need to merge treat (and covar) back into response data 
  
  #sort by treatvar and patt var, couldnt get treatvar working with patt so try sorting previously on treatvar   
  finaldatSS <-sts4Dpatt[order(patt),] 
  
  ndx = order(finaldatSS[,treatvar])
  finaldatS <- finaldatSS[ndx,]
  
  # consistent with Stata to get right order, drop id and treat cols
  drops <-c(idvar,treatvar)
  
  
  #also get the dummy pattern vectors by treatment
  #dummypatt<-unique(xSTSdummy)
  pattmat<-unique(STSdummy[,1:ncol(STSdummy)])
  
  # sort of hard coded for now may have to revisit this!
  pows2d <- sapply(1:ncol(pattmat),function(i) pattmat[,i]*2^(i-1)) 
  patt<-rowSums(pows2d)
  #want  join this to ex from mimix_group table
  #merge(uniqdat,sts4,by=idvar)
  pattmatd<-cbind(pattmat,patt)
  
  # this best for  count, implicit sort by treat
  # may have to generalise using function argument  later!
  ex1 <- sts4Dpatt %>%
    dplyr::group_by(sts4Dpatt$treat,sts4Dpatt$patt) %>%
    dplyr::summarise(X1 =n())
  ex1 <- dplyr::rename(ex1,treat='sts4Dpatt$treat',patt='sts4Dpatt$patt')
  # and now find cumX1
  ex1$X1cum <- cumsum(ex1$X1) 
  #want  join this to ex from mimix_group table
  # create index number to preserve ordr after merge
  ex1$exid <- 1:nrow(ex1)
  ex1id <-merge(ex1,pattmatd,by="patt")
  ex1s<-ex1id[order(ex1id$exid),]
  
 
  for (val in t(ntreat)) {
    print(paste0("prenormMC",val))
    #assign(paste0("prenorm2",val),subset(finaldat,treat==val))
  } 
  
  return(list(sts4Dpatt,finaldatS,finaldat,ntreat,ex1s,ex1,ex1id,pattmat,patt,ntime,M,refer,meth))
}
# Question is should the wide data be sorted for analysis or just ordered in fact, selecting out treatmet and rbinding? to provide analysis data set.   
# answer is should be sorted by patt as in Stata!
# therefore just select on treat when looping thru ntreat 


mi_impute <-function(idvar,timevar,depvar,covar) {
  #preprodata(depvar,treatvar,idvar,timevar,covar)
  tst<-pivot_wider(mxdata,id_cols=c(idvar),names_from=timevar,names_prefix=depvar,values_from=depvar)
  # add ont covariates and add on to resp list
  respvars<-c(names(tst[,-1]),covar)
  return(respvars)
}


# analse function to be used after main run to output  for summary stats
analse <- function(meth,no)  {
  assign( paste0("mata_all_new_rmna",meth), na.omit(mata_all_new))
  assign( paste0("mata_all_new_rmna",meth,no) , filter(get(paste0("mata_all_new_rmna",meth)),SNO == no))
  print(paste0("method= ",meth,"SNO = ",SNO))
  t(round(stat.desc(get(paste0("mata_all_new_rmna",meth,no))[,c("fev2","fev4","fev8","fev12")]),3)[c(1,9,13,4,8,5),])
}

# mata_S_miss something like [2 3 4] ,so is cc
CIR_loop <- function(c_mata_miss,mata_Means,MeansC)
{    
  miss_count <- length(c_mata_miss)
  mata_means <- mata_Means
  
  for (b in 1:miss_count)  {
    
    # if 1st col missing then no value before so need to check for that
    # note mimix counter just no rows in each patt
    # so main looping is over missin fields, ie miss_count
    # if 1st col must be ref  
    if (c_mata_miss[b] ==1) {   
      #print will cause objectto return 
      #print(paste0(miss_count))
      #for each patt_count (from mimix_group) 
      #test purposes
      #count<-2
      
      mata_means[b] <- MeansC[b] 
    } else  {
      mata_means[c_mata_miss[b]] = mata_means[(c_mata_miss[b]-1)]+ MeansC[(c_mata_miss[b])]- MeansC[(c_mata_miss[b])-1]
    } 
  }
  return(mata_means)
}

pttestf<- function(n1,n2,mn1,sd1,mn2,sd2) {
  pttest = pt((((mn1 - mn2) -0) /sqrt(sd1^2/n1+sd2^2/n2)),(n1+n2-2))
  return(pttest)
}

LMCF_loop <- function(c_mata_miss,mata_Means)
{
  miss_count <- length(c_mata_miss)
  mata_means <- mata_Means
  for (b in 1:miss_count)  {
    if (c_mata_miss[b] > 1) {
      mata_means[c_mata_miss[b]] <- mata_means[(c_mata_miss[b]-1)]      
    } 
  }
  return(mata_means)
}


#works
testread <-function(pathdat) {
  #mxdata <-read.csv("./asthma.csv")
  txdata <-read.csv(pathdat)
}
 txdata <-testread("asthma.csv")


