# program to read in data sets to Rmimix as in Stata


# from asthma data or other raw data wit repeated measures, 

# need to get to this 
preimputemvn1 <- read_dta("preimputemvn1.dta")


# create function for the 1st bit
#mimix depvar treatvar, id(varname) time(varname) [options]

# at present prepare data outside the function
library(tidyr)
library(dplyr)
library(plyr)
#read data
# need haven
library(haven)
mxdata <- read_dta("asthma.dta")  
#subset wont work in a user defined function ,subset needs to be evakuated with envir = data  so give up for now
# trying to write a subroutine for this part of data preprocessing


# now put all the functions together, possibly just run consecutively? 26/10/19
# this section to convert from long to wide data format
preprodata<- function(depvar,treatvar,idvar,timevar,covar)  {
  #extract relevant vars
  fevdata<-select(mxdata,idvar,depvar,timevar)
  # only want to widen the dose var, ie fev, so take the other variables
  #subdata<-select(mxdata,idvar,depvar,timevar)
  #reshape to wide and assign new names 
  # prefix must be supplied from input argument rather than hard coded, this canbe hard coded
  sts4<-pivot_wider(fevdata,id_cols=c(idvar),names_from=timevar,names_prefix="fev",values_from=depvar)
  #sts4 is just response data, can join later treat and covar cols from finaldat 
  #assumes no NA for these vars should add in checking routines !!
  uniqdat<-unique(mxdata[c(idvar,treatvar,covar)])
  
 #uniqdat<- uniqdat[order(treatvar),]
  #uniq to merge onto the wide data set sql joinon idjoin R
  finaldat<- merge(sts4,uniqdat,by=idvar)
  # now try and sort on treatvar, doesnt work so instead of sorting , select on treat
  # and rbind together
  #finaldat<- finaldat[order(treatvar),]
  
  # need find no. ntreat  to loop over
  #ntreat<-unique(unique(mxdata$treatvar))
  ntreatcol<-(select(mxdata,treatvar))
  ntreat <- unique(ntreatcol)
  ntimecol<-(select(mxdata,timevar))
  ntime<-unique(ntimecol)
  
       
  
  
  #in order to aggregate by pattern need create dummy vars
  #depvar <-  ifelse(!is.na(data$fev1),0,1)  
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
  # testing return(finaldat3)
  # then sort by patt( last col) and treat and find cumulative sequence  AND split on treatment var!
  # so need to merge treat (and covar) back into response data 
  
  #sort by treatvar and patt var, couldnt get treatvar working with patt so try sorting previously on treatvar   
  finaldatSS <-sts4Dpatt[order(patt),] 
 #not work
  #finaldatS <-finaldatS[with(finaldatS,order(treatvar)),] 
  ndx = order(finaldatSS[,treatvar])
  finaldatS <- finaldatSS[ndx,]
  
  # consistent with Stata to get right orderdrop id and treat cols
  drops <-c(idvar,treatvar)
  finaldatS<-finaldatS[,!(names(mata_raw) %in% drops)] 
  mata_Obs[order(mata_Obs$treat,mata_Obs$patt),]
  #finaldat3[with(finaldat3,order(treatvar)),]
  #return(finaldatS)
  
  #also get the dummy pattern vectors by treatment
  dummypatt<-unique(xSTSdummy)
  pattmat<-unique(STSdummy[,1:ncol(STSdummy)])
  #couldnt get this working
  #ival=0
  #for (val in (ntime)) {
  #  kmvar=get(paste0("pattmat[,\"fev",val,"\"]"))
   # patti <- patti+kmvar*2^ival
  #  ival=ival+1
  #}  
  # sort of hard coded for now will have to revisit this!
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
 
  # would lke a look up table for patt and missing vectors
  # make sure ids
  #test[[1]]$id <- as.numeric(test[[1]]$id)
  #test[[2]]$id <- as.numeric(test[[2]]$id)
  #finaldat2<-left_join(test[[1]],test[[2]], by=idvar)
  #finaldat2<-join(sts4,finaldat, by=idvar)
  # sts4Dpows2[order(sts4Dpows2$treat,] 
  
  for (val in t(ntreat)) {
    print(paste0("prenormMC",val))
    #assign(paste0("prenorm2",val),subset(finaldat,treat==val))
  } 
  
  return(list(sts4Dpatt,finaldatS,finaldat,ntreat,ex1s,ex1,ex1id,pattmat,patt,ntime))
}
# Question is should the wide data be sorted for analysis or just ordered in fact, selecting out treatmet and rbinding? to provide analysis data set.   
# answer is should be sorted by patt as in Stata!
# therefore just select on treat when looping thru ntreat 
testb<-preprodata()
test<-preprodata("fev","treat","id","time","base")

# vital to get the mata_obs correctly sorted! so corresponds with mimix_group lookup 
# to be consistent with Stata move the base col after the fevs!
mata_Obs <- test[[2]]
mata_OBS <- mata_Obs[order(mata_Obs$treat,mata_Obs$patt),]  

# consistent with Stata is fev1-4 ,base so drop id (,1) and treat cols

#So test[[5]] is the lookup table equiv to !


test[[1]]$id <- as.numeric(test[[1]]$id)
finaldat2<-join(test[[1]],test[[2]], by=id)
test<-preprodata("fev","treat","id","time","base")
head(test[[2]])
ntreat <-test[[2]]

# next is to find the mimix_group subset over which can loop and which provides the missing data pattern 
# now need find totals in each patt category as well as cum position  to act as counter later in program
# need to aggregate , collapse on treat and patt!
# firstly for each depvar create dummy vars representing missing/not


A <- function(x) {  ifelse(!is.na(x),0,1) }
sts4<-(test[[2]])
# col 3 if treat var exists
STSdummy<-apply(sts4[,3:ncol(sts4)],MARGIN=2,FUN=A)
# append to names
colnames(STSdummy) <- paste(colnames(STSdummy),'.1')
# merge back on to data 
sts4D<-cbind(sts4,STSdummy)
#create the patt var by powers of 2 ,easiest on STSdummy

#outer(row(STSdummy),col(STSdummy),FUN=function(x,c) x*2^(c-1) )
#B <- function(x) {}
#tryagg$patttest <- tryagg$Group.2+tryagg$Group.3*2^1+tryagg$Group.4*2^2+tryagg$Group.5*2^3 

# create powrs of 2
pows2 <- sapply(2:ncol(STSdummy),function(i) STSdummy[,i]*2^(i-2))
#need to add up to find patt
patt <- rowSums(pows2)
sts4Dpatt<-cbind(sts4D,patt)
#assign name to pattern var
#names(sts4Dpows2)[ncol(sts4Dpows2)]<-"patt"

# then sort by patt( last col) and treat and find cumulative sequence  AND split on treatment var! 
# 
sts4Dpows2[order(sts4Dpows2$treat,]
sts4Dpatt[with(sts4Dpows2,order(treat,patt)),]
# this best for  count
ex1 <- sts4Dpatt %>%
       group_by(sts4Dpatt$treat,sts4Dpatt$patt) %>%
         summarise(X1 =n())
ex1 <- rename(ex1,treat='sts4Dpatt$treat',patt='sts4Dpatt$patt')
# and now find cumX1
ex1$X1cum <- cumsum(ex1$X1) 


# now need find totals in each patt category as well as cum position  to act as counter later in program
# need to aggregate , collapse on treat and patt!
#tryagg <- aggregate(sort.datafev,by =list(sort.datafev$treat,sort.datafev$fev1,sort.datafev$fev2,sort.datafev$fev3,sort.datafev$fev4),FUN="sum")
#co
#mimix_group<-agregate[sts4Dpatt,by = list(treat,patt),FUN="sum")
#mimix_group<-summarise[sts4Dpatt,by = list(treat,patt),FUN=n())

test<-preprodata("fev","treat","id","time","base")
head(test[[2]])
ntreat <-test[[2]]
# so can loop over ntreat values
# create the ntreat data sets

for (val in ntreat) {
 # 
     dat_treat <- subset(test[[1]], treat==val,)
}

# no need for this as put in above
testfun <- function(treatvar) {
  #ntreat<-unique(unique(mxdata$treatvar))
  ntreatcol<-(select(mxdata,treatvar))
  ntreat <- unique(ntreatcol)
  return(ntreat)
}

# need find no. ntreat  to loop over
# create ntreat var

#head(mxdata$treat)
uniqfun<- function(treatvar)
ntreatx<-unique(unique(mxdata$treat))
for (val in t(ntreat)) {
  print(paste0("prenormMC",val))
  assign(paste0("prenorm2",val),subset(uniq_fevdata,treat==val)) 
}




}

# works
tstfev<-dcast(melt(fevmxdata,id.vars = c("id","time")),id~variable+time)


tryfun3 <- function(idvar,timevar) {
  sts3<-pivot_wider(fevmxdata,id_cols=c(idvar),names_from=timevar,names_prefix="fev",values_from=timevar)
  #tstfev3<-dcast(melt(fevmxdata,id.vars = c("idvar","timevar")),"idvar"~variable+"timevar")
}
  
tryfun1 <- function(depvar,treatvar,idvar,timevar,covar ) {
  #extract relevant vars for wide transform
  fevmxdata <- subset(mxdata, select = c(idvar,depvar,timevar))
  #reshape to wide
  sts3<-pivot_wider(fevmxdata,id_cols=c(idvar),names_from=timevar,names_prefix="fev",values_from=timevar)
 # fun2dat<-select(mxdata,idvar,depvar,timevar)
  # only want to widen the dose var, ie fev, so take the other variables
  submxdata <- subset(mxdata, select = c(idvar,treatvar,covar))
  #take rows with unique id out
  uniq<-unique(submxdata[c(idvar,treatvar,covar)])
  uniq_fevdata<-join(uniq,sts3,by=idvar)
  
  return(list(uniq_fevdata))
}

# note need to call as tryfun2("fev","id","time" )  
tryfun2 <- function(depvar,idvar,timevar) {
   fun2dat<-select(mxdata,idvar,depvar,timevar)
  return(fun2dat)
}
fun2dat <-tryfun2("fev","id","time")




normalize <- function(avec,annum) {
    norvec <- (avec/annum)^0.5
    return(norvec) 
}
normalvec<-normalize(50,5)


tryfun <- 
{
  #extract relevant vars
  #fevmxdata <- subset(mxdata, select = c(idvar,depvar,timevar))
  #fevmxdata <- c(mxdata$idvar,mxdata$depvar,mxdata$timevar)
  fevmxdata <- c(mxdata[idvar,depvar,timevar])
}


# need find no. ntreat  to loop over

# create ntreat var

#head(mxdata$treat)
ntreat<-unique(unique(mxdata$treat))
for (val in ntreat) {
      print(paste0("prenormMC",val))
  assign(paste0("prenorm2",val),subset(uniq_fevdata,treat==val)) 
}

# note return makes function element global!



# then loop thru treatments selectng each treatment subset but also have to sort by the data pattern!  




# ie a data frame with patient id and covariaes selected
#each row corresponding to data for each individual.

# stata just uses reshape
# mimix fev treat, id(id) time(time) method(j2r) refgroup(2) covariates(base) clear m(1000)  seed(101)

#syntax varlist, time(varname) id(varname) [ COVariates(varlist) methodvar(varname) METHOD(string) refgroupvar(varname) REFgroup(str) SAVing(string) clear m(integer 5) BURNin(integer 100) BURNBetween(integer 100) seed(integer 0) interim(string) iref(str) mixed regress]


# read asthma data set


#suggest 

reshapeData <- function(data, var.col, id.col="id")  {
    
} 

library(reshape2)
tst<-dcast(melt(mxdata,id.vars = c("id","time")),id~variable+time)

library(magrittr)
requireNamespace("tidyr")
mxdata %>% pivot_wider( names_from=c(time,values_from=fev))

# only want to widen the dose var, ie fev, so take the other variables
submxdata <- subset(mxdata, select = c(id,treat,base))
# need check no missing values in covariates and treatment
fevmxdata <- subset(mxdata, select = c(id,fev,time))
tstfev<-dcast(melt(fevmxdata,id.vars = c("id","time")),id~variable+time)

#must be more efficient to do quick  sub analysis on treat and base vars to make sure all there
#take rows with unique id out 
dis<-distinct(submxdata, id,treat,base)
uniq<-unique(submxdata[c("id","treat","base")])
#uniq to merge onto the wide data set sql joinon idjoin R

uniq_fevdata<-join(uniq,tstfev,by="id")


    
