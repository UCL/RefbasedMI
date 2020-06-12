
# to install if not already installed
#if(!require(norm2)) install.packages('norm2')
#if(!require(mice)) install.packages('mice')
#for select  function
#install.packages("dplyr")
#library(dplyr)

# for pivot_wider function
#install.packages("tidyr")
#library("tidyr")

# for emNorm
#install.packages("norm2")
#library(norm2)

# for cholsolve
#install.packages("sparseinv")
#library(sparseinv)

#  for analysis (stat.desc)
#install.packages("pastecs")
#library(pastecs)

# for Amelia
#install.packages("amelia")
#library(amelia)


#' @title LMCF_loop
#' @description process LMCF method
#' @details This is based on Suzie Cro's Stata program
#' @param c_mata_miss vector of col locaton of missing values , eg 5 6  
#' @param mata_Means vector of means after mcmc draws eg 17 1 16.8 15.5 14.6 13.2 
#' @return mata_means

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

#' @title CIR_loop
#' @description process CIR method
#' @details This is based on Suzie Cro's Stata program
#' @param c_mata_miss vector of col locaton of missing values , eg 5 6  
#' @param mata_Means vector of means after mcmc draws eg 17 1 16.8 15.5 14.6 13.2   
#' @param MeansC vector of means after mcmc draws using variance from reference group 
#' @return mata_means



CIR_loop <- function(c_mata_miss,mata_Means,MeansC)
  # mata_S_miss something like [2 3 4] ,so is cc
{
  # browser()
  miss_count <- length(c_mata_miss)
  mata_means <- as.data.frame(mata_Means)

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
      # MeansC is list so use [[1]]
      mata_means[b] <- MeansC[[1]][b]
    } else  {
      #filling in column at a time
      mata_means[c_mata_miss[b]] = mata_means[(c_mata_miss[b]-1)]+ MeansC[[1]][(c_mata_miss[b])]- MeansC[[1]][(c_mata_miss[b])-1]
    }
  }
  return(mata_means)
}




#' @title Causal_loop
#' @description process Causal method
#' @details This is based on "White,Royes,Best" paper
#' @param c_mata_miss vector of col locaton of missing values , eg 5 6  
#' @param mata_Means vector of means after mcmc draws eg 17 1 16.8 15.5 14.6 13.2   
#' @param MeansC vector of means after mcmc draws using variance from reference group 
#' @param K0 Causal constant for use with Causal method     
#' @param K1 exponential decaying Causal constant for use with Causal method  0<k1<1
#' @return mata_means


Causal_loop<- function(c_mata_miss,mata_Means,MeansC,K0,K1)
{

  miss_count <- length(c_mata_miss)
  mata_means <- as.data.frame(mata_Means)

  for (b in 1:miss_count)  {

    # if 1st col missing then no value before so need to check for that
    # note mimix counter just no rows in each patt
    # so main looping is over missin fields, ie miss_count
    # if 1st col must be refn
    if (c_mata_miss[b] ==1) {
      #print will cause objectto return
      #print(paste0(miss_count))
      #for each patt_count (from mimix_group)
      #test purposes
      #count<-2
      # MeansC is list so use [[1]]
      mata_means[b] <- MeansC[[1]][b]
    } else  {
      #filling in column at a time
     # mata_means[c_mata_miss[b]] = mata_means[(c_mata_miss[b]-1)]+ Kd*( MeansC[[1]][(c_mata_miss[b])]- MeansC[[1]][(c_mata_miss[b])-1])

    # assuming  mata_means[t]- MeansC[[1]][t] is depature from overall mean , ie MyCov in sas macro
    # find time of discontinuation (ignoring interim cases?)

      # establish time of last visit ( asunming no interims for now!)
      # eg c_mata_miss = (2,3,4) shows missing cols
      #lastvisit is time t in Ian's paper
      lastvisit <-min(c_mata_miss)-1
      # departure from overall mean at time t (lastvist) , active mean - ref mean at last visit
      ActRef_diff <-mata_means[lastvisit]-MeansC[[1]][lastvisit]

      # 1/5/20 need to compere v CIR , J2R and doesnt use terms in  formula 7
  # this was eqn 5 soln   2/6/20
  #mata_means[c_mata_miss[b]] <- MeansC[[1]][c_mata_miss[b]]+  ( Kd*ActRef_diff )
      #mata_means[c_mata_miss[b]] <-   *(mata_means[(c_mata_miss[b]-1)]- MeansC[[1]][(c_mata_miss[b])-1]) + (mata_means[(c_mata_miss[b])]-MeansC[[1]][(c_mata_miss[b])])

      # try eq 6 Above ok! so dont interfere
      # to implement eq 6 need calc current vist - lastvisit and use it to exponentiate 
      # the unit visit  diferences are  c_mata_miss[b]-lastvisit
  
      # 2/6/20 eqn7 but specify k0,k1
       mata_means[c_mata_miss[b]] <- MeansC[[1]][c_mata_miss[b]]+( K0*(K1^(c_mata_miss[b]-lastvisit))*ActRef_diff )
      
    }
  }
  return(mata_means)
}




#' @title regressimp
#' @description run regression on M imputed data set, combining as Rubin's rules
#' @details This is approach followed from  norm2 user manual
#' @export regressimp
#' @param dataf data-frame
#' @param regmodel regression model specfication
#' @return estimates of regression coefficients
#' @example
#' \dontrun{
#' regressimp(impdataset,"fev.12~treat+base")
#' }


regressimp <- function(dataf,regmodel)  {
  # to get the list
  implist1x <- split(dataf,dataf[,"II"])
  # so has M elements in list
  # can obtain a list of coefficients and their se's from a regression
  # declare list for estimates
  est.list <- as.list(NULL)
  # declare lists for se's
  std.err.list <- as.list( NULL )
  M<- utils::tail(dataf[,"II"],1)
  for( m in 1:M ){
    #mod<-lm(fev12~as.factor(treat)+base,data=kmlist1x[[m]] )
    #mod<-lm(head12~head_base+sex,data=implist1x[[m]] )
    #mod<-lm(HAMD17.TOTAL7~basval+HAMD17.TOTAL6,data=implist1x[[m]] )
    mod<-stats::lm(regmodel,data=implist1x[[m]] )
    est.list[[m]] <- stats::coef(summary(mod))[,1]
    std.err.list[[m]] <- stats::coef(summary(mod))[,2]
    }
  ## combine the results by rules of Barnard and Rubin (1999)
  ## with df.complete = 27, because a paired t-test in a dataset with
  ## N=28 cases has 27 degrees of freedom

  miResult <- norm2::miInference(est.list, std.err.list, df.complete=801)
  print(miResult)

  # trying mcerror
  #return(list(est.list,std.err.list))
}



#' @title analyselist
#' @description find descriptive stats on the  M imputed data set
#' @details select on patient id and find their means etc
#' @export analyselist
#' @param id patient identifier
#' @param datlist imputed dataset of M imputations
#' @param varlist list of derived variables ,varlist <- c("fev.2","fev.4","fev.8","fev.12","base")
#' @return printout of descriptve stats
#' @example
#' \dontrun{
#' varlist <- c("fev.2","fev.4","fev.8","fev.12","base")
#' analyselist(5099,impdataset,varlist)
#' }


analyselist <-function(id,datlist,varlist) {
  datano <- subset(datlist,id==datlist$.id)
   # numbers denote the descriptive stats to display
   t(round(pastecs::stat.desc(datano)[,varlist],3)[c(1,9,13,4,8,5),])
 }

#' @title AddDelta
#' @description add delta's to imputed values
#' @details adding delta values after wthdrawal
#' @param vec_tst  vector of visit names
#' @param ncovar number covariates
#' @param mata_imp the imputed values (as well as the complete) 
#' @param delta vector of delta values (a values in Roger's) for each vist time
#' @param dlag vector of dlag values (b values in Roger's) for each vist time
#' @return mata_imp the adjusted imputed vaues (and unadjusted non-missing)



#define function 19/04 to add delta's to imputed values
AddDelta<-function(vec_tst,ncovar,mata_imp,delta,dlag)  {
 #browser()
  # need read dlag in mimix argument, so temp here
  
  # create vector of 1 and 0s
  #browser()  no space before .miss 12/5/20, stat at 3rd col skipping GI II
  #onezero<-sapply(vec_tst[3:(2+length(vec_tst)-ncovar)], function(x) return(mata_imp[1,paste0(x,".miss")]))
  onezero<- sapply(vec_tst, function(x) return(mata_imp[1,paste0(x,".miss")]))

  #drop covariates (which should be complete anyway)
  dropcovar <- ncovar+1
  onezero<-onezero[c(dropcovar:length(onezero))]
  # then 1st and last  ,set last0 as last of complete  before the missing values start
  # but as to go in if because min(0,0,..) cause warnings
  #lastVisit <- min(which(unlist(onezero)==1))

  # check not interim , ie no gaps in onezero seq
  # in which case just leave wout delta adjustment but print warnnig msg
  # determine whether interim by noting if a zero appears after any one, 
  # ie if 1st occurence of 1 to left of any 0  
  if (max(which(unlist(onezero)==0)) >  min(which(unlist(onezero)==1))) {
    # do nothing as interim 
  } 
  else
  { 
  
  # if length = max then no gaps
  # if all 0's the no missing so no adjustment required
  #if (sum(which(unlist(onezero)!=0) ) ) {
    # check no gaps, ie interims, not quite right!
    #if ( length(which(unlist(onezero)==1)) == max(which(unlist(onezero)==1)) ){
    #  lastVisit <- min(which(unlist(onezero)==1)) }
  lastvisit <- min(which(unlist(onezero)==1))-1
      # so add appropriate delta to imputed values after last visit
  
  #lastvisit is p in JamesRogers paper 
   super_delta<-0
   v<-0

   
   
  # for (v in lastVisit:length(onezero)) {
    # when dlag used super_delta<- super_delta + delta[v]* dlag[v-1]
    # we only increment delta when missing, so skip if non imssing
    # 1st row should be same value for all rows in the same pattern group, gives warning othewise

    # range of values to be changed in mata_imp given by the cols start:end
    start<- 2+ncovar+1
    end <- 2+length(vec_tst)
  #mata_imp[,start:end]
    #loop ove selected range in mata_imp

  for (j in (start+lastvisit):end) {

    #4:8
    #if (mata_imp[1,start+v]==1) {
      # adjust from last visit (ignore interims for now)
      #for (v in lastVisit:length(onezero)) {
      # count from start to end
      v<- v+1
       super_delta <- super_delta + delta[j-start+1]*dlag[v]
       mata_imp[j] <- mata_imp[j] + super_delta
  }
       #browser()
    # mata_new values start in col 3 so must add 2 to index
    # needs adjusting for when interim case non missing! test dummy missing vars
   # jump<-length(vec_tst)
   # mata_imp[2+v]<- ifelse(mata_imp[2+v+1+jump]==1,mata_imp[2+v] + super_delta,mata_imp[2+v])
   #{
    #print(paste0("interim missing, check delta adjustment ",onezero))
   # mata_imp[2+v] <- mata_imp[2+v] + super_delta
    
} 
  #if not interim
    return(mata_imp)
}


