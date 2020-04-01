

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
#' @param c_mata_miss vector of missing
#' @param mata_Means fill-in here
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
#' @param c_mata_miss vector of missing
#' @param mata_Means fill-in here
#' @param MeansC fill-in
#' @return mata_means


# mata_S_miss something like [2 3 4] ,so is cc
CIR_loop <- function(c_mata_miss,mata_Means,MeansC)
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

#' @title pttestf
#' @description significance test between mean of 2 groups
#' @details This is based on ..
#' @param n1 sample size  group1
#' @param n2 samplesize group
#' @param mn1 mean value in sample 1
#' @param mn2 mean value in sample 2
#' @param sd1 standard deviation in sample 1
#' @param sd2 stadard deviation sample 2
#' @return p value
#' @examples
#' pttestf(200,200,1.302,0.623,1.367,0.573)

#'
pttestf<- function(n1,n2,mn1,sd1,mn2,sd2) {
  pttest = pt((((mn1 - mn2) -0) /sqrt(sd1^2/n1+sd2^2/n2)),(n1+n2-2))
  return(pttest)
}


#' @title regressimp
#' @description run regression on M imputed data set, combinging as Rubin's rules
#' @details This is approach followed from  norm2 user manual
#' @param dataf data-frame
#' @param regmodel regression model specfication
#' @return estimates of regression coefficients
#' @example
#' regressimp(impdataset,"fev12~treat+base")


regressimp <- function(dataf,regmodel)  {
  # to get the list
  implist1x <- split(dataf,dataf[,"II"])
  # so has M elements in list
  # can obtain a list of coefficients and their se's from a regression
  # declare list for estimates
  est.list <- as.list(NULL)
  # declare lists for se's
  std.err.list <- as.list( NULL )
  M<- tail(dataf[,"II"],1)
  for( m in 1:M ){
    #mod<-lm(fev12~as.factor(treat)+base,data=kmlist1x[[m]] )
    #mod<-lm(head12~head_base+sex,data=implist1x[[m]] )
    #mod<-lm(HAMD17.TOTAL7~basval+HAMD17.TOTAL6,data=implist1x[[m]] )
    mod<-lm(regmodel,data=implist1x[[m]] )
    est.list[[m]] <- coef(summary(mod))[,1]
    std.err.list[[m]] <- coef(summary(mod))[,2] }
  ## combine the results by rules of Barnard and Rubin (1999)
  ## with df.complete = 27, because a paired t-test in a dataset with
  ## N=28 cases has 27 degrees of freedom
  miResult <- miInference(est.list, std.err.list, df.complete=801)
  print(miResult)
}
