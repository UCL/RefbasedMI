

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
