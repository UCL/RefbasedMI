utils::globalVariables(c("tst2", "Imp_Interims", "M", ".imp", "tmptreat", "methodvar", "referencevar", "patt", "treat", "classtreatvar"))
# to install if not already installed


## @import  mice


# if(!require(norm2)) install_url(‘https://cran.r-project.org/src/contrib/Archive/norm2/norm2_2.0.3.tar.gz’)
## @import norm2

# if(!require(data.table)) install.packages('data.table')
# library('data.table')

# for cholsolve
# if(!require(sparseinv)) install.packages('sparseinv')
# library(sparseinv)

#  for analysis (stat.desc)
# if(!require(pastecs)) install.packages('pastecs')
# library(pastecs)

#  for summarize
# if(!require(Hmisc)) install.packages('Hmisc' , repos = "http://cran.us.r-project.org")
# library(Hmisc)


# note pastecs and data.table both have first and last functions
# so avoid conflict by specifying specfic function in pastec



LMCF_loop <- function(c_mata_miss, mata_Means) {
  miss_count <- length(c_mata_miss)
  mata_means <- mata_Means
  for (b in 1:miss_count) {
    if (c_mata_miss[b] > 1) {
      mata_means[c_mata_miss[b]] <- mata_means[(c_mata_miss[b] - 1)]
    }
  }
  return(mata_means)
}




CIR_loop <- function(c_mata_miss, mata_Means, MeansC)
# c_mata_miss something like [2 3 4] ,so is cc
{
  # browser(text="1102")
  miss_count <- length(c_mata_miss)
  # this not right?
  # mata_means <- as.data.frame(mata_Means)
  mata_means <- (mata_Means)

  # diagnostics
  # cat(paste("\nc_mata_miss="))
  # print(c_mata_miss)
  # cat(paste("\nmata_means="))
  # print(mata_means)
  # cat(paste("\nMeansC="))
  # print(MeansC)
  # browser(text="2801")
  for (b in 1:miss_count) {

    # if 1st col missing then no value before so need to check for that
    # note mimix counter just no rows in each patt
    # so main looping is over missin fields, ie miss_count
    # if 1st col must be ref
    if (c_mata_miss[b] == 1) {
      # test purposes
      # count<-2
      # MeansC is list so use [[1]]
      mata_means[b] <- MeansC[[1]][b]
    } else {
      # filling in column at a time
      # mata_means causes error when indiv specific use so try unlist
      mata_means[c_mata_miss[b]] <- unlist(mata_means)[(c_mata_miss[b] - 1)] + MeansC[[1]][(c_mata_miss[b])] - MeansC[[1]][(c_mata_miss[b]) - 1]
    }
  }

  # diagnostics
  # cat(paste("\nafterloop mata_means="))
  # print(mata_means)

  return(mata_means)
}



Causal_loop <- function(c_mata_miss, mata_Means, MeansC, K0, K1) {


  miss_count <- length(c_mata_miss)
  mata_means <- as.data.frame(mata_Means)

  lastvisit <- min(c_mata_miss) - 1
  # note lastvisit cannot be zero,(though is time 0) so recode to 1
  lastvisit <- ifelse(lastvisit == 0, 1, lastvisit)
  for (b in 1:miss_count) {

    # if 1st col missing then no value before so need to check for that
    # note mimix counter just no rows in each patt
    # so main looping is over missin fields, ie miss_count
    # if 1st col must be refn
    if (c_mata_miss[b] == 1) {
      # MeansC is list so use [[1]]
      mata_means[b] <- MeansC[[1]][b]
    } else {
      # filling in column at a time

      # assuming  mata_means[t]- MeansC[[1]][t] is depature from overall mean , ie MyCov in sas macro
      # find time of discontinuation (ignoring interim cases?)

      # establish time of last visit ( asunming no interims for now!)
      # eg c_mata_miss = (2,3,4) shows missing cols
      # lastvisit is time t in Ian's paper  which is discontinuation at visit t , ie the last visit before all missing va;ues
      # this assumes missing values occur after lastvisit
      # but what happens when xxxo is the pattern?  ie c_mata_miss = 2,3,4 so baseline is lastvist but 5 is sa well?!
      # in asthma data largely montone so assume min true , BUT some wont be, eg X0XXO  245 has 2 lastvist values
      # noting the interims -  case will be interim if  final visit , ie at  tmax non-missing!

      # if interim instead try something like

      #  consider 2 4 5
      if (max(c_mata_miss) < length(mata_Means)) {
        # MAR so just use the treatment arm mata_means
        # lastvisit <-max(c_mata_miss)+1, ie interims, assign to ref man or Arm mean (just defaults ?)
        mata_means[c_mata_miss[b]] <- MeansC[[1]][c_mata_miss[b]]
      }


      # departure from overall mean at time t (lastvist) , active mean - ref mean at last visit
      #  ActRef_diff <-mata_means[lastvisit]-MeansC[[1]][lastvisit]

      # need to compere v CIR , J2R and doesnt use terms in  formula 7
      # this was eqn 5 soln   2/6/20

      # try eq 6 Above ok! so dont interfere
      # to implement eq 6 need calc current vist - lastvisit and use it to exponentiate
      # the unit visit  diferences are  c_mata_miss[b]-lastvisit

      # but specify k0,k1
      # mata_means[c_mata_miss[b]] <- MeansC[[1]][c_mata_miss[b]]+( K0*(K1^(c_mata_miss[b]-lastvisit))*ActRef_diff )

      # try Ians  formula
      # reference mean MeansC[[1]][c_mata_miss[b]] at time u
      # active mean time t  (t lastvist)  mata_means[lastvisit)]
      # refernce mean time t (t lasvist)   MeansC[[1]][lastvisit]

      v_u <- as.numeric(sub(".*\\.", "", colnames(mata_means[c_mata_miss]))[b])

      # test whether index lastvisit has a "*.number" format (ie whether base covariate)
      # if (grep('.*\\.','',colnames(mata_means[lastvisit]))==0L  )
      # the base covariate (or the final covariate if many)  is the lastvisit which means every value is missing ! , ie v_t is set to 0
      # 21/10
      #  if c_mata_miss = 2 4 5 then lastvisit has 2 values, ie 1 and 3 so need to check if 2nd last value appropriate
      #  interim missing if
      # define c_mata_miss[0] = 0 as should always be complete (but check knock-on effects as length(c_mata_miss) changed)
      # this doesnt work   c_mata_miss[0] = 0 so try making exception when b=1,


      # note following only has effect such as X0XX becasue only occurence for 2 lastvisits in asthma
      if (b > 1) {
        if ((c_mata_miss[b - 1] + 1) != c_mata_miss[b]) {
          lastvisit <- c_mata_miss[b] - 1
        }
      }


      if (suppressWarnings(is.na(as.numeric(sub(".*\\.", "", colnames(mata_means[lastvisit])))))) {
        v_t <- 0
      } else {
        v_t <- as.numeric(sub(".*\\.", "", colnames(mata_means[lastvisit])))
      }

      # outcome at visit u after discontinuation at visit t , tmax is maximum no. visits
      # Yt(R)  =  MeansC[[1]][lastvisit])  , ie reference group mean at last vist
      # so need calc E[Y(A)]t for treat grp A up to lastvisit t , then from t+1 to tmax is mean for ref group
      # so only issue is to  entwine mata_means from A up to and including lastvist, then MeansC[[1]] from lastvisit+1 to tmax
      # for 1st term on RHS depends,if missing before last visit then use mata_means, if missing after then MeansC[[1]]


      # but if lastvisit = tmax then no reference based treatment
      #   if (lastvisit == length(mata_Means)) {
      #     mata_means[c_mata_miss[b]] <- MeansC[[1]][c_mata_miss[b]]+K0*(K1^(v_u-v_t))*(mata_means[c_mata_miss[b]-1]-MeansC[[1]][c_mata_miss[b]]-1)

      # note brackets round (c_mata_miss[b])-1]
      # also assigmnr <- or =??

      # mata_means[c_mata_miss[b]] <- mata_means+MeansC[[1]][c_mata_miss[b]]+K0*(K1^(v_u-v_t))*(mata_means[lastvisit]-MeansC[[1]][lastvisit])


      mata_means[c_mata_miss[b]] <- MeansC[[1]][c_mata_miss[b]] + K0 * (K1^(v_u - v_t)) * (mata_means[lastvisit] - MeansC[[1]][lastvisit])
      # mata_means[c_mata_miss[b]] =  MeansC[[1]][c_mata_miss[b]] + K0*(K1^(v_u-v_t))*(mata_means[(c_mata_miss[b]-1)]-MeansC[[1]][(c_mata_miss[b])-1])
      # CIR  mata_means[c_mata_miss[b]] = mata_means[(c_mata_miss[b]-1)]+ MeansC[[1]][(c_mata_miss[b])]- MeansC[[1]][(c_mata_miss[b])-1]
    } # if not interim
  }

  return(mata_means)
}



# dont include rthis in the final package ,  for test purposes only!



analyselist <- function(id, datlist, varlist) {

  datano <- subset(datlist, id == datlist$.id)
  cat(paste0("\ncase = ", id))
  cat(paste0("\n treatarm = ", subset(datano$treat, datano$.imp == 0), "\n"))
  # numbers denote the descriptive stats to display
  t(round(pastecs::stat.desc(datano)[, varlist], 8)[c(1, 9, 13, 4, 8, 5), ])
}



AddDelta <- function(vec_tst, covar, mata_imp, delta, dlag) {
  # infinite loop ?

  # create vector of 1 and 0s
  # no space before .miss 12/5/20, stat at 3rd col skipping GI II
  # onezero<-sapply(vec_tst[3:(2+length(vec_tst)-ncovar)], function(x) return(mata_imp[1,paste0(x,".miss")]))

  # simpler way to eliminate covars, need preserve vec_tsts for later
  # vec_tst<- setdiff(vec_tst,covar)
  onezero <- sapply(setdiff(vec_tst, covar), function(x) {
    return(mata_imp[1, paste0(x, ".miss")])
  })

  # drop covariates (which should be complete anyway - should che ck!!)
  # order now is depvars followed by covars so drop last terms
  # hence this method  no longer valid
  # dropcovar <- ncovar+1
  # onezero<-onezero[c(dropcovar:length(onezero))]
  # onezero<-onezero[c(1:(length(onezero)-ncovar))]

  # then 1st and last  ,set last0 as last of complete  before the missing values start
  # but as to go in if because min(0,0,..) cause warnings
  # lastVisit <- min(which(unlist(onezero)==1))

  # check not interim , ie no gaps in onezero seq
  # in which case just leave wout delta adjustment but print warnnig msg
  # determine whether interim by noting if a zero appears after any one,
  # ie if 1st occurence of 1 to left of any 0

  # set lastvisit before discontinuation
  lastvisit <- max(which(unlist(onezero) == 0))
  # check whether any missing after this
  if (lastvisit > max(which(unlist(onezero) == 1))) {
    # if (max(which(unlist(onezero)==0)) >  min(which(unlist(onezero)==1))) {
    # do nothing as no discontinutions
  } else
  # do for discontinuations
  {

    # lastvisit is p in JamesRogers paper
    super_delta <- 0
    v <- 0

    # we only increment delta when missing, so skip if non-missing
    # 1st row should be same value for all rows in the same pattern group, gives warning othewise
    # needs be inserted here as previous passed by argument
    ncovar <- length(covar)
    # range of values to be changed in mata_imp given by the cols start:end
    # the covars are now after the depvars so wrong to include them
    # start<- 2+ncovar+1
    start <- 2 + 1
    # end <- 2+length(vec_tst)
    end <- 2 + length(setdiff(vec_tst, covar))


    # loop over selected range in mata_imp

    for (j in (start + lastvisit):end) {


      # adjust from last visit (ignore interims for now)
      # for (v in lastVisit:length(onezero)) {
      # count from start to end
      v <- v + 1
      super_delta <- super_delta + delta[j - start + 1] * dlag[v]
      mata_imp[j] <- mata_imp[j] + super_delta
    }

    # mata_new values start in col 3 so must add 2 to index
    # needs adjusting for when interim case non missing! test dummy missing vars
    # jump<-length(vec_tst)
    # mata_imp[2+v]<- ifelse(mata_imp[2+v+1+jump]==1,mata_imp[2+v] + super_delta,mata_imp[2+v])
    # {
    # print(paste0("interim missing, check delta adjustment ",onezero))
    # mata_imp[2+v] <- mata_imp[2+v] + super_delta
  } # if not interim


  # need report after final call
  if (sum(is.na(mata_imp)) != 0) {
    cat(paste0("\nWARNING!!! unimputed data values, possibly due to mis-specified delta"))
  }

  return(mata_imp)
}
