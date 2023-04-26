# RefBasedMI: utility functions

utils::globalVariables(c("tst2", "Imp_Interims", "M", ".imp", "tmptreat", "methodvar", "referencevar", "patt", "treat", "classtreatvar"))

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




CIR_loop <- function(c_mata_miss, mata_Means, MeansC) {
  miss_count <- length(c_mata_miss)
  mata_means <- (mata_Means)

  for (b in 1:miss_count) {
    # if 1st col missing then no value before so need to check for that
    # note mimix counter just no rows in each patt
    # so main looping is over missin fields, ie miss_count
    # if 1st col must be ref
    if (c_mata_miss[b] == 1) {
      # MeansC is list so use [[1]]
      mata_means[b] <- MeansC[[1]][b]
    } 
    else {
      # filling in column at a time
      mata_means[c_mata_miss[b]] <- unlist(mata_means)[(c_mata_miss[b] - 1)] + MeansC[[1]][(c_mata_miss[b])] - MeansC[[1]][(c_mata_miss[b]) - 1]
    }
  }
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
    } 
    else {
      # filling in column at a time
      if (max(c_mata_miss) < length(mata_Means)) {
        # MAR so just use the treatment arm mata_means
        mata_means[c_mata_miss[b]] <- MeansC[[1]][c_mata_miss[b]]
      }
      v_u <- as.numeric(sub(".*\\.", "", colnames(mata_means[c_mata_miss]))[b])

      # note following only has effect such as X0XX becasue only occurence for 2 lastvisits in asthma
      if (b > 1) {
        if ((c_mata_miss[b - 1] + 1) != c_mata_miss[b]) {
          lastvisit <- c_mata_miss[b] - 1
        }
      }
      if (suppressWarnings(is.na(as.numeric(sub(".*\\.", "", colnames(mata_means[lastvisit])))))) {
        v_t <- 0
      } 
      else {
        v_t <- as.numeric(sub(".*\\.", "", colnames(mata_means[lastvisit])))
      }
      mata_means[c_mata_miss[b]] <- MeansC[[1]][c_mata_miss[b]] + K0 * (K1^(v_u - v_t)) * (mata_means[lastvisit] - MeansC[[1]][lastvisit])
    } # if not interim
  }

  return(mata_means)
}



analyselist <- function(id, datlist, varlist) {
  datano <- subset(datlist, id == datlist$.id)
  message(paste0("\ncase = ", id))
  message(paste0("\n treatarm = ", subset(datano$treat, datano$.imp == 0), "\n"))
  # numbers denote the descriptive stats to display
  t(round(pastecs::stat.desc(datano)[, varlist], 8)[c(1, 9, 13, 4, 8, 5), ])
}



AddDelta <- function(vec_tst, covar, mata_imp, delta, dlag) {
 
  onezero <- sapply(setdiff(vec_tst, covar), function(x) {
    return(mata_imp[1, paste0(x, ".miss")])
  })
  
  # set lastvisit before discontinuation
  lastvisit <- max(which(unlist(onezero) == 0))
  # check whether any missing after this
  if (lastvisit > max(which(unlist(onezero) == 1))) {
    # do nothing as no discontinutions
  } 
  else {
    # do for discontinuations
    # lastvisit is p in JamesRogers paper
    super_delta <- 0
    v <- 0

    # we only increment delta when missing, so skip if non-missing
    # 1st row should be same value for all rows in the same pattern group, gives warning othewise
    # needs be inserted here as previous passed by argument
    ncovar <- length(covar)
    # range of values to be changed in mata_imp given by the cols start:end
    # the covars are now after the depvars so wrong to include them
    start <- 2 + 1
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
  } # if not interim


  # need report after final call
  if (sum(is.na(mata_imp)) != 0) {
    warning(paste0("\nUnimputed data values, possibly due to mis-specified delta"))
  }

  return(mata_imp)
}
