# RefBasedMI: functions to pre-process data

# Function to print data using message
print_in_message <- function(x)
{
  message(paste(capture.output(print(x)), collapse = "\n"))
}

# preprocess data for group method
preprodata <-
  function(data,
           covar,
           depvar,
           treatvar,
           tmptreat,
           idvar,
           timevar,
           M,
           reference,
           method = NULL, initial_levels_treat) {

    # change order put covar last

    fevdata <- get("data")[c(idvar, depvar, timevar, treatvar, covar)]
    # extract covar cols 1 row per id to merge onto the wide data
    uniqdat <- unique(get("data")[c(idvar, covar, treatvar)])
    ntreatcol <- get("data")[c(treatvar)]
    ntimecol <- get("data")[c(timevar)]


    # reshape from long to wide longitudinal data with as many depvars as time points

    # no covar
    fevdata <- get("data")[c(idvar, depvar, timevar, treatvar)]
    uniqcovar <- unique(get("data")[c(idvar, covar)])

    sts4 <-
      stats::reshape(
        as.data.frame(fevdata),
        v.names = depvar,
        timevar = timevar,
        idvar = idvar,
        direction = "wide"
      )

    # assumes no NA for these vars should add in checking routines !!
    # check how many covars used


    finaldat <- merge(sts4, uniqcovar, by = idvar)


    # need find no. ntreat  to loop over
    ntreat <- unique(ntreatcol)
    ntime <- unique(ntimecol)


    # in order to aggregate by pattern need create dummy vars
    # for sts4 array drop 1st col and create newvars for rest

    testfevdata <- get("data")[c(idvar, depvar, timevar)]
    sts4dummy <-
      stats::reshape(
        as.data.frame(testfevdata),
        v.names = depvar,
        timevar = timevar,
        idvar = idvar,
        direction = "wide"
      )
    STSdummy <-
      apply(sts4dummy[, grepl(depvar, names(sts4dummy))], MARGIN = 2, function(x)
        ifelse(!is.na(x), 0, 1))

    # append to names  try paste0 as otherwise space before miss
    colnames(STSdummy) <- paste0(colnames(STSdummy), ".miss")
    # merge back on to data
    sts4D <- cbind(finaldat, STSdummy)
    # create powrs of 2
    pows2 <-
      sapply(1:ncol(STSdummy), function(i)
        STSdummy[, i] * 2 ^ (i - 1))
    # need to add up to find patt
    patt <- rowSums(pows2)
    # sts4Dpatt<-cbind(sts4D,patt)
    sts4Dpatt <- cbind(sts4D, patt)

    # but also add on covars (should be no missing so dont contribute to patt)
    # ie should be  cols of 0's

    # have to cope when no covars specified
    if (length(covar) != 0) {
      tmp_covpatt <-
        apply(as.data.frame(sts4Dpatt[, covar]), MARGIN = 2, function(x)
          ifelse(!is.na(x), 0, 1))
      # add names
      colnames(tmp_covpatt) <- paste0(c(covar), ".miss")
      # than combine below the dummies onto the finaldat
      #  move covar to  col following depvar
      sts4Dpatt <- cbind(finaldat, STSdummy, tmp_covpatt, patt)
    } 
    else {
      sts4Dpatt <- cbind(finaldat, STSdummy, patt)
    }


    # In order  to get patt with all missing patts

    Overall_patt <-
      unique(sts4Dpatt[grepl(".miss", colnames(sts4Dpatt))])
    patt <- unique(sts4Dpatt[, "patt"])

    # then combine depvar and covar patt!
    all_patt <- cbind(Overall_patt, patt)


    # then sort by patt( last col) and treat and find cumulative sequence  AND split on treatment var!
    # so need to merge treat (and covar) back into response data

    # sort by treatvar and patt var, couldnt get treatvar working with patt so try sorting previously on treatvar
    finaldatSS <- sts4Dpatt[order(patt),]

    ndx <- order(finaldatSS[, treatvar])
    finaldatS <- finaldatSS[ndx,]


    # need sort by patt and align with mg (lookup) ,to get lookup right merge using patt
    finaldatS <-
      sts4Dpatt[order(sts4Dpatt[, treatvar], sts4Dpatt[, "patt"]),]

    # consistent with Stata to get right order, drop id and treat cols
    drops <- c(idvar, treatvar)


    # also get the dummy pattern vectors by treatment
    pattmat <- unique(STSdummy[, 1:ncol(STSdummy)])


    pows2d <-
      sapply(1:ncol(pattmat), function(i)
        pattmat[, i] * 2 ^ (i - 1))
    patt <- rowSums(pows2d)
    # want  join this to ex from mimix_group table
    pattmatd <- cbind(pattmat, patt)


    # and now find cumX1 cumulative no. cases in each pattern/treatment group

    sts4Dpatt$X1 <- 1

    ex1 <-
      Hmisc::summarize(sts4Dpatt$X1,
                       by = Hmisc::llist(sts4Dpatt[, treatvar], sts4Dpatt$patt),
                       FUN = sum)
    # to rename
    newnames <- c(treatvar, "patt", "X1")
    names(ex1) <- newnames
    # and now find cumX1
    ex1$X1cum <- cumsum(ex1$X1)

    ex1$X1cum <- cumsum(ex1$X1)
    # want  join this to ex from mimix_group table
    # create index number to preserve ordr after merge
    ex1$exid <- 1:nrow(ex1)
    ex1id <- merge(ex1, pattmatd, by = "patt")
    ex1s <- ex1id[order(ex1id$exid),]
    names(ex1)[names(ex1) == "X1"] <- "cases"
    names(ex1)[names(ex1) == "X1cum"] <- "cumcases"
    test_ex1 <-
      merge(ex1, all_patt, by = "patt")[order(merge(ex1, all_patt, by = "patt")$exid),]


    stopifnot(reference %in% t(ntreat))
    # rename to more user-friendly
    names(ex1s)[names(ex1s) == "X1"] <- "patients"

    # dont want to display so
    ex1s$X1cum <- NULL
    ex1s$exid <- NULL
    # prefer to call pattern
    names(ex1s)[names(ex1s) == "patt"] <- "pattern"

    message(paste0("   ", "\n\nSummary of missing data pattern by ", treatvar, ":\n\n"))


    # setting row names NULL automatically produces sequential index

    rownames(ex1s) <- NULL
    # put labels to put back original treat  levels (when not orig 1,2..)
    ex1s[, treatvar] <-
      ordered(ex1s[, treatvar], labels = initial_levels_treat)

    # tmptreat is factor which could  lead to wrong order so convert to numeric
    ex1s[, treatvar] <-
      sort(as.numeric(as.character(ex1s[, treatvar])))
    # above replaces below

    print_in_message(ex1s)


    # so finaldatS is sorted by treat,patt  data
    # test_ex1 is mg lookup table for missing dummies, all_patt is missing pattern table
    return(list(
      finaldatS,
      ntreat,
      test_ex1,
      all_patt,
      ntime,
      M,
      reference,
      method
    ))
  }


# preprocess data for individual method
preproIndivdata <-
  function(data,
           covar,
           depvar,
           treatvar,
           idvar,
           timevar,
           M,
           reference = NULL,
           method = NULL,
           methodvar,
           referencevar) {

    # check covars complete
    stopifnot(sum(is.na(get("data")[, covar])) == 0)

    methodL <- unique(get("data")[, methodvar])

    # change order for covar didnt work so take out treatvar covar and use uniqdat
    fevdata <-
      get("data")[c(idvar, depvar, timevar, methodvar, referencevar)]
    # now covar added to data list so need need for unique?
    uniqdat <- unique(get("data")[c(idvar, covar, treatvar)])
    ntreatcol <- get("data")[c(treatvar)]
    ntimecol <- get("data")[c(timevar)]



    # generate names depvar#time
    sts4 <-
      stats::reshape(
        fevdata,
        v.names = depvar,
        timevar = timevar,
        idvar = idvar,
        direction = "wide"
      )

    # assumes the covars all non-missing

    # merge on the covariates and treatment to names
    finaldat <- merge(sts4, uniqdat, by = idvar)

    # in order to aggregate by pattern need create dummy vars


    testfevdata <- get("data")[c(idvar, depvar, timevar)]
    sts4dummy <-
      stats::reshape(
        testfevdata,
        v.names = depvar,
        timevar = timevar,
        idvar = idvar,
        direction = "wide"
      )
    STSdummy <-
      apply(sts4dummy[, grepl(depvar, names(sts4dummy))], MARGIN = 2, function(x)
        ifelse(!is.na(x), 0, 1))

    # append to names,
    colnames(STSdummy) <- paste0(colnames(STSdummy), ".miss")

    # merge back on to data  ,ie fevdata then testfevdata
    sts4D <- cbind(finaldat, STSdummy)

    # create powrs of 2


    pows2 <-
      sapply(1:ncol(STSdummy), function(i)
        STSdummy[, i] * 2 ^ (i - 1))
    # need to add up to find patt
    patt <- rowSums(pows2)
    sts4Dpatt <- cbind(sts4D, patt)

    # if covars specified...
    if (length(covar) != 0) {
      tmp_covpatt <-
        apply(as.data.frame(sts4Dpatt[, covar]), MARGIN = 2, function(x)
          ifelse(!is.na(x), 0, 1))
      # add names
      colnames(tmp_covpatt) <- paste0(c(covar), ".miss")
      #  move covar to  col following depvar
      sts4Dpatt <- cbind(finaldat, STSdummy, tmp_covpatt, patt)
    } 
    else {
      sts4Dpatt <- cbind(finaldat, STSdummy, patt)
    }

    # In order  to get patt with all missing patts

    Overall_patt <-
      unique(sts4Dpatt[grepl(".miss", colnames(sts4Dpatt))])
    patt <- unique(sts4Dpatt[, "patt"])

    # then combine depvar and covar patt!
    all_patt <- cbind(Overall_patt, patt)

    # patt has just been created so use $ notation wheras treatvar from input argument and need to by methodindiv[1] as well

    # has to be also sorted into methodvar 1 and 2 grps
    finaldatS <-
      sts4Dpatt[order(sts4Dpatt[, treatvar], sts4Dpatt[, methodvar], sts4Dpatt[, referencevar], sts4Dpatt[, "patt"]),]

    # consistent with Stata to get right order, drop id and treat cols
    drops <- c(idvar, treatvar)


    # also get the dummy pattern vectors by treatment

    pattmat <- unique(STSdummy[, 1:ncol(STSdummy)])


    pows2d <-
      sapply(1:ncol(pattmat), function(i)
        pattmat[, i] * 2 ^ (i - 1))
    patt <- rowSums(pows2d)
    # want  join this to ex from mimix_group table

    pattmatd <- cbind(pattmat, patt)


    # X1 is a count variable of 1's'
    sts4Dpatt$X1 <- 1

    ex1 <-
      Hmisc::summarize(
        sts4Dpatt$X1,
        by = Hmisc::llist(sts4Dpatt[, treatvar], sts4Dpatt[, methodvar], sts4Dpatt[, referencevar], sts4Dpatt$patt),
        FUN = sum
      )
    # to rename

    newnames <- c(treatvar, methodvar, referencevar, "patt", "X1")
    names(ex1) <- newnames
    # and now find cumX1
    ex1$X1cum <- cumsum(ex1$X1)

    names(ex1)[names(ex1) == "X1"] <- "cases"
    names(ex1)[names(ex1) == "X1cum"] <- "cumcases"
    # create index number to preserve ordr after merge  , #want  join this to ex from mimix_group table
    ex1$exid <- 1:nrow(ex1)
    ex1id <- merge(ex1, pattmatd, by = "patt")
    ex1s <- ex1id[order(ex1id$exid),]

    rownames(ex1s) <- NULL
    names(ex1s)[names(ex1s) == "X1"] <- "patients"
    # dont want to display cumcases

    ex1s$exid <- NULL

    message(paste0("Summary missing pattern:\n\n"))


    # changing to patients has knock on effect in mg so leave as cases

    # dont want to display so
    ex1s$cumcases <- NULL
    ex1s$exid <- NULL
    # prefer to call pattern
    names(ex1s)[names(ex1s) == "patt"] <- "pattern"

    # prints out the summay table
    print_in_message(ex1s)
    # dont need finaldat,exlid

    # pattmat is just the missing dummies
    # patt vector of patterns need find no. ntreat  to loop over

    ntreatcol <- get("data")[c(treatvar)]
    ntreat <- unique(ntreatcol)
    ntimecol <- get("data")[c(timevar)]
    ntime <- unique(ntimecol)

    test_ex1 <-
      merge(ex1, all_patt, by = "patt")[order(merge(ex1, all_patt, by = "patt")$exid),]

    # error chk

    # find unique values for referencevar to check against ntreat values
    refencevars <- unique(get("data")[, referencevar])
    stopifnot(refencevars %in% t(ntreat))

    # main outputs have to be test_ex1 and finaldatS which should correspond
    return(
      list(
        finaldatS,
        ntreat,
        test_ex1,
        all_patt,
        pattmat,
        patt,
        ntime,
        M,
        methodvar,
        referencevar
      )
    )
}



# alternative logic for individual method
ifmethodindiv <-
  function(methodvar,
           referencevar,
           mg,
           m,
           M,
           paramBiglist,
           i,
           treatvar,
           c_mata_nonmiss,
           c_mata_miss,
           mata_miss,
           mata_nonmiss,
           K0,
           K1) {

    if (!is.na(methodvar)) {
      # methodvar needs editing this bit
      trtgp <- mg[i, treatvar]
      # the refernce group comes from the indvidual colunm!
      refergp <- mg[i, referencevar]
    }

    # without this (mg[i,methodvar[1]]) will be cir etc
    methindiv <- mg[i, methodvar]

    methindiv <-
      ifelse((methindiv == "j2r" | methindiv == "J2R" | methindiv == "j2R" | methindiv == "J2r"), 3,
        ifelse((methindiv == "CR" | methindiv == "cr" | methindiv == "Cr" | methindiv == "cR"), 2,
          ifelse((methindiv == "MAR" | methindiv == "mar" | methindiv == "Mar" |
          methindiv == "MAr" | methindiv == "Mr" | methindiv == "MR"), 1,
            ifelse((methindiv == "CIR" | methindiv == "cir" | methindiv == "CIr" | methindiv == "cliR"), 4,
              ifelse((toupper(methindiv) == "CAUSAL" | toupper(methindiv) == "CASUAL" | toupper(methindiv) == "CUASAL"), 6,
                ifelse((methindiv == "LMCF" | methindiv == "lmcf" | methindiv == "Last" | methindiv == "last"), 5,
                  9
                )
              )
            )
          )
        )
      )

    # only done methods 3,4 so far, so Mar need correcting
    # MAR
    if (methindiv == 1) {
      mata_means <- paramBiglist[[M * (trtgp - 1) + m]][1]
      mata_means <- (mata_means[[1]])

      Sigmatrt <- paramBiglist[[M * (trtgp - 1) + m]][2]
      Sigma <- Sigmatrt

      S11 <- Sigmatrt[[1]][c_mata_nonmiss, c_mata_nonmiss]

      # to ensure col pos same as stata
      S12 <-
        matrix(Sigmatrt[[1]][c_mata_nonmiss, c_mata_miss], nrow = length(c_mata_nonmiss))
      S22 <- Sigmatrt[[1]][c_mata_miss, c_mata_miss]
    }
    # 'J2R'
    else if (methindiv == 3) {
      # changed saving the result into  just the param file, list of 2 so can use list index here
      # treatments are 1.. M then M+1 ..2M .. etc

      # 6/3/20 need editing after changing suffixes in paramBiglist
      mata_means_trt <- paramBiglist[[M * (trtgp - 1) + m]][1]
      mata_means_ref <- paramBiglist[[M * (refergp - 1) + m]][1]

      # below causes error after using >1 covars and mata_nonmiss has covar.1, not proper covar names
      mata_means_t <-
        lapply(
          mata_means_trt,
          FUN = function(x)
            x * mata_nonmiss
        )

      mata_means_r <-
        lapply(
          mata_means_ref,
          FUN = function(x)
            x * mata_miss
        )

      # so when all missing  1,1,1, ... then all contributing comes from reference means
      mata_means <- unlist(mata_means_r) + unlist(mata_means_t)
      mata_means <- (as.matrix(t(mata_means)))
      # and preserve names

      ############# SIGMA is from paramsigma  the reference group ################

      SigmaRefer <- paramBiglist[[M * (refergp - 1) + m]][2]
      Sigma <- SigmaRefer

      # note use of [[1]] as is matrix rather than list,

      S11 <- SigmaRefer[[1]][c_mata_nonmiss, c_mata_nonmiss]

      # to ensure rows and cols as should reflect their stucture use matrix
      S12 <-
        matrix(SigmaRefer[[1]][c_mata_nonmiss, c_mata_miss], nrow = length(c_mata_nonmiss))
      S22 <- SigmaRefer[[1]][c_mata_miss, c_mata_miss]

    }
    else if (methindiv == 2) {
      # no need to use Sigmatrt here
      mata_means <- paramBiglist[[M * (refergp - 1) + m]][1]
      # convert from list to matrix
      mata_means <- (mata_means[[1]])

      SigmaRefer <- paramBiglist[[M * (refergp - 1) + m]][2]

      Sigma <- SigmaRefer
      S11 <- SigmaRefer[[1]][c_mata_nonmiss, c_mata_nonmiss]
      S12 <-
        matrix(SigmaRefer[[1]][c_mata_nonmiss, c_mata_miss], nrow = length(c_mata_nonmiss))
      S22 <- SigmaRefer[[1]][c_mata_miss, c_mata_miss]
    }
    else if (methindiv == 4) {
      # need to use Sigmatrt as in j2r
      # pre-deviating use mean of trt gp up to last obs time bfore deviating, post-deviating use mean from ref grp

      mata_Means <- unlist(paramBiglist[[M * (trtgp - 1) + m]][1])

      # convert from list to matrix
      MeansC <- paramBiglist[[M * (refergp - 1) + m]][1]

      # might be better to copy mimix algorithm
      mata_means <- CIR_loop(c_mata_miss, mata_Means, MeansC)

      SigmaRefer <- paramBiglist[[M * (refergp - 1) + m]][2]
      Sigma <- SigmaRefer
      # when reading in Stata sigmas
      # needs to to ths as tibble will fail in cholesky

      S11 <- SigmaRefer[[1]][c_mata_nonmiss, c_mata_nonmiss]
      S12 <-
        matrix(SigmaRefer[[1]][c_mata_nonmiss, c_mata_miss], nrow = length(c_mata_nonmiss))
      S22 <- SigmaRefer[[1]][c_mata_miss, c_mata_miss]
    }
    else if (methindiv == 5) {
      mata_Means <- paramBiglist[[M * (trtgp - 1) + m]][1]

      mata_means <- LMCF_loop(c_mata_miss, mata_Means)
      mata_means <- (mata_means[[1]])

      Sigmatrt <- paramBiglist[[M * (trtgp - 1) + m]][2]
      Sigma <- Sigmatrt
      S11 <- Sigmatrt[[1]][c_mata_nonmiss, c_mata_nonmiss]
      S12 <-
        matrix(Sigmatrt[[1]][c_mata_nonmiss, c_mata_miss], nrow = length(c_mata_nonmiss))
      S22 <- Sigmatrt[[1]][c_mata_miss, c_mata_miss]
    }
    # need to account for interims though
    else if (methindiv == 6) {
      mata_Means <- paramBiglist[[M * (trtgp - 1) + m]][1]
      # convert from list to matrix
      mata_Means <- mata_Means[[1]]

      MeansC <- paramBiglist[[M * (refergp - 1) + m]][1]

      mata_means <- Causal_loop(c_mata_miss, mata_Means, MeansC, K0, K1)

      SigmaRefer <- paramBiglist[[M * (refergp - 1) + m]][2]
      # when reading in Stata sigmas

      S11 <- SigmaRefer[[1]][c_mata_nonmiss, c_mata_nonmiss]
      S12 <-
        matrix(SigmaRefer[[1]][c_mata_nonmiss, c_mata_miss], nrow = length(c_mata_nonmiss))
      S22 <- SigmaRefer[[1]][c_mata_miss, c_mata_miss]

      Sigma <- SigmaRefer
    }
    return(list(mata_means, Sigma, S11, S12, S22))
  }
