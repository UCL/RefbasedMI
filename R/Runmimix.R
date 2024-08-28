#' @title Reference-based multiple imputation of longitudinal data
#' @description Performs reference-based multiple imputation of longitudinal data where data are missing after treatment discontinuation. Methods available are missing at random, jump to reference, copy reference, copy increments in reference, last mean carried forward, the causal model, and delta-adjustment.
#' @details The program works through the following steps:
#'  \enumerate{
#'  \item Set up a summary table based on treatment arm and missing data pattern
#'       (i.e. which timepoints are unobserved)
#'  \item Fit a multivariate normal distribution to each treatment arm using MCMC methods in package norm2
#'  \item Impute all interim missing values under a MAR assumption, looping over treatments and patterns
#'  \item Impute all post-discontinuation missing values under the user-specified assumption,
#'        looping over treatments and patterns (and over methodvar and referencevar if specified)
#'  \item Perform delta-adjustment if specified
#'  \item Repeat steps 2-5 M times and form into a single data frame
#' }
#' @details The baseline value of the outcome could be handed as an outcome, but this would allow a treatment effect at baseline. We instead recommend handling it as a covariate.
#' @details The program is based on Suzie Cro's Stata program mimix
#' @details The user can use the as.mids() function in the mice package to convert the output data to mids data type and then perform analysis using Rubin's rules.
#' @export RefBasedMI
#' @import mice
#' @param data Dataset in long format
#' @param covar Baseline covariate(s): must be complete (no missing values)
#' @param depvar Outcome variable
#' @param treatvar Treatment group variable: can be numeric or character
#' @param idvar Participant identifiervariable
#' @param timevar  Variable indicating time point for repeated measures
#' @param method Reference-based imputation method: must be "J2R", "CR", "CIR", "MAR", "Causal" or "LMCF"
#' @param reference  Reference group for "J2R", "CIR", "CR" methods, or control group for causal method: can be numeric or string
#' @param methodvar Variable in dataset specifying individual method
#' @param referencevar Variable in dataset specifying reference group for individual method
#' @param K0 Causal constant for use with Causal method
#' @param K1 Exponential decaying causal constant for use with Causal method
#' @param delta Optional vector of delta values to add onto imputed values (non-mandatory) (a's in Five_Macros user guide), length equal to number of time points
#' @param dlag Optional vector of delta values to add onto imputed values (non-mandatory) (b's in Five_Macros user guide), length equal to number of time points
#' @param M Number of imputations to be created
#' @param seed  Seed value: specify this so that a new run of the command will give the same imputed values
#' @param prior  Prior when fitting multivariate normal distributions: can be one of "jeffreys" (default), "uniform" or "ridge"
#' @param burnin  Number of burn-in iterations when fitting multivariate normal distributions
#' @param bbetween  Number of iterations between imputed data sets when fitting multivariate normal distributions
#' @param mle Use with extreme caution: do improper imputation by drawing from the model using the maximum likelihood estimates. This does not allow for uncertainty in the MLEs and invalidates interval estimates from Rubin's rules.
#' @return A data frame containing the original data stacked above the M imputed data sets. The original ID variable (idvar) is renamed as .id. A new variable .imp indicates the original data (.imp=0) or the imputed data sets (.imp=1,...,M).
#' @examples
#' # Perform jump to reference imputation on asthma trial data, with reference arm 1 
#' asthmaJ2R <- RefBasedMI(data=asthma, depvar=fev, treatvar=treat, 
#'  idvar=id, timevar=time, method="J2R", reference=1, M=5, seed=54321)

#' # Fit regression model to each imputed data set by treating output data frame as mids object
#' library(mice)
#' fit <- with(data = as.mids(asthmaJ2R), lm(fev ~ factor(treat), subset=(time==12)))

#' # Find pooled treatment effects using Rubin's rules 
#' summary(pool(fit))

# v0.2.0
# @param mle logical option to Use maximum likelihood parameter estimates instead of MCMC draw parameters
# mimix<- function(data,covar=NULL,depvar,treatvar,idvar,timevar,M=1,reference=NULL,method=NULL,seed=101,prior="jeffreys",burnin=1000,bbetween=NULL,methodvar=NULL,referencevar=NULL,delta=NULL,dlag=NULL,K0=1,K1=1,mle=FALSE) {

RefBasedMI<- function(data,covar=NULL,depvar,treatvar,idvar,timevar,method=NULL,reference=NULL,methodvar=NULL,referencevar=NULL,
                 K0=NULL,K1=NULL,delta=NULL,dlag=NULL,M=1,seed=101,prior="jeffreys",burnin=1000,bbetween=NULL,mle=FALSE)
  {

  # test if "data set does not exist!!"
  assertthat::assert_that( length(get("data"))>0 )

  if  (!any((class(get("data"))) == "data.frame")) {stop("data must be type dataframe")}

  # if testinterims then want method to 1stly be MAR this forces interims to be estimated as MAR by default
  # deparse needed as argument in quotes
  depvar <- deparse(substitute(depvar))
  if (length(grep(paste0("^", depvar, "$"), names(get("data")))) == 0) {
    stop(paste(depvar, "(depvar) not in data"))
  }
  treatvar <- deparse(substitute(treatvar))
  if (length(grep(paste0("^", treatvar, "$"), names(get("data")))) == 0) {
    stop(paste(treatvar, "(treatvar) not in data"))
  }
  idvar <- deparse(substitute(idvar))
  if (length(grep(paste0("^", idvar, "$"), names(get("data")))) == 0) {
    stop(paste(idvar, "(idvar) not in data"))
  }
  timevar <- deparse(substitute(timevar))
  if (length(grep(paste0("^", timevar, "$"), names(get("data")))) == 0) {
    stop(paste(timevar, "(timevar) not in data"))
  }
  # check that reference is category of treatment var
  # and is not null (because method not needed if indiv specific cols requested)
  if (!is.null(reference) ) {
      if (!any(as.character(as.matrix(get("data")[,(substitute(treatvar))]))==reference)) { stop("reference must be a category of treatment") }
  }

  # must be null when using methodvar
  #if (!is.null(method)) {method<-deparse(substitute(method))}

  # try read covar without quotes
  scovar<-substitute(covar)

  # check if covar null
  if (length(scovar) > 0 ) {
    # term on rhs copied & pasted!  "‘c’"
    # this  unicode ("\U2018","c","\U2019") for Left & right single quotation mark
    # package wont accept non-ascii char so replace
    # if (sQuote(scovar)[[1]]== "‘c’")
    if (sQuote(scovar)[[1]]== paste0("\U2018","c","\U2019") ) {
       covarname <- vector(mode = "list", length = (length(scovar) - 1))
       for (i in 2:length(scovar)) {
         # check if covar exists 
         if (length(grep(paste0("^", scovar[[i]], "$"), names(get("data")))) == 0) {
           stop(paste(scovar[[i]], "(covar) not in data"))
         }
         covarname[[i - 1]] <-
           names(get("data"))[[grep(paste0("^", scovar[[i]], "$"), names(get("data")))]]
       }
       covarvector <- as.vector(unlist(covarname))
    } 
    else {
       covarvector <-
         names(get("data"))[[grep(paste0("^", scovar, "$"), names(get("data")))]]
    }
    covar <- covarvector
  }
  #else covar just NULL

  # print a warning if reference = NULL



  reference<-substitute(reference)

  # make sure reference converted to numeric by using lookup between unique treatvar and unique tmptreat

  # check reference consistent with treatvar


  if (!is.null(method) & (toupper(method) != "MAR") & (toupper(method) != "LMCF")  ) {
     if (is.null(reference)) {stop("\nStopped !! reference value NULL, required for \"J2R\",\"CIR\",\"CR\",\"Causal\" ")}
  }

  # check treatvar in sorted order
  if (is.unsorted(do.call("order",data.frame(get("data")[,idvar]))) ) {
    stop("\nStopped - warning !! ", idvar,"\n in input data requires to be in sorted order ")
  }

  # try recoding treat, eg 2,3 into 1,2,...
  # should work whether treatvar numeric or char


  # first save class of original treatvar - for use at end
  classtreatvar <- class(unlist(get("data")[, treatvar]))

  # ensure not a tibble  - a when using  readr to read csv data

  tmptreat <- factor(unlist(as.data.frame(get("data"))[, treatvar]))
  initial_levels_treat <- levels(tmptreat)
  levels(tmptreat) <- 1:(nlevels(tmptreat))


  data[,treatvar]<-as.numeric(as.character(tmptreat))
 
  if (!is.null(reference)) {
    reference<-which(initial_levels_treat==reference)
  }

  testinterim<-1
  if (!inherits(mle, "logical") & mle !=0 & mle !=1 ) { stop("mle must be logical value") }
  
  if  (!any((class(get("data"))) == "data.frame")) {stop("data must be type dataframe")}

  # to put quotes in method
  if (!is.null(substitute(method))) {method<-paste0("",(substitute(method)),"")}
  if (!is.null(substitute(methodvar))) {methodvar<-paste0("",(substitute(methodvar)),"")}

  # this won't work no quoted etc
  if (!xor(is.null(method), is.null(methodvar)) ) {stop("Either method or methodvar must be specified but NOT both") }

  # establish whether specifying individual or group by creating a flag var
  if (!is.null(method)) {
    flag_indiv <-0
    # Causal constant must be number and check it exists if Causal specified

    if (toupper(method)=="CAUSAL" |
        toupper(method)== "CASUAL" |
        toupper(method)== "CUASAL") {
      if (missing(K0))  {
        stop("K0 Causal constant not specified")
      }
      if (missing(K1) & (K0 != 0))  {
        stop("K1 Causal constant not specified")
      }
      if (K0 < 0) {
        warning("K0 Causal constant negative.. ")
      }
      if (K0 > 1) {
        warning("K0 Greater than 1.. ")
      }
      if (!(K1 >= 0 & K1 <= 1)) {
        stop("K1 Causal constant not in range 0..1 ")
      }
    } 
    # Causal constant must be number
    K0<- as.numeric(K0)
    K1<- as.numeric(K1)
  }
  if (is.null(method) & !is.null(methodvar)) {
    flag_indiv <- 1
  }

  # find number of covars
  ncovar_i = length(covar)
  
  # if covars exist
  # check that covars are complete AND integer/factor
  if (length(covar)!=0) {
    if (sum(!stats::complete.cases(get("data")[,covar]))!=0) {stop("covariates not complete!!")}
    stopifnot( (sapply((get("data")[,covar]), is.factor)) | (sapply((get("data")[,covar]), is.numeric)) )
  }

  # note if no covar then treat the first depvar level as a covar , eg covar = fev.0.

  # build in checks specifically  for delta
  ntimecol<- get("data")[c(timevar)]
  ntime<-nrow(unique(ntimecol))
  # only run if delta specified
  if (!is.null(delta)) {
    stopifnot(length(delta) == ntime)

    #set dlag to default if 1 1 1 ...if NULL
    if (is.null(dlag)) {
      dlag <- rep(1, length(delta))
    } 
    else  {
      stopifnot(length(dlag) == ntime)
    }
  }


  set.seed(seed)


  # assign all characters in meth to UPPER case
  # check first if meth exists. If statement needed otherwise meth change form Null to character(0)
  if (!is.null(method) | !length(method)==0) {
    method <- toupper(method)
    stopifnot(
      (
        method == "MAR" |
          method == "J2R" | method == "J2" | method == "JR" |
          method == "CIR" | method == "CLIR" |
          method == "CR" |
          method == "LMCF" | method == "LAST"  |
          method == "CAUSAL" | method == "CASUAL"
      ),
      is.numeric(M),
      is.character(depvar),
      is.character(idvar),
      is.character(timevar)
    )
  }


  if (!is.null(method) ) {
    testlist<- preprodata(data,covar,depvar,treatvar,tmptreat,idvar,timevar,M,reference,method, initial_levels_treat)
        reference <- testlist[[7]]
    # need to recode meth
    method<-  ifelse( (  method=="J2R" |method=="J2"|method=="JR" ),3,
                ifelse( ( method=="CR"  ),2,
                  ifelse( ( method=="MAR" | method=="MR" ),1,
                    ifelse( ( method=="CIR" |method=="CLIR" ),4,
                      ifelse( ( method=="LMCF" | method=="LAST" ),5,
                        ifelse( ( method=="CAUSAL" | method=="CASUAL" | method=="CUASAL"),6,9))))))

    # for user specified
  } 
  else if (!is.null(methodvar) ) {
    testlist =  preproIndivdata(data,covar,depvar,treatvar,idvar,timevar,M,reference,method,methodvar,referencevar)
  }


  ntreat<-sort(unlist(testlist[[2]]))
  #sort it !!


  finaldatS<-testlist[[1]]

  mg<-testlist[[3]]


  # vital to get the mata_obs correctly sorted! so corresponds with mimix_group lookup
  # to be consistent with Stata move the base col after the fevs!

  mata_Obs <- testlist[[1]]

  # move treatvar to end by deleting and merging back in
  Obs_treat<-mata_Obs[,c(treatvar,methodvar,referencevar)]
  mata_ObsX<- mata_Obs[,!(names(mata_Obs) %in% c(treatvar,methodvar,referencevar))]
  # combine back the extracted cols
  mata_Obs <- cbind(mata_ObsX,Obs_treat)
  # set name for treatvar otherwise defaults to Obs_treat
  # this ok when methodvar null and tested ok when not null  names just become methodvar and referncevar
  names(mata_Obs)[names(mata_Obs)=="Obs_treat"]<-treatvar


  tst<-stats::reshape(as.data.frame(get("data")[,c(idvar,depvar,timevar)]),v.names = depvar,timevar = timevar,idvar=idvar,direction="wide")


  # change order to put covar last as in Stata  2603

  tst2<-c(names(tst[,-1]),covar)


  # put commas in
  tst3<-paste(tst2,collapse = ",")

  #create input data sets for each treatment from which to model

  for (val in 1:length(ntreat)) {
    assign(paste0("prenormdat",val),subset(finaldatS,finaldatS[,treatvar]==val))
  }


  # CREATE AN EMPTY MATRIX FOR COLLECTING IMPUTED DATA
  # just an empty row but take dimensions and col types fron mata_obs

  GI<-c(0)
  II<-c(0)

  SNO <-mata_Obs[1,1]
  names(SNO)<-idvar
  mata_ObsX <- mata_Obs[,c(2:(ncol(mata_Obs)),1)]
  mata_all_new <- cbind(GI,II,mata_ObsX)
  mata_all_newlist <- vector('list',M*nrow(mg))


  # run the mcmc simulations over the treatment groups

  # create a matrix for param files, a beta and sigma matrices
  # create emptylist for each treat and multiple m's

  paramBiglist <- vector('list',length(ntreat)*M)
  for (val in 1:length(ntreat)) {
    assign(paste0("paramBiglist",val),vector('list',M))
  }
  # create a matrix ntreat by M dimension to store paramBiglists
  paramMatrix<-matrix(1:(length(ntreat)*M),nrow=length(ntreat),ncol=M)



  cumiter<-0

  message(paste0("\nFitting multivariate normal model by ",treatvar,":\n ") )

  for (val in 1:length(ntreat)) {

    kmvar=get(paste0("prenormdat",val))
    prnormobj<-assign(paste0("prnormobj",val), subset(kmvar, select=c(tst2)))
    #create empty list for each treat

    # want to use respective  treat level

    message(paste0("\n",treatvar," = ", (initial_levels_treat)[val],"\nperforming mcmcNorm for m = 1 to ",M,"\n") )
    for(m in 1:M) {

      # suppress warnings regarding solution  near boundary, see norm2 user guide, also mimix about this problem

      # need to error check when ridge or invwish used, the accompanying parameter values supplied.

      if ( prior[1] == "ridge" ) {
        # find sd of depvar over all times for ridge option
        sd_depvar<- stats::sd((get("data")[,depvar]),na.rm=TRUE)
        if ( is.na(prior[2])) { prior[2]<-(sd_depvar*0.1) }
      }
      # invwish not implemented!

      # mle false or 0, true or 1
      if (mle==FALSE) {
        # WARN if not enough data

        # doesn't suppress msgs capture_condition(emResultT<-(emNorm(prnormobj,prior = priorvar[1],prior.df=priorvar[2])) )
        # if error then want to print otherwise dont show

        invisible(capture.output(emResultT<-(emNorm(prnormobj,prior = prior[1],prior.df=prior[2])) ))

        # now test whether emResult created - if not need to see the error msg
        if (is.null(emResultT)) {emResultT<-(emNorm(prnormobj,prior = prior[1],prior.df=prior[2])) }

        if (length(grep("negative definite",emResultT$msg ))>0) {
          message((emResultT$msg))
          message("please disregard UNDECLARED() message - not the error!")
          # UNDECLARED()
        }
        mcmcResultT <- (
            mcmcNorm(
              emResultT,
              iter = burnin,
              multicycle = bbetween,
              prior = prior[1],
              prior.df = prior[2]
            )
          )
        # try for when using mle!

      } 
      else {
        invisible(capture.output(emResultT <-
                                   (
                                     emNorm(prnormobj, prior = prior[1], prior.df = prior[2])
                                   )))
        # for mle
        mcmcResultT <- emResultT
      }

      # cumiter needs to be used as greater than M after 1st treatment
      cumiter<-cumiter+1

      paramBiglist[[ cumiter]] <- mcmcResultT$param
      assign(paste0("paramBiglist",val,"_",m), mcmcResultT$param)
    }
    message(paste0("\nmcmcNorm Loop finished.\n"))
  }

  #store paraBiglist in a single structure


  # can repeat interactively from here

  ############################################ start big loop #########################################
  # now loop over the lookup table mg, looping over every pattern - make sure mata_Obs sorted same way!

  message(paste0("\n\nNumber of original missing values = ", sum(is.na(mata_Obs)), "\n"))
  # declare iterate for saving data

  # not for indiv-specific
  if (flag_indiv==0) {
   message("\nImputing interim missing values under MAR:\n\n")
  } 
  else {
    message("\nImputing missing values using individual-specific method:\n\n")
  }


  # initialise interim
  interim<-0
  # construct structure to save interim ids but this get reinitialised too many times!
  # this may have been causing errors after .id replaced by idvar
  interim_id<- mata_Obs[c(mg[1,1]),idvar]


  m_mg_iter<-0
  for (i in 1:nrow(mg)) {

    # define mata_miss as vector of 1's denoting missing using col names ending i ".missing"
    # this section to be amended to cope with multiple covariates
    mata_miss <- mg[i,grep("*..miss",colnames(mg)),drop=F]

    # find no rows to create covar vector to cbind with mg
    numrows<- nrow(mata_miss)
    mata_nonmiss <- (ifelse(mata_miss==0,1,0))  #define mata-nonmiss from miss

    # need transform nonmiss,miss to c lists - ie. index the
    c_mata_miss<-which(mata_miss==1)
    c_mata_nonmiss<-which(mata_nonmiss==1)
    # eg no missing is c(1,2,3,4,5)

    # count of pattern by treatment
    cnt<- mg$cases[i]

    # treatment grp
    trtgp<- mg[i,treatvar]

    pattern <- mg$patt[i]
    #message("\ntrtgp = ", trtgp)

    if  (!is.null(method) ) {
      # only print imputed case, ie interims
      #message("\n",treatvar ," = ", trtgp,"patt = ",pattern,"number cases = ", cnt)
    } 
    else if(!is.null(methodvar) ) {
      message(treatvar ," = ",trtgp,methodvar," = ",sprintf("%-10s",as.character(mg[i,methodvar])))
      message(referencevar," = ",as.character(mg[i,referencevar]),"pattern = ",pattern,"number patients = ", cnt,"\n")
    }

    trtgpindex<-which(trtgp==ntreat)
    referindex<-which(reference==ntreat)

    # multiple  simulations start here within the pattern loop #########
    for ( m in  1:M)  {
      #FOR INDIVIDUALS WITH NO MISSING DATA COPY COMPLETE DATA INTO THE NEW DATA MATRIX mata_all_new `m' TIMES
      m_mg_iter<-m_mg_iter+1

      # if no missing values
      if(length(c_mata_miss)==0 ) {
        # start and end row positions
        st<-mg[i,"cumcases"]-mg[i,"cases"]+1
        en <-mg[i,"cumcases"]
        # id (SNO) is 1st col try changed 0812
        SNO<-mata_Obs[c(st:en),idvar]

        mata_new<-mata_Obs[c(st:en),]
        # then just move  id col to last
        mata_new<-(mata_new[,c(2:(ncol(mata_new)),1)])

        GI <- array(data=mg[i,treatvar],dim=c(mg[i,"cases"],1))
        #II  no imputations
        II <- array(data=m,dim=c(mg[i,"cases"],1))

        SNO<-mata_Obs[c(st:en),(names(mata_Obs) %in% c(idvar))]
        names(SNO)<-idvar

        mata_new=cbind(GI,II,mata_new)

        mata_all_newlist[[m_mg_iter]]=mata_new
      } 
      else {
        # need to distinguish between meth and methodindiv

        if (flag_indiv==0 ) {

          referindex<- reference
          # FOR INDIVIDUALS WITH  MISSING DATA  `m' TIMES
          # dependent on method chosen
          # 'MAR'

          # if testinterim calc then want to process interims as MAR regardless of value of method
          if (method== 1 | testinterim==1)  {

            mata_means <- paramBiglist[[M*(trtgpindex-1)+m]][1]

            # convert from list element to matrix
            mata_means <- mata_means[[1]]

            Sigmatrt <- paramBiglist[[M*(trtgpindex-1)+m]][2]
            S11 <-Sigmatrt[[1]][c_mata_nonmiss,c_mata_nonmiss]

            # to ensure col pos same as stata
            S12 <-matrix(Sigmatrt[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss))
            S22 <-Sigmatrt[[1]][c_mata_miss,c_mata_miss]

            if (testinterim==1) {
              # need to set method to MAR for 1st pass not works here try above

              # need to find the interim cols so as to set values as MAR
              miss_count <- length(c_mata_miss)

              # now covar has moved to end need adjust by length of covar
              if ((miss_count == 1) & (c_mata_miss[1] < length(mata_means)-length(covar)) ) { # 1 missing and not at end point 5115
                interim<-1
              } 
              else if (miss_count >1) {
                for (b in 2:(miss_count)) {    # last miss_count end pt check separately below
                  if ((c_mata_miss[b-1]+1) != c_mata_miss[b]) { # so next entry after c_mata_miss[b] is non-missing
                    interim<-1
                    if (m==1) {
                      # in case more than 1 interim in patt group
                      for (it_interim in 1:mg[i,"cases"]) {
                        }
                    }
                      # construct vector to save interims ids
  
                  } 
                  else {
                    #note that if all missing then wont be interims and c_mata_nonmiss is integer(0) ie NUL
                    # now covar has moved to end need adjust by excluding positions of covars at the end
          
                    c_mata_nonmiss_nocov <-c_mata_nonmiss[c(which( c_mata_nonmiss <= (length(mata_means) - length(covar)))) ]
          
                    if (length(c_mata_nonmiss[c(which( c_mata_nonmiss <=   (length(mata_means) - length(covar) )     ))]) != 0 )
                      {
                        if ( (length(c_mata_nonmiss)!=0) & (c_mata_miss[b-1]+1 == c_mata_miss[b]) & ( c_mata_miss[b-1] < max(c_mata_nonmiss_nocov)))   
                          #{print("tf")}
                          { 
                            interim<-1 
                          } #need to include outside the for loop  when condition b=miss_count
                    }
                  } #if
                } #for
              }
          
              # now covar has moved to end need adjust by length of covar

              # only works for last missing so instead

              deplen<- length(mata_means)-length(covar)
              if ( length(setdiff(c(c_mata_miss[1]:deplen),c_mata_miss)) != 0 ) {
                interim<-1
                #note there is another message line
                if (m==1) {
                  message(treatvar," = ",initial_levels_treat[trtgp],"pattern = ",pattern,"number patients = ", cnt,"\n")
                }

              } #if

            } #testinterim  processed save  interim_ids

          }

          # causal method uses same matrices as CIR with K parameter

          ############# individual analysis #########################
          #  dont need to do MAR on interims because this could be set by useer within the dataset if required
          # so just need to do in  pass2 
        } 
        else if (flag_indiv==1) {
          # call function for  indiv
          indparamlist  <- ifmethodindiv(methodvar,referencevar,mg,m,M,paramBiglist,i,treatvar,c_mata_nonmiss,c_mata_miss,mata_miss,mata_nonmiss,K0,K1)
  
          mata_means<- indparamlist[[1]]
          Sigma <- indparamlist[[2]]
          S11 <- indparamlist[[3]]
          S12 <- indparamlist[[4]]
          S22 <- indparamlist[[5]]
        }


        # loop still open for row(mg)

        ###################### MNAR IMPUTATION ################
        # need insert routine when ALL missing values
        # else
        ###################### MNAR IMPUTATION ################


        #make sure these are single row vectors! as mistake in LMCF but have to be duplicate rows so add ,s
        #and move after dup fun


        mata_means<- do.call("rbind",replicate(mg[i,"cases"],mata_means,simplify=FALSE))


        if (!is.null(nrow(mata_means)) )  {
          m1 <- mata_means[,c_mata_nonmiss]
          m2 <- mata_means[,c_mata_miss]
        } 
        else {
          m1 <- mata_means[c_mata_nonmiss]
          m2 <- mata_means[c_mata_miss]
        }

        #need a counter to accumulate j  - easiest way is to create a cumulative col in mimix_group
        j <- mg[i,"cases"]

        k <- mg[i,"cumcases"]
        startrow <-(k-j+1)
        stoprow  <-(k)

        #WARNING make sure no id col in 1st
        preraw <-(mata_Obs[c(startrow:stoprow),!(names(mata_Obs) %in% c(idvar))])
        raw1 <- preraw[,c_mata_nonmiss]

        # when all missing data, the length of c_mata_nonmiss must exclude the no. of covariates
        # note  no observed data includes no base line as well

        # covars always non missing

        if (length(c_mata_nonmiss)-length(covar)==0)  {
          # change  because error list obj cannot be coerced to double
          # replaced this with S22 because S11 S12 have 0 values when so try
          U <- try(chol(S22), silent=T)
          if (inherits(U,"try-error")) stop("Error: the covariance matrix for drawing the imputations is not positive definite.")
          
          # generate inverse normal, same as used below
          miss_count<-sum(mata_miss)
          Z <- stats::qnorm(matrix(stats::runif( mg[i,"cases"]* miss_count,0,1),mg[i,"cases"],miss_count))

          # raw and m1 null fields so hut use m2
          meanval = as.matrix(m2)

          # falls over so change
          mata_y1 <-  m2 +Z%*%(U)

          #set dimensions mata_new to mata_y1
          #define new matrix from observed,  id column  (the last)
          mata_new <- preraw
          # mata_new has to be already defined
          mata_new[,c_mata_miss] <- (mata_y1)
          GI <- array(data=mg[i,treatvar],dim=c(mg[i,"cases"],1))
          #II  no imputations
          II <- array(data=m,dim=c(mg[i,"cases"],1))
          # SNO just id col
          SNO <- mata_Obs[c(startrow:stoprow),1]

          #needed as 1 column
          SNO<-as.data.frame(SNO)
          colnames(SNO)[which(colnames(SNO)=="SNO")]<-idvar
          mata_new=cbind(GI,II,mata_new,SNO)
          mata_all_newlist[[m_mg_iter]]=mata_new
        } 
        else {
          t_mimix=solve(S11,S12)
          conds <-  S22-t(S12)%*%t_mimix
  
  
          # below for CR but need checks work with J2R
  
  
          if (mg[i,"cases"] == 1) {
              meanval = (m2) + as.matrix(raw1 - m1)%*%as.matrix(t_mimix)
          }  
          else {
            meanval = as.matrix(m2) + (as.matrix(raw1 - m1)%*%as.matrix(t_mimix))
          }
  
          U <- try(chol(conds), silent=T)
          if (inherits(U,"try-error")) stop("Error: the covariance matrix for drawing the imputations is not positive definite.")
          miss_count=rowSums(mata_miss)
  
          # generate inverse normal
          Z<-stats::qnorm(matrix(stats::runif( mg[i,"cases"]* miss_count,0,1),mg[i,"cases"],miss_count))
  
          mata_y1 = meanval+Z%*%(U)
  
          #define new matrix from observed,  id column  (the last)
          mata_new <- preraw
  
          # assigning the columns where the missing values
  
          if(length(c_mata_miss)==0 ) { 
            mata_new[,c(1:length(tst2))] <- preraw[,c(1:length(tst2))]
          }
          else {
            mata_new[,c_mata_miss] <- mata_y1
          }
  
          # SNO just id col
          SNO <- mata_Obs[c(startrow:stoprow),1]
          names(SNO)<-names(mata_Obs[1])
          # GI treatment grp column 1 (here),II imputation number col, mata_new matrix then SNO is id col.
          GI <- array(data=mg[i,1],dim=c(mg[i,"cases"],1))
          #II  no imputations
          II <- array(data=m,dim=c(mg[i,"cases"],1))
  
          # this works but better to pre-initialise data structure outside loop
          mata_new<-cbind(GI,II,mata_new,SNO)
          names(mata_new)[[ncol(mata_new)]] <-idvar
  
          mata_all_newlist[[m_mg_iter]]=mata_new
  
        }


      } # ( m in  1:M) so insert delta module just before this
    } # for row[mg]

    # check catching interims correct;
    if (flag_indiv==0) {
      if ( (interim==1) & (length(c_mata_miss)!=0)) {
        # save interim ids
        # this line only catches last interim case in the pattern group but if more than 1 case will omit the previous
        # hence need introduce a interim counter or catch all cases
        for (iter_interim in 1:mg[i,"cases"]) {
          interim_id<-  rbind(interim_id,mata_Obs[c(mg[i,"cumcases"]-mg[i,"cases"]+iter_interim),idvar] )
        }
        #re-set interim flag
        interim<-0
      }
    }

  } # for M STOP HERE!!


  impdataset<-getimpdatasets(list(mata_all_newlist,mg,M,method,idvar))


  # only perform following if not individual method  as only 1 pass for that

  if (flag_indiv==0) {
    if (testinterim==1){ # 01/12 try

      interim_id<-as.data.frame(interim_id)
      colnames(interim_id)<-idvar

      rawplusinterim <- fillinterims(impdataset,interim_id,M,idvar,covar)

      #check .id in correrct ccols for mata_Obs
      Imp_Interims<-rawplusinterim[[2]]

      # 1sty obtain rawplusinterim_1
      test1611impD<-rawplusinterim[[2]]
      # check if no interims
      if(nrow(interim_id) !=0) {
        test1611imp1<-subset(as.matrix(test1611impD[test1611impD$.imp==1,]))
      }
      # need match with mata_Obs
      impMarint0<-rawplusinterim[[1]]

      # all cols  need be in same pose, .id need brough to 1st col from last col

      if(nrow(interim_id) !=0) {
        impMarint0[impMarint0[,idvar] %in% test1611imp1[,idvar], ] <- test1611imp1
      }

      # no interim ids to match with so simply

      impMarint1 <-impMarint0
      # need to drop the old patt!!
      impMarint1nopatt<-as.data.frame(impMarint1)[,!(names(as.data.frame(impMarint1)) %in% c("patt"))]

      STSdummy<- apply(as.data.frame(impMarint1nopatt)
                       [,grepl(paste0(depvar,".","[0-9]"),names(as.data.frame(impMarint1nopatt)))],
                       MARGIN=2,function(x) ifelse(!is.na(x),0,1))
      #careful using grep because if same phrases in covar then will be duplicated like head_base and head so must use paste as above

      colnames(STSdummy) <- paste0(colnames(STSdummy),'.miss')

      #change because covar now last col
      if (length(covar)!=0 ) {
        covar_miss<-data.frame(matrix(0,ncol=length(covar),nrow=nrow(STSdummy)))
        names(covar_miss)<-paste0(covar,".miss")
        STSdummy<-cbind(STSdummy,covar_miss)
      }

      sts4D<-(cbind(impMarint1nopatt,STSdummy))

      pows2 <- sapply(1:ncol(STSdummy),function(i) STSdummy[,i]*2^(i-1))
      #need to add up to find patt
      patt <- rowSums(pows2)
      sts4Dpatt<-cbind(sts4D,patt)
      # in case zero covars
      
      # and now find cumX1 cumulative no. cases in each pattern/treatment group
      
      sts4Dpatt$X1<-1
      finaldatSS <-sts4Dpatt[order(sts4Dpatt[,treatvar],sts4Dpatt$patt),]
      
      ex1<-Hmisc::summarize(finaldatSS$X1, by=Hmisc::llist(finaldatSS[,treatvar],finaldatSS$patt),FUN=sum)
      newnames <- c( treatvar,"patt","X1")
      names(ex1)<-newnames
      ex1$X1cum <- cumsum(ex1$X1)
      ex1$exid <- 1:nrow(ex1)
      names(ex1)[names(ex1)=="X1"]<-"cases"
      names(ex1)[names(ex1)=="X1cum"]<-"cumcases"
      #Now find mg table
      # In order  to get patt with all missing patts
      Overall_patt<-unique(sts4Dpatt[grepl(".miss",colnames(sts4Dpatt))])
      patt<- unique(sts4Dpatt[,"patt"])
      
      #then combine depvar and covar patt!
      all_patt<-cbind(Overall_patt,patt)
      
      test_ex1<-merge(ex1,all_patt,by="patt")[order(merge(ex1,all_patt,by="patt")$exid),]
      #  test_ex1 exactly like mg so only need adjust finnaldatSS by taking out differnt names
      finaldatSS<- finaldatSS[,!(names(finaldatSS)) %in%c(".imp","X1")]
      
      # now have achieved the mg table and msata_Obs so no need to call preprodata 2nd time!
      }  # testinterim
    }
  else {
    # also run report on na's
    
    message(paste0("\nNumber of final na values = ", sum(is.na(subset(impdataset,impdataset$.imp>0)))))
    
    # return long data set
    varyingnames<-names((impdataset)[,grepl(paste0(depvar,"\\."),colnames(impdataset))])
    varyingtimes <- gsub(".*\\.","",(varyingnames))
    impdatalong<- reshape(impdataset,varying = varyingnames ,direction="long",sep=".",v.names = depvar,timevar=timevar,times=varyingtimes)
    # re-set time levels to original as could have been transformed to 1,2..)
    # timevar must be numeric (integer?)
    impdatalong[,timevar]<-as.numeric(impdatalong[,timevar])
    
    # merge onto original data set
    impdatamerge<-(merge(get("data"),impdatalong,by.x = c(idvar,timevar),by.y = c(idvar,timevar),all.y=TRUE))
    # can delete all .x's
    impdatamerge<-(impdatamerge[,-c(grep(("\\.x"),colnames(impdatamerge)))])
    # and remove all .y suffixes
    names(impdatamerge)<-gsub("\\.y","",names(impdatamerge))
    # finally re-order
    impdatamergeord<-(impdatamerge[order(impdatamerge[,".imp"],impdatamerge[,timevar]),])
    
    # if idvar ="id" then wil be duplicste id cols name so  delete the 2nd occurence which is the last col
    # this when individual specific cols dataset has idvar equal to id
    if (idvar =="id") {
      impdatamergeord[,ncol(impdatamergeord)]<-NULL
    }
    else {
      # need .id variable  to enable  use of mice in after-analysis
      # overwrite values inid col
      impdatamergeord[,"id"]<- impdatamergeord[,idvar]
    }
    names(impdatamergeord)[names(impdatamergeord)=="id"]<-".id"
    impdatamergeord<- within(impdatamergeord,rm(patt))
    return(impdatamergeord)
   }


  #  finaldatSS only created if interim
  if (nrow(interim_id) !=0) {
    ntreat <- unique(finaldatSS[c(treatvar)])
    mg <- test_ex1
    mata_Obs<- finaldatSS
  }

  colx<-grep(treatvar,colnames(mata_Obs))
  # check if not already last
  if (colx<length(mata_Obs)) {
    mata_Obs.reorder<-mata_Obs[,c(1:(colx-1),(colx+1):length(mata_Obs),colx)]
    mata_Obs<-mata_Obs.reorder
  }
  idcol<- grep(paste0(idvar),colnames(mata_Obs))
  # id to 1st
  # if not already 1st!
  if (idcol !=1) {
    mata_Obs.reorder_id<-mata_Obs[,c(idcol,1:(idcol-1),(idcol+1):length(mata_Obs))]
    # move id to 1st col
    mata_Obs<-mata_Obs.reorder_id
  }


  #*CREATE AN EMPTY MATRIX FOR COLLECTING the after mar  IMPUTED DATA
  GI<-c(0)
  II<-c(0)
  SNO <-mata_Obs[1,1]
  dropid<-c("id")
  mata_ObsX<-(mata_Obs[,c(2:(ncol(mata_Obs)),1)])
  mata_all_new <- cbind(GI,II,mata_ObsX)

  # this just initialises (maybe can do in lass2Loop?)
  mata_all_newlist <- vector('list',M*nrow(mg))

  if (flag_indiv==1) {
    impdataset$treat<-factor(impdataset$treat)
    levels(impdataset$treat)<-initial_levels_treat
    return(impdataset)
  } 
  else {
    # Imp_interims  consists of unimputed record followed by imputed record for each m
    testpass2impdatset<- pass2Loop(Imp_Interims,method,mg,ntreat,depvar,covar,treatvar,
                                   tmptreat,classtreatvar,reference,trtgp,mata_Obs,
                                   mata_all_newlist,paramBiglist,idvar,flag_indiv,M,
                                   delta,dlag,K0,K1,timevar,data, tst2, initial_levels_treat)
  }

}

############ END OF DEFINING FUNCTION REFBASEDMI ###############

getimpdatasets <- function(varlist){

  # to obtain M imputed data sets
  # dimension of data set, nrows in pattern times no imputations,
  # note sub data sets wi have different cols if completely missing so
  mata_all_newlist<-  varlist[1]
  mg<-(varlist[2])
  M<- unlist(varlist[3])
  method<- unlist(varlist[4])
  idvar<-unlist(varlist[5])

  # extract from nested list
  # combine into data set containing M imputed datasets
  mata_all_newData1x <- do.call(rbind,mata_all_newlist[[1]])
  # then sort (by imputation and patient id) into M data sets and split into M lists
  impdatasets <- mata_all_newData1x[order(mata_all_newData1x$II,mata_all_newData1x[,idvar]),]

  #############################################
  # now recreate orig data set by selecting 1st imputed data set and setting NAs using .miss dummies
  imp1st<-impdatasets[impdatasets$II=="1",]
  col_miss<-(imp1st[,grepl(".miss",colnames(imp1st))])
  # get header of missing
  hd_miss <- colnames(imp1st[,grepl(".miss",colnames(imp1st))])

  # strp .miss to gt orig vars names
  hd_new <- gsub(pattern=".miss",replacement = "", hd_miss)
  # converting the 1's in .miss columns to NAs then
  # replace dummy 0,1 by 0,na
  imp1st_NA <- apply(col_miss, MARGIN =2, function(x) ifelse(x==1,NA,x) )

  #assign  0 imputation number to unimputed dataset
  imp1st$II <-0
  # overwrite .miss cols with 0,NA  instead of 0,1
  imp1st[,hd_new]<-(imp1st[,hd_new] +imp1st_NA)

  # now combine recreated original with impute data
  impdatasetsmiss<-rbind(imp1st,impdatasets)

  # no reason to keep the dummies
  # and they cause a warning on setting  as.mids() function
  impdatasets <- impdatasetsmiss[,-grep(".miss",colnames(impdatasetsmiss))]

  # change II, SNo to be consistent with mids format

  names( impdatasets)[names(impdatasets)=="II"]<-".imp"
  #1212 just drop SNO?? needed because are pass2 duplicate id cols would be generated!
  impdatasets$.id <-NULL

  # change row names to e sequential
  rownames(impdatasets)<-NULL
  #drop GI
  impdatasets$GI <-NULL

  # report and check number na's

  if (sum(is.na(subset(impdatasets,impdatasets$.imp>0))) !=0 ) { warning(paste0("\nUnimputed data values")) }
  # write which model processed
  # but not when indiv method used
  if (length(method) !=0 ) {
  if (method==3) {model<-"J2R"} 
    else if (method==2) {model<-"CR"} 
    else if (method==1) {model<-"MAR"} 
    else if (method==4) {model<-"CIR"} 
    else if (method==5) {model<-"LMCF"} 
    else if (method==6) {model<-"CAUSAL"}
  } # length

  return(impdatasets)
}

############ END OF DEFINING FUNCTION getimpdatasets ###############

pass2Loop <- function(Imp_Interims,method,mg,ntreat,depvar,covar,treatvar,tmptreat,classtreatvar,reference,trtgp,mata_Obs,
                     mata_all_newlist, paramBiglist,idvar,flag_indiv,M,delta,dlag,K0,K1,timevar,data, tst2, initial_levels_treat)
{
  # this doesn't call proprocess as data already in wide format
  # for reporting purposes try here rather than runmimix
  if  (method==3) {model<-"J2R"  } 
  else if (method==2) {model<-"CR"  } 
  else if (method==1) {model <-"MAR"  } 
  else if (method==4) {model<-"CIR"  } 
  else if (method==5) {model<-"LMCF"  } 
  else if (method==6) {model<-"CAUSAL" }

  # NOTE check output interims correctly?
  # interims imputed
  No_ints_Impd <- 0
  # need nrow(nrow in case of data frame 0 obs. of  1 variable:
  if (length(nrow(nrow(Imp_Interims))>0)) {
    No_ints_Impd <- sum(is.na(subset(Imp_Interims,.imp==0)))-sum(is.na(subset(Imp_Interims,.imp==1)))
  }
  message(paste0("\nNumber of post-discontinuation missing values = ",sum(is.na(mata_Obs)),"\n"))
  message(paste0("\nImputing post-discontinuation missing values under ",model,":\n\n"))

  m_mg_iter<-0
  #  create a dup col to preserve original numeric levels

  mg[,"orig_treat"] <- mg[,treatvar]

  for (i in 1:nrow(mg)) {
    # define mata_miss as vector of 1's denoting missing using col names ending i ".missing"
    mata_miss <- mg[i,grep("*..miss",colnames(mg)),drop=F]

    #find no rows to create covar vector to cbind with mg
    numrows<- nrow(mata_miss)

    mata_nonmiss <- (ifelse(mata_miss==0,1,0))  #define mata-nonmiss from miss

    # need transform nonmiss,miss to c lists - ie. index the
    c_mata_miss<-which(mata_miss==1)
    c_mata_nonmiss<-which(mata_nonmiss==1)
    # eg no missing is c(1,2,3,4,5)

    # count of pattern by treatment
    cnt<- mg$cases[i]

    # treatment grp

    # assign before recode to preservefor trtgpindex a few lines down
    # instead of  trtgp<- mg[i,treatvar]
    trtgp<- mg[i,"orig_treat"]

    # but note in case need to change to the character values
    mg[,treatvar] <- ordered(mg[,treatvar], labels=as.character(unique(tmptreat)))
    levels(mg[,treatvar]) <- initial_levels_treat

    pattern <- mg$patt[i]

    if  (!is.null(method) ) {

      #try this converting factor to numeric to ensure correct ordering
      mg[,treatvar]<-sort(as.numeric(as.character(mg[,treatvar])))
      message(paste(treatvar," = ", mg[i,treatvar],"pattern = ",pattern,"number patients = ",cnt,"\n")) 
    } 
    else if(!is.null(methodvar) ) {
      message(treatvar," = ",trtgp,methodvar," = ",as.character(mg[i,methodvar]),
          referencevar," = ",as.character(mg[i,referencevar]),
          "pattern = ",pattern,"number patients = ", cnt,"\n\n")
    }

    #necessary!
    trtgpindex<-which(trtgp==ntreat)

    # try thIs when lmcf or mar? ie no or NULL  reference
    # make sure reference not applicable for MAR or LMCF

    if (length(reference) !=0)  {
       referindex<-which(reference==ntreat)
    }

    # multiple  simulations start here within the pattern loop #########
    for ( m in  1:M)  {
      # FOR INDIVIDUALS WITH NO MISSING DATA COPY COMPLETE DATA INTO THE NEW DATA MATRIX mata_all_new `m' TIMES
      m_mg_iter<-m_mg_iter+1

      # at every m we need to replace the already imputed interims estiated from fillinterims subroutune

      # only need create these once, not every iteration/imputation, so test  if the last imputation has already been created
      # because if it has then so have all the others!

      # but not if no interims  !
      if (nrow(Imp_Interims) !=0 )    {
        if (!exists((paste0("Imp_Interims_",M)))  )  {
          assign(paste0("Imp_Interims_",m),subset(as.matrix(Imp_Interims[Imp_Interims$.imp==m,])))
        }
        # also create 0 for use in adddelta
  
        Imp_Interims_0<- subset(as.matrix(Imp_Interims[Imp_Interims$.imp==0,]))
  
        ImpInters <-    get(paste0("Imp_Interims_",m))
  
        colnames(ImpInters)[colnames(ImpInters)==treatvar]<-"treat"
        test_Imp <- subset(ImpInters,select=-c(.imp,patt,treat))
        #  interim for aNTIDEP data")
        # the last col (id) has to be moved to the 1st!
        test_Imp_r<-test_Imp[,c(ncol(test_Imp),1:(ncol(test_Imp)-1))]
  
        # just use as a look up table ! so
        df1<-as.data.frame(mata_Obs)
        # move id to 1st col prob should check cols agree
  
        test_Imp<- as.data.frame(test_Imp)[c(idvar,setdiff(colnames(test_Imp),idvar))]
  
        df2<-as.data.frame(test_Imp)
  
        # find depvar vars
        depcols<-setdiff( grep(paste0(depvar),names(df1)) , grep('.miss',names(df1)))
        # assume these corresponds with interim lookup ! prob need a check here!
  
        # actually faster using match than fmatch
  
        # needed to edit for antidepressant
        depvarnames<-colnames(ImpInters)[grepl(paste0(depvar,"\\."),colnames(ImpInters))]
        matchseq<-match(ImpInters[,idvar],mata_Obs[,idvar])
        for (jj in 1:length(matchseq) ) {
          mata_Obs[,depvarnames][matchseq[jj],]<-as.data.frame(ImpInters)[,depvarnames][jj,]
        }
      }
      
      # if no missing values
      if(length(c_mata_miss)==0 ) {
        # start and end row positions
        st<-mg[i,"cumcases"]-mg[i,"cases"]+1
        en <-mg[i,"cumcases"]
        # id (SNO) is 1st col
        # either reorganise mata_Obs so id is 1st col  or use id as col name 0812

        SNO<-mata_Obs[c(st:en),idvar]
        SNO<-as.data.frame(SNO)
        colnames(SNO)[which(colnames(SNO)=="SNO")]<-idvar

        mata_new<- (mata_Obs[c(st:en),])[,-grep(idvar,colnames(mata_Obs))]

        #treat defined within fun
        GI <- array(data=mg[i,treatvar],dim=c(mg[i,"cases"],1))
        #II  no imputations
        II <- array(data=m,dim=c(mg[i,"cases"],1))
        mata_new=cbind(GI,II,mata_new,SNO)

        # change back name from SNO
        names(mata_new)[[ncol(mata_new)]] <-idvar

        mata_all_newlist[[m_mg_iter]]=mata_new
      } 
      else {
        # need to distinguish between meth and methodindiv

        if (flag_indiv==0 ) {

          referindex<- reference
          #FOR INDIVIDUALS WITH  MISSING DATA  `m' TIMES
          # dependent on method chosen
          
          # 'MAR'
          if (method== 1)  {

            #  checking methd used for interims in J2R works as MAR


            mata_means <- paramBiglist[[M*(trtgpindex-1)+m]][1]

            # convert from list element to matrix
            mata_means <- mata_means[[1]]


            Sigmatrt <- paramBiglist[[M*(trtgpindex-1)+m]][2]
            S11 <-Sigmatrt[[1]][c_mata_nonmiss,c_mata_nonmiss]
            # to ensure col pos same as stata
            S12 <-matrix(Sigmatrt[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss))
            S22 <-Sigmatrt[[1]][c_mata_miss,c_mata_miss]
          } 
          
          # 'J2R'
          else if (method == 3 ) {
            # changed saving the result into  just the param file, list of 2 so can use list index here
            #treatmnets are 1.. M then M+1 ..2M .. etc
            mata_means_trt <- paramBiglist[[M*(trtgpindex-1)+m]][1]
            mata_means_ref <- paramBiglist[[M*(referindex-1)+m]][1]

            mata_means_t <-lapply(mata_means_trt,FUN = function(x) x*mata_nonmiss)


            mata_means_r <-lapply(mata_means_ref,FUN = function(x) x*mata_miss)
            # so when all missing  1,1,1, ... then all contributing comes from reference means
            mata_means <- unlist(mata_means_r)+unlist(mata_means_t)
            mata_means<-(as.matrix(t(mata_means)))


            ############# SIGMA is from paramsigma  the reference group ################


            SigmaRefer <- paramBiglist[[M*(referindex-1)+m]][2]
            Sigmatrt <- paramBiglist[[M*(trtgpindex-1)+m]][2]
            # note use of [[1]] as is matrix rathe than list,
            S11 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_nonmiss]
            # causes non-def error in conds
            #to ensure rows and cols as should reflect their stucture use matrix
            S12 <-matrix(SigmaRefer[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss))
            S22 <-SigmaRefer[[1]][c_mata_miss,c_mata_miss]
          } 

          # 'CR'
          else if (method==2) {
            # no need to use Sigmatrt here
            mata_means <- paramBiglist[[M*(referindex-1)+m]][1]
            # convert from list to matrix
            mata_means <- mata_means[[1]]

            SigmaRefer <- paramBiglist[[M*(referindex-1)+m]][2]
            S11 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_nonmiss]
            S12 <-matrix(SigmaRefer[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss) )
            S22 <-SigmaRefer[[1]][c_mata_miss,c_mata_miss]
          }

          # 'CIR'
          else if (method==4) {
            # need to use Sigmatrt as in j2r
            # pre-deviating use mean of trt gp up to last obs time bfore deviating, post-deviating use mean from ref grp

            # put equiv to mimix
            mata_Means <- paramBiglist[[M*(trtgpindex-1)+m]][1]
            # convert from list to matrix
            mata_Means <- mata_Means[[1]]
            MeansC <-  paramBiglist[[M*(referindex-1)+m]][1]

            # might be better to copy mimix algoRITHM

            mata_means<-CIR_loop(c_mata_miss,mata_Means,MeansC)
            #returns mata_means as single row
            # then duplicate over patt rows

            SigmaRefer <- paramBiglist[[M*(referindex-1)+m]][2]
            # when reading in Stata sigmas

            S11 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_nonmiss]
            S12 <-matrix(SigmaRefer[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss))
            S22 <-SigmaRefer[[1]][c_mata_miss,c_mata_miss]
          }

          # 'LMCF'
          else if (method==5) {
            mata_Means <- paramBiglist[[M*(trtgpindex-1)+m]][1]
            # convert from list to matrix
            mata_Means <- mata_Means[[1]]
            # no ref MeansC <- mata_means_ref
            mata_means<-LMCF_loop(c_mata_miss,mata_Means)


            Sigmatrt <- paramBiglist[[M*(trtgpindex-1)+m]][2]
            # when reading in Stata sigmas
            S11 <-Sigmatrt[[1]][c_mata_nonmiss,c_mata_nonmiss]
            S12 <-matrix(Sigmatrt[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss))
            S22 <-Sigmatrt[[1]][c_mata_miss,c_mata_miss]
          }
          
          # causal method - uses same matrices as CIR with K parameter
          else if (method==6) {
            mata_Means <- paramBiglist[[M*(trtgpindex-1)+m]][1]
            # convert from list to matrix
            mata_Means <- mata_Means[[1]]
            MeansC <-  paramBiglist[[M*(referindex-1)+m]][1]

            mata_means<-Causal_loop(c_mata_miss,mata_Means,MeansC,K0,K1)

            SigmaRefer <- paramBiglist[[M*(referindex-1)+m]][2]
            # when reading in Stata sigmas
            S11 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_nonmiss]
            S12 <-matrix(SigmaRefer[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss))
            S22 <-SigmaRefer[[1]][c_mata_miss,c_mata_miss]
          }
          ############# individual analysis #########################
        
        } 
        else if (flag_indiv==1) {
          # call function for  indiv
          indparamlist  <- ifmethodindiv(methodvar,referencevar,mg,m,M,paramBiglist,i,treatvar,
                                         c_mata_nonmiss,c_mata_miss,mata_miss,mata_nonmiss,K0,K1)
          mata_means<- indparamlist[[1]]
          Sigma <- indparamlist[[2]]
        }


        # loop still open for row(mg)

        ###################### MNAR IMPUTATION ################
        # need insert routine when ALL missing values
        # else
        ###################### MNAR IMPUTATION ################

        #make sure these are single row vectors! as mistake in LMCF but have to be duplicate rows so add ,s
        #and move after dup fun


        # need to replicate mata_means to same number rows as data pattern group

        mata_means<-mata_means[rep(seq_len(nrow(mata_means)),each=mg$cases[i]),]

        if (!is.null(nrow(mata_means)) )  {
          m1 <- mata_means[,c_mata_nonmiss]
          m2 <- mata_means[,c_mata_miss]
        } else {
          m1 <- mata_means[c_mata_nonmiss]
          m2 <- mata_means[c_mata_miss]
        }

        #then mata_obs is the sequential selection of rows according to the mimix_group variable X1 values

        #need a counter to accumulate j  - easiest way is to ceate a cumulstive col in mimix_group

        j <- mg[i,"cases"]

        k <- mg[i,"cumcases"]
        startrow <-(k-j+1)
        stoprow  <-(k)


        # drop id  cols from raw data

        preraw<- (mata_Obs[c(startrow:stoprow),])[,-1]
        raw1 <- preraw[,c_mata_nonmiss]

        if (length(c_mata_nonmiss)-length(covar)==0)  {
          ## routine copied from mimix line 1229

          U <- try(chol(S22), silent=T)
          if (inherits(U,"try-error")) stop("Error: the covariance matrix for drawing the imputations is not positive definite.")
          
          # generate inverse normal, same as used below
          miss_count<-sum(mata_miss)
          Z <- stats::qnorm(matrix(stats::runif( mg[i,"cases"]* miss_count,0,1),mg[i,"cases"],miss_count))

          # raw and m1 null fields so hut use m2
          meanval = as.matrix(m2)

          mata_y1 = m2 + Z%*%(U)

          #set dimensions mata_new to mata_y1
          #define new matrix from observed,  id column  (the last)
          mata_new <- preraw
          # mata_new has to be already defined
          mata_new[,c_mata_miss] <- (mata_y1)
          GI <- array(data=mg[i,treatvar],dim=c(mg[i,"cases"],1))
          #II  no imputations
          II <- array(data=m,dim=c(mg[i,"cases"],1))
          # SNO just id col

          SNO <- mata_Obs[c(startrow:stoprow),idvar]
          #make sure col name consistent, ie idvar
          SNO<-as.data.frame(SNO)
          colnames(SNO)[which(colnames(SNO)=="SNO")]<-idvar
          mata_new=cbind(GI,II,mata_new,SNO)
          mata_all_newlist[[m_mg_iter]]=mata_new
        } 
        else {
          #so S12 must be declare as a matrix ! as number otherwise is class number
          t_mimix=solve(S11,S12)
          conds <-  S22-t(S12)%*%t_mimix

          if (mg[i,"cases"] == 1) {
            meanval = (m2) + as.matrix(raw1 - m1)%*%as.matrix(t_mimix)
          }  
          else {
            meanval = as.matrix(m2) + (as.matrix(raw1 - m1)%*%as.matrix(t_mimix))
          }

          U <- try(chol(conds), silent=T)
          if (inherits(U,"try-error")) stop("Error: the covariance matrix for drawing the imputations is not positive definite.")
          miss_count=rowSums(mata_miss)

          # generate inverse normal
          Z<-stats::qnorm(matrix(stats::runif( mg[i,"cases"]* miss_count,0,1),mg[i,"cases"],miss_count))
          # check same input parameters for inverse norm gen as in stata

          mata_y1 = meanval+Z%*%(U)

          #define new matrix from observed, id column  (the last)
          mata_new <- preraw

          # assigning the columns where the missing values

          # so if no missing then just copy full values into mata_new columns
          if(length(c_mata_miss)==0 ) { mata_new[,c(1:length(tst2))] <- preraw[,c(1:length(tst2))]
          } 
          else {
            mata_new[,c_mata_miss] <- mata_y1
          }

          # need to be consistent for mata_all_newlist!  so try change SNO to patient name
          SNO <- mata_Obs[c(startrow:stoprow),idvar]
          #SNO <- mata_ObsX[,ncol(mata_Obs)]
          # GI treatment grp column 1 (here),II imputation number col, mata_new matrix then SNO is id col.
          GI <- array(data=mg[i,1],dim=c(mg[i,"cases"],1))
          #II  no imputations
          II <- array(data=m,dim=c(mg[i,"cases"],1))

          #this works but betteR to pre-initialise data structure outsidE loop
          mata_new<-cbind(GI,II,mata_new,SNO)
          # change back name from SNO
          names(mata_new)[[ncol(mata_new)]] <-idvar

          #assume delta to be used if specified in input argument
          if (length(delta != 0) ) {
            ncovar_i<-length(covar)
            if (is.null(dlag)) {
              dlag <- rep(1,length(delta))
            }
            mata_new <-  AddDelta(tst2, covar,mata_new,delta,dlag)
          }
          mata_all_newlist[[m_mg_iter]]=mata_new
        }
      } # ( m in  1:M) so insert delta module just before this
    } #for row[mg]
  } #for M STOP HERE!!

  impdataset<-getimpdatasets(list(mata_all_newlist,mg,M,method,idvar))
  
  # but need to adjust orig data set to set interims back to missing
  # only needed if there are interims!
  if (nrow(Imp_Interims)!=0) {
    assign(paste0("Imp_Interims_",0),subset(as.matrix(Imp_Interims[Imp_Interims$.imp==0,])))
    # find depvar vars  cols
    ImpInters <-    get(paste0("Imp_Interims_",0))
    # note keep .imp so can use as match with impdaset and also presreve same col pos
  
    colnames(ImpInters)[colnames(ImpInters)==treatvar]<-"treat"
    test_Imp <- subset(ImpInters,select=-c(patt,treat))
  
    test_Imp<-as.data.frame(test_Imp)
  
    depcolsf<- grep(paste0(depvar),names(impdataset))
    # assume  corresponds with interim lookup ! prob need a check here!!
  
    for ( pos in 1:length(depcolsf)) {
      # match by .imp and id
      impdataset[,depcolsf[pos]][match(paste(test_Imp[,idvar],test_Imp$.imp),paste(impdataset[,idvar],impdataset$.imp))]<-test_Imp[,depcolsf[pos]]
    }
    # moved from getimpdatasets fun

  } 
  #  if nrow(Imp_interims)
  message(paste0("\nNumber of final missing values = ", sum(is.na(subset(impdataset,impdataset$.imp>0))), "\nEnd of RefBasedMI\n"))

  # does it have to be ordered!? yes purpose to put on labels
  impdataset[,ncol(impdataset)-1] <- ordered(impdataset[,ncol(impdataset)-1],labels=levels(tmptreat))

  # in order to output in long format with original data set
  varyingnames<-names((impdataset)[,grepl(paste0(depvar,"\\."),colnames(impdataset))])
  varyingtimes <- gsub(".*\\.","",(varyingnames))

  # reshape has reserved word id so make sure rename if id
  if (idvar=="id") { 
    names(impdataset)[names(impdataset)=="id"]<-"SNOx"  
  }
  impdatalong<- reshape(impdataset,varying = varyingnames ,direction="long",sep=".",
                        v.names = depvar,timevar=timevar,times=varyingtimes)
  # this would do if we didn't want all original data cols
  # but we prob do so need merge orig data
  # re-set time levels to original as could have been transformed to 1,2 as in antidepressant data)

  impdatalong[,timevar]<-as.numeric(impdatalong[,timevar])

  #  now need check if idvar=id
  # if so then delete the id col (because now called SNOx)
  if (idvar=="id") { impdatalong$id<-NULL }
  # now need to rename SNO back to id for the merge with the original input data
  names(impdatalong)[names(impdatalong)=="SNOx"]<-"id"

  impdatamerge<-(merge(get("data"),impdatalong,by.x = c(idvar,timevar),by.y = c(idvar,timevar),all.y=TRUE))
  # can delete all .x's
  impdatamerge<-(impdatamerge[,-c(grep(("\\.x"),colnames(impdatamerge)))])
  # and remove all .y suffixes
  names(impdatamerge)<-gsub("\\.y","",names(impdatamerge))
  # sort
  # for sorted output ( as in input data)
  impdatamergeord<-(impdatamerge[order(impdatamerge[,".imp"],impdatamerge[,idvar],impdatamerge[,timevar]),])

  # Make sure levels are same as original data
  impdatamergeord[,treatvar]<-factor(impdatamergeord[,treatvar])
  levels(impdatamergeord[,treatvar])<-initial_levels_treat
  
  # copy  class of treatvar same as in input data

  class(impdatamergeord[,treatvar])<-classtreatvar

  impdatamergeord[,treatvar] <- levels(impdatamergeord[,treatvar])[(impdatamergeord[,treatvar])]
  # need to repeat this in case has changed to character
  class(impdatamergeord[,treatvar])<-classtreatvar



  # drop patt
  impdatamergeord$patt<-NULL

  return(impdatamergeord)   # pass2loop end
}

############ END OF DEFINING FUNCTION pass2Loop ###############

fillinterims<- function(impdata,interims,Mimp=M,idvar,covar ) {

  # make sure impdata sorted by patient number then no worries about sorting by set key

  impMarint_dt <- data.table::as.data.table(impdata)
  interims_dt <- data.table::as.data.table(interims)

  # no need merge empty interims
  if (nrow(interims) !=0 ) {
    tmpdata<-merge(impdata, interims ,by.x=  idvar,by.y= idvar)
    # but to be the same as impMarint_dt[interims_dt]) need to move id to last col and sort
    tmpdata<-tmpdata[,c(2:(ncol(tmpdata)),1)]
    impboth<-tmpdata[order(tmpdata[,idvar],tmpdata$.imp),]
  } 
  else  {
    tmpdata<- impdata
    impboth<-tmpdata[order(tmpdata[,idvar],tmpdata$.imp),]
  }

  test10<-sapply(impboth,function(x) ifelse(is.na(x) ,1,0) )
  # subtract end cols patt treat id as well as covars  , 3 because patt, treat and idvar cols
  # needed in this form to account when no interims
  test10x<-as.data.frame(test10)[,c(1:(ncol(as.data.frame(test10))-3-length(covar)))]

  # find max non-missing
  lastvalid<-apply(test10x,1, function(x) max(which(x==0))  )
  #merge back

  test1611<-cbind(impboth,lastvalid)

  # over many imps, find mean value by interims for imp>0 then combineback to
  # have .imp0 and mean row

  test1611imp0<- subset(test1611,.imp==0)
  #declare void matrix for rbinding
  test1611impD <- test1611imp0
  for (val in 1:Mimp) {
    rbind(assign(paste0("test1611imp",val),subset(test1611,.imp==val | .imp==0)),test1611imp0)

    test1611impx<-as.matrix(get(paste0("test1611imp",val )))
    for (r in seq(from = 1, to = nrow(test1611impx) - 1, by = 2)) {
      for (j in 2:(ncol(test1611impx) - 3)) {
        # if element na and must be bfore lastvalid non missing
        if (is.na(test1611impx[r, j]) &
            (j < test1611impx[r, ncol(test1611impx)])) {
          # then shift value from below row
          test1611impx[r, j] = test1611impx[r + 1, j]
        }
      }
    }
    # keep imp0's and recode non imp0 to imp (m) number
    test1611impxm <-  subset(test1611impx,test1611impx[,".imp"]==0)
    test1611impxm[,".imp"]<-val

    # build up over M imps
    # impxm is m'th imputed data set
    test1611impD <- rbind(test1611impD,test1611impxm)
  }

  # so need to insert into original unimputed data set for each imputation and process each data set into the 2nd pass
  
  test1611impD$lastvalid <-NULL
  #extract m'th imp
  test1611impD02<-as.matrix(test1611impD[test1611impD$.imp==2,])

  # is this necessary??
  impMarint <- as.matrix(impMarint_dt)

  impMarint0<-as.matrix(impMarint_dt[impMarint_dt$.imp==0,])
  # this works well (but MUST Be MATRICES)

  if (nrow(interims) !=0) {
    return(list(impMarint0, test1611impD))
  } 
  else {
    return(list(impMarint0, interims))
  }
}

############ END OF DEFINING FUNCTION fillinterims ###############
