#' ###################################################################################
#' Rmimix.R                                                                          #
#' R program to mimic stata program mimix                                            #
#' ie reference based imputation                                                     #
#' Note 1st part is to set up a summary table based on missing data pattern-         #
#' mg  mimix_group                                                                   #
#' reflects the pattern and treatment group configuration of the raw data            #
#' then acts as a looping mechanism                                                  #
#' norm2 is used as MCMC multivariate normal                                         #
#' this version 22/5/2020                                                            #
#' calls functions listed in functions.R file                                        #
#' function preprodata prepares data and and finds mg (the mimix group)              #
#' function Runmimix  performs major analysis                                        #
#' the required packages as listed in utilities file                                 #
#' v0.0.3                                                                            #
#' Author : Kevin McGrath                                                            #
#' ###################################################################################


#' @title mimix
#' @description main wrapper for running mimix
#' @details This is based on Suzie Cro's Stata program
#' @export mimix
#' @param data  datset in wide format
#' @param covar covariates and base depvar must be complete (no missing vaules)
#' @param depvar dependent variable
#' @param treatvar treatment group , recoded to 1,2,..
#' @param idvar patient id
#' @param M number of imputations
#' @param timevar time point for repeated measure
#' @param refer  reference group for j2r,cir,cr mrthods
#' @param meth RBI method
#' @param seedval  seed value to obtain same outputs
#' @param priorvar  prior,defualt Jeffries, also uniform or ridge
#' @param burnin  burnin value
#' @param bbetween  value between iterations in mcmc
#' @param methodindiv  2 element vector designating variables in data specifying individual method and reference group
#' @param delta vector of delta values to add onto imputed values (non-mandatory)
#' @param Kd Causal constant for use with Causal method
#' @return impdataset the m impute data-sets appended to the "missing values" data-set in wide format
#' @example
#' \dontrun{
#' mimix("asthma",c("base"),"fev","treat","id","time",10,1,"J2R",101,"jeffreys",1000,NULL,NULL,c(0.5,0.5,1,1),0.6)
#' }

mimix<- function(data,covar=NULL,depvar,treatvar,idvar,timevar,M=1,refer=NULL,meth=NULL,seedval=101,priorvar="jeffreys",burnin=1000,bbetween=NULL,methodindiv=NULL,delta=NULL,Kd=NULL) {

  # insert error checks  HERE

  #check not both meth and methodundiv specified.
  stopifnot(meth=="NULL" | methodindiv=="NULL")
  # establish whether specifying individual or group by creating a flag var
  if (!is.null(meth)) {
    flag_indiv <-0
  # Causal constant must be number and check it exists if Causal specified
  #stopifnot(meth=="Causal" & !missing(Kd))
  if (toupper(meth)=="CAUSAL" | toupper(meth)== "CASUAL" | toupper(meth)== "CUASAL") {
    if (missing(Kd))  {stop("Kd Causal constant not specified")}
  } #Causal constant must be number
    Kd<- as.numeric(Kd)
  }
  if (is.null(meth) & !is.null(methodindiv[1]) ) {flag_indiv<-1 }

  #find no covars
  ncovar_i = length(covar)
  # if covars exist
  # check that covars are complete AND integer !!
  # if covar=NULL  17/05/20   #discuss.analyticsvidyha.com
  if (length(covar)!=0) {
    stopifnot(sum(!stats::complete.cases(get(data)[,covar]))==0)
    stopifnot(sapply((get(data)[,covar]),is.numeric))
    #stopifnot(is.numeric((get(data)[,covar])))
    #stopifnot(apply((get(data)[,covar]),2,is.numeric) )
    # really should be integer
    # if need to convert
    # test<-as.data.frame(apply((get(data)[,covar]),2,as.integer) )
  }





  # note if no covar then treat the first depvar level as a covar , eg covar = fev.0.
  # sO must preform small routne tO transform datA set INTO baseline covars.




  set.seed(seedval)

#21/05
  # assign all characters in meth to UPPER case
  # check i1st if meth exists! if statement needed otherwise meth change form Null to character(0)
  if (!is.null(meth) | !length(meth)==0 ) {
  meth <- toupper(meth)
  stopifnot( (meth == "MAR" | meth=="MR"|
              meth=="J2R" | meth=="J2" | meth=="JR"|
              meth=="CIR" | meth=="CLIR" |
              meth=="CR" |
              meth=="LMCF" | meth=="LAST"  |
              meth=="CAUSAL" | meth== "CASUAL" | meth== "CUASAL"),

             is.numeric(M),
             is.character(depvar),
             # treatvar could be char or numeric ? then refer must be same type
             #is.character(treatvar),
             is.character(idvar),
             is.character(timevar)
             #is.character(covar[[2]],
             # only if covar exists  so if NULL noneed to errr check but check if covar Null then...
             #is.character(covar)
            )
  }


  if (!is.null(meth) ) {
    testlist = do.call( preprodata,list(data,covar,depvar,treatvar,idvar,timevar,M,refer,meth))
    refer <- testlist[[7]]



    # stopifnot(refer %in% ntreat)
    #no need for this?
    #meth <- testlist[[10]]
    #need to recode meth
    meth<-  ifelse( (  meth=="J2R" |meth=="J2"|meth=="JR" ),3,
                    ifelse( ( meth=="CR"  ),2,
                            ifelse( ( meth=="MAR" | meth=="MR" ),1,
                                    ifelse( ( meth=="CIR" |meth=="CLIR" ),4,
                                            ifelse( ( meth=="LMCF" | meth=="LAST" ),5,
                                              ifelse( ( meth=="CAUSAL" | meth=="CASUAL" | meth=="CUASAL"),6,9))))))


  }
  else if (!is.null(methodindiv[1]) ) {
    testlist = do.call( preproIndivdata,list(data,covar,depvar,treatvar,idvar,timevar,M,refer,meth,methodindiv))
    # need to re-set meth for individual
    #meth <- testlist[[8]][1]

  }




  ntreat<-sort(unlist(testlist[[2]]))
  #sort it !!



  #ntreat <-testlist$ntreat
  #stopifnot()
  finaldatS<-testlist[[1]]
  #finaldatS <- testlist$finaldatS
  mg<-testlist[[3]]

  # vital to get the mata_obs correctly sorted! so corresponds with mimix_group lookup
  # to be consistent with Stata move the base col after the fevs!
  #mata_Obs <- testlist$finaldatS
  mata_Obs <- testlist[[1]]

  # move treatvar to end by deleting and merging back in

  # extract treatvar and methodindiv vars
  Obs_treat<-subset(mata_Obs,select=c(treatvar,methodindiv))
  mata_ObsX<- mata_Obs[,!(names(mata_Obs) %in% c(treatvar,methodindiv))]
  # combine back the extracted cols
  mata_Obs <- cbind(mata_ObsX,Obs_treat)





  tst<-stats::reshape(get(data)[,c(idvar,depvar,timevar)],v.names = depvar,timevar = timevar,idvar=idvar,direction="wide")

  tst2<-c(covar,names(tst[,-1]))



  # put commas in
  tst3<-paste(tst2,collapse = ",")

  #create input data sets for each tment from which to model
  #browser()
  for (val in 1:length(ntreat)) {
    print(paste0("prenormdat",val))
    # this fails as treatvar not substantiated
    # assign(paste0("prenormdat",val),subset(finaldatS,treatvar==val))
    assign(paste0("prenormdat",val),subset(finaldatS,finaldatS[,treatvar]==val))
  }



  #*CREATE AN EMPTY MATRIX FOR COLLECTING IMPUTED DATA
  # just an empty row but take dimensions and col types fron mata_obs
  #mata_all_new<-array(data=1,dim=c(1,2+ncol(mata_Obs)))
  GI<-c(0)
  II<-c(0)
  SNO <-mata_Obs[1,1]
  dropid<-c("id")
  #31/03
 # browser()
  mata_ObsX<-mata_Obs[,!(names(mata_Obs) %in% dropid)]
  mata_all_new <- cbind(GI,II,mata_ObsX,SNO)

  # browser()
  # is this line needed?
  #mata_all_new[ mata_all_new>=0] <-NA
  mata_all_newlist <- vector('list',M*nrow(mg))
  #Warning message:
  #In Ops.factor(left, right) : ‘>=’ not meaningful for factors



  # run the mcmc simulations over the treatment grps
  # create a matrix for param files, a beta  and sigma matrices
  #paramMatrix<-matrix(1:(nrow(ntreat)*M),nrow=M,dimnames=list(c(1:M),c(1:nrow(ntreat))))



  #create  emptylist for each treat and multiple m's

  paramBiglist <- vector('list',length(ntreat)*M)
  for (val in 1:length(ntreat)) {
    assign(paste0("paramBiglist",val),vector('list',M))
  }

  #create a matrix ntreat by M dimension to store paramBiglists
  paramMatrix<-matrix(1:(length(ntreat)*M),nrow=ntreat,ncol=M)


  #paramMatrixT<-matrix(,nrow=2,ncol=M)
  #start_time <- proc.time()
  #instead iof start _time try system_time
  #paramBetaMatrixT <- matrix(1,nrow=1,ncol=5)
  #paramSigmaMatrixT <- matrix(,nrow=5,ncol=5)
  cumiter<-0
  #system.time(



  for (val in 1:length(ntreat)) {

    kmvar=get(paste0("prenormdat",val))
    prnormobj<-assign(paste0("prnormobj",val), subset(kmvar, select=c(tst2)))
    #create  emptylist for each treat

    paste0("Looping for treatment = ",val," performing mcmcNorm for m = 1 to ",M)
    for(m in 1:M) {

      # supppress warnings regarding solution  near boundary, see norm2 user guide, also mimix about this problem

      #need to error check when ridge or invwish used, the accompanying parameter values supplied.



      if ( priorvar[1] == "ridge" ) {
        # find sd of depvar over all times for ridge option
        sd_depvar<- stats::sd((get(data)[,depvar]),na.rm=TRUE)
        if ( is.na(priorvar[2])) { priorvar[2]<-(sd_depvar*0.1) }}
      # invwish not implemented!
      #if ( priorvar[1] == "invwish" ) { stopifnot(priorvar[2]>0 & priorvar[3]>0 ) }


      # doesnt suppress msgs capture_condition(emResultT<-(norm2::emNorm(prnormobj,prior = priorvar[1],prior.df=priorvar[2])) )
      emResultT<-(norm2::emNorm(prnormobj,prior = priorvar[1],prior.df=priorvar[2]))
      #mcmcResultT<- (mcmcNorm(emResultT,iter=1000,multicycle = NULL,prior = priorvar[1],prior.df = priorvar[2]))
      mcmcResultT<- (norm2::mcmcNorm(emResultT,iter=burnin,multicycle = bbetween,prior = priorvar[1],prior.df = priorvar[2]))

      # msg from emNorm
      #Note: Finite-differencing procedure strayed outside
      #parameter space; solution at or near boundary
      #OCCURRED IN: estimate_worst_frac in MOD norm_engine


      cumiter<-cumiter+1
      #browser()
      # to see if can save directly to val indice
      paramBiglist[[cumiter]] <- mcmcResultT$param
      assign(paste0("paramBiglist",val,"_",m), mcmcResultT$param)
    }

  }

  #store paraBiglist in a single structure

  # ) # system.time
  print(paste0("mcmcNorm Loop finished, m = ",M))




  # can repeat interactively from here

  ############################################ start big loop #########################################
  # now loop over the lookup table mg, looping over every pattern - make sure mata_Obs sorted same way!

  # declare iterate for saving data
  m_mg_iter<-0
  for (i in 1:nrow(mg))
  {

    # define mata_miss as vector of 1's denoting missing using col names ending i ".missing"
    # this section to be amended to cope with multiple covariates

    mata_miss <- mg[i,grep("*..miss",colnames(mg)),drop=F]
    #mata_miss <- mimix_group[i,c(2,3,4,5)]     #defimne mata_miss
    #assumes covariate non missing
    # now >1 covariate code must be amended

    #find no rows to create covar vector to cbind with mg
    numrows<- nrow(mata_miss)


    mata_nonmiss <- (ifelse(mata_miss==0,1,0))  #define mata-nonmiss from miss

    # need transform nonmiss,miss to c lists - ie. index the
    c_mata_miss<-which(mata_miss==1)
    c_mata_nonmiss<-which(mata_nonmiss==1)
    # eg no missing is c(1,2,3,4,5)


    # count of pattern by treatment
    cnt<- mg$X1[i]

    # treatment grp
    trtgp<- mg[i,treatvar]

    pattern <- mg$patt[i]
    #cat("\ntrtgp = ", trtgp)


    if  (is.null(methodindiv[1]) ) {
      cat("\ntrtgp = ", trtgp,"patt = ",pattern,"no patients = ", cnt)
    }else if(!is.null(methodindiv[1]) ) {
      cat("\ntrtgp = ", trtgp,"method= ",as.character(mg[i,methodindiv[1]]),"refgp=",as.character(mg[i,methodindiv[2]]),"patt = ",pattern,"no patients = ", cnt)
    }

    #need to convert (relate) treatment group to position in ntreat (create Pindex vector)
    #unneceassary now recoded

    trtgpindex<-which(trtgp==ntreat)
    referindex<-which(refer==ntreat)

    # multiple  simulations start here within the pattern loop #########
    for ( m in  1:M)  {
      #*FOR INDIVIDUALS WITH NO MISSING DATA COPY COMPLETE DATA INTO THE NEW DATA MATRIX mata_all_new `m' TIMES
      #if `pat' == 0{
      m_mg_iter<-m_mg_iter+1
      # if no missing values
      if(length(c_mata_miss)==0 ) {
        # start and end row positions
        st<-mg[i,"X1cum"]-mg[i,"X1"]+1
        en <-mg[i,"X1cum"]
        # id (SNO) is 1st col
        SNO<-mata_Obs[c(st:en),1]
        #9/5/20
        #browser()
        # this doesnt delete treat
        #mata_new <- mata_Obs[c(st:en),2:ncol(mata_Obs)]
        mata_new<-mata_Obs[c(st:en),!(names(mata_Obs) %in% c(idvar))]
        #browser() # treat defined within fun
        GI <- array(data=mg[i,treatvar],dim=c(mg[i,"X1"],1))
        #II  no imputations
        II <- array(data=m,dim=c(mg[i,"X1"],1))
        mata_new=cbind(GI,II,mata_new,SNO)
        #names(mata_all_new)<-names(mata_new)
        mata_all_newlist[[m_mg_iter]]=mata_new
        # mata_all_new=rbind(mata_all_new,mata_new)
      } else {
        # need to distinguish between meth and methodindiv

        if (flag_indiv==0 ) {

          referindex<- refer
          #FOR INDIVIDUALS WITH  MISSING DATA  `m' TIMES
          # dependent on method chosen
          # 'MAR'
          if (meth== 1)  {

            #mata_means <- get(paste0("param",trtgp,m))[1]
            mata_means <- paramBiglist[[M*(trtgpindex-1)+m]][1]
            # mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]
            # convert from list element to matrix
            mata_means <- mata_means[[1]]


            #Sigmatrt <- get(paste0("param",trtgp,m))[2]
            Sigmatrt <- paramBiglist[[M*(trtgpindex-1)+m]][2]
            S11 <-Sigmatrt[[1]][c_mata_nonmiss,c_mata_nonmiss]
            #S12 <-Sigmatrt[[1]][c_mata_nonmiss,c_mata_miss]
            # to ensure col pos same as stata
            S12 <-matrix(Sigmatrt[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss))
            S22 <-Sigmatrt[[1]][c_mata_miss,c_mata_miss]

            Sigma<-Sigmatrt
          }
          # 'J2R'
          else if (meth == 3 ) {

            # changed saving the result into  just the param file, list of 2 so can use list index here
            #treatmnets are 1.. M then M+1 ..2M .. etc
            mata_means_trt <- paramBiglist[[M*(trtgpindex-1)+m]][1]
            mata_means_ref <- paramBiglist[[M*(referindex-1)+m]][1]


            #  browser()
            # below causes error after using >1 covars and mata_nonmiss has covar.1, not proper covar names

            # try not unisting because error only fr J2r merod when on patien in patt
            #

            mata_means_t <-lapply(mata_means_trt,FUN = function(x) x*mata_nonmiss)
            #mata_means_t <- unlist(mata_means_trt)*mata_nonmiss
            # print(paste0("mata_means_trt, mata_nonmiss= ",mata_means_trt,mata_miss))

            mata_means_r <-lapply(mata_means_ref,FUN = function(x) x*mata_miss)
           #mata_means_r <- unlist(mata_means_ref)*mata_miss
            # so when all missing  1,1,1, ... then all contributing comes from reference means
            mata_means <- unlist(mata_means_r)+unlist(mata_means_t)
            #try this 11/04
            mata_means<-(as.matrix(t(mata_means)))
            # and preserve names
           # 11/04 colnames(mata_means) <-  colnames(mata_means_t[[1]])
            #replicate to number of rows defined by X1
            #mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]


            ############# SIGMA is from paramsigma  the reference group ################


            # do we ever  use SigmaTrt !!?? in j2r?
            # answer is ??, need to use SigmaTrt for the predeviation observations, ie up to where they go missing
            # only after they go missing (trailing missing) need to use the SigmaRef

            #9/3/20
            SigmaRefer <- paramBiglist[[M*(referindex-1)+m]][2]
            #SigmaRefer <- paramBiglist[[M*(referindex-1)+m]][2]
            #SigmaRefer <- get(paste0("paramBiglist",refer,m))[2]
            Sigmatrt <- paramBiglist[[M*(trtgpindex-1)+m]][2]


            # note use of [[1]] as is matrix rathe than list,

            S11 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_nonmiss]
            # causes non-def error in conds

            #to ensure rows and cols as should reflect their stucture use matrix
            S12 <-matrix(SigmaRefer[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss))
            S22 <-SigmaRefer[[1]][c_mata_miss,c_mata_miss]

            #edit due to passing  Sigma after loop, ie save sigma
            #  Sigma <- SigmaRefer[[1]]

            Sigma<-SigmaRefer
          }
          # 'CR'
          else if (meth==2) {
            # no need to use Sigmatrt here
            mata_means <- paramBiglist[[M*(referindex-1)+m]][1]
            #mata_means <- get(paste0("param",refer,m))[1]
            # convert from list to matrix
            mata_means <- mata_means[[1]]

            SigmaRefer <- paramBiglist[[M*(referindex-1)+m]][2]
            S11 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_nonmiss]
            S12 <-matrix(SigmaRefer[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss) )
            S22 <-SigmaRefer[[1]][c_mata_miss,c_mata_miss]

            Sigma<-SigmaRefer
          }
          # 'CIR'
          else if (meth==4)

          {
            # need to use Sigmatrt as in j2r
            # pre-deviating use mean of trt gp up to last obs time bfore deviating, post-deviating use mean from ref grp



            # put equiv to mimix
            #mata_Means <- get(paste0("param",trtgp,m))[1]
            mata_Means <- paramBiglist[[M*(trtgpindex-1)+m]][1]
            # convert from list to matrix
            mata_Means <- mata_Means[[1]]
            #mata_Means <-  get(paste0("parambeta",trtgp,m))
            #MeansC <-  get(paste0("param",refer,m))[1]
            MeansC <-  paramBiglist[[M*(referindex-1)+m]][1]

            #might be better to copy mimix algol

            mata_means<-CIR_loop(c_mata_miss,mata_Means,MeansC)
            #returns mata_means as single row
            # then duplicate over patt rows
            #replicate to number of rows defined by X1
            # mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]

            #SigmaRefer <- get(paste0("paramsigma",refer,m))

            #SigmaRefer <- get(paste0("param",refer,m))[2]
            SigmaRefer <- paramBiglist[[M*(referindex-1)+m]][2]
            # when reading in Stata sigmas


            S11 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_nonmiss]
            S12 <-matrix(SigmaRefer[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss))
            S22 <-SigmaRefer[[1]][c_mata_miss,c_mata_miss]

            Sigma<- SigmaRefer
          }
          # 'LMCF'
          else if (meth==5) {

            #mata_Means <-  get(paste0("param",trtgp,m))[1]
            mata_Means <- paramBiglist[[M*(trtgpindex-1)+m]][1]
            # convert from list to matrix
            mata_Means <- mata_Means[[1]]
            # no ref MeansC <- mata_means_ref
            mata_means<-LMCF_loop(c_mata_miss,mata_Means)
            #mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]


            #Sigmatrt <- get(paste0("param",trtgp,m))[2]
            Sigmatrt <- paramBiglist[[M*(trtgpindex-1)+m]][2]
            # when reading in Stata sigmas
            S11 <-Sigmatrt[[1]][c_mata_nonmiss,c_mata_nonmiss]
            S12 <-matrix(Sigmatrt[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss))
            S22 <-Sigmatrt[[1]][c_mata_miss,c_mata_miss]

            Sigma<- Sigmatrt
          }  #if meth=5
          # causal method uses same matrices as CIR with K parameter
         else if (meth==6)  {
           mata_Means <- paramBiglist[[M*(trtgpindex-1)+m]][1]
           # convert from list to matrix
           mata_Means <- mata_Means[[1]]
           #mata_Means <-  get(paste0("parambeta",trtgp,m))
           #MeansC <-  get(paste0("param",refer,m))[1]
           MeansC <-  paramBiglist[[M*(referindex-1)+m]][1]
           #put Kd tempval
           #Kd =0 eq0iv to J2R?
           #Kd =1 equiv to CIR
           #Kd<-0.8

           mata_means<-Causal_loop(c_mata_miss,mata_Means,MeansC,Kd)

          #this temporary  fpr test purposes until algo decided upon
           SigmaRefer <- paramBiglist[[M*(referindex-1)+m]][2]
           # when reading in Stata sigmas


           S11 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_nonmiss]
           S12 <-matrix(SigmaRefer[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss))
           S22 <-SigmaRefer[[1]][c_mata_miss,c_mata_miss]

           Sigma<- SigmaRefer

         }
          ############# individual analysis #########################
        }
        #else if(!is.null(methodindiv[1]))
        else if (flag_indiv==1) {

          # call function for  indiv

        indparamlist  <- ifmethodindiv(methodindiv,mg,m,M,paramBiglist,i,treatvar,c_mata_nonmiss,c_mata_miss,mata_miss,mata_nonmiss)
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

        S11 <-Sigma[[1]][c_mata_nonmiss,c_mata_nonmiss]
        S12 <-matrix(Sigma[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss))
        S22 <-Sigma[[1]][c_mata_miss,c_mata_miss]

        # need to repictae mata_means to same numbe rrows as data pattern group

        # causing problems replace with simpler
        mata_means<-mata_means[rep(seq_len(nrow(mata_means)),each=mg$X1[i]),]



        #perhaps put zeros in ? to make sure same length??
        if (!is.null(nrow(mata_means)) )  {
          m1 <- mata_means[,c_mata_nonmiss]
          m2 <- mata_means[,c_mata_miss]
        } else {
          m1 <- mata_means[c_mata_nonmiss]
          m2 <- mata_means[c_mata_miss]
        }


        #then mata_obs is the sequential selection of rows according to the mimix_group variable X1 values

        #need a counter to accumulate j  - easiest way is to ceate a cumulstive col in mimix_group



        j <- mg[i,"X1"]

        #for debug
        #print(paste0(" count in patt = ", j))

        k <- mg[i,"X1cum"]
        startrow <-(k-j+1)
        stoprow  <-(k)



        #raw1 <- mata_obs[, c_mata_nonmiss]
       # preraw <- mata_Obs[c(startrow:stoprow),2:ncol(mata_Obs)]

        # drop id  cols from raw data
        preraw <-(mata_Obs[c(startrow:stoprow),!(names(mata_Obs) %in% c(idvar))])
        raw1 <- preraw[,c_mata_nonmiss]

        ##### try inserting routine for all missing values here NOt really necessary!! ### 20/1/20 #####
        ## when all missing data, the length of c_mata_nonmiss must exclude the no. of covariates
        ## temporary fix set as 1, ie no. of covars 20/01/20
        # nOTE think no observed data includes no base line as well?


        # all missing ? didnt think we did ths scenario, ie base depvar always complete?
        if (length(c_mata_nonmiss)==0)  {
          ## routine copied from mimix line 1229

          # change 9/5/20 becaus error list obj cannot be coerced to double
          U <- chol(Sigma[[1]])

          # generate inverse normal, same as used below
          #for debug 20/01
          #   print(paste0('mg[i,X1] =',mg[i,"X1"],' miss_count= ',miss_count))
          miss_count<-sum(mata_miss)
          Z <- stats::qnorm(matrix(stats::runif( mg[i,"X1"]* miss_count,0,1),mg[i,"X1"],miss_count))

          # raw and m1 null fields so hut use m2
          meanval = as.matrix(m2)

          mata_y1 = meanval+Z%*%(U)
          #set dimensions mata_new to mata_y1
          #define new matrix from observed,  id column  (the last)
          mata_new <- preraw
          # mata_new has to be already defined
          mata_new[,c_mata_miss] <- (mata_y1)
          GI <- array(data=mg[i,treatvar],dim=c(mg[i,"X1"],1))
          #II  no imputations
          II <- array(data=m,dim=c(mg[i,"X1"],1))
          # SNO just id col
          SNO <- mata_Obs[c(startrow:stoprow),1]
          mata_new=cbind(GI,II,mata_new,SNO)
          mata_all_newlist[[m_mg_iter]]=mata_new

          # 9/5/20 inser else here

        } else {


        #so S12 must be declare as a matrix ! as number otherwise is class number
        #try transpose as was failing compatible error
               #note y must have same no rows as Q
        #solve Qx=y
        #t_mimix =cholsolve(Q=mS11,y=(mS12))
        #tryig solve instead ax=b
        t_mimix=solve(S11,S12)
        conds <-  S22-t(S12)%*%t_mimix
        #
        #meanval = as.matrix(m2) + as.matrix(raw1 - m1)%*%as.matrix(t_mimix)
        # below for CR but need checks work with J2R
        # need edit m2 as errormsg "data frame with 0 columns and 1 row"
        # 9/3/20
        #browser() the else condition not works ,m2 needs matrix

        ## important calculation to CHECK !!!
        #depends whether meth or methodindiv used
        #really this needs to be independent of which merh used
        # if (meth=='J2R')
        #  if  (mg[i,methodindiv[1]] == 3) {
        #{  1*1 MATRIX ?


        # } else if (mg[i,methodindiv[1]] == 4)   {
        # m2 careful as matrix needs to be vertical , test change 17/03
        # this works for accupuncture data but not asthma, because when one patient raw daa is 1 by n matrix so m2 needs to be horizontal, ie not a matrix

        #19/03  11/04 works for CR but then not for J2R for 5333 patt=7 one patient
        if (mg[i,"X1"] == 1) {
            meanval = (m2) + as.matrix(raw1 - m1)%*%as.matrix(t_mimix)
        }  else {
        meanval = as.matrix(m2) + (as.matrix(raw1 - m1)%*%as.matrix(t_mimix))
         }


        #24/02browser()
        U <- chol(conds)
        # mg[i,X1] is equiv to Stata counter, miss_count is no. of missing, so
        miss_count=rowSums(mata_miss)

        # gen erate inverse normal
        Z<-stats::qnorm(matrix(stats::runif( mg[i,"X1"]* miss_count,0,1),mg[i,"X1"],miss_count))
        # check same input parameters for inverse norm gen as in stata

        #for debug 20/01
        #  print(paste0('mg[i,X1] =',mg[i,"X1"],' miss_count= ',miss_count))

        #mata_y1 = meanval+Z%*%t(U)  14/01/20 try without transpose because Stata Cholesky lower tri, R upper (or vice versa)
        mata_y1 = meanval+Z%*%(U)

        #define new matrix from observed,  id column  (the last)
        mata_new <- preraw

        # assigning the columns where the missing values

        # was nct (as in stata) = no ntimes + ncovar
        # so if no missing then just copy full values into mata_new columns
        #if(length(c_mata_miss)==0 ) { mata_new[,c(1:length(tst2))] <- mata_Obs[,c(2:length(tst2))]
        if(length(c_mata_miss)==0 ) { mata_new[,c(1:length(tst2))] <- preraw[,c(1:length(tst2))]
        }else
          {
          mata_new[,c_mata_miss] <- mata_y1
        }
        # if(length(c_mata_miss)!=0 )

        #save U,Z for imp(m) and patt(i)
        # this worked for Z



        #assuming this  from Stata if "`interim'"==""{
        # SNO just id col
        SNO <- mata_Obs[c(startrow:stoprow),1]
        #SNO <- mata_ObsX[,ncol(mata_Obs)]
        # GI treatment grp column 1 (here),II imputation number col, mata_new matrix then SNO is id col.
        GI <- array(data=mg[i,1],dim=c(mg[i,"X1"],1))
        #II  no imputations
        II <- array(data=m,dim=c(mg[i,"X1"],1))

        #doesnt need SNO, as id already in
        #mata_new<-cbind(GI,II,mata_new,SNO)

        #this works but bette to pre initialise data structure outsidr loop
        mata_new<-cbind(GI,II,mata_new,SNO)


        #assume delta to be used if specified in input argument
        if (length(delta != 0) ) {
          #browser()
        mata_new <-  AddDelta(tst2,  ncovar_i,mata_new,delta)
        }

        mata_all_newlist[[m_mg_iter]]=mata_new
        #mata_all_new<-rbind(mata_all_new,mata_new)
        #mata_all_new<- na.omit(mata_all_new)
        #stata equilv \ is row bind NOt cbind!!?

        #9/5/20 end else c_mata_nonmiss ==0 here
        }


      } # ( m in  1:M) so insert delta module just before this
    } #for row[mg]

  } #for M StOP HERE!!


  # 23/03 addon to save as implist
  #  browser()
  #  dimlist <- (nrow(mg[[1]])*M)

  # extract from nested list
  # combine into data set containing M imputed datasets
  #  mata_all_newData1x <- do.call(rbind,mata_all_newlist[[1]])
  # then sort (by imputation and patient id) into M data sets and split into M lists
  #  impdatasets <- mata_all_newData1x[order(mata_all_newData1x$II,mata_all_newData1x$SNO),]
  #  impdatasets2 <- mata_all_newlist[order(mata_all_newlist$II,mata_all_newlist$SNO),]
  impdataset<-getimpdatasets(list(mata_all_newlist,mg,M,meth))

  #return(list(mata_all_newlist,mg,M,meth))
  return(impdataset)
} # for mimix test

#' @title getimpdatasets
#' @description to obtain the M imputed data set from the output list into one dataset
#' @details This combines the imputations found from the M pattern groups
#' @param varlist  list of data containing imputed values from the M pattern groups
#' @return impdataset




getimpdatasets <- function(varlist){
  #9/5/20
 #browser()
  # to obtain M imputed data sets
  # dimension of data set, nrows in pattern times no imputations,
  # note sub data sets wi have different cols if completely missing so
  mata_all_newlist<-  varlist[1]
  mg<-(varlist[2])
  M<- unlist(varlist[3])
  meth<- unlist(varlist[4])

  #dimlist <- (nrow(mg[[1]])*M)

  # extract from nested list
  # combine into data set containing M imputed datasets
  mata_all_newData1x <- do.call(rbind,mata_all_newlist[[1]])
  # then sort (by imputation and patient id) into M data sets and split into M lists
  impdatasets <- mata_all_newData1x[order(mata_all_newData1x$II,mata_all_newData1x$SNO),]


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


  #impdatasets <- as.mids(tmpdatasets, .id="SNO",.imp="II")


  return(impdatasets)
}



# so no almost reday for mice!
#fit<-with(data= as.mids(impanticausalun, .id="SNO",.imp="II"), expr = lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+POOLED.INVESTIGATOR+PATIENT.SEX))
#summary(pool(fit))

