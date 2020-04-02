

#' @title Runmimix
#' @description main wrapper for running mimix
#' @details This is based on Suzie Cro's Stata program
#' @export Runmimix
#' @param data  datset in wide format
#' @param covar covariates and base depvar must be complete
#' @param depvar dependent variable
#' @param treatvar treatment group , recoded to 1,2,..
#' @param idvar patient id
#' @param timevar time point for repeated measure
#' @param M  number of imputations
#' @param refer  reference group for j2r,cir,cr mrthods
#' @param meth RBI method
#' @param seedval  seed value to obtain same outputs
#' @param priorvar  prior,defualt Jeffries, also uniform or ridge
#' @param burnin  burnin value
#' @param bbetween  value between iterations in mcmc
#' @param methodindiv  columns in dat for individual method and reference group
#' @return impdataset
#' @examples
#' impdatasets <- Runmimix("acupuncture",c("head_base","sex","age"),"head","treat","id","time",10,1,"J2R",201,c("jeffreys"),1000,NULL,NULL)
#' impdataset<- Runmimix("asthma",c("base"),"fev","treat","id","time",100,1,"J2R",101,"jeffreys",1000,NULL,NULL)

Runmimix<- function(data,covar,depvar,treatvar,idvar,timevar,M=1,refer,meth,seedval=101,priorvar,burnin=1000,bbetween=NULL,methodindiv) {
  #browser()
  # establish whether indivual or group by cretung a flag var
  if (!is.null(meth)) { flag_indiv <<-0 }
  if (is.null(meth) & !is.null(methodindiv[1]) ) {flag_indiv<-1 }


  #
  #find no covars
  ncovar_i = length(covar)
  #browser(1)
  print(paste0("ncovar_i = ",ncovar_i))
  print(paste0("covar= ",covar[1]))


  #if prior =ridge then need fing default ridge parameter ,try 0.1sd(depvar)
  #browser()
  #recode refer to number within treatvar if not already numeric/integer no need for this now treat recoded previously
  # if ((class(mxdata[,treatvar])!="integer" & class(mxdata[,treatvar])!="numeric") )  {refer = which(unique(mxdata[,treatvar])==refer) }
  #recode treatvar to number, doesnt like this
  #for  (i in 1:length(unique(mxdata[,treatvar])) ) {mxdata[,treatvar]=i  }
  # danger to use transform
  #transform(mxdata,treatvar=as.numeric(treatvar) )
  # if char convert to numeric, if numeric then leave
  #  if (class(mxdata[,treatvar])=="char") {mxdata[,treatvar]<-as.numeric(mxdata[,treatvar])}
  #treatvar <- as.numeric(treatvar)
  # treatvar and refer convert to numeric and check refer in treatvar but as.nemeric not appropriate here
  # need to recode instead
  #treatvar<- as.numeric(treatvar)
  #refer <- as.numeric(refer)
  #set.seed(101)

  # to recode refer, need use a lookup table



  set.seed(seedval)
  #set.seed(unlist(tail(kmargs,n=1)))

  #mxdata<- readdata(data)


  #do not read in last elemnt LAST otherwise get an unused argument error msg

  #this outputs the summary pattern table
  #testlist = do.call( preprodata,kmargs[-length(kmargs)])
  #seed, prior not used
  # testlit = do.call( preprodata,kmargs[c(-(length(kmargs)-1),-(length(kmargs)))])
  # 15/03browser()

  if (!is.null(meth) ) {
    testlist = do.call( preprodata,list(data,covar,depvar,treatvar,idvar,timevar,M,refer,meth))
    refer <- testlist[[9]]

    # stopifnot(refer %in% ntreat)
    #no need for this?
    #meth <- testlist[[10]]
    #need to recode meth
    meth<-  ifelse( ( meth=="j2r" | meth=="J2R" |meth=="j2R"|meth=="J2r" ),3,
                    ifelse( ( meth=="CR" | meth=="cr" |meth=="Cr"|meth=="cR" ),2,
                            ifelse( ( meth=="MAR" | meth=="mar" |meth=="Mar"|meth=="MAr"|meth=="Mr"|meth=="MR" ),1,
                                    ifelse( ( meth=="CIR" | meth=="cir" |meth=="CIr"|meth=="cliR" ),4,
                                            ifelse( ( meth=="LMCF" | meth=="lmcf" |meth=="Last"|meth=="last" ),5,9)))))
    #ifelse(mxdata$methodvar=="j2r",3,ifelse(mxdata$methodvar=="cir",4,9))

  }
  else if (!is.null(methodindiv[1]) ) {
    testlist = do.call( preproIndivdata,list(data,covar,depvar,treatvar,idvar,timevar,M,refer,meth,methodindiv))
    # need to re-set meth for individual
    meth <- testlist[[9]][1]

  }

  #browser()
  # find sd depvar


  # 10/02/20
  #testlist = preprodata(kmargs)
  #browser()
  # returns list from preprodata function
  ntreat<-sort(unlist(testlist[[3]]))
  #sort it !!


  #ntreat <-testlist$ntreat
  #stopifnot()
  finaldatS<-testlist[[2]]
  #finaldatS <- testlist$finaldatS
  mg<-testlist[[4]]
  #mg <- testlist$ex1s
  # vital to get the mata_obs correctly sorted! so corresponds with mimix_group lookup
  # to be consistent with Stata move the base col after the fevs!
  #mata_Obs <- testlist$finaldatS
  mata_Obs <- testlist[[2]]
 # M <- testlist[[8]]


  # so tsts2 is th equiv to mi_impute so try to put in norm2
  # est tis tomorrow 17/01
  #tst2 <- mi_impute("id","time","fev","base")
  #list("fev","treat","id","time","base",10,2,"J2R",101)
  #pre move, tst2 <- mi_impute(kmargs[[3]],kmargs[[4]],kmargs[[1]],kmargs[[5]])
  #browser()
  #tst2 <- mi_impute(kmargs[[4]],kmargs[[5]],kmargs[[2]],kmargs[[1]])
  tst2 <- mi_impute(data,idvar,timevar,depvar,covar)
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

  #paramMatrix<-matrix(,nrow=(nrow(ntreat)*M),ncol=2)

  #create  emptylist for each treat and multiple m's
#browser()
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

  #try ser up as many Result data files as treatments instead of one big file?.
  # browser()

  for (val in 1:length(ntreat)) {

    kmvar=get(paste0("prenormdat",val))
    prnormobj<-assign(paste0("prnormobj",val), subset(kmvar, select=c(tst2)))
    #create  emptylist for each treat

    paste0("Looping for treatment = ",val," performing mcmcNorm for m = 1 to ",M)
    for(m in 1:M) {

      # supppress warnings regarding solution  near boundary, see norm2 user guide, also mimix about this problem

      #need to error check when ridge or invwish used, the accompanying parameter values supplied.

      # find sd of depvar over all times?
      sd_depvar<- sd((get(data)[,depvar]),na.rm=TRUE)

      if ( priorvar[1] == "ridge" ) { if ( is.na(priorvar[2])) { priorvar[2]<-(sd_depvar*0.1) }}
      # invwish not implemented!
      if ( priorvar[1] == "invwish" ) { stopifnot(priorvar[2]>0 & priorvar[3]>0 ) }

      # emResultT<-(emNorm(prnormobj,prior = priorvar[1],prior.df=priorvar[2],prior.sscp=priorvar[3]))
      #  mcmcResultT<- (mcmcNorm(emResultT,iter=1000,multicycle = NULL,prior = priorvar[1],prior.df = priorvar[2],prior.sscp=priorvar[3]))
      # browser()
      #if prior =ridge then need fing default ridge parameter ,try 0.1sd(depvar)


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

  # try and use lappy instead of loop for M


  # can repeat interactively from here
  # now loop over the lookup table mg, looping over every pattern - make sure mata_Obs sorted same way!

  # declare iterate for saving data
  m_mg_iter<-0
  for (i in 1:nrow(mg))
  {
    # define mata_miss as vector of 1's denoting missing using col names ending i ".missing"
    # this section to be amended to cope with multiple covariates

    mata_miss <- mg[i,grep("*..missing",colnames(mg)),drop=F]
    #mata_miss <- mimix_group[i,c(2,3,4,5)]     #defimne mata_miss
    #assumes covariate non missing
    # now >1 covariate code must be amended

    #find no rows to create covar vector to cbind with mg
    numrows<- nrow(mata_miss)

    #mata_miss$covar.1 <-0                       #assuming cov col is non missing
    #browser(2)
    # so create matrix with names as covariates.1 and 0 row to signify no missing values
    covarsdf<- provideDimnames(matrix(0:0,ncol=length(covar)),base=list(paste0(""),paste0(covar,".missing")) )
    mata_miss <- cbind(mata_miss,covarsdf)


    mata_nonmiss <- (ifelse(mata_miss==0,1,0))  #define mata-nonmiss from miss

    # need transform nonmiss,miss to c lists - ie. index the
    c_mata_miss<-which(mata_miss==1)
    c_mata_nonmiss<-which(mata_nonmiss==1)
    # eg no missing is c(1,2,3,4,5)


    # count of pattern by treatment
    cnt<- mg$X1[i]

    # treatment grp
    trtgp<- mg$treat[i]

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
        mata_new <- mata_Obs[c(st:en),2:ncol(mata_Obs)]
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

            Sigma<<-Sigmatrt
          }
          # 'J2R'
          else if (meth == 3 ) {

            # changed saving the result into  just the param file, list of 2 so can use list index here
            #treatmnets are 1.. M then M+1 ..2M .. etc
            mata_means_trt <- paramBiglist[[M*(trtgpindex-1)+m]][1]
            mata_means_ref <- paramBiglist[[M*(referindex-1)+m]][1]


            #  browser()
            # below causes error after using >1 covars and mata_nonmiss has covar.1, not proper covar names
            mata_means_t <- unlist(mata_means_trt)*mata_nonmiss
            # print(paste0("mata_means_trt, mata_nonmiss= ",mata_means_trt,mata_miss))

            mata_means_r <- unlist(mata_means_ref)*mata_miss
            # so when all missing  1,1,1, ... then all contributing comes from reference means
            mata_means <- mata_means_r+mata_means_t
            # and preserve names
            colnames(mata_means) <- colnames(mata_means_trt)
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

            Sigma<<-SigmaRefer
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

            Sigma<<-SigmaRefer
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

            Sigma<<- SigmaRefer
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

            Sigma<<- Sigmatrt
          }  #if meth

          ############# individual analysis #########################
        }
        #else if(!is.null(methodindiv[1]))
        else if (flag_indiv==1) {

          # call function for  indiv

          mata_means <- ifmethodindiv(methodindiv,mg,m,M,paramBiglist,i,treatvar,c_mata_nonmiss,c_mata_miss,mata_miss,mata_nonmiss)

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
        preraw <- mata_Obs[c(startrow:stoprow),2:ncol(mata_Obs)]
        raw1 <- preraw[,c_mata_nonmiss]

        ##### try inserting routine for all missing values here NOt really necessary!! ### 20/1/20 #####
        ## when all missing data, the length of c_mata_nonmiss must exclude the no. of covariates
        ## temporary fix set as 1, ie no. of covars 20/01/20
        # nOTE think no observed data includes no base line as well?
        if (length(c_mata_nonmiss)==0)  {
          ## routine copied from mimix line 1229
          U <- chol(Sigma)

          # generate inverse normal, same as used below
          #for debug 20/01
          #   print(paste0('mg[i,X1] =',mg[i,"X1"],' miss_count= ',miss_count))
          Z <- qnorm(matrix(runif( mg[i,"X1"]* miss_count,0,1),mg[i,"X1"],miss_count))

          mata_y1 = meanval+Z%*%(U)
          #set dimensions mata_new to mata_y1
          mata_new <- dim(mata_y1)
          GI <- array(data=mg[i,"treat"],dim=c(mg[i,"X1"],1))
          #II  no imputations
          II <- array(data=m,dim=c(mg[i,"X1"],1))
          mata_new=cbind(GI,II,mata_new,SNO)
          mata_all_newlist[[m_mg_iter]]=mata_new

        }


        #so S12 must be declare as a matrix ! as number otherwise is class number
        #try transpose as was failing compatible error
        #tS12<-t(S12)
        mS11<-as.matrix(S11)
        mS12 <-as.matrix(S12)
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
        #17/03
        #   browser()
        # try change 16/03
        #meanval <- as.matrix(m2)+as.matrix(raw1 - m1)%*%as.matrix(t_mimix)

        # } else if (mg[i,methodindiv[1]] == 4)   {
        # m2 careful as matrix needs to be vertical , test change 17/03
        # this works for accupuncture data but not asthma, because when one patient raw daa is 1 by n matrix so m2 needs to be horizontal, ie not a matrix

        #19/03    if (mg[i,"X1"] == 1) {
        #    meanval = (m2) + as.matrix(raw1 - m1)%*%as.matrix(t_mimix)
        # }  else {
        meanval = as.matrix(m2) + (as.matrix(raw1 - m1)%*%as.matrix(t_mimix))
        #  }


        #24/02browser()
        U <- chol(conds)
        # mg[i,X1] is equiv to Stata counter, miss_count is no. of missing, so
        miss_count=rowSums(mata_miss)
        #  print("miss_count = ")
        #  print (miss_count)
        # gen erate inverse normal
        Z<-qnorm(matrix(runif( mg[i,"X1"]* miss_count,0,1),mg[i,"X1"],miss_count))
        # check same input parameters for inverse norm gen as in stata

        #for debug 20/01
        #  print(paste0('mg[i,X1] =',mg[i,"X1"],' miss_count= ',miss_count))

        #mata_y1 = meanval+Z%*%t(U)  14/01/20 try without transpose because Stata Cholesky lower tri, R upper (or vice versa)
        mata_y1 = meanval+Z%*%(U)

        #define new matrix from observed,  id column  (the last)
        mata_new <- preraw
        #10/03/20
        #browser()
        # assigning the columns where the missing values
        #17/03  browser()
        if(length(c_mata_miss)==0 ) { mata_new <- mata_Obs[,c(2:nct+1)]
        }else {
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

        mata_all_newlist[[m_mg_iter]]=mata_new
        #mata_all_new<-rbind(mata_all_new,mata_new)
        #mata_all_new<- na.omit(mata_all_new)
        #stata equilv \ is row bind NOt cbind!!?

        #}
      } #if patt ==0r if all  missing
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
} # for runmimix test

#' @title getimpdatasets
#' @description to obtain the M imputed data set from the output list into one dataset
#' @details This combines the imputations found from the M pattern groups
#' @param varlist  list of data containing imputed values from the M pattern groups
#' @return impdataset




getimpdatasets <- function(varlist){
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
  # to get the list
  # implist1x <- split(impdatasets,impdatasets[,"II"])
  return(impdatasets)
}


#' @title mi_impute
#' @description  create a vector defining imputation variables
#' @details similar to mi_impute in mimix Stata programx
#' @param data data set
#' @param idvar patient id
#' @param timevar time point for repeated measure
#' @param depvar dependent variable
#' @param covar covariates and base depvar must be complete
#' @return list of variables




mi_impute <-function(data,idvar,timevar,depvar,covar) {
  #preprodata(depvar,treatvar,idvar,timevar,covar)
  tst<-reshape(get(data)[,c(idvar,depvar,timevar)],v.names = depvar,timevar = timevar,idvar=idvar,direction="wide")
  #tst<-tidyr::pivot_wider(get(data),id_cols=c(idvar),names_from=timevar,names_prefix=depvar,values_from=depvar)
  # add ont covariates and add on to resp list
  respvars<-c(names(tst[,-1]),covar)
  return(respvars)
}
