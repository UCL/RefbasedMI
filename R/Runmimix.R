#' @title mimix
#' @description main wrapper for running mimix
#' @details This is based on Suzie Cro's Stata program
#' @details sets up a summary table based on missing data pattern- mg  mimix_group                                                                   
#' @details reflects the pattern and treatment group configuration of the raw data
#' @details then acts as a looping mechanism, norm2 is used as MCMC multivariate normal 
#' @export mimix
#' @param data  Dataset in wide format
#' @param covar Covariates - may include the baseline value of depvar. Must be complete (no missing values).
#' @param depvar Dependent (outcome) variable
#' @param treatvar Treatment group, coded 1,2,..
#' @param idvar Participant id
#' @param timevar Time point for repeated measure
#' @param M Number of imputations to be created
#' @param reference  Reference group for J2R, CIR, CR methods
#' @param method Reference-based imputation method: must be ...
#' @param seed  Seed value. Specify this so that a new run of the command will give the same imputed values.
#' @param prior  Prior for the variance-covariance matrix when fitting multivariate normal distributions. Jeffreys (default), uniform or ridge
#' @param burnin  Number of burn-in iterations when fitting multivariate normal distributions.
#' @param bbetween  Number of iterations between imputed data sets when fitting multivariate normal distributions.
#' @param methodvar   vector designating variable in data-set specifying individual method 
#' @param referencevar   vector designating variable in data-set specifying reference group corresponding to individual method
#' @param delta vector of delta values to add onto imputed values (non-mandatory) (a's in Rogers paper),length as number of time points
#' @param dlag vector of delta values to add onto imputed values (non-mandatory) (b's in Rogers paper),length as number of time points
#' @param K0 Causal constant for use with Causal method
#' @param K1 exponential decaying Causal constant for use with Causal method
#' @param mle logical option to Use maximum likelihood parameter estimates instead of MCMC draw parameters
#' @return impdataset the M imputed data-sets appended to the "missing values" data-set in wide format
#' @examples
#' \dontrun{
#' mimix("asthma",c("base"),"fev","treat","id","time",5,1,"J2R",,,,,,,,,1,0.5,)
#' fit<-with(data= as.mids(impdataCausal), expr = lm(fev.12~treat+base))
#' summary(pool(fit))
#' }

mimix<- function(data,covar=NULL,depvar,treatvar,idvar,timevar,M=1,reference=NULL,method=NULL,seed=101,prior="jeffreys",burnin=1000,bbetween=NULL,methodvar=NULL,referencevar=NULL,delta=NULL,dlag=NULL,K0=1,K1=1,mle=FALSE) {
  # 6/11 try account for interims J2R MAR
  # if testinterims then want method to 1stly be MAR
  # this forces interims to be estimated as MAR by default
  testinterim<-1
  #browser(text="2912")
  if (class(mle) !="logical" & mle !=0 & mle !=1 ) { stop("mle must be logical value") } 
    
  
  # 29/10 change mor user friendly
  # stopifnot(class(get(data)) == "data.frame")
   if (class(get(data)) != "data.frame") {stop("data must be type dataframe")}
  
  # insert error checks  HERE
  #check not both meth and methodundiv specified.
  # stopifnot(method=="NULL" | methodvar=="NULL")
  if (!(is.null(method) | is.null(methodvar))) {stop("Either method or methodvar must be specified but NOT both") }
  # establish whether specifying individual or group by creating a flag var
  if (!is.null(method)) {
    flag_indiv <-0
  # Causal constant must be number and check it exists if Causal specified
  #stopifnot(meth=="Causal" & !missing(Kd))
    
  if (toupper(method)=="CAUSAL" | toupper(method)== "CASUAL" | toupper(method)== "CUASAL") {
    if (missing(K0))  {stop("K0 Causal constant not specified")}
    if (missing(K1) & (K0!=0) )  {stop("K1 Causal constant not specified")}
    #if (!(K0>=0 & K0<=1)) {warning("K1 Causal constant not in range 0..1 ")}
    if (K0<0) {warning("K0 Causal constant negative.. ")}
    if (K0>1) {warning("K0 Greater than 1.. ")}
    if (!(K1>=0 & K1<=1)) {stop("K1 Causal constant not in range 0..1 ")}
  } #Causal constant must be number
    K0<- as.numeric(K0)
    K1<- as.numeric(K1)
  }
  if (is.null(method) & !is.null(methodvar) ) {flag_indiv<-1 }

  #find no covars
  ncovar_i = length(covar)
  # if covars exist
  # check that covars are complete AND integer !!
  # NO  complete but could be factors 
  # if covar=NULL  17/05/20   #discuss.analyticsvidyha.com
  #browser()
  if (length(covar)!=0) {
   # stopifnot(sum(!stats::complete.cases(get(data)[,covar]))==0)
    if (sum(!stats::complete.cases(get(data)[,covar]))!=0) {stop("covariates not complete!!")}
    stopifnot( (sapply((get(data)[,covar]), is.factor)) | (sapply((get(data)[,covar]), is.numeric)) )
    
    # really should be integer
    # if need to convert
    # test<-as.data.frame(apply((get(data)[,covar]),2,as.integer) )
  }


  # note if no covar then treat the first depvar level as a covar , eg covar = fev.0.
  # sO must preform small routne tO transform datA set INTO baseline covars.

 # build in checks specifically  for delta
 
  ntimecol<- get(data)[c(timevar)]
  ntime<-nrow(unique(ntimecol))
  # only run if delta specified
  if (!is.null(delta)) {
   
     stopifnot(length(delta)==ntime)
 
     #set dlag to default if 1 1 1 ...if NULL
     if (is.null(dlag)) {
       dlag <- rep(1,length(delta))
       } else  {
       stopifnot(length(dlag)==ntime)
       }   
  }
    

  
  set.seed(seed)

#21/05
  # assign all characters in meth to UPPER case
  # check i1st if meth exists! if statement needed otherwise meth change form Null to character(0)
  if (!is.null(method) | !length(method)==0 ) {
  method <- toupper(method)
  stopifnot( (method == "MAR" | method=="MR"|
              method=="J2R" | method=="J2" | method=="JR"|
              method=="CIR" | method=="CLIR" |
              method=="CR" |
              method=="LMCF" | method=="LAST"  |
              method=="CAUSAL" | method== "CASUAL" | method== "CUASAL"),

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


  if (!is.null(method) ) {
    testlist = do.call( preprodata,list(data,covar,depvar,treatvar,idvar,timevar,M,reference,method))
#browser(text="0720")
        reference <- testlist[[7]]



    # stopifnot(refer %in% ntreat)
    #no need for this?
    #meth <- testlist[[10]]
    #need to recode meth
    method<-  ifelse( (  method=="J2R" |method=="J2"|method=="JR" ),3,
                    ifelse( ( method=="CR"  ),2,
                            ifelse( ( method=="MAR" | method=="MR" ),1,
                                    ifelse( ( method=="CIR" |method=="CLIR" ),4,
                                            ifelse( ( method=="LMCF" | method=="LAST" ),5,
                                              ifelse( ( method=="CAUSAL" | method=="CASUAL" | method=="CUASAL"),6,9))))))


  }
  # for user specified
  else if (!is.null(methodvar) ) {
    testlist = do.call( preproIndivdata,list(data,covar,depvar,treatvar,idvar,timevar,M,reference,method,methodvar,referencevar))
    # need to re-set meth for individual
    #meth <- testlist[[8]][1]
    
  }


#browser(text="0312")

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

  #browser()
  # extract treatvar and methodindiv vars
  # fix 12/06 
  # Obs_treat<-subset(mata_Obs,select=c(treatvar,methodvar))
  Obs_treat<-mata_Obs[,c(treatvar,methodvar,referencevar)]
  mata_ObsX<- mata_Obs[,!(names(mata_Obs) %in% c(treatvar,methodvar,referencevar))]
  # combine back the extracted cols
  mata_Obs <- cbind(mata_ObsX,Obs_treat)
  #set name for treatvar otherwise defaults to Obs_treat 
  # this ok when methodvar null and tested ok when not null 12/06 names just become methodva and referncevar 
  names(mata_Obs)[names(mata_Obs)=="Obs_treat"]<-treatvar




  tst<-stats::reshape(get(data)[,c(idvar,depvar,timevar)],v.names = depvar,timevar = timevar,idvar=idvar,direction="wide")

  #1612 make readable in pass2
  tst2<<-c(covar,names(tst[,-1]))



  # put commas in
  tst3<-paste(tst2,collapse = ",")

  #create input data sets for each tment from which to model
  #browser()
  for (val in 1:length(ntreat)) {
    #print(paste0("prenormdat",val))
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

    cat(paste0("\n treatment = ",val,"\n performing mcmcNorm for m = 1 to ",M) )
    for(m in 1:M) {

      # supppress warnings regarding solution  near boundary, see norm2 user guide, also mimix about this problem

      #need to error check when ridge or invwish used, the accompanying parameter values supplied.

  #browser()

      if ( prior[1] == "ridge" ) {
        # find sd of depvar over all times for ridge option
        sd_depvar<- stats::sd((get(data)[,depvar]),na.rm=TRUE)
        if ( is.na(prior[2])) { prior[2]<-(sd_depvar*0.1) }}
      # invwish not implemented!
      #if ( priorvar[1] == "invwish" ) { stopifnot(priorvar[2]>0 & priorvar[3]>0 ) }

      # statsgeek suggestion prof of concept ,will have to change emNorm line
 # browser(text="2701")
 # test jomo 2601

      # mle false or 0, true or 1     
  
      if (mle==FALSE) { 
      # doesnt suppress msgs capture_condition(emResultT<-(norm2::emNorm(prnormobj,prior = priorvar[1],prior.df=priorvar[2])) )
     invisible(capture.output(emResultT<-(norm2::emNorm(prnormobj,prior = prior[1],prior.df=prior[2])) ))
      #mcmcResultT<- (mcmcNorm(emResultT,iter=1000,multicycle = NULL,prior = priorvar[1],prior.df = priorvar[2]))
     mcmcResultT<- (norm2::mcmcNorm(emResultT,iter=burnin,multicycle = bbetween,prior = prior[1],prior.df = prior[2]))
        # try for when using mle!  
        # mcmcResultT<- emResultT 
        # if using jomo
     #setnburn=1000
    #invisible(capture.output(testimp<- jomo::jomo.MCMCchain(prnormobj,nburn=burnin,meth=common, output=0))) 
      
     }
      else {
      invisible(capture.output(emResultT<-(norm2::emNorm(prnormobj,prior = prior[1],prior.df=prior[2])) ))  
        # for mle
      mcmcResultT<- emResultT  
     # mcmcResultT <- norm2::impNorm(emResultT, method="predict")
      }
      
      
      
      # msg from emNorm
      #Note: Finite-differencing procedure strayed outside
      #parameter space; solution at or near boundary
      #OCCURRED IN: estimate_worst_frac in MOD norm_engine

      # cumiter needs to be used as greater than M after 1st treatment 
      cumiter<-cumiter+1
      #browser(text="1412")
      # to see if can save directly to val indice
      # introduce here now using jomo  replace setnburn by burnin
      #mcmcResultT<- list(testimp[[2]][,,setnburn],testimp[[3]][,,setnburn])
      #mcmcResultT$param<- list(testimp[[2]][,,setnburn],testimp[[3]][,,setnburn])
      
      paramBiglist[[ cumiter]] <- mcmcResultT$param
      assign(paste0("paramBiglist",val,"_",m), mcmcResultT$param)
    }
    cat(paste0("\n mcmcNorm Loop finished, m = ",M,"\n"))
  }

  #store paraBiglist in a single structure

  # ) # system.time
  #print(paste0("mcmcNorm Loop finished, m = ",M))




  # can repeat interactively from here

  ############################################ start big loop #########################################
  # now loop over the lookup table mg, looping over every pattern - make sure mata_Obs sorted same way!

  # declare iterate for saving data
  cat("\n Starting imputation")
  #initialise interim 
  interim<-0
  # construct structure to savde interim ids but this get reinitialised to many times!
  interim_id<- mata_Obs[c(mg[1,1]),"id"] 
  #interim_pos <- c(0,0,0,0,0)
  
  #18/11 see if rawplusinterim works here!
  #browser()
  
  m_mg_iter<-0
  for (i in 1:nrow(mg))
  {

    # define mata_miss as vector of 1's denoting missing using col names ending i ".missing"
    # this section to be amended to cope with multiple covariates

    #browser(text="0101")
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
    cnt<- mg$cases[i]
#browser()
    # treatment grp
    trtgp<- mg[i,treatvar]

    pattern <- mg$patt[i]
    #cat("\ntrtgp = ", trtgp)
    

 #browser()
    if  (!is.null(method) ) {
      cat("\n",treatvar ," = ", trtgp,"patt = ",pattern,"number cases = ", cnt) }
        else if(!is.null(methodvar) ) {
      cat("\n",treatvar ," = ", trtgp,methodvar," = ",as.character(mg[i,methodvar]),referencevar," = ",as.character(mg[i,referencevar]),"patt = ",pattern,"number cases = ", cnt)
    }
# browser(text="0501")
    #need to convert (relate) treatment group to position in ntreat (create Pindex vector)
    #unnecessary now recoded
 #browser()
    trtgpindex<-which(trtgp==ntreat)
    referindex<-which(reference==ntreat)
  
        # multiple  simulations start here within the pattern loop #########
    for ( m in  1:M)  {
      #*FOR INDIVIDUALS WITH NO MISSING DATA COPY COMPLETE DATA INTO THE NEW DATA MATRIX mata_all_new `m' TIMES
      #if `pat' == 0{
      m_mg_iter<-m_mg_iter+1
      
      # construct structure to savde interim ids but this get reinitialised to many times!
      #11/11 interim_id<- mata_Obs[c(mg[i,1]),"id"] 
      
      
      # if no missing values
      if(length(c_mata_miss)==0 ) {
        # start and end row positions
        st<-mg[i,"cumcases"]-mg[i,"cases"]+1
        en <-mg[i,"cumcases"]
        # id (SNO) is 1st col try changed 0812
         SNO<-mata_Obs[c(st:en),1]
        #SNO <- mata_Obs[c(st:en),".id"]
        #21/07/20
        #browser()
        # this doesnt delete treat
        #mata_new <- mata_Obs[c(st:en),2:ncol(mata_Obs)]
         # duplicate col headbasemiss.1?
     #    browser(text="0301")
        mata_new<-mata_Obs[c(st:en),!(names(mata_Obs) %in% c(idvar))]
        #browser() # treat defined within fun
        GI <- array(data=mg[i,treatvar],dim=c(mg[i,"cases"],1))
        #II  no imputations
        II <- array(data=m,dim=c(mg[i,"cases"],1))
        mata_new=cbind(GI,II,mata_new,SNO)
        #names(mata_all_new)<-names(mata_new)
       # browser(text="0201")
        mata_all_newlist[[m_mg_iter]]=mata_new
        # mata_all_new=rbind(mata_all_new,mata_new)
        # else if there are missing values
      } else {
        # need to distinguish between meth and methodindiv

        if (flag_indiv==0 ) {

          referindex<- reference
          #FOR INDIVIDUALS WITH  MISSING DATA  `m' TIMES
          # dependent on method chosen
          # 'MAR'
          #browser()
          
          #2611 if testinterim cale then want o proces interims as MAR regardles of value of method
          if (method== 1 | testinterim==1)  {
            
            # 11/11 checking methd used for interims in J2R works as MAR 

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

            
            #11/11 flag interims put in as option in command
           # testinterim<-1
            
            #initialise interim 
           # interim<-0
            
            # is this necessary??
            if (testinterim==1) {
            # need to set method to MAR for 1st pass not works here try above
            #2611method<-1  
            # need to find the interim cols so as to set values as MAR
            miss_count <- length(c_mata_miss)
    
            # cndition if miss_count =1 then only 1 missing  and if Not at end eg c_mata_miss = 5 then by defn is an interim!
            #browser() 
         #   browser(text="0101") # interim at ... this case prog to fail
            #tryng to put sigma from MAR in ,works for 533 but not 5051 but once did before putting in j2r sigma at bottom, 
            # MAR sigma gives correc tinterim but incorrect subsequents  
            # 1st phase to find interims and apply MAR sigmas
         #   browser(text="0101") 
        #  browser(text="0401")
            if ((miss_count == 1) & (c_mata_miss[1] < length(mata_means)) ) # 1 missing and not at end point 5115 
            {
              interim<-1
                    #     cat(paste0("miss+_count=1 paramBiglist= ?"))
            } else if (miss_count >1)  # if more than 1 missing as in acupunt data 
             
            {
              #assuming 1st col non-missing so start at  5333 is failing xxx0 234  being treated as a J2R!
              for (b in 2:(miss_count)) {    # last miss_count end pt check separately below 
                #              # c_mata_miss vector of missing locations so if next col not missing then will be gap in sequence 
                # c_mata_miss[b] is interim so overwrite with MAR value
                # adjust when only last col missing (patt=8) ie c_mata_miss = 5  
                #  245 => b = 1,2,3
                if ((c_mata_miss[b-1]+1) != c_mata_miss[b])  # so next entry after c_mata_miss[b] is non-missing
                  # so interim found but what if next one missing, need to check 
                { 
                  #browser() MAR sigma
                  interim<-1
                 
                  #browser()
                  #cat(paste0("interim at ",m_mg_iter ,"SNO=",SNO, "b=",b  ))
                  if (m==1) {
                  #cat(paste0("\ninterim at ",m_mg_iter , "b=",b , "id= ", mata_Obs[c(mg[i,"cumcases"]),"id"] ))
                  # in case more than 1 interim in patt group
                    for (it_interim in 1:mg[i,"cases"]) {
                      cat(paste0("\ninterim at id= ", mata_Obs[c(mg[i,"cumcases"]-mg[i,"cases"]+it_interim),idvar] ))      
                      #cat(paste0("\ninterim at id= ", mata_Obs[c(mg[i,"cumcases"]),"id"] ))
                      }
                  }  
                    # construct vector to save interims ids
               #11/11   interim_id <- rbind(interim_id, mata_Obs[c(mg[i,"cumcases"]),"id"] )  
                
                   
                  
                } else 
                # note that if all missing then wont be interims and c_mata_nonmiss is integer(0) ie NUL
          if (length(c_mata_nonmiss)!=0) {(c_mata_miss[b-1]+1 == c_mata_miss[b]) & ( c_mata_miss[b-1] < max(c_mata_nonmiss)) }   # need to check there is a non-missing to the right, eg 23 5 , 234 all interims!  
                  #  cat(paste0("check nonmissing to right")) trying to catch 5333 7/ this s not going to affect when b=miss_count so need add for this condition   
                { interim<-1
                } #need to include outside the for loop  when condition b=miss_count
              } #if
            } #for
            #  else if doesnt work need to process the miss_count element because not processed in for loop  
            #else if ( (b==miss_count) & ( c_mata_miss[b] < length(mata_means)) )
            if  ((c_mata_miss[miss_count]) < length(mata_means) )
            {
              interim<-1
              
             #note there is another cat line 
              
              if (m==1) { 
               # browser(text="1801")
            # seems to be the correct position !     
               #cat(paste0("\ninterim at id= ", mata_Obs[c(mg[i,"cumcases"]),idvar] ))
                for (it_interim in 1:mg[i,"cases"]) {
                  cat(paste0("\ninterim at id= ", mata_Obs[c(mg[i,"cumcases"]-mg[i,"cases"]+it_interim),idvar] ))      
                }             
              }
              
           #11/11   interim_id<-  rbind(interim_id, mata_Obs[c(mg[i,"cumcases"]),"id"])  
              
              # end sigma test successful                 
            } #if  
           # browser()
          } #testinterim  processed save  interim_ids 
            
            
           # Sigma<-Sigmatrt
          }
          # 'J2R'
          
          # can delete in 1st pass just MAR??
          
          # 'CR'
          
          # 'CIR'
          
          # 'LMCF'
          
          # causal method uses same matrices as CIR with K parameter
         
          ############# individual analysis #########################
        }
        #else if(!is.null(methodindiv[1]))
        else if (flag_indiv==1) {

          # call function for  indiv
     #     browser(text="0501")
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

      #10/08/20 not necessary/incorrect?  
        #S11 <-Sigma[[1]][c_mata_nonmiss,c_mata_nonmiss]
        #S12 <-matrix(Sigma[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss))
        #S22 <-Sigma[[1]][c_mata_miss,c_mata_miss]

        # need to repictae mata_means to same numbe rrows as data pattern group
#      browser(text="2701")
        # causing problems replace with simpler
       mata_means<-mata_means[rep(seq_len(nrow(mata_means)),each=mg$cases[i]),]
       # mata_means<- do.call("rbind",replicate(mg[i,"cases"],mata_means,simplify=FALSE))
        
 #9/11
        # browser()

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



        j <- mg[i,"cases"]

        #for debug
        #print(paste0(" count in patt = ", j))

        k <- mg[i,"cumcases"]
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

 # check when no covar
       # browser(text="0401")
        # all missing ? didnt think we did ths scenario, ie base depvar always complete?
        # but what if no covar , ie no base depvar! need find suitable Sigma
        # ie when no covar selected and All missing then no Sigmas have been calculated 
        if (length(c_mata_nonmiss)==0)  {
          ## routine copied from mimix line 1229

          # change 9/5/20 becaus error list obj cannot be coerced to double
          #U <- chol(Sigma[[1]])
          # replaced this with S22 because S11 S12 have 0 values when so try
          U <- chol(S22)

          # generate inverse normal, same as used below
          #for debug 20/01
          #   print(paste0('mg[i,X1] =',mg[i,"X1"],' miss_count= ',miss_count))
          miss_count<-sum(mata_miss)
          Z <- stats::qnorm(matrix(stats::runif( mg[i,"cases"]* miss_count,0,1),mg[i,"cases"],miss_count))

          # raw and m1 null fields so hut use m2
          meanval = as.matrix(m2)

          mata_y1 = meanval+Z%*%(U)
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
        if (mg[i,"cases"] == 1) {
            meanval = (m2) + as.matrix(raw1 - m1)%*%as.matrix(t_mimix)
        }  else {
        meanval = as.matrix(m2) + (as.matrix(raw1 - m1)%*%as.matrix(t_mimix))
         }


        #24/02browser()
        U <- chol(conds)
        # mg[i,cases] is equiv to Stata counter, miss_count is no. of missing, so
        miss_count=rowSums(mata_miss)

        # gen erate inverse normal
        Z<-stats::qnorm(matrix(stats::runif( mg[i,"cases"]* miss_count,0,1),mg[i,"cases"],miss_count))
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
        GI <- array(data=mg[i,1],dim=c(mg[i,"cases"],1))
        #II  no imputations
        II <- array(data=m,dim=c(mg[i,"cases"],1))

        #doesnt need SNO, as id already in
        #mata_new<-cbind(GI,II,mata_new,SNO)

        #this works but bette to pre initialise data structure outsidr loop
        mata_new<-cbind(GI,II,mata_new,SNO)


    # this no longer applicable after 1st pass , only after 2nd plass
        #assume delta to be used if specified in input argument
    #    if (length(delta != 0) ) {
         # browser()
        #  browser(text="ts2") 
    #    mata_new <-  AddDelta(tst2,  ncovar_i,mata_new,delta,dlag)
    #    }

        mata_all_newlist[[m_mg_iter]]=mata_new
        #mata_all_new<-rbind(mata_all_new,mata_new)
        #mata_all_new<- na.omit(mata_all_new)
        #stata equilv \ is row bind NOt cbind!!?

        #9/5/20 end else c_mata_nonmiss ==0 here
        }


      } # ( m in  1:M) so insert delta module just before this
    } #for row[mg]
# check for interim 11/11 
  # browser(text="0201")
    if ( (interim==1) & (length(c_mata_miss)!=0))
      {
      #save interim ids
     # browser(text="1112 interimid")
      # this line only cacthes last interim case in th pattern group bu if more than 1 case will omit the previous
      # hence need introduce a interim counter or catc all cases 
     #interim_id<-  rbind(interim_id, mata_Obs[c(mg[i,"cumcases"]),idvar])
      for (iter_interim in 1:mg[i,"cases"]) {
        interim_id<-  rbind(interim_id,mata_Obs[c(mg[i,"cumcases"]-mg[i,"cases"]+iter_interim),idvar] ) 
      }
      #re-set interim flag
      interim<-0  
    }
    
    
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

  # 11/11 if interim option then save results also interim ids
 # if (testinterim==1){
#    impdatasetMAR<-getimpdatasets(list(mata_all_newlist,mg,M,method)) 
#  }
  # reset testinterim to run other methods, eg J2R 
  
  # save te MAR data-set and need the .imp=0's to replace in the final data output after pass2
 # browser(text="0912 check mata_Obs id col")
 # browser(text="1801")
  impdataset<-getimpdatasets(list(mata_all_newlist,mg,M,method))

  # if regression requested
  #If (regress = TRUE) {
  #     impdatamids <- as.mids(impdataset,  .id="SNO",.imp="II")
  #     fit<-with(data= impdatamids, expr = lm(impdatamids[,length(mpdatamids)-4]~impdatamids[,treatvar])
  #     paste0(summary(pool(fit)))
  #}
  #return(list(mata_all_newlist,mg,M,meth))

  
      #11/11 return(impdataset)
  # if runnin interims
  # for each cycle ,ie m value, method MAR, insert MAR nterims into .imp=0, then rerun
  # using interims as orig data ,   
  #browser(text="interims")
 
  # only perform following if not individual method  as only 1 pass for that
  # 0501
  if (flag_indiv!=1) {   
    if  (testinterim==1){  
     
     #2311 stat from after MAR this not works
     #2811 need to change this as fillointerims works by taking mean of interims over M imputaions so usng same interims over all M data sets
     # in 2nd pass but his incorrect, need us eall M imputed interims so pass M diif dat sets to 2ns pass 
     #29/11
     #browser(text="m?")
     #impdatset is M imputed datasets  based on MAR so for interim cases need to replace .imp 0's interims with imputed values.  
     # and save for M data sets
    # browser(text="M? 0112")
     #212 return list of raw data plus the imputed interims 
     #rawplusinterim <- fillinterims(impdataset,interim_id,Mimp=M)
    
     # browser(text="0202")
     rawplusinterim <- fillinterims(impdataset,interim_id,M)
     Imp_Interims<<-rawplusinterim[[2]] 
     #0312
    # browser(text="0312")
     # 1sty obtain rawplusinterim_1
     test1611impD<-rawplusinterim[[2]]
     test1611imp1<-subset(as.matrix(test1611impD[test1611impD$.imp==1,]))
     # need match with mata_Obs
     impMarint0<-rawplusinterim[[1]]
     
     impMarint0[impMarint0[,'.id'] %in% test1611imp1[,'.id'], ] <- test1611imp1
     impMarint1 <-impMarint0 
     # need to drop the old patt!!
     impMarint1nopatt<-as.data.frame(impMarint1)[,!(names(as.data.frame(impMarint1)) %in% c("patt"))]
  
     STSdummy<- apply(as.data.frame(impMarint1nopatt)[,grepl(paste0(depvar,".","[0-9]"),names(as.data.frame(impMarint1nopatt)))],MARGIN=2,function(x) ifelse(!is.na(x),0,1))
  #careful using grep because if same phrses in covar then will be duplicated like head_base and head so must use paste as above
     
        colnames(STSdummy) <- paste0(colnames(STSdummy),'.miss')
     #also if covar then have to include these
 #    browser(text="0301")
     if (length(covar)!=0 )
      { 
       #assume covars complete hence .miss just col of 0's 
       impMarint1nopatt[,paste0(covar,".miss")]<- impMarint1nopatt[,covar]
       impMarint1nopatt[,paste0(covar,".miss")]<- 0
      } 
     
     sts4D<-(cbind(impMarint1nopatt,STSdummy))
        
     
    
     # now can respahe from wide to long BUT DO WE HAVE TO?? NO!
  #  patt <- rowSums(pows2)
     pows2 <- sapply(1:ncol(STSdummy),function(i) STSdummy[,i]*2^(i-1))
     #need to add up to find patt
     patt <- rowSums(pows2)
     sts4Dpatt<-cbind(sts4D,patt)
     # in case zero covars
     # is thks really neceaasary ,gives duplicae covar?
    #if (length(covar) !=0) {
    #   tmp_covpatt<-apply(as.data.frame(sts4Dpatt[,covar]),MARGIN=2,function(x) ifelse(!is.na(x),0,1))
       #add names
    #   colnames(tmp_covpatt) <- paste0(c(covar),".miss")
       # than combine below the dummies onto the finaldat
    #   sts4Dpatt<-cbind(impMarint1nopatt,tmp_covpatt,STSdummy,patt)
    # }   else {
    #   sts4Dpatt<-cbind(impMarint1nopatt,STSdummy,patt)
    # }
     
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
     #cbind(impMarint1,STSdummy,patt))
     #Now find mg table 
     # In order  to get patt with all missing patts
     Overall_patt<-unique(sts4Dpatt[grepl(".miss",colnames(sts4Dpatt))])
     patt<- unique(sts4Dpatt[,"patt"])
   
     #then combine depvar and covar patt!
     all_patt<-cbind(Overall_patt,patt)
     # this is thr mg table!
     # ensure no duplicate covars as in acupuncture data head_base.miss.1, gives error if no miss.1
     #all_patt<-all_patt[,-grep("*.miss.1",names(all_patt))]
  #   browser(text="0301") 
     test_ex1<-merge(ex1,all_patt,by="patt")[order(merge(ex1,all_patt,by="patt")$exid),]
     # I think test_ex1 exactly like mg so only need adjust finnaldatSS by taking out differnt ames
     finaldatSS<- finaldatSS[,!(names(finaldatSS)) %in%c(".imp","X1")]
     # to here 0312
     # now have achieved the mg table and msata_Obs so no need to call preprodata 2nd time! 

   
   }else {
     return(impdataset)
   }   

      #load("rawint8.Rda")
 
  # process the data with filled interims to find mg etc
  # 2411 seems treat is in wrong position in the output?
  #browser(text="call proprocess on 2nd pass")
 
   # superseded 0312 
  #test2211<-do.call(preprodata,list(data2,covar,depvar,treatvar,idvar,timevar,M,reference,method))
  #instead of calling preprodata must now adjust for test2211 elements
  #browser(text="processing")
  # for reporting purposes movd to pass2loop
  #if  (method==3) {model<-"J2R"}
  #else if ( method==2 ) {model<-"CR"}
  #else if  ( method==1 ) {model <-"MAR"}
  #else if  ( method==4) {model<-"CIR"}
  #else if  ( method==5) {model<-"LMCF"}
  #else if (method==6)   {model<-"CAUSAL" }
  
  cat(paste0("\n\nending interim imputations "))
  
  # from above after 1st call preprodata
 # ntreat<-sort(unlist(test2211[[2]]))
  # trea pos seems to cause problems so make sure last col
  
  #finaldatS<-test2211[[1]]
  finaldataS<-finaldatSS
  ntreat <- unique(finaldataS[c(treatvar)])
    #finaldatS <- testlist$finaldatS
  #mg<-test2211[[3]]
  #browser(text="0101")
  mg <- test_ex1
    # vital to get the mata_obs correctly sorted! so corresponds with mimix_group lookup
  # to be consistent with Stata move the base col after the fevs!
  #mata_Obs <- testlist$finaldatS
 # mata_Obs <- test2211[[1]]
 # browser(text="0301")
  mata_Obs<- finaldatSS
  #trweat pos need to be last?
  colx<-grep(treatvar,colnames(mata_Obs))
  mata_Obs.reorder<-mata_Obs[,c(1:(colx-1),(colx+1):length(mata_Obs),colx)]
  mata_Obs<-mata_Obs.reorder
  # 1012 try also putting here the change of position for id also   
 # browser(text="1112") # is just .id  works here not idvar!! 
  #idcol<- grep(paste0(".",idvar),colnames(mata_Obs))
  idcol<- grep(paste0(".id"),colnames(mata_Obs))
  # beginning mata_Obs.reorder<-mata_Obs[,c(idcol,2:(idcol-1),(idcol+1):length(mata_Obs))]
  #this only has to be done once !! better to do it at beinning of call
  mata_Obs.reorder<-mata_Obs[,c(1:(idcol-1),(idcol+1):length(mata_Obs),idcol)]
  mata_Obs<-mata_Obs.reorder
  
  
  #*CREATE AN EMPTY MATRIX FOR COLLECTING the after mar  IMPUTED DATA
  GI<-c(0)
  II<-c(0)
  SNO <-mata_Obs[1,1]
  dropid<-c("id")
  #31/03
  # browser(text="mata_Obsx")
  mata_ObsX<-mata_Obs[,!(names(mata_Obs) %in% dropid)]
  mata_all_new <- cbind(GI,II,mata_ObsX,SNO)
  #mata_all_new[ mata_all_new>=0] <-NA
  # this just initialises (maybe can do in lass2Loop?) 
  mata_all_newlist <- vector('list',M*nrow(mg))
  # after draws
  #return(list(mata_means,S11,S12,S22) )
  
# something not quite right with final data set as still has NAs in it for .imp=1 ,I think treat in wrong  col. may be not trwting fev.12 
  # error occuring in mata_all_newlist 
  # seems treat is gettingin the way between base and fev?!
  # passing in imputed interims by m
#browser(text="0312")
  #Imp_Interims<-rawplusinterim[[2]]    
  
  #0912 id in wrong position in mata_Obs?? 
 
  # acupuncture data mg has reprsted head_base 
  # note Imp_Interims has interims repeated  M times as well as .imp=0 
  # the fault is in mg head head_base.miss and head_base.miss.1
  # in mata_Obs  but have to be in for asthma!!
  
  
 # browser(text="0702")
 testpass2impdatset<- pass2Loop(Imp_Interims,method,mg,ntreat,depvar,covar,treatvar,reference,trtgp,mata_Obs,mata_all_newlist,paramBiglist,idvar,flag_indiv,M,delta,dlag,K0,K1)
 #browser(text="0702")
 # this doesnt call proprocess as data already in wide format  

 # this leagcay prob delete    
 # Method3(m,M,trtgpindex,referindex,paramBiglist,mata_nonmiss,mata_miss,c_mata_nonmiss,c_mata_miss)

  } #0501 
  
  
} #  2311 ORIGINAL END for mimix test


# rawplusinterim dataset ,now produce new pattern
#rawplusinterim<-fillinterims(impdatasetMARnoDnoBase)
# having found data ,now produce new pattern table using rawplusinterim, .imp==0 as input data but
# not raw input data, only at final dat stage (ie output from preprocess  

#try writing wrapper for mimix
#mimix<- function(data,covar=NULL,depvar,treatvar,idvar,timevar,M=1,reference=NULL,method=NULL,seed=101,prior="jeffreys",burnin=1000,bbetween=NULL,methodvar=NULL,referencevar=NULL,delta=NULL,dlag=NULL,K0=1,K1=1) {
  





#21/11
#rawplusinterim <- fillinterims(impdataset,interim_id)
#rawinterim0<-(subset((as.data.frame(rawplusinterim)),.imp==0))

# below not appropriate because need long format but no need to do anyway , instead  
# need to substitute back into mata_Obs 
# do.call( preprodata,list( rawplusinterim ,covar,depvar,treatvar,idvar,timevar,M,reference,method))


#' @title getimpdatasets
#' @description to obtain the M imputed data set from the output list into one dataset
#' @details This combines the imputations found from the M pattern groups
#' @param varlist  list of data containing imputed values from the M pattern groups
#' @return impdatasets


getimpdatasets <- function(varlist){
  #12129/5/20
 #browser(text="1212")
  # to obtain M imputed data sets
  # dimension of data set, nrows in pattern times no imputations,
  # note sub data sets wi have different cols if completely missing so
  mata_all_newlist<-  varlist[1]
  mg<-(varlist[2])
  M<- unlist(varlist[3])
  method<- unlist(varlist[4])

  #dimlist <- (nrow(mg[[1]])*M)

  # extract from nested list
  # combine into data set containing M imputed datasets
  mata_all_newData1x <- do.call(rbind,mata_all_newlist[[1]])
  # then sort (by imputation and patient id) into M data sets and split into M lists
  impdatasets <- mata_all_newData1x[order(mata_all_newData1x$II,mata_all_newData1x$SNO),]

 #browser()
  #############################################
  # now recreate orig data set by selecting 1st imputed data set and setting NAs using .miss dummies
 # browser(text="0201")
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
  # change II, SNo to be consistent with mids format

  names( impdatasets)[names(impdatasets)=="II"]<-".imp"
  #1212 just drop SNO?? needed because are pass2 duplicate id cols would be generated!
  impdatasets$.id <-NULL
  names( impdatasets)[names(impdatasets)=="SNO"]<-".id"
  
  # change row names to e sequential
  rownames(impdatasets)<-NULL
  #drop GI
  impdatasets$GI <-NULL
  
  # report and check number na's
 # browser(text="1801")
  cat(paste0("\n\nnumber of original na values = ", sum(is.na(subset(impdatasets,impdatasets$.imp==0)))))
 # cat(paste0("\nnumber of final na values = ", sum(is.na(subset(impdatasets,impdatasets$.imp>0)))))
  #browser()
  if (sum(is.na(subset(impdatasets,impdatasets$.imp>0))) !=0 ) { cat(paste0("\nWARNING! unimputed data values")) }
  # write which model processed
  # but not when indiv method used
  if (length(method) !=0 ) { 
  if  (method==3) {model<-"J2R"}
  else if ( method==2 ) {model<-"CR"}
  else if  ( method==1 ) {model <-"MAR"}
  else if  ( method==4) {model<-"CIR"}
  else if  ( method==5) {model<-"LMCF"}
  else if (method==6)   {model<-"CAUSAL" }
  cat(paste0("\n\nImputed data based on model ", model))
  
  } # length
  
  return(impdatasets)
}


#' mimix: A package for Reference-based imputation for longitudinal clinical trials with protocol deviation
#' 
#' similar to the Stata mimix function
#' 
#' The mimix package contains the functions preprodata and preproIndivdata to 
#'  process long longitudinal data into wide data format
#'  
#'  pass2Loop performs 2nd pass after interims found by MAR 
#'  
#'  Also the function Addelta to add delta adjustment to the imputed estimates 
#' @docType package
#' @name mimix
NULL     

#' @title pass2Loop
#' @description 2nd pass for specifued model after 1st pass MAR ran
#' @details reads the summary table based on missing data pattern- mg  mimix_group                                                                   
#' @details reflects the pattern and treatment group configuration of the raw data
#' @details then acts as a looping mechanism, norm2 is used as MCMC multivariate normal 
#' @param Imp_Interims Interim cases
#' @param method - Specified model to run Reference-based imputation method
#' @param mg  the summary table based on missing data patern
#' @param ntreat vector of treatment groups
#' @param depvar response variable 
#' @param covar covariate variable(s)  
#' @param treatvar Treatment group, coded 1,2,..
#' @param reference  Reference group for J2R, CIR, CR methods
#' @param trtgp treatmet grp
#' @param mata_Obs raw data with interims imputed 
#' @param mata_all_newlist   raw data with interims imputed in list
#' @param paramBiglist list of MCMC beta and Sigma parameters 
#' @param idvar Participant id
#' @param flag_indiv flag whether specified individual column in data  
#' @param M number of imputations
#' @param delta vector of delta values to add onto imputed values (non-mandatory) (a's in Rogers paper),length as number of time points
#' @param dlag vector of delta values to add onto imputed values (non-mandatory) (b's in Rogers paper),length as number of time points
#' @param K0 Causal constant for use with Causal method
#' @param K1 exponential decaying Causal constant for use with Causal method
#' @return impdataset the M imputed data-sets appended to the "missing values" data-set in wide format
#' @examples
#' \dontrun{
#' testpass2impdatset<- pass2Loop(Imp_Interims,method,mg,ntreat,depvar,treatvar,reference,trtgp,mata_Obs,mata_all_newlist,paramBiglist,idvar,flag_indiv,M,delta,K0,K1)
#'}

# error in mata_miss and mata_nonmiss have duplicate hesd_base cols

# make sure mg has the covariate.miss!

pass2Loop<- function(Imp_Interims,method,mg,ntreat,depvar,covar,treatvar,reference,trtgp,mata_Obs,mata_all_newlist, paramBiglist,idvar,flag_indiv,M,delta,dlag,K0,K1)
{  
  #browser(text="0702")
  # mle option? 
  # this doesnt call proprocess as data already in wide format 
  # for reporting purposes try here rather than runmimix
  if  (method==3) {model<-"J2R"}
  else if ( method==2 ) {model<-"CR"}
  else if  ( method==1 ) {model <-"MAR"}
  else if  ( method==4) {model<-"CIR"}
  else if  ( method==5) {model<-"LMCF"}
  else if (method==6)   {model<-"CAUSAL" }
  cat(paste0("\n\nbegining processing  ",model))
  
  # 2411 mata_all_newlist <- vector('list',M*nrow(mg))
  # browser(text="0312")
  # chreck cols mg alin with mata_Obs 
  # moved to runmimix 1012
  # idcol<- grep(paste0(".",idvar),colnames(mata_Obs))
  # beginning mata_Obs.reorder<-mata_Obs[,c(idcol,2:(idcol-1),(idcol+1):length(mata_Obs))]
  #this only has to be done once !! better to do it at beinning of call
  #   mata_Obs.reorder<-mata_Obs[,c(1:(idcol-1),(idcol+1):length(mata_Obs),idcol)]
  #   mata_Obs<-mata_Obs.reorder
  
  m_mg_iter<-0
  for (i in 1:nrow(mg))
  {
  
    # define mata_miss as vector of 1's denoting missing using col names ending i ".missing"
    # this section to be amended to cope with multiple covariates
    #browser(text="0101")
    #mg comes in on acupunture data with repeated base_head  
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
    cnt<- mg$cases[i]
    #browser()
    # treatment grp
    trtgp<- mg[i,treatvar]
    
    pattern <- mg$patt[i]
    #cat("\ntrtgp = ", trtgp)
    
    
    #browser()
    if  (!is.null(method) ) {
      cat("\n",treatvar ," = ", trtgp,"patt = ",pattern,"number cases = ", cnt) }
    else if(!is.null(methodvar) ) {
      cat("\n",treatvar ," = ", trtgp,methodvar," = ",as.character(mg[i,methodvar]),referencevar," = ",as.character(mg[i,referencevar]),"patt = ",pattern,"number cases = ", cnt)
    }
    
    #need to convert (relate) treatment group to position in ntreat (create Pindex vector)
    #unneceassary now recoded
    #browser()
    trtgpindex<-which(trtgp==ntreat)
     
    # try ths when lmcf or mar? . ie no or NULL  refernce !! 0702
    # make sure reference not applicable for MAR or LMCF, this seems to take care of the problem!
    if (length(reference) !=0)  {
       referindex<-which(reference==ntreat)  
     }
    
    # multiple  simulations start here within the pattern loop #########
    for ( m in  1:M)  {
      #*FOR INDIVIDUALS WITH NO MISSING DATA COPY COMPLETE DATA INTO THE NEW DATA MATRIX mata_all_new `m' TIMES
      #if `pat' == 0{
      m_mg_iter<-m_mg_iter+1
      
      # 0312 at every m we need to replace the already imputed interims estiated from fillinterims subroutune
      #impMarint0[impMarint0[,'.id'] %in% test1611impD02[,'.id'], ] <- test1611impD02
      
      # only need create these once, not every iteration/imputation, so test  if the last imputation has already been created 
      # because if it has then so have all the others!     
      if (!exists((paste0("Imp_Interims_",M)))  )  { 
        assign(paste0("Imp_Interims_",m),subset(as.matrix(Imp_Interims[Imp_Interims$.imp==m,])))
      } 
      # also create 0 for use in adddelta ?
      Imp_Interims_0<- subset(as.matrix(Imp_Interims[Imp_Interims$.imp==0,]))
      
      # browser(text="0412")
      # prblems wit hte matching, has to be ame cols but ordre affected so create index
      # tet_mata_Obs<-mata_Obs
      #  tet_mata_Obs$index<-1:nrow(tet_mata_Obs)
      #  tet_mata_Obs_drop<- subset(tet_mata_Obs,select= c(.id,patt,treat,index))
      
      #  tet_mata_Obs_x <-  subset(tet_mata_Obs,select= -c(patt,treat,index))
      #save the .miss cols
      #  tet_mata_Obs_miss<-tet_mata_Obs_x[,grep(".miss",colnames(tet_mata_Obs_x))]
      #drop the .miss cols
      #  tet_mata_Obs_x<- tet_mata_Obs_x[,-grep(".miss",colnames(tet_mata_Obs_x))]
      
      # have have same cols in interim file 
      ImpInters <-    get(paste0("Imp_Interims_",m))
 #     browser(text="1112") #NOTE depvar = fev in antidep data!!??r
       # set treatname ,may have to do others sometime?
      colnames(ImpInters)[colnames(ImpInters)==treatvar]<-"treat"
      test_Imp <- subset(ImpInters,select=-c(.imp,patt,treat))
     # browser(text="11121 interim for aNTIDEP data")
      # the last col (id) has to be moved to the 1st!
      test_Imp_r<-test_Imp[,c(ncol(test_Imp),1:(ncol(test_Imp)-1))]
      
      # now merge Imp_Interims from m'th iteration
      #mata_Obs[mata_Obs[,'.id'] %in% get(paste0("Imp_Interims_",m))[,'.id'],] <- get(paste0("Imp_Interims_",m))
      
      ## it is essential that the 2 files (have the same cols?) and in the same col order
      # note setkey wil sort by key
      #  mata_Obsdt2<-data.table(mata_Obs)
      #  test_Impdt<-data.table(test_Imp)
      #  setkey(test_Impdt,'.id')
      
      #0812 try
      # setting row seq to preserve order after merge , failng to ovewrite imp >1 values? 
      #       mata_Obs$index<-1:nrow(mata_Obs)
      #        Result<-(merge(mata_Obs,test_Imp,all.x=T))
      #        Result<-Result[order(Result$index),]
      #        mata_Obs<-subset(Result,select= -c(index))
      
      # just use as a look upi table ! so
      
      df1<-as.data.frame(mata_Obs)
      df2<-as.data.frame(test_Imp)
      
      # repeat for all depvar cols
      #0812 seems this good methd!?
      #df1$fev.2[match(df2$.id,df1$.id)]<-df2$fev.2)
      # find depvar vars 
      depcols<-setdiff( grep(paste0(depvar),names(df1)) , grep('.miss',names(df1)))
      # assume thes coorespnd with interim lookup ! prob ned a check here!!
    
   
      # actually faster using match than fmatch
    
      for ( pos in 1:length(depcols)) {
        df1[,depcols[pos]][match(df2$.id,df1$.id)]<-df2[,depcols[pos]]
      }      
    
      
      mata_Obs <- df1 
      # check that id col at end ratherthan 1012 beginning  0912?
      #   idcol<- grep(paste0(".",idvar),colnames(mata_Obs))
      # beginning mata_Obs.reorder<-mata_Obs[,c(idcol,2:(idcol-1),(idcol+1):length(mata_Obs))]
      #this only has to be done once !! better to do it at beinning of call
      #    mata_Obs.reorder<-mata_Obs[,c(1:(idcol-1),(idcol+1):length(mata_Obs),idcol)]
      #    mata_Obs<-mata_Obs.reorder
      
      
      
      
      # this alt methd can delete once above proven to work
      #      commonNames <- names(df1)[which(colnames(df1) %in% colnames(df2))]
      #      commonNames <- commonNames[commonNames != ".id"]
      #     dfmerge<- merge(df1,df2,by=".id",all=T)
      #   for(i in commonNames){
      #      left <- paste(i, ".x", sep="")
      #      right <- paste(i, ".y", sep="")
      #      dfmerge[is.na(dfmerge[left]),left] <- dfmerge[is.na(dfmerge[left]),right]
      #      dfmerge[right]<- NULL
      #      colnames(dfmerge)[colnames(dfmerge) == left] <- i
      #    }
      
      #mata_Obsdt2[mata_Obsdt2[,'.id'] %in% test_Impdt[,'.id]',] <- test_Impdt
      # but id pos changes from m =1 to m=2
      #tet_mata_Obs_x[tet_mata_Obs_x[,'.id'] %in% test_Imp[,'.id'],] <- test_Imp
      
      # have to rejoin and sort back to original order
      #   tet_mata_Obs_a<-merge(tet_mata_Obs_x,tet_mata_Obs_drop,by=".id")
      #    tet_mata_Obs_s<-tet_mata_Obs_a[order(tet_mata_Obs_a$index),]
      #    tet_mata_Obs_af<- cbind(tet_mata_Obs_s,tet_mata_Obs_miss)
      # now reset mata_Obs
      # make sure index dropped but also check posiyion of id!
      #    tet_mata_Obs_af <-tet_mata_Obs_af[,-c(which(colnames(tet_mata_Obs_af)=="index"))]
      #    mata_Obs <- tet_mata_Obs_af
      
      #  the row  orders are ok but the col order may be wrong ,have to check!
      
      # construct structure to savde interim ids but this get reinitialised to many times!
      #11/11 interim_id<- mata_Obs[c(mg[i,1]),"id"] 
      
      #browser(text="no missing values")
      # if no missing values
      if(length(c_mata_miss)==0 ) {
        # start and end row positions
        st<-mg[i,"cumcases"]-mg[i,"cases"]+1
        en <-mg[i,"cumcases"]
        # id (SNO) is 1st col
        # either reorganise mata_Obs so id is 1st col  or use id as col name 0812    
        #SNO<-mata_Obs[c(st:en),1]
        SNO<-mata_Obs[c(st:en),".id"]
        
        #21/07/20
        #browser()
        # this doesnt delete treat
        #mata_new <- mata_Obs[c(st:en),2:ncol(mata_Obs)]
        mata_new<-mata_Obs[c(st:en),!(names(mata_Obs) %in% c(idvar))]
        #browser() # treat defined within fun
        GI <- array(data=mg[i,treatvar],dim=c(mg[i,"cases"],1))
        #II  no imputations
        II <- array(data=m,dim=c(mg[i,"cases"],1))
        mata_new=cbind(GI,II,mata_new,SNO)
        #names(mata_all_new)<-names(mata_new)
    #    browser(text="0201")
        mata_all_newlist[[m_mg_iter]]=mata_new
        # mata_all_new=rbind(mata_all_new,mata_new)
        # else if there are missing values
      } else {
        # need to distinguish between meth and methodindiv
        
        if (flag_indiv==0 ) {
          
          referindex<- reference
          #FOR INDIVIDUALS WITH  MISSING DATA  `m' TIMES
          # dependent on method chosen
          # 'MAR'
          #browser(text="2611")
          
          
          if (method== 1)  {
            
            # 11/11 checking methd used for interims in J2R works as MAR 
            
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
            
            
            #11/11 flag interims put in as option in command
            # testinterim<-1
            
            #initialise interim 
            # interim<-0
            
            
            
            
            # Sigma<-Sigmatrt
          }
          # 'J2R'
          else if (method == 3 ) {
            
            # changed saving the result into  just the param file, list of 2 so can use list index here
            #treatmnets are 1.. M then M+1 ..2M .. etc
            mata_means_trt <- paramBiglist[[M*(trtgpindex-1)+m]][1]
            mata_means_ref <- paramBiglist[[M*(referindex-1)+m]][1]
            
            
            #browser(text="pass2method")
            # below causes error after using >1 covars and mata_nonmiss has covar.1, not proper covar names
            
            # try not unisting because error only fr J2r merod when on patien in patt
            #
            # 2/11
            # browser()
        #    browser(text="3112")
            # mata_nonmiss and mata_miss toolage in acupuncture data
            # repetition of head_base ,the covariate !
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
            
            
            
            
            
            ################test here just to see if get same tment effect inserting J2Rsigma##########
            #9/3/20 7/11/20
            ##### yes this has the effect of obtaining same treastment effect as J2R , ie 0.1166 so hopefully run this after running MAr    
            #      SigmaRefer <- paramBiglist[[M*(referindex-1)+m]][2]
            #      #SigmaRefer <- paramBiglist[[M*(referindex-1)+m]][2]
            #SigmaRefer <- get(paste0("paramBiglist",refer,m))[2]
            #        Sigmatrt <- paramBiglist[[M*(trtgpindex-1)+m]][2]
            # note use of [[1]] as is matrix rathe than list,
            #        S11 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_nonmiss]
            # causes non-def error in conds
            #to ensure rows and cols as should reflect their stucture use matrix
            #        S12 <-matrix(SigmaRefer[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss))
            #        S22 <-SigmaRefer[[1]][c_mata_miss,c_mata_miss] 
            #################### delete above !!! afte test ######################################
            #10/11
            #    browser()
            
          } #method
          # 'CR'
          else if (method==2) {
            # no need to use Sigmatrt here
            mata_means <- paramBiglist[[M*(referindex-1)+m]][1]
            #mata_means <- get(paste0("param",refer,m))[1]
            # convert from list to matrix
            mata_means <- mata_means[[1]]
            
            SigmaRefer <- paramBiglist[[M*(referindex-1)+m]][2]
            S11 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_nonmiss]
            S12 <-matrix(SigmaRefer[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss) )
            S22 <-SigmaRefer[[1]][c_mata_miss,c_mata_miss]
            
            #Sigma<-SigmaRefer
          }
          # 'CIR'
          else if (method==4)
            
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
            
            # Sigma<- SigmaRefer
          }
          # 'LMCF'
          else if (method==5) {
            #browser(text="0602")
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
            
            #   Sigma<- Sigmatrt
          }  #if meth=5
          # causal method uses same matrices as CIR with K parameter
          else if (method==6)  {
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
            #10/08/20          
            #browser(text="1301")
            
            mata_means<-Causal_loop(c_mata_miss,mata_Means,MeansC,K0,K1)
            #  browser()
            #this temporary  for test purposes until algo decided upon
            SigmaRefer <- paramBiglist[[M*(referindex-1)+m]][2]
            # when reading in Stata sigmas
            
            
            S11 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_nonmiss]
            S12 <-matrix(SigmaRefer[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss))
            S22 <-SigmaRefer[[1]][c_mata_miss,c_mata_miss]
            #10/08/20  
            #Sigma<- SigmaRefer
            
          }
          ############# individual analysis #########################
        }
        #else if(!is.null(methodindiv[1]))
        else if (flag_indiv==1) {
          
          # call function for  indiv
          
          indparamlist  <- ifmethodindiv(methodvar,referencevar,mg,m,M,paramBiglist,i,treatvar,c_mata_nonmiss,c_mata_miss,mata_miss,mata_nonmiss,K0,K1)
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
        
        #10/08/20 not necessary/incorrect?  
        #S11 <-Sigma[[1]][c_mata_nonmiss,c_mata_nonmiss]
        #S12 <-matrix(Sigma[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss))
        #S22 <-Sigma[[1]][c_mata_miss,c_mata_miss]
        
        # need to repictae mata_means to same numbe rrows as data pattern group
        
        # causing problems replace with simpler
        mata_means<-mata_means[rep(seq_len(nrow(mata_means)),each=mg$cases[i]),]
       # mata_means<- do.call("rbind",replicate(mg[i,"cases"],mata_means,simplify=FALSE))
        
        #9/11
        # browser()
        
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
        
        
        
        j <- mg[i,"cases"]
        
        #for debug
        #print(paste0(" count in patt = ", j))
        
        k <- mg[i,"cumcases"]
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
         #0401 U <- chol(Sigma[[1]])
          # replaced this with S22 because S11 S12 have 0 values when so try
          U <- chol(S22)
          
          # generate inverse normal, same as used below
          #for debug 20/01
          #   print(paste0('mg[i,X1] =',mg[i,"X1"],' miss_count= ',miss_count))
          miss_count<-sum(mata_miss)
          Z <- stats::qnorm(matrix(stats::runif( mg[i,"cases"]* miss_count,0,1),mg[i,"cases"],miss_count))
          
          # raw and m1 null fields so hut use m2
          meanval = as.matrix(m2)
          
          mata_y1 = meanval+Z%*%(U)
          #set dimensions mata_new to mata_y1
          #define new matrix from observed,  id column  (the last)
          mata_new <- preraw
          # mata_new has to be already defined
          mata_new[,c_mata_miss] <- (mata_y1)
          GI <- array(data=mg[i,treatvar],dim=c(mg[i,"cases"],1))
          #II  no imputations
          II <- array(data=m,dim=c(mg[i,"cases"],1))
          # SNO just id col
          # change 1012
          #SNO <- mata_Obs[c(startrow:stoprow),1]
          SNO <- mata_Obs[c(startrow:stoprow),".id"]
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
          if (mg[i,"cases"] == 1) {
            meanval = (m2) + as.matrix(raw1 - m1)%*%as.matrix(t_mimix)
          }  else {
            meanval = as.matrix(m2) + (as.matrix(raw1 - m1)%*%as.matrix(t_mimix))
          }
          
          
          #24/02browser()
          U <- chol(conds)
          # mg[i,cases] is equiv to Stata counter, miss_count is no. of missing, so
          miss_count=rowSums(mata_miss)
          
          # gen erate inverse normal
          Z<-stats::qnorm(matrix(stats::runif( mg[i,"cases"]* miss_count,0,1),mg[i,"cases"],miss_count))
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
          # chg 1012SNO <- mata_Obs[c(startrow:stoprow),1]
          SNO <- mata_Obs[c(startrow:stoprow),".id"]
          #SNO <- mata_ObsX[,ncol(mata_Obs)]
          # GI treatment grp column 1 (here),II imputation number col, mata_new matrix then SNO is id col.
          GI <- array(data=mg[i,1],dim=c(mg[i,"cases"],1))
          #II  no imputations
          II <- array(data=m,dim=c(mg[i,"cases"],1))
          
          #doesnt need SNO, as id already in
          #mata_new<-cbind(GI,II,mata_new,SNO)
          
          #this works but bette to pre initialise data structure outsidr loop
          mata_new<-cbind(GI,II,mata_new,SNO)
          
          
          #assume delta to be used if specified in input argument
          if (length(delta != 0) ) {
            # browser()
            
            # reset ncovar_i and create dlag if null] browser(text="1612")
           # browser(text="1712")
            ncovar_i<-length(covar)
            if (is.null(dlag)) {
              dlag <- rep(1,length(delta))
            }  
            mata_new <-  AddDelta(tst2,  ncovar_i,mata_new,delta,dlag)
           # mata_new <-  AddDelta(tst2,  ncovar_i,mata_new,delta,dlag,Imp_Interims_0)
          }
          
          mata_all_newlist[[m_mg_iter]]=mata_new
          #mata_all_new<-rbind(mata_all_new,mata_new)
          #mata_all_new<- na.omit(mata_all_new)
          #stata equilv \ is row bind NOt cbind!!?
          
          #9/5/20 end else c_mata_nonmiss ==0 here
        }
        
        
      } # ( m in  1:M) so insert delta module just before this
    } #for row[mg]
    # check for interim 11/11 
    # browser()
    
    
  } #for M StOP HERE!!
  # browser(text="passtoloop")
  #browser(text="1801")
  impdataset<-getimpdatasets(list(mata_all_newlist,mg,M,method))
  # but need to adjust orig data set to set interims back to missing
  # get .imp=0 's
  assign(paste0("Imp_Interims_",0),subset(as.matrix(Imp_Interims[Imp_Interims$.imp==0,])))
  # find depvar vars  cols 
  ImpInters <-    get(paste0("Imp_Interims_",0))
  # note keep .imp so can use as match with impdaset and also presreve same col pos
 # browser(text="1112 treat? ")
  colnames(ImpInters)[colnames(ImpInters)==treatvar]<-"treat"
  test_Imp <- subset(ImpInters,select=-c(patt,treat))
  
  #df1<-as.data.frame(mata_Obs)
  test_Imp<-as.data.frame(test_Imp)
  
  
  depcolsf<- grep(paste0(depvar),names(impdataset)) 
  # assume thes coorespnd with interim lookup ! prob ned a check here!!
  
  for ( pos in 1:length(depcolsf)) {
    #if (impdataset$.imp==0 ) 
    #impdataset[,depcolsf[pos]][match(test_Imp$.id,impdataset$.id)]<-test_Imp[,depcolsf[pos]]
    impdataset[,depcolsf[pos]][match(paste(test_Imp$.id,test_Imp$.imp),paste(impdataset$.id,impdataset$.imp))]<-test_Imp[,depcolsf[pos]]
  }  
  # moved from getimpdatasets fun
  cat(paste0("\nnumber of final na values = ", sum(is.na(subset(impdataset,impdataset$.imp>0)))))
  cat(paste("\ntest pass2 in runmimx"))
  return(impdataset) 
}     

#' @title fillinterims
#' @description fills missing interims distinguising from post-discontinuation
#' @details checks methodindiv not null
#' @param impdata the data with missing values 
#' @param interims the interim cases with estimated MAR values
#' @param Mimp the number of imputations specified , ie M total imputsations
#' @return list of 1st data set with interims imputed plus M interim cases of each interim case to be matche in 2nd pass  


fillinterims<- function(impdata,interims,Mimp=M ) {
  #2811
  #0312 browser(text="not means")
  #2711browser(text="find estimate over all imps")
  #convert to data.table
  #2611 browser(text="fillinterims")
  impMarint_dt <- data.table::as.data.table(impdata)
  interims_dt <- data.table::as.data.table(interims)
  
  #browser(text="setkey 1112")
  
  data.table::setkey(impMarint_dt,.id)
  data.table::setkey(interims_dt,V1)
  #merge
  (impMarint_dt[interims_dt])
  # convert to 0, 1 =missing
  sapply(impMarint_dt[interims_dt],function(x) ifelse(is.na(x) ,1,0) )
  # exclude last non-response cols
  test10<-sapply(impMarint_dt[interims_dt],function(x) ifelse(is.na(x) ,1,0) )
  test10x<-test10[,c(1:(ncol(test10)-3))]
  # find max non-missing
  lastvalid<-apply(test10x,1, function(x) max(which(x==0))  )
  #merge back
  test1611<-cbind(impMarint_dt[interims_dt],lastvalid)
  #if (test1611$.id == shift(test1611$.id,-1) ) {
  
  # over many imps, find mean value by interims for imp>0 then combineback to
  # have .imp0 and mean row
  
  # 2811 no need for means now just need to send  M'th values of fev 2-fev12
  #test2611 <-setDT(as.data.frame(subset(test1611,.imp>0)))[,lapply(.SD,mean),by= .id,]
  #2911 test2611 <- setDT(as.data.frame(subset(test1611,.imp>0)))[,lapply(.SD,mean),by= .(id,imp)]
  # find.imp0's  test3611[order(test3611$.id),][order(assign(paste0("test1611imp",val),subset(test1611,.imp==val))$.id),]
  test1611imp0<- subset(test1611,.imp==0)
  #declare void matrix for rbinding 
  test1611impD <- test1611imp0
  for (val in 1:Mimp) {
    rbind(assign(paste0("test1611imp",val),subset(test1611,.imp==val | .imp==0)),test1611imp0)
    #test1611impX<- rbind(test1611imp0,test1611)
    #assign(paste0("test1611imp",val),subset(test1611,.imp==val))
    # get(paste0("test1611imp",val))<-as.matrix(get(paste0("test1611imp",val)))
    #}
    # now tes1611imp2 etc alrready sorted with .imp0 and .impm records fror each interim record so jus tned to proces in loop  
    # no need sort 
    #test1611imp2[order(test1611imp2$.id,test1611imp2$.imp),]f  (is.na(test3611s[k,i]) & (i<test3611s[k,ncol(test3611s)] )) {
    
    test1611impx<-as.matrix(get(paste0("test1611imp",val )))
    for  (r in seq(from =1 , to= nrow(test1611impx)-1,by=2)) {
      for (j in 2:(ncol(test1611impx)-3)) {
        # if element na and must be bfore lastvalid non missing
        if (is.na(test1611impx[r,j]) & (j<test1611impx[r,ncol(test1611impx)]) ) {
          # then shift value from below row
          test1611impx[r,j] = test1611impx[r+1,j]
        }
      }
    }
    # keep imp0's and recode non imp0 to imp (m) number
    test1611impxm <-  subset(test1611impx,test1611impx[,".imp"]==0)
    test1611impxm[,".imp"]<-val
    
    # build up ovder M imps
    test1611impD <- rbind(test1611impD,test1611impxm)
    
  }     
  #browser(text="1801")   not right
  #cat(paste0("\nnumber of interims final na values = ", sum(is.na(subset(test1611impD, test1611impD$.imp==0)))))
# browser(text="0501")
  # test1611impD has the interim cases imputed for each imputation number m
  # so need to insert into original unimputed data set for each imputstion and process each dat set into the 2nd pass     
  
  
  # take lastvalid off
  #test1611impDl<-test1611impD[,-ncol(test1611impD)]
  test1611impD$lastvalid <-NULL
  #extract m'th imp
  test1611impD02<-as.matrix(test1611impD[test1611impD$.imp==2,])
  #test1611impD02<-as.matrix(test1611impDl[test1611impDl[,1]==2,])
  
  impMarint <- as.matrix(impMarint_dt)
  impMarint0<-as.matrix(impMarint_dt[impMarint_dt$.imp==0,])
  # this works well (but MUST Be MATRICES)   
  #212 impMarint0[impMarint0[,'.id'] %in% test1611impD02[,'.id'], ] <- test1611impD02
  # so .imp=0 now has interims filled and can be treated as original data  
  
  return(list(impMarint0,test1611impD))
}


