

# trying to write a subroutine for this part of data preprocessing
A <- function(x) {  ifelse(!is.na(x),0,1) }


  
#readdata
readdata <-function(data) {
  # Specify full path to data, e.g /directory/path/to/NIRData.csv
  #path_to_my_data <- file.path("directory","path","to","NIRData.csv")
   # Check if file exists at that path
  stopifnot(file.exists(paste0("./",data)))
  mxdata <-read.csv(paste0("./",data),fileEncoding = "UTF-8-BOM")
  return(mxdata)
}


# this section to find mg, the missing  pattern group and also converts from long to wide data format
preprodata<- function(depvar,treatvar,idvar,timevar,covar,M,refer,meth)  {
  #extract relevant vars
 
  #more informative in error msg to use this explicit and  
  #and put in one statement 
  stopifnot( (meth == "MAR" | meth=="J2R" | meth=="CIR" | meth=="CR" | meth=="LMCF"),
             #is.numeric(refer),
             is.numeric(M),
             is.character(depvar),
             is.character(treatvar),
             is.character(idvar),
             is.character(timevar),
             is.character(covar)    )
  
  #stopifnot(any(meth_valid_values==meth))
        #stopifnot(is.numeric(refer))
        #stopifnot(is.numeric(M))
  #mxdata <- do.call( readdata,kmargs[1])
 
  fevdata<-select(mxdata,idvar,depvar,timevar)
  # only want to widen the dose var, ie fev, so take the other variables
  #reshape to wide and assign new names 
  # prefix must be supplied from input argument rather than hard coded, this canbe hard coded
  sts4<-pivot_wider(fevdata,id_cols=c(idvar),names_from=timevar,names_prefix=depvar,values_from=depvar)
  #sts4 is just response data, can join later treat and covar cols from finaldat 
  #assumes no NA for these vars should add in checking routines !!
  uniqdat<-unique(mxdata[c(idvar,covar,treatvar)])
  
  
  finaldat<- merge(sts4,uniqdat,by=idvar)
  # now try and sort on treatvar, doesnt work so instead of sorting , select on treat
 
  
  # need find no. ntreat  to loop over
  #ntreat<-unique(unique(mxdata$treatvar))
  ntreatcol<-(select(mxdata,treatvar))
  ntreat <- unique(ntreatcol)
  ntimecol<-(select(mxdata,timevar))
  ntime<-unique(ntimecol)
  
  
  #in order to aggregate by pattern need create dummy vars
  
  # for sts4 array drop 1st col and create newvars for rest
  STSdummy<-apply(sts4[,2:ncol(sts4)],MARGIN=2,FUN=A)
  # append to names
  colnames(STSdummy) <- paste(colnames(STSdummy),'.1')
  # merge back on to data 
  sts4D<-cbind(finaldat,STSdummy)
  # create powrs of 2
  pows2 <- sapply(1:ncol(STSdummy),function(i) STSdummy[,i]*2^(i-1))
  #need to add up to find patt
  patt <- rowSums(pows2)
  #sts4Dpatt<-cbind(sts4D,patt)
  sts4Dpatt<-cbind(sts4D,patt)
  
  # then sort by patt( last col) and treat and find cumulative sequence  AND split on treatment var!
  # so need to merge treat (and covar) back into response data 
  
  #sort by treatvar and patt var, couldnt get treatvar working with patt so try sorting previously on treatvar   
  finaldatSS <-sts4Dpatt[order(patt),] 
  
  ndx = order(finaldatSS[,treatvar])
  finaldatS <- finaldatSS[ndx,]
  
  # consistent with Stata to get right order, drop id and treat cols
  drops <-c(idvar,treatvar)
  
  
  #also get the dummy pattern vectors by treatment
  #dummypatt<-unique(xSTSdummy)
  pattmat<-unique(STSdummy[,1:ncol(STSdummy)])
  
  # sort of hard coded for now may have to revisit this!
  pows2d <- sapply(1:ncol(pattmat),function(i) pattmat[,i]*2^(i-1)) 
  patt<-rowSums(pows2d)
  #want  join this to ex from mimix_group table
  #merge(uniqdat,sts4,by=idvar)
  pattmatd<-cbind(pattmat,patt)
  
  # this best for  count, implicit sort by treat
  # may have to generalise using function argument  later!
  # this bit groups and aggregates but can try using base (aggregate) instead! .
  #  
  ex1 <- sts4Dpatt %>%
    # test chg 17/01 treat needs made generic
    # substantiate var?
    
    dplyr::group_by(sts4Dpatt[,c(treatvar)],sts4Dpatt$patt) %>%
    #dplyr::group_by(sts4Dpatt$treat,sts4Dpatt$patt) %>%
    dplyr::summarise(X1 =n())
  ex1 <- dplyr::rename(ex1,treat="sts4Dpatt[, c(treatvar)]",patt='sts4Dpatt$patt')
  #ex1 <- dplyr::rename(ex1,treat='sts4Dpatt$treat',patt='sts4Dpatt$patt')
  # and now find cumX1
  ex1$X1cum <- cumsum(ex1$X1) 
  #want  join this to ex from mimix_group table
  # create index number to preserve ordr after merge
  ex1$exid <- 1:nrow(ex1)
  ex1id <-merge(ex1,pattmatd,by="patt")
  ex1s<-ex1id[order(ex1id$exid),]
  
  #try using aggregating instead 
  #aggregate(x=sts4Dpatt,by=list(treat,patt),FUN="sum")
  
  
  
 
 # for (val in t(ntreat)) {
 #    print(paste0("prenormMC",val))
    #assign(paste0("prenorm2",val),subset(finaldat,treat==val))
 #  } 
  
  
  #error chk
  stopifnot(refer %in% t(ntreat))
  print("summary missing pattern")
  print(ex1s)
  return(list(sts4Dpatt,finaldatS,finaldat,ntreat,ex1s,ex1,ex1id,pattmat,patt,ntime,M,refer,meth))
}
# Question is should the wide data be sorted for analysis or just ordered in fact, selecting out treatmet and rbinding? to provide analysis data set.   
# answer is should be sorted by patt as in Stata!
# therefore just select on treat when looping thru ntreat 


# Main function 
Runmimix<- function(depvar,treatvar,idvar,timevar,covar,M=1,refer=1,meth=NULL,seedval=101) {
  # 
  
  #set.seed(101)
  
  set.seed(seedval)
  #set.seed(unlist(tail(kmargs,n=1)))
  
  #mxdata<- readdata(data)
  
  
  #do not read in last elemnt LAST otherwise get an unused argument error msg
  testlist = do.call( preprodata,kmargs[-length(kmargs)])
  
  
  # returns list from preprodata function
  ntreat<-unlist(testlist[[4]])
  #stopifnot()
  finaldatS<-testlist[[2]]
  mg<-testlist[[5]]
  # vital to get the mata_obs correctly sorted! so corresponds with mimix_group lookup 
  # to be consistent with Stata move the base col after the fevs!
  mata_Obs <- testlist[[2]]
  M <- testlist[[11]]
  refer <- testlist[[12]]
  stopifnot(refer %in% ntreat)
  meth <- testlist[[13]]
  
  # so tsts2 is th equiv to mi_impute so try to put in norm2
  # est tis tomorrow 17/01 
  #tst2 <- mi_impute("id","time","fev","base")
  #list("fev","treat","id","time","base",10,2,"J2R",101)
  tst2 <- mi_impute(kmargs[[3]],kmargs[[4]],kmargs[[1]],kmargs[[5]])
  # put commas in
  tst3<-paste(tst2,collapse = ",")
  
  #create input data sets for each tment from which to model 
  for (val in t(ntreat)) {
    print(paste0("prenormdat",val))
    assign(paste0("prenormdat",val),subset(finaldatS,treat==val))
  }  
  
  
  
  #*CREATE AN EMPTY MATRIX FOR COLLECTING IMPUTED DATA 	
  # just an empty row but take dimensions and col types fron mata_obs
  #mata_all_new<-array(data=1,dim=c(1,2+ncol(mata_Obs)))
  GI<-c(0)
  II<-c(0)
  SNO <-mata_Obs[1,1]
  dropid<-c("id")
  mata_ObsX<-mata_Obs[,!(names(mata_Obs) %in% dropid)] 
  mata_all_new <- cbind(GI,II,mata_ObsX,SNO)
  mata_all_new[ mata_all_new>=0] <-NA 
  mata_all_newlist <- vector('list',M*nrow(mg))
  #Warning message:
  #In Ops.factor(left, right) : ‘>=’ not meaningful for factors
  
  
  
  # run the mcmc simulations over the treatment grps 
  # create a matrix for param files, a beta  and sigma matrices
  #paramMatrix<-matrix(1:(nrow(ntreat)*M),nrow=M,dimnames=list(c(1:M),c(1:nrow(ntreat))))
  
  #paramMatrix<-matrix(,nrow=(nrow(ntreat)*M),ncol=2)
  
  #create  emptylist for each treat and multiple m's
  paramBiglist <- vector('list',length(ntreat)*M)
  #create 
  
  #paramMatrixT<-matrix(,nrow=2,ncol=M)
  start_time <- proc.time()
  #instead iof start _time try system_time
  #paramBetaMatrixT <- matrix(1,nrow=1,ncol=5)
  #paramSigmaMatrixT <- matrix(,nrow=5,ncol=5)
  iter<-0
  system.time(
    
    #try ser up as many Result data files as treatments instead of one big file?.
    
    
    for (val in t(ntreat)) {
      #prnormobj <- subset(prenormdat2, select=c(tst2))
      #prnormobj <- assign(paste0("prenormdat",val),subset(finaldatS,treat==val))
      kmvar=get(paste0("prenormdat",val))
      prnormobj<-assign(paste0("prnormobj",val), subset(kmvar, select=c(tst2)))
      #create  emptylist for each treat 
      # assign(paste0("paramTRTlist",val),vector('list',M))
      print(paste0("Looping for treatment = ",val," performing mcmcNorm for m = 1 to ",M))
      for(m in 1:M) {  
        #emResultT<-emNorm(testprnormobj,prior = "ridge",prior.df=0.5)
        # supppress warnings regrading solution  near boundary, see norm2 user guide, also mimix about this problem
        # may have to change prior '
        emResultT<-suppressWarnings(emNorm(prnormobj,prior = "jeffreys"))
        #print(paste0("running emNorm"))
        mcmcResultT<- suppressWarnings(mcmcNorm(emResultT,iter=1000,multicycle = NULL,prior = "jeffreys"))
        #print(paste0("running mcmcNorm"))
        #try saving parm files  to a matrix instead of indiv parm files
        
        # assign doesnty make much differenc in tghis loop , perhaps efficiency comes from not calling them later?
        
        # keep fo rnow to test againt biglist 
        # teszted ok for j2r
        #assign(paste0("param",val,m),mcmcResultT$param)
        
        iter<-iter+1    
        paramBiglist[[iter]] <- mcmcResultT$param
        
        
        }  
      
    }
  ) # system.time 
  print(paste0("mcmcNorm Loop finished, m = ",M))
  
  # try and use lappy instead of loop for M
  
  
  
  
  # can repeat interactively from here
  # now loop over the lookup table mg, looping over every pattern
  # declare iterate for saving data
  m_mg_iter<-0
  for (i in 1:nrow(mg))
  {
    #define mata_miss as vector of 1's denoting missing using col names ending i ".1"
    mata_miss <- mg[i,grep("*..1",colnames(mg)),drop=F]
    #mata_miss <- mimix_group[i,c(2,3,4,5)]     #define mata_miss
    #assumes covariate non missing
    mata_miss$covar.1 <-0                       #assuming cov col is non missing
    mata_nonmiss <- (ifelse(mata_miss==0,1,0))  #define mata-nonmiss from miss
    
    # need transform nonmiss,miss to c lists - ie. index the  
    c_mata_miss<-which(mata_miss==1)
    c_mata_nonmiss<-which(mata_nonmiss==1)
    # eg no missing is c(1,2,3,4,5)
    
    #debug
    #print("mata_miss")
    #print(mata_miss)
    
    # count of pattern by treatment
    cnt<- mg$X1[i]
    
    # treatment grp
    trtgp<- mg$treat[i]
    
    pattern <- mg$patt[i]
    #cat("\ntrtgp = ", trtgp)
    cat("\ntrtgp = ", trtgp,"patt = ",pattern,"no patients = ", cnt)
    
    #need to convert (relate) treatment group to position in ntreat (create Pindex vector) 
    #something like
    #position(trtgp,"ntreat")
    # in fact this should do it
    trtgpindex<-which(trtgp==ntreat)
    referindex<-which(refer==ntreat) 
    
    # multiple  simulations start here within the pattern loop #########  
    for ( m in  1:M)  { 
      #*FOR INDIVIDUALS WITH NO MISSING DATA COPY COMPLETE DATA INTO THE NEW DATA MATRIX mata_all_new `m' TIMES
      #if `pat' == 0{ 
      m_mg_iter<-m_mg_iter+1
      if(length(c_mata_miss)==0 ) {
        st<-mg[i,"X1cum"]-mg[i,"X1"]+1
        en <-mg[i,"X1cum"]
        SNO<-mata_Obs[c(st:en),1]
        mata_new <- mata_Obs[c(st:en),2:ncol(mata_Obs)]
        GI <- array(data=mg[i,"treat"],dim=c(mg[i,"X1"],1))
        #II  no imputations 
        II <- array(data=m,dim=c(mg[i,"X1"],1))
        mata_new=cbind(GI,II,mata_new,SNO)
        #names(mata_all_new)<-names(mata_new)
        mata_all_newlist[[m_mg_iter]]=mata_new
        # mata_all_new=rbind(mata_all_new,mata_new)
      } else {
        
        #FOR INDIVIDUALS WITH  MISSING DATA  `m' TIMES  
        # dependent on method chosen 
        
        if (meth== 'MAR')  {
          
          #mata_means <- get(paste0("param",trtgp,m))[1]  
          mata_means <- paramBiglist[[M*(trtgpindex-1)+m]][1]
          # mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]
          # convert from list element to matrix
          mata_means <- mata_means[[1]]
          
          
          #Sigmatrt <- get(paste0("param",trtgp,m))[2]
          Sigmatrt <- paramBiglist[[M*(trtgpindex-1)+m]][2]
          S11 <-Sigmatrt[[1]][c_mata_nonmiss,c_mata_nonmiss]
          S12 <-Sigmatrt[[1]][c_mata_nonmiss,c_mata_miss]
          S22 <-Sigmatrt[[1]][c_mata_miss,c_mata_miss]
          
          
        }
        else if (meth == 'J2R') { 
          
          # changed saving the result into  just the param file, list of 2 so can use list index here
          #treatmnets are 1.. M then M+1 ..2M .. etc
          mata_means_trt <- paramBiglist[[M*(trtgpindex-1)+m]][1]
          mata_means_ref <- paramBiglist[[M*(referindex-1)+m]][1]
          
          #mata_means_trt <- get(paste0("param",trtgp,m))[1]
          #mata_means_ref <- get(paste0("param",refer,m))[1]
          #mata_means_trt <- get(paste0("mcmcResultT",trtgp,m,"$param$beta")) 
          #mata_means_ref <- get(paste0("mcmcResultT",refer,m,"$param$beta"))
          
          #mata_means <- get(paste0("parambeta",trtgp,m))
          #one way is to element multiply (because 1,0) then add 
          mata_means_t <- unlist(mata_means_trt)*mata_nonmiss
          mata_means_r <- unlist(mata_means_ref)*mata_miss
          # so when all missing  1,1,1, ... then all contributing comes from reference means      
          mata_means <- mata_means_r+mata_means_t
          # and preserve names   
          colnames(mata_means) <- colnames(mata_means_trt)
          #replicate to number of rows defined by X1 
          #mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]
          
          
          ############# SIGMA is from paramsigma  the reference group ################
          
          
          # do we ever  use SigmaTrt !!?? in j2r? 
          
          SigmaRefer <- paramBiglist[[M*(referindex-1)+m]][2]
          #SigmaRefer <- get(paste0("param",refer,m))[2] 
          
          # get(paste0("mcmcResultT",refer,m,"$param$sigma"))
          # when reading in Stata sigmas
          # needs to to ths as tibble will fail in cholesky
          #SigmaRefer <- as.matrix(get(paste0("paramsigmaStata",refer,m)))
          
          #print(paste0("paramsigma",refer,m))
          # print(paste0("refer,m = ",refer," ",m))   
          #SigmaRefer <- Reduce(rbind,res1sigma.list[m])
          
          # note use of [[1]] as is matrix rathe than list, 
        
          S11 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_nonmiss]
          #to ensure rows and cols as should reflect their stucture use matrix 
          S12 <-matrix(SigmaRefer[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss))
          S22 <-SigmaRefer[[1]][c_mata_miss,c_mata_miss]
         
          #edit due to passing  Sigma after loop, ie save sigma
          Sigma <- SigmaRefer[[1]]
          #print("J2R Sigma = ")
          #print(Sigma)
          #print ("S11 = " )
          #print(S11)
          #print("S12 = ")
          #print(S12)
        }
        else if (meth=='CR') {
          mata_means <- paramBiglist[[M*(referindex-1)+m]][1]
          #mata_means <- get(paste0("param",refer,m))[1]
          # convert from list to matrix
          mata_means <- mata_means[[1]]
          #mata_means_r <- unlist(mata_means_ref)*1
          #mata_means_r <- unlist(mata_means_ref)*mata_miss
          #mata_means <- mata_means_r
          # and preserve names   
          #colnames(mata_means) <- colnames(mata_means_trt)
          #replicate to number of rows defined by X1 
          #mata_means <- mata_means_r
          #colnames(mata_means) <- colnames(mata_means_ref)
          
          #mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]
          
          #SigmaRefer <- get(paste0("param",refer,m))[2]
          SigmaRefer <- paramBiglist[[M*(referindex-1)+m]][2]
          S11 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_nonmiss]
          S12 <-matrix(SigmaRefer[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss) )
          S22 <-SigmaRefer[[1]][c_mata_miss,c_mata_miss]
        }
        else if (meth=='CIR') {
          # pre-deviating use mean of trt gp up to last obs time bfore deviating, post-deviating use mean from ref grp 
          
          #mata_means_trt <- get(paste0("parambeta",trtgp,m)) 
          #mata_means_ref <- get(paste0("parambeta",refer,m))
          
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
          # needs to to ths as tibble will fail in cholesky
          #SigmaRefer <- as.matrix(get(paste0("paramsigmaStata",refer,m)))
          
          #print(paste0("paramsigma",refer,m))
          #SigmaRefer <- Reduce(rbind,res1sigma.list[m])
          
          S11 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_nonmiss]
          S12 <-matrix(SigmaRefer[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss))
          S22 <-SigmaRefer[[1]][c_mata_miss,c_mata_miss]
        }  
        else if (meth=='LMCF') { 
          
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
        }  #if meth
        
      #debug
      #  print(SigmaRefer[[1]])
      #  print(paste0("c_mata_nonmiss,",c_mata_nonmiss))
      #  print(c_mata_nonmiss)
      #  print(paste0("c_mata_miss,",c_mata_miss))
      #  print(c_mata_miss)
        
        # loop still open for row(mg)
        
        ###################### MNAR IMPUTATION ################
        # need insert routine when ALL missing values
        # else
        ###################### MNAR IMPUTATION ################
        
        
        #make sure these are single row vectors! as mistake in LMCF but have to be duplicate rows so add ,s
        #and move after dup fun
        #Error in mata_means[, c_mata_nonmiss] : incorrect number of dimensions
        #m1 <- mata_means[c_mata_nonmiss]
        #m2 <- mata_means[c_mata_miss]
        mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]
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
        
        # for debug
        #  print("startrow stoprow = ")
        #  print(startrow)
        #  print(stoprow)
        
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
           # print("U = ")
          #  print(U)
          #  print("miss_count = ")
          #  print (miss_count)
            # generate inverse normal, same as used below
            #for debug 20/01  
         #   print(paste0('mg[i,X1] =',mg[i,"X1"],' miss_count= ',miss_count)) 
            Z <- qnorm(matrix(runif( mg[i,"X1"]* miss_count,0,1),mg[i,"X1"],miss_count))
          #  print("Z")
          #  print(Z)
            mata_y1 = meanval+Z%*%(U) 
            #set dimensions mata_new to mata_y1 
            mata_new <- dim(mata_y1)
            GI <- array(data=mg[i,"treat"],dim=c(mg[i,"X1"],1))
            #II  no imputations 
            II <- array(data=m,dim=c(mg[i,"X1"],1))
            mata_new=cbind(GI,II,mata_new,SNO)
            mata_all_newlist[[m_mg_iter]]=mata_new
            
        } 
        
        # for debug    
        #print("mata_raw = ")
       # print("S11 S12")
      #  print(S11)
      #  print(S12)
      #  print(paste0("typeS12= ",class(S12)))
      #  print( dim(S12))
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
        if (meth=='J2R') {
          meanval = as.matrix(m2) + as.matrix(raw1 - m1)%*%as.matrix(t_mimix)
        } else  {
          meanval = (m2) + as.matrix(raw1 - m1)%*%as.matrix(t_mimix)
        } 
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
        
        # assigning the columns where the missing values  
        if(length(c_mata_miss)==0 ) { mata_new <- mata_Obs[,c(2:nct+1)]
        }else {
          mata_new[,c_mata_miss] <- mata_y1
        }
        # if(length(c_mata_miss)!=0 )
        
        #save U,Z for imp(m) and patt(i)
        # this worked for Z
        
        #what i this below? 
        #presumablywas  for testing?
        #assign(paste0("Z11",i,"_imp",m),subset(Z))
        #assign(paste0("U11",i,"_imp",m),subset(U))
        # assign(paste0("S12",i,"_imp",m),S12)
        # assign(paste0("S22",i,"_imp",m),S22)
        # assign(paste0("mata_y1_11",i,"_imp",m),mata_y1)
        # assign(paste0("conds11",i,"_imp",m),conds) 
        # no idea why below cuases error!
        # assign(paste0("t_mimix11",i,"_imp",m),t_mimix)
        
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
  
  return(list(mata_all_newlist,mg,M,meth))
  
} # for runmimix test


mi_impute <-function(idvar,timevar,depvar,covar) {
  #preprodata(depvar,treatvar,idvar,timevar,covar)
  tst<-pivot_wider(mxdata,id_cols=c(idvar),names_from=timevar,names_prefix=depvar,values_from=depvar)
  # add ont covariates and add on to resp list
  respvars<-c(names(tst[,-1]),covar)
  return(respvars)
}


# analse function to be used after main run to output  for summary stats
analse <- function(meth,no)  {
  assign( paste0("mata_all_new_rmna",meth), na.omit(mata_all_new))
  assign( paste0("mata_all_new_rmna",meth,no) , filter(get(paste0("mata_all_new_rmna",meth)),SNO == no))
  print(paste0("method= ",meth,"SNO = ",no))
  t(round(stat.desc(get(paste0("mata_all_new_rmna",meth,no))[,c("fev2","fev4","fev8","fev12")]),3)[c(1,9,13,4,8,5),])
}


#analselist slow so try optimse by pre-declaring saved data structure (as matrix) instead of rbind
analyselist <-function(meth,no) {
   subSNOx <- head(mata_all_newlist[[1]][[1]],1) 
   subSNOx[subSNOx>=0] <-NA
  # subSNOxmata_newlist <- vector('list',M*nrow(mg))
   #subSNOx_newMatrix <- matrix(, nrow=(M*nrow(mg)), ncol=ncol(subSNOx) )
for (i in 1:(nrow(mg[[1]])*M) )  {
  subSNO <- subset((mata_all_newlist[[1]][[i]]),SNO== no)
  #for (j in 1:ncol(subSNO)) {
  #subSNOx_newMatrix[i,j]=subSNO[i,j]
  #}
  subSNOx<- rbind(subSNOx,subSNO)
  subSNOx<-na.omit(subSNOx)
  
}
print(paste0("meth = ",meth))  
#t(round(stat.desc(subSNOx)[,c("head3","head12","head_base")],3)[c(1,9,13,4,8,5),])
t(round(stat.desc(subSNOx)[,c("fev2","fev4","fev8","fev12")],3)[c(1,9,13,4,8,5),])
}




# now changed to mata_all_newlist
#for (i in 1:M)


# mata_S_miss something like [2 3 4] ,so is cc
CIR_loop <- function(c_mata_miss,mata_Means,MeansC)
{    
  miss_count <- length(c_mata_miss)
  mata_means <- mata_Means
  
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

pttestf<- function(n1,n2,mn1,sd1,mn2,sd2) {
  pttest = pt((((mn1 - mn2) -0) /sqrt(sd1^2/n1+sd2^2/n2)),(n1+n2-2))
  return(pttest)
}

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


#works
testread <-function(pathdat) {
  #mxdata <-read.csv("./asthma.csv")
  txdata <-read.csv(pathdat)
}
 #txdata <-testread("asthma.csv")
 

 
  


