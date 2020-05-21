

#' @title preprodata
#' @description process data into wide format for group method
#' @details checks method finds missingness pattern
#' @param data  data in long format
#' @param covar covariates and base depvar complete
#' @param depvar dependent variable
#' @param treatvar treatment group
#' @param idvar patient id
#' @param timevar time variable for repeated visit
#' @param M number imputations
#' @param refer refernce group
#' @param meth RBI method
#'
#'
#' @return list of outputs


#preprocess data for group method
preprodata<- function(data,covar,depvar,treatvar,idvar,timevar,M,refer,meth=NULL)  {
  #extract relevant vars

  browser()
  # if covar null (or 1st depvar complete) then ceate baseval covar

  #if (length(covar) ==0) {

  if (length(covar) ==999) {

      #get 1st row each id group


    rowgp<-get(data)[,head(.SD,1), by=idvar]
    rowgpx<- as.data.frame(rowgp)[,c(idvar,depvar)]
    names(rowgpx)[names(rowgpx)==(depvar)]<-"baseval"
    #drop 1st row each id group ,assuming 1st time 0

    depetedata<- get(data)[!(get(data)[,get(timevar)]==0),]
    # now need to merge 1st row back onto depleted df creating new covar var
    data<-merge(depetedata,rowgpx, by=(idvar))
    # now create covar
    covar <- c("baseval")
    # depetedasthma<- testdfasthma[!(testdfasthma$time==0),]
    #move from below
    fevdata<- as.data.frame(data)[,c(idvar,depvar,timevar)]
    uniqdat<-  unique(as.data.frame(data)[c(idvar,covar,treatvar)])
    ntreatcol<- as.data.frame(data)[c(treatvar)]
    ntimecol<- as.data.frame(data)[c(timevar)]
    # this incomplete for now - need to complete later.
 }# else {
    # fevdata<-dplyr::select(get(data),idvar,depvar,timevar)
    fevdata<- get(data)[c(idvar,covar,depvar,timevar,treatvar)]
    # now covar added to data list so need need for unique?
    uniqdat<-unique(get(data)[c(idvar,covar,treatvar)])
    ntreatcol<- get(data)[c(treatvar)]
    ntimecol<- get(data)[c(timevar)]
 # }

  # browser() for checking covar Null has to be done for indiv case as well below!
  #  (is.na(covar))
  ### checking covars complete, ie non- missing
  #tryCatch(stopifnot(sum(is.na(mxdata[,covar]))==0,error=stop("Error: not all covariates are complete !!")))

  # this errr catcher invalid when cova NULL
  # stopifnot(sum(is.na(get(data)[,covar]))==0)

  #more informative in error msg to use this explicit and
  #and put in one statement
  #put in Runmimix program

#

 # fevdata<-dplyr::select(get(data),idvar,depvar,timevar)
#  fevdata<- get(data)[c(idvar,depvar,timevar)]
  # only want to widen the dose var, ie fev, so take the other variables
  #reshape to wide and assign new names
  # prefix must be supplied from input argument rather than hard coded, this canbe hard coded

  #replace tidyr function
  #sts4<-tidyr::pivot_wider(fevdata,id_cols=c(idvar),names_from=timevar,names_prefix=depvar,values_from=depvar)
  sts4<-stats::reshape(fevdata,v.names = depvar,timevar = timevar,idvar=idvar,direction="wide")
  #sts4 is just response data, can join later treat and covar cols from finaldat

  #assumes no NA for these vars should add in checking routines !!
  #check how many covars used

  #nocovar=covar[[1]]

  #tcovar<-unlist(covar)
  #uniqdat<-unique(mxdata[c(idvar,unlist(tstcovar),treatvar)])
  print(paste0("covar=",covar))
#  uniqdat<-unique(get(data)[c(idvar,covar,treatvar)])
  print(utils::head(uniqdat))

  finaldatOld<- merge(sts4,uniqdat,by=idvar)
  # now no need to merge because covar ,treat already specified in fevdata
  finaldat<- sts4


  print(utils::head(finaldat))
  # now try and sort on treatvar, doesnt work so instead of sorting , select on treat


  # need find no. ntreat  to loop over
  #ntreat<-unique(unique(mxdata$treatvar))
 # ntreatcol<-(dplyr::select(get(data),treatvar))
#  ntreatcol<- get(data)[c(treatvar)]
  ntreat <- unique(ntreatcol)
 # ntimecol<-(dplyr::select(get(data),timevar))
#  ntimecol<- get(data)[c(timevar)]
  ntime<-unique(ntimecol)


  #in order to aggregate by pattern need create dummy vars

  # for sts4 array drop 1st col and create newvars for rest
  #STSdummy<-apply(sts4[,2:ncol(sts4)],MARGIN=2,FUN=A)
  #STSdummy<-apply(sts4[,2:ncol(sts4)],MARGIN=2,function(x) ifelse(!is.na(x),0,1))
  #select out the depvar.time variables to calc patt

    # instead of below (covars could have name of depvar)
  #STSdummy<- apply(sts4[,grepl(paste0(depvar,"."),names(sts4))],MARGIN=2,function(x) ifelse(!is.na(x),0,1))

  testfevdata<- get(data)[c(idvar,depvar,timevar)]
  sts4dummy<-stats::reshape(testfevdata,v.names = depvar,timevar = timevar,idvar=idvar,direction="wide")
  STSdummy<- apply(sts4dummy[,grepl(depvar,names(sts4dummy))],MARGIN=2,function(x) ifelse(!is.na(x),0,1))

#9/5/20
#browser()

  # append to names  try paste0 as otherwise space before miss
  colnames(STSdummy) <- paste0(colnames(STSdummy),'.miss')
  # merge back on to data
  sts4D<-cbind(finaldat,STSdummy)
  # create powrs of 2
  pows2 <- sapply(1:ncol(STSdummy),function(i) STSdummy[,i]*2^(i-1))
  #need to add up to find patt
  patt <- rowSums(pows2)
  #sts4Dpatt<-cbind(sts4D,patt)
  sts4Dpatt<-cbind(sts4D,patt)

  # but also add on covars (should be no missing so dont contribute to patt)
  # ie should be  cols of 0's
  if (length(covar) !=0) {
    tmp_covpatt<-apply(as.data.frame(sts4Dpatt[,covar]),MARGIN=2,function(x) ifelse(!is.na(x),0,1))
    #add names
    colnames(tmp_covpatt) <- paste0(c(covar),".miss")
    # than combine below the dummies onto the finaldat
  #8/5/20sts4Dpatt<-cbind(sts4D,tmp_covpatt,patt)

    sts4Dpatt<-cbind(finaldat,tmp_covpatt,STSdummy,patt)
  }   else {
    sts4Dpatt<-cbind(finaldat,STSdummy,patt)
  }


  # In order  to get patt with all missing patts

  Overall_patt<-unique(sts4Dpatt[grepl(".miss",colnames(sts4Dpatt))])
  patt<- unique(sts4Dpatt[,"patt"])

  #then combine depvar and covar patt!
  all_patt<-cbind(Overall_patt,patt)


  # then sort by patt( last col) and treat and find cumulative sequence  AND split on treatment var!
  # so need to merge treat (and covar) back into response data

  #sort by treatvar and patt var, couldnt get treatvar working with patt so try sorting previously on treatvar
  finaldatSS <-sts4Dpatt[order(patt),]

  ndx = order(finaldatSS[,treatvar])
  finaldatS <- finaldatSS[ndx,]


  #need sort by patt and align with mg (lookup) ,t oget looku right merge using patt
  finaldatS <- sts4Dpatt[order(sts4Dpatt[,treatvar],sts4Dpatt[,"patt"]),]

  # consistent with Stata to get right order, drop id and treat cols
  drops <-c(idvar,treatvar)


  #also get the dummy pattern vectors by treatment
  #dummypatt<-unique(xSTSdummy)
  pattmat<-unique(STSdummy[,1:ncol(STSdummy)])

  #
  pows2d <- sapply(1:ncol(pattmat),function(i) pattmat[,i]*2^(i-1))
  patt<-rowSums(pows2d)
  #want  join this to ex from mimix_group table
  #merge(uniqdat,sts4,by=idvar)
  pattmatd<-cbind(pattmat,patt)


  # and now find cumX1

  sts4Dpatt$X1<-1

  ex1<-Hmisc::summarize(sts4Dpatt$X1, by=Hmisc::llist(sts4Dpatt[,treatvar],sts4Dpatt$patt),FUN=sum)
  #to rename
  newnames <- c( treatvar,"patt","X1")
  names(ex1)<-newnames
  # rename(ex1,treat="sts4Dpatt[, c(treatvar)]",patt='sts4Dpatt$patt',methodvar=`sts4Dpatt[, c(methodindiv[1])]`,referencevar=`sts4Dpatt[, c(methodindiv[2])]`)
  # and now find cumX1
  ex1$X1cum <- cumsum(ex1$X1)

   ex1$X1cum <- cumsum(ex1$X1)
  #want  join this to ex from mimix_group table
  # create index number to preserve ordr after merge
  ex1$exid <- 1:nrow(ex1)
  ex1id <-merge(ex1,pattmatd,by="patt")
  ex1s<-ex1id[order(ex1id$exid),]

#8/5/20
  test_ex1<-merge(ex1,all_patt,by="patt")[order(merge(ex1,all_patt,by="patt")$exid),]

  stopifnot(refer %in% t(ntreat))
  print("summary missing pattern")
  print(ex1s)
  #  return(list(sts4Dpatt,finaldatS,finaldat,ntreat,ex1s,ex1,ex1id,pattmat,patt,ntime,M,refer,meth))

#  return(list(sts4Dpatt,finaldatS,ntreat,ex1s,ex1,pattmat,patt,ntime,M,refer,meth))

 # return(list(sts4Dpatt,finaldatS,ntreat,ex1s,pattmat,patt,ntime,M,refer,meth))
 #8/5/20

  #so finaldatS is sorted by treat,patt  data
  # test_ex1 is mg lookup table for missing dummies
  # all_patt is missing pattern table
    return(list(finaldatS,ntreat,test_ex1,all_patt,ntime,M,refer,meth))

}

#' @title preproIndivdata
#' @description process data into wide format for individual-specified method
#' @details checks methodindiv finds missingness pattern
#' @param data  data in long format
#' @param covar covariates and base depvar complete
#' @param depvar dependent variable
#' @param treatvar treatment group
#' @param idvar patient id
#' @param timevar time variable for repeated visit
#' @param M number imputations
#' @param refer reference group must be NULL
#' @param meth RBI method  must be NULL
#' @param methodindiv list of c(method, reference) column location in data specifying individual RBI methods and references
#' @return list of outputs




#preprocess data for individual method
preproIndivdata<- function(data,covar,depvar,treatvar,idvar,timevar,M,refer=NULL,meth=NULL,methodindiv)  {
  #extract relevant vars
# 11/05/20
#browser()
  #check covars complete
  stopifnot(sum(is.na(get(data)[,covar]))==0)
  #tryCatch(stopifnot(sum(is.na(mxdata[,covar]))!=0,error=stop("Error: not all covariates are complete !!")))
  #more informative in error msg to use this explicit and
  #and put in one statement
  methodL <- unique(get(data)[,methodindiv[1]])
  # stopifnot( (methodL == "MAR" | methodL=="j2r" | methodL=="cir" | methodL=="CR" | methodL=="LMCF"| methodL=="null"),
  #is.numeric(refer),
  #     is.numeric(M),
  #     is.character(depvar),
  # treatvar could be char or numeric ? then refer must be same type
  #is.character(treatvar),
  #     is.character(idvar),
  #     is.character(timevar),
  #is.character(covar[[2]],
  #     is.character(covar) )
  # convert to numic should be done outside function in main.

  #11/05/20
  fevdata<- get(data)[c(idvar,covar,depvar,timevar,treatvar,methodindiv)]
  # now covar added to data list so need need for unique?
  uniqdat<-unique(get(data)[c(idvar,covar,treatvar)])
  ntreatcol<- get(data)[c(treatvar)]
  ntimecol<- get(data)[c(timevar)]


  #select out the id and depvar to generate names for depvar#time
  #fevdata<-dplyr::select(get(data),idvar,depvar,timevar)
##  fevdata<-get(data)[c(idvar,depvar,timevar)]

  #generate names depvar#time
  #sts4<-tidyr::pivot_wider(fevdata,id_cols=c(idvar),names_from=timevar,names_prefix=depvar,values_from=depvar)
  sts4<-stats::reshape(fevdata,v.names = depvar,timevar = timevar,idvar=idvar,direction="wide")

  # assumes the covars all non-missing
## uniqdat<-unique(get(data)[c(idvar,covar,treatvar,methodindiv)])

  # merge on the covariates and treatment to names
##  finaldat<- merge(sts4,uniqdat,by=idvar)
  finaldat<-sts4

  #in order to aggregate by pattern need create dummy vars

  # for sts4 array drop 1st col and create newvars for rest thes are the dummy missing patterns
##  STSdummy<-apply(sts4[,2:ncol(sts4)],MARGIN=2,function(x) ifelse(!is.na(x),0,1))
  # append to names
##  colnames(STSdummy) <- paste(colnames(STSdummy),'.missing')
  # merge back on to data
##  print(utils::head(STSdummy))
##  sts4D<-cbind(finaldat,STSdummy)
  # create powrs of 2

  #11/05/20
  testfevdata<- get(data)[c(idvar,depvar,timevar)]
  sts4dummy<-stats::reshape(testfevdata,v.names = depvar,timevar = timevar,idvar=idvar,direction="wide")
  STSdummy<- apply(sts4dummy[,grepl(depvar,names(sts4dummy))],MARGIN=2,function(x) ifelse(!is.na(x),0,1))

  # append to names, try paste0 12/5/20
  colnames(STSdummy) <- paste0(colnames(STSdummy),'.miss')

  # merge back on to data  ,ie fevdata then testfevdata
  sts4D<-cbind(finaldat,STSdummy)

  # create powrs of 2


  pows2 <- sapply(1:ncol(STSdummy),function(i) STSdummy[,i]*2^(i-1))
  #need to add up to find patt
  patt <- rowSums(pows2)
  sts4Dpatt<-cbind(sts4D,patt)
  print(utils::head(sts4Dpatt))

  # if covars specified...
  if (length(covar) !=0) {
    tmp_covpatt<-apply(as.data.frame(sts4Dpatt[,covar]),MARGIN=2,function(x) ifelse(!is.na(x),0,1))
    #add names
    colnames(tmp_covpatt) <- paste0(c(covar),".miss")

    sts4Dpatt<-cbind(finaldat,tmp_covpatt,STSdummy,patt)
  }   else {
    sts4Dpatt<-cbind(finaldat,STSdummy,patt)
  }

  # In order  to get patt with all missing patts

  Overall_patt<-unique(sts4Dpatt[grepl(".miss",colnames(sts4Dpatt))])
  patt<- unique(sts4Dpatt[,"patt"])

  #then combine depvar and covar patt!
  all_patt<-cbind(Overall_patt,patt)






  #patt has just been created so use $ notation wheras treatvar from input argument and need to by methodindiv[1] as well

  finaldatSS <-sts4Dpatt[order(patt),]

  ndx = order(finaldatSS[,treatvar],finaldatSS[,methodindiv[1]],finaldatSS[,methodindiv[2]])
  finaldatS <- finaldatSS[ndx,]


  #12/05/20
  #try this
  #need sort by patt and align with mg (lookup) ,t oget looku right merge using patt
  #finaldatS <- sts4Dpatt[order(sts4Dpatt[,treatvar],sts4Dpatt[,"patt"]),]
  # has to be also sorted into methodindiv 1 and 2 grps
  finaldatS <- sts4Dpatt[order(sts4Dpatt[,treatvar],sts4Dpatt[,c(methodindiv[1])],sts4Dpatt[,c(methodindiv[2])],sts4Dpatt[,"patt"]),]

  # consistent with Stata to get right order, drop id and treat cols
  drops <-c(idvar,treatvar)


  #also get the dummy pattern vectors by treatment
  #dummypatt<-unique(xSTSdummy)
  pattmat<-unique(STSdummy[,1:ncol(STSdummy)])


  pows2d <- sapply(1:ncol(pattmat),function(i) pattmat[,i]*2^(i-1))
  patt<-rowSums(pows2d)
  #want  join this to ex from mimix_group table
  #merge(uniqdat,sts4,by=idvar)
  pattmatd<-cbind(pattmat,patt)
  print(pattmatd)



  #
  #X1 is a count variable of 1's'
  sts4Dpatt$X1<-1

  ex1<-Hmisc::summarize(sts4Dpatt$X1, by=Hmisc::llist(sts4Dpatt[,treatvar],sts4Dpatt[,c(methodindiv[1])],sts4Dpatt[,c(methodindiv[2])],sts4Dpatt$patt),FUN=sum)
 #to rename
  newnames <- c( treatvar,"methodvar","referencevar","patt","X1")
  names(ex1)<-newnames
  # rename(ex1,treat="sts4Dpatt[, c(treatvar)]",patt='sts4Dpatt$patt',methodvar=`sts4Dpatt[, c(methodindiv[1])]`,referencevar=`sts4Dpatt[, c(methodindiv[2])]`)
  # and now find cumX1
  ex1$X1cum <- cumsum(ex1$X1)

  # create index number to preserve ordr after merge  , #want  join this to ex from mimix_group table
  ex1$exid <- 1:nrow(ex1)
  ex1id <-merge(ex1,pattmatd,by="patt")
  ex1s<-ex1id[order(ex1id$exid),]
  print(ex1)
  print(ex1s)
  # dont need finaldat,exlid

  # pattmat is just the missing dummies
  # patt vector of patterns
  # need find no. ntreat  to loop over
  #ntreat<-unique(unique(mxdata$treatvar))
 # ntreatcol<-(dplyr::select(get(data),treatvar))
  ntreatcol<-get(data)[c(treatvar)]
  ntreat <- unique(ntreatcol)
 # ntimecol<-(dplyr::select(get(data),timevar))
  ntimecol<-get(data)[c(timevar)]
  ntime<-unique(ntimecol)


  #11/05/20
  test_ex1<-merge(ex1,all_patt,by="patt")[order(merge(ex1,all_patt,by="patt")$exid),]

  #error chk

  # find unique values for referencevar to check against ntreat values
  refencevars <- unique(get(data)[,methodindiv[2]])
  stopifnot(refencevars %in% t(ntreat))
  #print("summary missing pattern")
#remove ex1 seems to cause problems!

  # main outputs have to be test_ex1 and finaldatS which should correspond
  return(list(finaldatS,ntreat,test_ex1,all_patt,pattmat,patt,ntime,M,methodindiv))
  # return(list(sts4Dpatt,finaldatS,ntreat,ex1s,pattmat,patt,ntime,M,methodindiv))
}



#' @title ifmethodindiv
#' @description alternative logic for individual method
#' @details checks methodindiv not null
#' @param methodindiv  c(method,reference) cols
#' @param mg  pattern lookup table
#' @param m where we are in the imputations
#' @param M number of total imputations.
#' @param paramBiglist  list of Beta and Sigma parameters from mcmc
#' @param i in loop throug mg rows
#' @param treatvar treatment group
#' @param c_mata_nonmiss  0,1 vector of nonmissing
#' @param c_mata_miss 0,1 vector of missing
#' @param mata_miss position of missing values in repeated time visits
#' @param mata_nonmiss positions of nonmissing values
#' @return mata_means



# alternative logic for individual method
ifmethodindiv <- function(methodindiv,mg,m,M,paramBiglist,i,treatvar, c_mata_nonmiss,c_mata_miss,mata_miss,mata_nonmiss)
{


  #17/03
  # browser()
  # assign paramBiglist to inividual
  #for j in 1:val
  #assign(paste0("paramBiglist",val,"_",m),paramBiglist) # check index correct!

  #} else
  #browser()
  if(!is.na(methodindiv[1])) {
    # methodindiv needs editing this bit
    trtgp <- mg[i,treatvar]
    #the refernce group comes from the indvidual colunm!
    refergp <- mg[i,methodindiv[2]]
  }

  # without this (mg[i,methodindiv[1]]) will be cir etc
  methindiv<- mg[i,methodindiv[1]]

  methindiv<-  ifelse( ( methindiv=="j2r" | methindiv=="J2R" |methindiv=="j2R"|methindiv=="J2r" ),3,
                       ifelse( ( methindiv=="CR" | methindiv=="cr" |methindiv=="Cr"|methindiv=="cR" ),2,
                               ifelse( ( methindiv=="MAR" | methindiv=="mar" |methindiv=="Mar"|methindiv=="MAr"|methindiv=="Mr"|methindiv=="MR" ),1,
                                       ifelse( ( methindiv=="CIR" | methindiv=="cir" |methindiv=="CIr"|methindiv=="cliR" ),4,
                                               ifelse( ( methindiv=="LMCF" | methindiv=="lmcf" |methindiv=="Last"|methindiv=="last" ),5,9)))))

#browser()
  # only done methods  3,4 so far, so Mar need correcting
  #MAR
  if  (methindiv== 1)  {
    #mata_means <- get(paste0("param",trtgp,m))[1]
    mata_means <- paramBiglist[[M*(trtgp-1)+m]][1]
    #  mata_means <- get(paste0("paramBiglist",trtgp,"_",m))[1]

    # mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]
    # convert from list element to matrix
    #    mata_means <- mata_means[[1]]


    #Sigmatrt <- get(paste0("param",trtgp,m))[2]
    Sigmatrt <- paramBiglist[[M*(trtgp-1)+m]][2]
    #Sigmatrt <- get(paste0("paramBiglist",trtgp,"_",m))[2]
    Sigma <- Sigmatrt


  }
  # 'J2R'
  else if ( methindiv == 3)  {

    # changed saving the result into  just the param file, list of 2 so can use list index here
    #treatmnets are 1.. M then M+1 ..2M .. etc

    #6/3/20 need editing after chnaging suffixes in paramBiglist
    mata_means_trt <- paramBiglist[[M*(trtgp-1)+m]][1]
    mata_means_ref <- paramBiglist[[M*(refergp-1)+m]][1]
    #mata_means_ref <- paramBiglist[[M*(refergp-1)+m]][1]

    # 9/3/20  need editing after chnaging suffixes in paramBiglist
    #mata_means_trt <- get(paste0("paramBiglist",trtgp,"_",m))[1]
    #mata_means_ref <- get(paste0("paramBiglist",refergp,"_",m))[1]


    # below causes error after using >1 covars and mata_nonmiss has covar.1, not proper covar names
    # edits 12/04 same as in runmimix
    mata_means_t <- lapply(mata_means_trt,FUN = function(x) x*mata_nonmiss)
    #mata_means_t <- unlist(mata_means_trt)*mata_nonmiss
    # print(paste0("mata_means_trt, mata_nonmiss= ",mata_means_trt,mata_miss))

    mata_means_r <- lapply(mata_means_ref, FUN = function(x) x*mata_miss)
    #mata_means_r <- unlist(mata_means_ref)*mata_miss
    # so when all missing  1,1,1, ... then all contributing comes from reference means
    mata_means <- unlist(mata_means_r)+unlist(mata_means_t)
    mata_means <- (as.matrix(t(mata_means)))
    # and preserve names
   # colnames(mata_means) <- colnames(mata_means_trt)
    #replicate to number of rows defined by X1
    #mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]


    ############# SIGMA is from paramsigma  the reference group ################


    # do we ever  use SigmaTrt !!?? in j2r?
    # answer is yes, need to use SigmaTrt for the predeviation observations, ie up to where they go missing
    # only after they go missing (trailing missing) need to use the SigmaRef

    #SigmaRefer <- get(paste0("paramBiglist",refergp,"_",m))[2]
    SigmaRefer <- paramBiglist[[M*(refergp-1)+m]][2]
    Sigma <- SigmaRefer
    #Sigmatrt <- get(paste0("paramBiglist",trtgp,"_",m))[2]


    # note use of [[1]] as is matrix rathe than list,

    S11 <-SigmaRefer[[1]][c_mata_nonmiss,c_mata_nonmiss]
    # causes non-def error in conds
    # S11 <-Sigmatrt[[1]][c_mata_nonmiss,c_mata_nonmiss]
    #to ensure rows and cols as should reflect their stucture use matrix
    S12 <-matrix(SigmaRefer[[1]][c_mata_nonmiss,c_mata_miss],nrow=length(c_mata_nonmiss))
    S22 <-SigmaRefer[[1]][c_mata_miss,c_mata_miss]
    # prob have to change this to loop to fill in matrix


  }
  #else if (meth=='CR') {
  else if ( methindiv == 2)
  {
    # no need to use Sigmatrt here
    mata_means <- paramBiglist[[M*(refergp-1)+m]][1]


    #mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]

    #SigmaRefer <- get(paste0("param",refer,m))[2]
    SigmaRefer <- paramBiglist[[M*(refergp-1)+m]][2]
    # SigmaRefer <- get(paste0("paramBiglist",refergp,"_",m))[2]
    Sigma <- SigmaRefer

  }
  #else if (meth=='CIR')
  else if ( methindiv == 4)
  {
    # need to use Sigmatrt as in j2r
    # pre-deviating use mean of trt gp up to last obs time bfore deviating, post-deviating use mean from ref grp

    #
    mata_Means<-paramBiglist[[M*(trtgp-1)+m]][1]

    # convert from list to matrix
    # mata_Means <- mata_Means[[1]]
    #mata_Means <-  get(paste0("parambeta",trtgp,m))
    #MeansC <-  get(paste0("param",refer,m))[1]
    # MeansC <-  get(paste0("paramBiglist",refergp,"_",m))[1]
    MeansC<- paramBiglist[[M*(refergp-1)+m]][1]

    #might be better to copy mimix algol

    mata_means<-CIR_loop(c_mata_miss,mata_Means,MeansC)
    #returns mata_means as single row
    # then duplicate over patt rows
    #replicate to number of rows defined by X1
    # mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]

    #SigmaRefer <- get(paste0("paramsigma",refer,m))

    #SigmaRefer <- get(paste0("param",refer,m))[2]
    #SigmaRefer <- get(paste0("paramBiglist",refergp,"_",m))[2]
    SigmaRefer<- paramBiglist[[M*(refergp-1)+m]][2]
    Sigma <- SigmaRefer
    # when reading in Stata sigmas
    # needs to to ths as tibble will fail in cholesky
    #SigmaRefer <- as.matrix(get(paste0("paramsigmaStata",refer,m)))

    #print(paste0("paramsigma",refer,m))
    #SigmaRefer <- Reduce(rbind,res1sigma.list[m])


  }
  # else if (meth=='LMCF')
  else if ( methindiv == 5) {

    #mata_Means <-  get(paste0("param",trtgp,m))[1]
    mata_Means <- paramBiglist[[M*(trtgp-1)+m]][1]
    #mata_Means <- get(paste0("paramBiglist",trtgp,"_",m))[1]
    # convert from list to matrix
    # mata_Means <- mata_Means
    # no ref MeansC <- mata_means_ref
    mata_means<-LMCF_loop(c_mata_miss,mata_Means)
    #mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]


    #Sigmatrt <- get(paste0("param",trtgp,m))[2]
    Sigmatrt <- paramBiglist[[M*(trtgp-1)+m]][2]
    #Sigmatrt <- get(paste0("paramBiglist",trtgp,"_",m))[2]
    Sigma <- Sigmatrt


  }

  return(list(mata_means,Sigma))
}



