
# try looping as in mimix
# best way to agree I'm doing same thing as mimix

# therefore looping thru matrices




if (meth == 'CIR') { 
  
  mata_means_trt <- get(paste0("parambeta",trtgp,m)) 
  mata_means_ref <- get(paste0("parambeta",refer,m))
  
  # put equiv to mimix 
  mata_Means <-  mata_means_trt
  MeansC <- mata_means_ref
  
  #mata:cc =mata_S_miss equiv to c_mata_miss  and miss_count is no.of missing fields
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
      
        mata_means[b] <- MeansC[b] 
     } else  {
        mata_means[c_mata_miss[b]] = mata_means[(c_mata_miss[b]-1)]+ MeansC[(c_mata_miss[b])]- MeansC[(c_mata_miss[b])-1]
     } 
  }
     return(mata_means)
  }
  # then duplicate over patt rows
  #replicate to number of rows defined by X1 
  mata_means<-mata_means[rep(seq(nrow(mata_Means)),each=mg$X1[i]),]
  
   
  #LMCF
  #Last Mean carried forward
  # put equiv to mimix 
  mata_Means <-  mata_means_trt
  # no ref MeansC <- mata_means_ref
  LMCF_loop <- function(c_mata_miss,mata_Means)
     {
    miss_count <- length(c_mata_miss)
    mata_means <- mata_Means
    for (b in 1:miss_count)  {
      if (c_mata_miss[b] > 1) {
        mata_means[c_mata_miss[b]] = mata_means[(c_mata_miss[b]-1)]      
        } 
    }
    return(mata_means)
  }
  
  
  #J2R  Jump to refernce
  mata_Means <-  mata_means_trt
  MeansC     <-  mata_means_ref
  J2R_loop <- function(c_mata_miss,mata_Means,MeansC)
  {
      
      miss_count <- length(c_mata_miss)
      mata_means <- mata_Means
        for (b in 1:miss_count)  {
          mata_means[c_mata_miss[b]] =  MeansC[(c_mata_miss[b])]
        } 
      return(mata_means)
  }
  
  
  
  #one way is to element multiply (because 1,0) then add 
  mata_means_t <- unlist(mata_means_trt)*mata_nonmiss
  mata_means_r <- unlist(mata_means_ref)*mata_miss
  # so when all missing  1,1,1, ... then all contributing comes from reference means      
  mata_means <- mata_means_r+mata_means_t
  # and preserve names   
  colnames(mata_means) <- colnames(mata_means_trt)
  #replicate to number of rows defined by X1 
  mata_means<-mata_means[rep(seq(nrow(mata_means)),each=mg$X1[i]),]
