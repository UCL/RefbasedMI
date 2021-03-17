#' Sample data: asthma trial
#' 
#' A data set containing  asthma trial data as used in the Stata mimix help file
#' The primary outcome variable is fev, measured at 2,4,8,12 weeks
#'
#'
#' @format A data frame containing 732 rows and 5 columns
#'   \describe{
#'   \item{id}{patient identifier}
#'   \item{time}{}
#'   \item{treat}{}
#'   \item{base}{covariate}
#'   \item{fev}{outcome variable}
#'   }
#' @examples
#' \dontrun{
#'  impJ2Rridge<-(mimix(data=asthma,covar=c("base"),depvar=fev,treatvar=treat,idvar=id,timevar=time,
#'     method="J2R",reference=1,delta=c(0.5,0.5,1,1 ),M=5,seed=101,prior=ridge,burnin=1000)
#'  library(mice) 
#'  fitJ2R<-with(data= as.mids(subset(impJ2Rridge[[2]],time==12)),lm(fev~treat+base))    
#'  summary(pool(fitJ2R))    
#' }     
"asthma"



#' Sample data: antidepressant trial
#' 
#' A data set containing antidepressant trial data as described in paper 
#'   by White,Royes,Best (2019)
#' The primary outcome is HAMD17.TOTAL measured at visit number 4,5,6,7.  
#'
#' @format  dataframe containing 688 rows and 14 columns
#'  \describe{
#'  \item{PATIENT.NUMBER}{}
#'  \item{HAMA.TOTAL}{}
#'  \item{PGI_IMPROVEMENT}{}
#'  \item{VISIT...VISIT.3.DATE}{}
#'  \item{VISIT.NUMBER}{}
#'  \item{TREATMENT.NAME}{}
#'  \item{PATIENT.SEX}{}
#'  \item{POOLED.INVESTIGATOR}{}
#'  \item{basval}{}
#'  \item{HAMD17.TOTAL}{outcome variable}
#'  \item{change}{}
#'  \item{miss_flag}{}
#'  \item{methodcol}{individual-specific method}
#'  \item{referencecol}{individual-specific reference arm}
#'  }
#' @examples
#' \dontrun{
#'  # Run with  covariates "basval" and "PATIENT.SEX" using columns within data to specify
#'  # method and reference 
#'  impIndiv <- mimix(data=antidepressant,covar=c("basval","PATIENT.SEX"),depvar=HAMD17.TOTAL,
#'         treatvar=TREATMENT.NAME,idvar=PATIENT.NUMBER,
#'          timevar=VISIT.NUMBER,methodvar="methodcol",referencevar="referencecol",M=5,seed=54321)
#'  library(mice)+basval+ression
#'  fit<-with(data= as.mids(impIndiv[[1]]),lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+PATIENT.SEX))
#'  summary(pool(fit))  
#'  impantdep <- mimix(data=antidepressant,covar=c("basval","PATIENT.SEX"),depvar=HAMD17.TOTAL,
#'          treatvar=TREATMENT.NAME,idvar=PATIENT.NUMBER,
#'          timevar=VISIT.NUMBER,method="J2R",reference=1,M=2,seed=54321)
#'  fitdep21<-with(data= as.mids(subset(impantdep[[2]],VISIT.NUMBER==7)),
#'               lm(HAMD17.TOTAL~TREATMENT.NAME))
#'  summary(pool(fitdep21)) 
#' }
"antidepressant"


#' Sample data: acupuncture trial 
#'
#' A data set containing results of a randomised, double-blind, parallel-group 
#'     comparing active treatment with placebo
#' The primary outcome is head, measured at time 3 and 12 
#'
#'
#' @format A data frame with 802 rows and 11 columns
#'  \describe{
#'  \item{id}{}
#'  \item{time}{}
#'  \item{age}{}
#'  \item{sex}{}
#'  \item{migraine}{}
#'  \item{chronicity}{}
#'  \item{practice_id}{}
#'  \item{treat}{}
#'  \item{head_base}{covariate}
#'  \item{head}{outcome variable}
#'  \item{withdrawal_reason}{}
#'  }
#' @examples 
#' \dontrun{
#'  impCausalref1 <- mimix(data=acupuncture,covar=c("head_base","sex"),depvar=head,treatvar=treat,
#'        idvar=id,timevar=time,
#'        method="Causal",reference=1,K0=1,K1=0.5,M=5,seed=54321)
#'  library(mice)         
#'  fitacup<-with(as.mids(subset(impCausalref1[[2]],time==12)), lm(head~treat+head_base+sex))
#'  summary(pool(fitacup))    
#' }
"acupuncture"

# eample output file 
# impdataCausal.
# 
#  A dataset containing imputed values obtained from mimix on asthma data 
#  
#   @format A data frame with 1098 rows and 9 columns
#   \describe{
#   \item{.imp}{imputation number}
#   \item{base}{covariate}
#   \item{fev.2}{response at time 2}
#   \item{fev.4}{response at time 4}
#   \item{fev.8}{response at time 8}
#   \item{fev.12}{response at time 12}
#   \item{patt}{missing data pattern}
#   \item{treat}{treatment group}
#   \item{.id}{case id}
#  }



# "mpdataCausal"   