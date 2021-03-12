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
#'   \item{fev}{dependent variable}
#'   }
#' @examples
#' \dontrun{
#'  impJ2Rridge<-(mimix(data=asthma,covar=c("base"),depvar=fev,treatvar=treat,idvar=id,timevar=time,
#'     method="J2R",reference=1,delta=c(0.5,0.5,1,1 ),M=5,seed=101,prior=ridge,burnin=1000)
#'  library(mice)     
#'  fit<-with(data= as.mids(impJ2Rridge),lm(fev.12~treat))
#'  summary(pool(fit))    
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
#'  \item{HAMD17.TOTAL}{dependent variable}
#'  \item{change}{}
#'  \item{miss_flag}{}
#'  \item{methodcol}{individual-specific method}
#'  \item{referencecol}{individual-specific reference arm}
#'  }
#' @examples
#' \dontrun{
#'  impIndiv <- mimix(data=antidepressant,covar=c("basval","PATIENT.SEX"),depvar=HAMD17.TOTAL,
#'  treatvar=TREATMENT.NAME,idvar=PATIENT.NUMBER,
#'  timevar=VISIT.NUMBER,methodvar="methodcol",referencevar="referencecol",M=5,seed=54321)
#'  library(mice)  
#'  fit<-with(data= as.mids(impantiIndivDt),lm(HAMD17.TOTAL.7~TREATMENT.NAME+basval+PATIENT.SEX))
#'  summary(pool(fit))    
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
#'  \item{head}{dependent variable}
#'  \item{withdrawal_reason}{}
#'  }
#' @examples 
#' \dontrun{
#'  impCausalref1 <- mimix(data=acupuncture,covar=c("head_base","sex"),depvar=head,treatvar=treat,
#'        idvar=id,timevar=time,
#'        method=Causal,reference=1,K0=1,K1=0.5,M=5,seed=54321)
#'  library(mice)         
#'  fit<-with(as.mids(impCausalref1), lm(head.12~treat+head_base+sex))  
#'  summary(pool(fit))    
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