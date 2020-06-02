<a href ="https://www.ctu.mrc.ac.uk/"><img src="MRCCTU_at_UCL_Logo.png" width="50%" /></a>

# mimix
 *0.0.4*

# an R package for Reference-based multiple imputation for sensitivity analysis of longitudinal trials with protocol deviation

We have ported the functionality of the Stata program **mimix**  into R. 

The purpose of mimix is as described in the paper

> Reference-based sensitivity analysis via multiple imputation for longitudinal trials with protocol deviation
by Suzie Cro, Tim P. Morris, Michael G. Kenward, and James R. Carpenter
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5796638/

and (within Stata) type  "ssc install mimix" to install code  and "help mimix" to read help documetation
 

The 5 methods (plus Causal) available for sensitivity analysis are
 
|  Method         | option cmd             | reference group required |
| --------------- | --------------- | --------------------  |
| Randomized-arm                | MAR |  n |
| Jump to reference	            | J2R |  y |
| Copy increments in reference	| CIR |  y |
| Copy reference	              | CR  |  y |
| Last mean carried forward	    | LMCF|  n |
| Causal	                       | Causal|  y |

For explanation of the Causal model refer to 

> https://www.tandfonline.com/doi/full/10.1080/10543406.2019.1684308


For an explanation of the Delta adjustment of imputed values see James Rogers SAS programs and user-guide at 
https://missingdata.lshtm.ac.uk/files/2017/04/Five_Macros20171010.zip 



# installation

within the R environment type 

    library(devtools) 

followed by 

    install_github("UCL/mimix")



# usage

covariates and basevalue response must be non-missing, also all covariates to be converted to numeric, for example 
antidepressant$PATIENT.SEX <- as.numeric(antidepressant$PATIENT.SEX)

otherwise you get the following warning 
*Warning messages:
1: In emNorm.default(prnormobj, prior = priorvar[1], prior.df = priorvar[2]) :
  Factors in argument "y" converted to mode "numeric".
2: In emNorm.default(prnormobj, prior = priorvar[1], prior.df = priorvar[2]) :
  Factors in argument "y" converted to mode "numeric".
3*


arguments in function mimix() showing default values
NOTE - either meth and methodIndiv to be specified but NOT both

mimix(data,covar,depvar,treatvar,idvar,timevar,M=1,refer,meth,seedval=101,priorvar,burnin=1000,bbetween=NULL,methodindiv,delta=NULL) 

# examples


### Jump to reference (J2R) with  placebo treatment group as reference, delta adjustment 

impdatasetJ2R<-mimix("asthma",c("base"),"fev","treat","id","time",2,1,"J2R",101,"jeffreys",1000,NULL,NULL,c(0.5,0.5,1,1 ) )   

check outputs for individual patient id

varlist <- c("fev.2","fev.4","fev.8","fev.12","base")

analyselist(5017,impdatasetJ2R,varlist)

run regression on imputed data-sets, combining using Rubin's rules

set file as a mids class, specifying id and imputation number columns  
impdata= as.mids(impdatasetJ2R, .id="SNO",.imp="II")

fit specified model to each imputed data set and pool results together (Rubin's rules), functions from mice package
fit<-with(impdata, lm(fev.12~treat+base))

summary(pool(fit))



### Using the Causal method

impdataCausal <- mimix("antidepressant",c("basval","POOLED.INVESTIGATOR","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"Causal",101,c("jeffreys"),1000,NULL,NULL)


### Individual specific method, with delta adjustment

impdataInd <- mimix("antidepressant",c("basval","POOLED.INVESTIGATOR","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,NULL,101,c("jeffreys"),1000,NULL,c("methodvar","referencevar"),c(0.5,0.5,1,1 ))


