<a href ="https://www.ctu.mrc.ac.uk/"><img src="MRCCTU_at_UCL_Logo.png" width="50%" /></a>

# mimix
 *0.0.6*

# An R package for Reference-based multiple imputation for sensitivity analysis of longitudinal trials with protocol deviation

We have ported the functionality of the Stata program **mimix**  into R. 

The purpose of mimix is as described in the paper

> Reference-based sensitivity analysis via multiple imputation for longitudinal trials with protocol deviation
by Suzie Cro, Tim P. Morris, Michael G. Kenward, and James R. Carpenter
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5796638/

and (within Stata) type  "ssc install mimix" to install code  and "help mimix" to read help documentation
 

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


For an explanation of the Delta adjustment of imputed values see James Roger's SAS programs and user-guide at 
https://missingdata.lshtm.ac.uk/files/2017/04/Five_Macros20171010.zip 

For details of the norm2 package which supplies the function mcmcNorm - the MCMC algorithm for incomplete multivariate normal data
https://rdrr.io/cran/norm2/src/R/norm2.R

# installation

within the R environment type 

    if(!require(devtools)) install.packages('devtools') 

    library(devtools) 

followed by 

    install_github("UCL/mimix")


# usage

mimix(data,covar,depvar,treatvar,idvar,timevar,..... options..... 


Arguments in function mimix() 

**data**	       dataset in wide (longitudinal data) format

**covar**       baseline covariates and/or baseline depvar, must be complete (no missing vaules) and treated as numeric

**depvar**	     dependent variable

**treatvar**   	treatment group , recoded to 1,2,..

**idvar**	      patient id

**timevar**	    time point for repeated measure

**M**	          number of imputations
 
**refer**	      reference group for j2r,cir,cr methods

**method**	       RBI method

**seed**	    seed value to obtain same outputs

**priorr**    prior tu use in mcmcNorm, default jeffreys, uniform  or ridge

**burnin**	     burnin value

**bbetween**	   value between iterations in mcmc

**methodvar**  2 element vector designating variables in data specifying individual method and reference group

**delta**       vector of delta values to add onto imputed values (a values in Roger's paper) (non-mandatory)

**dlag          vector of dlag values (b values in Roger's oaper)

**K0**	         Causal constant for use with Causal method

**K1**	         exponential decaying Causal constant for use with Causal method



# examples


### Jump to reference (J2R) with  placebo treatment group as reference, delta adjustment 

impdatasetJ2R<-mimix("asthma",c("base"),"fev","treat","id","time",2,1,"J2R",101,"jeffreys",1000,NULL,NULL,c(0.5,0.5,1,1),(1,1,1,1) )   

run regression on imputed data-sets, combining using Rubin's rules

set file as a mids class, specifying id and imputation number columns  

fit specified model to each imputed data set (assigned as mids class) and pool results together (Rubin's rules),
functions from mice package

library(mice)
fit<-with(as.mids(impdatasetJ2R), lm(fev.12~treat+base))
summary(pool(fit))



### Using the Causal method

Example
impdata <- mimix("antidepressant",c("basval","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",
                                                                                 100,1,"Causal",101,c("jeffreys"),1000,NULL,NULL,,1,1)


### Individual specific method, with delta adjustment
NOTE - either meth and methodIndiv to be specified but NOT both

impdataInd <- mimix("antidepressant",c("basval","POOLED.INVESTIGATOR","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,NULL,101,c("ridge"),1000,NULL,c("methodvar","referencevar"),c(0.5,0.5,1,1 ),c(1,1,2,2))



