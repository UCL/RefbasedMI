<a href ="https://www.ctu.mrc.ac.uk/"><img src="MRCCTU_at_UCL_Logo.png" width="50%" /></a>

# RefBasedMI (mimix)
 *0.0.22*

# An R package for Reference-based multiple imputation for sensitivity analysis of longitudinal trials with protocol deviation

# Warning: this package is still under development

We have ported the functionality of the Stata program **mimix**  into R plus added some extra functionaity including options for the Causal model, plus Delta adjustment  . 

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

For an explanation of the **Causal model** see the paper by White,Royes and Best https://arxiv.org/abs/1705.04506 where the maintained treatment effect after treatment discontinuation is either constant (specified by K0) or decays exponentially by being multiplied by specifying a K1 value ( 0<= K1 <=1 ) in every period. Values of K0=1 , K1=0 are equivalent to method J2R, whilst values K0=1,K1=1 are equivalent to method CIR.       

https://www.tandfonline.com/doi/full/10.1080/10543406.2019.1684308


For an explanation of the **Delta adjustment** of imputed values see James Roger's SAS programs and user-guide under "Reference-based MI via Multivariate Normal RM (the "five macros" and MIWithD)" at  
https://www.lshtm.ac.uk/research/centres-projects-groups/missing-data#dia-working-group


  Specifying delta and dlag allows imputations to differ sytematically from RBI methods. 
  They provide an increment which is added on to all values imputed after 
  treatment discontinuation, but not to interim (intermediate) missing values. 
  Values of delta are cumulated after treatment discontinuation.
  For example, for an individual who discontinued treatment at the 2nd time point, we take 
  the vector of delta's starting at the 3rd time point and add their cumulative sums to the imputed values. 
  Specifying dlag modifies this behaviour, so that the vector of delta's starting at the 3rd time point is 
  multipled elementwise by the vector dlag.
  The formula for the increment at time k for an individual who discontinued after time p is
   b_1xa_{p+1} + b_2xa_{p+2} + ... + b_{k-p}xa_k
   where delta=(a_1,a_2,...) and 
         dlag=(b_1,b_2,...). 
  A common increment of 3 at all time points after treatment discontinuation is achieved 
  by setting  delta=c(3,3,3,...) and dlag=c(1,0,0,...), both vectors having the length of the number of time points.


For details of the norm2 package which supplies the function mcmcNorm - the MCMC algorithm for incomplete multivariate normal data
https://rdrr.io/cran/norm2/src/R/norm2.R

# installation

within the R environment type 

    if(!require(devtools)) install.packages('devtools') 

    library(devtools) 

followed by 

    devtools::install_github("UCL/RefBasedMI")


# usage

**mimix**(data,covar,depvar,treatvar,idvar,timevar,..... options..... 


Arguments in function mimix() 

**data**	       dataset in wide (longitudinal data) format

**covar**       baseline covariates and/or baseline depvar, must be complete (no missing vaules) and treated as numeric

**depvar**	     dependent variable

**treatvar**   	treatment group , recoded to 1,2,..

**idvar**	      patient id

**timevar**	    time point for repeated measure

**M**	          number of imputations
 
**reference**	      reference group for j2r,cir,cr,Causal methods

**method**	       RBI method

**seed**	    seed value to obtain same outputs

**prior**     prior to use in mcmcNorm, default jeffreys, uniform  or ridge

**burnin**	     burnin value

**bbetween**	   value between iterations in mcmc

**methodvar**   variable in data specifying individual method 

**referencevar**   variable in data specifying individual reference group

**delta**       vector of delta values to add onto imputed values (a values in Roger's paper) (non-mandatory)

**dlag**        vector of dlag values (b values in Rogers paper)

**K0**	         Causal constant for use with Causal method

**K1**	         exponential decaying Causal constant for use with Causal method



# examples

## Sample data: asthma trial

### J2R analysis with control as reference
asthmaJ2R <- mimix(data = asthma, covar = base, depvar = fev, treatvar = treat,	idvar = id, timevar = time, method = "J2R", reference = 2, M = 5, seed = 101, 
	prior = "ridge", burnin = 1000)
 
### Analysis

fit specified model to each imputed data set (assigned as mids class) and pool results together (Rubin's rules),
functions from mice package

library(mice)

fit<-with(as.mids(subset(asthmaJ2R,time==12)), lm(fev~treat+base))

summary(pool(fit))


### Delta method - all values are 1 unit lower than expected under J2R
impJ2Rridge <- mimix(data = asthma, covar = c(base), depvar = fev, treatvar = treat,	idvar = id, timevar = time, method = "J2R", reference = 2, 
	delta = c(-1, 0, 0, 0), M = 5, seed = 101, prior = "ridge")      

fit<-with(as.mids(subset(impJ2Rridge,time==12)), lm(fev~treat+base))

summary(pool(fit))



## Sample data: antidepressant trial

### Mixed methods

methodcol and referencecol are variables in the data set 

NOTE - either method or methodvar to be specified but NOT both

antidepIndiv <- mimix(data = antidepressant, covar = c(basval, PATIENT.SEX),depvar = HAMD17.TOTAL, treatvar = TREATMENT.NAME, idvar = PATIENT.NUMBER, 
    timevar = VISIT.NUMBER, methodvar = methodcol, referencevar = referencecol, 
	  M = 5, seed = 54321)        

### Analysis

antidepIndiv <- with(data =  as.mids(subset(antidepIndiv, VISIT.NUMBER == 7)),	lm(HAMD17.TOTAL ~ TREATMENT.NAME + basval + PATIENT.SEX))
summary(pool(antidepIndiv))  

## Sample data: acupuncture trial 

### Causal model: treatment effect halves every 1 time unit 
### after treatment discontinuation
Note K0=1,K1=0 equivalent to J2R,  K0=1,K1=1 equivalent to CIR 

acuCausal <- mimix(data = acupuncture, covar = c(head_base), depvar = head,	treatvar = treat, idvar = id, timevar = time, method = "Causal", 
	reference = 1, K0 = 1, K1 = 0.5, M = 5, seed = 54321)
 
### Analysis

acufit <- with(as.mids(subset(impCausalref, time == 12)),	lm(head ~ treat + head_base ))
summary(pool(acufit))    








