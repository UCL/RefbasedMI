# mimix
<h2>Multiple imputation sensitivity analysis for longitudinal trials using R</h2> 

We aim to port the functionality of the Stata program **mimix**  into R 

The purpose of mimix is as described in the paper

Reference-based sensitivity analysis via multiple imputation for longitudinal trials with protocol deviation
by Suzie Cro, Tim P. Morris, Michael G. Kenward, and James R. Carpenter
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5796638/

The 5 methods available for sensitivity analysis are
 
|  Method         | option cmd             | reference group required |
| --------------- | --------------- | --------------------  |
| Randomized-arm                | MAR |  n |
| Jump to reference	            | J2R |  y |
| Copy increments in reference	| CIR |  y |
| Copy reference	              | CR  |  y |
| Last mean carried forward	    | LMCF|  n |
| Causal	                       | Causal|  y |

and an option for Delta adjustmnent to the imputed values

The R program does not provide the interim option available in Stata (where the individual has data observed later ) .

It does have the methodvar option (where  different imputation methods are specific for different individuals). 

# installation

within the R environment tyoe 

library(devtools) 

followed by 

install_github("UCL/mimix")


# usage

ensure all covariates are numeric, for example
antidepressant$PATIENT.SEX <- as.numeric(antidepressant$PATIENT.SEX)

otherwise you get the following warning 
Warning messages:
1: In emNorm.default(prnormobj, prior = priorvar[1], prior.df = priorvar[2]) :
  Factors in argument "y" converted to mode "numeric".
2: In emNorm.default(prnormobj, prior = priorvar[1], prior.df = priorvar[2]) :
  Factors in argument "y" converted to mode "numeric".
3

to run mimix()
ensure data set loaded, eg antidepressant

run mimix to save the imputed data-sets

Causal method
impdataCausal <- mimix("antidepressant",c("basval","POOLED.INVESTIGATOR","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,"Causal",101,c("jeffreys"),1000,NULL,NULL)

Individual specific, includind delta adjustment
impdataInd <- mimix("antidepressant",c("basval","POOLED.INVESTIGATOR","PATIENT.SEX"),"HAMD17.TOTAL","TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",100,1,NULL,101,c("jeffreys"),1000,NULL,c("methodvar","referencevar"),c(0.5,0.5,1,1 ))



