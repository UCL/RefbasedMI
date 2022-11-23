# mimixhelpfile.R
# mimix help file examples
# IW 23/11/2022

setwd("C:/ado/ian/RefBasedMI")
load("data/asthma.RData")

# Start with CC to check correct dataset loading
summary(lm(data=asthma, formula=fev~as.factor(treat)+base, subset=(time==12)))
M=500

#    Multiple imputation assuming the response variable fev is MAR
#use "data\asthma.dta", clear
#mimix fev treat, id(id) time(time) method(mar) covariates(base) clear m(500) seed(101)
#mi estimate, mcerror: reg fev i.treat base if time==12
impMAR <- RefBasedMI(data=asthma, covar=base, depvar=fev, treatvar=treat, 
                      idvar=id, timevar=time, M=M, method="MAR", seed=101, 
                      prior="jeffreys", burnin=1000, bbetween=NULL, methodvar=NULL)
fit <- with(as.mids(impMAR), lm(fev~as.factor(treat)+base), subset=(time==12))
impMARres <- summary(pool(fit))
impMARres

#    Multiple imputation and regression analysis assuming last mean carried forward for the response variable fev
#use "data\asthma.dta", clear
#mimix fev treat, id(id) time(time) method(lmcf) covariates(base) clear m(500) regress seed(101)
#mi estimate, mcerror: reg fev i.treat base if time==12
impLMCF <- RefBasedMI(data=asthma, covar=base, depvar=fev, treatvar=treat, 
                      idvar=id, timevar=time, M=M, method="LMCF", reference=1, seed=101, 
                      prior="jeffreys", burnin=1000, bbetween=NULL,methodvar=NULL)
fit <- with(as.mids(impLMCF), lm(fev~as.factor(treat)+base), subset=(time==12))
impLMCFres <- summary(pool(fit))
impLMCFres

#    Multiple imputation and regression analysis assuming jump to reference for the response variable fev, with placebo=2 as the reference
#use "data\asthma.dta", clear
#mimix fev treat, id(id) time(time) method(j2r) refgroup(2) covariates(base) clear m(500) regress seed(101)
#mi estimate, mcerror: reg fev i.treat base if time==12
impJ2R1 <- RefBasedMI(data=asthma, covar=base, depvar=fev, treatvar=treat, 
                      idvar=id, timevar=time, M=M, method="J2R", reference=1, seed=101, 
                      prior="jeffreys", burnin=1000, bbetween=NULL, methodvar=NULL)
fit <- with(as.mids(impJ2R1), lm(fev~as.factor(treat)+base), subset=(time==12))
impJ2R1res <- summary(pool(fit))
impJ2R1res

#    Saving the imputed dataset with filename mimix_example, assuming copy increments in reference for the response variable fev, with placebo=2 as the reference
#use "data\asthma.dta", clear
#mimix fev treat, id(id) time(time) method(cir) refgroup(2) covariates(base) m(500) seed(101)
#mi estimate, mcerror: reg fev i.treat base if time==12
impCIR1 <- RefBasedMI(data=asthma, covar=base, depvar=fev, treatvar=treat, 
                      idvar=id, timevar=time, M=M, method="CIR", reference=1, seed=101, 
                      prior="jeffreys", burnin=1000, bbetween=NULL, methodvar=NULL)
fit <- with(as.mids(impCIR1), lm(fev~as.factor(treat)+base), subset=(time==12))
impCIR1res <- summary(pool(fit))
impCIR1res
