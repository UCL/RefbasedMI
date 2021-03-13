asthmaJ2R <- mimix(data=asthma, covar='base', depvar=fev, treatvar=treat, idvar=id, timevar=time,
                   method='J2R', reference=1, M=5, seed=101, prior=ridge, burnin=1000)
library(mice)     
fit <- with(data= as.mids(asthmaJ2R), lm(fev.12 ~ treat + base))
summary(pool(fit))