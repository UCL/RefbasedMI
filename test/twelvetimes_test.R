# Test RefBasedMI with artificial data with 12 times, 4 groups and no covariate
# IW 22nov2021, updated 6apr2023
# twelvetimes20.csv  has  20 individuals / group
# twelvetimes200.csv has 200 individuals / group

# Problem: with 20 / group, norm2 reports "MCMC procedure aborted" in each imputation, 
#   yet RefBasedMI reports a final "mcmcNorm Loop finished" and it soldiers on
#   Can we make it abort more nicely?

# With 200 / group, norm2 succeeds

library(tidyverse)

# n = 200 / group

t12n200 <- read_csv("twelvetimes200.csv")
summary(t12n200)
table(t12n200$group,t12n200$time)

impJ2R1 <- RefBasedMI(data=t12n200,
                      depvar=yvar,
                      treatvar=group,
                      idvar=myid,
                      timevar=time,
                      M=2,
                      reference=1,
                      method="J2R",
                      seed=101,
                      prior="jeffreys",
                      burnin=1000,
                      bbetween=NULL,
                      methodvar=NULL
)


# n = 20 / group

t12n20 <- read_csv("twelvetimes20.csv")
summary(t12n20)
table(t12n20$group,t12n20$time)

impJ2R1 <- RefBasedMI(data=t12n20,
                      depvar=yvar,
                      treatvar=group,
                      idvar=myid,
                      timevar=time,
                      M=2,
                      reference=1,
                      method="J2R",
                      seed=101,
                      prior="jeffreys",
                      burnin=1000,
                      bbetween=NULL,
                      methodvar=NULL
)


