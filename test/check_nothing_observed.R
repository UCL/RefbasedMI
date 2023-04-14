# check_nothing_observed.R
# IW 4/1/2022, updated 6/4/2023
#
# explore how individuals with no outcomes are handled
# do this by adding a large number to all in treat=2, 
#   making one individual in treat=2 wholly missing, 
#   and viewing their imputed values with ref=1
#
# Finding: only MAR and LMCF impute 10000 larger, consistent with only these using the own-arm mean


# select the individual
idc = min(asthma[asthma$treat==2,"id"])


# modify the data
asthma2 <- asthma
asthma2[asthma2$treat==2,"fev"] <- asthma2[asthma2$treat==2,"fev"] + 10000
asthma2[asthma2$treat==2 & asthma2$id==idc,"fev"] <- NA

asthma2 %>% filter(id==idc)

asthma2 %>% filter(!is.na(fev)) %>% 
  group_by(treat, time) %>% 
  summarise(n=n(), fevmean=mean(fev), fevsd=sd(fev))

# impute various ways
impnooutcomes2mar <- RefBasedMI(data=asthma2,
                                depvar=fev,
                                treatvar=treat,
                                idvar=id,
                                timevar=time,
                                M=2,
                                method="mar",
                                seed=101,
                                prior="jeffreys",
                                burnin=1000,
                                bbetween=NULL,
                                methodvar=NULL
)

impnooutcomes2j2r <- RefBasedMI(data=asthma2,
                                depvar=fev,
                                treatvar=treat,
                                idvar=id,
                                timevar=time,
                                M=2,
                                method="j2r",
                                reference=1,
                                seed=101,
                                prior="jeffreys",
                                burnin=1000,
                                bbetween=NULL,
                                methodvar=NULL,
                                delta=c(1,2,3,4), dlag=c(.1,.1,.1,.05)
)

impnooutcomes2cir <- RefBasedMI(data=asthma2,
                                depvar=fev,
                                treatvar=treat,
                                idvar=id,
                                timevar=time,
                                M=2,
                                method="cir",
                                reference=1,
                                seed=101,
                                prior="jeffreys",
                                burnin=1000,
                                bbetween=NULL,
                                methodvar=NULL,
                                delta=c(1,2,3,4), dlag=c(.1,.1,.1,.05)
)

impnooutcomes2lmcf <- RefBasedMI(data=asthma2,
                             depvar=fev,
                             treatvar=treat,
                             idvar=id,
                             timevar=time,
                             M=2,
                             method="lmcf",
                             reference=1,
                             seed=101,
                             prior="jeffreys",
                             burnin=1000,
                             bbetween=NULL,
                             methodvar=NULL,
                             delta=c(1,2,3,4), dlag=c(.1,.1,.1,.05)
)

impnooutcomes2causal <- RefBasedMI(data=asthma2,
                             depvar=fev,
                             treatvar=treat,
                             idvar=id,
                             timevar=time,
                             M=2,
                             method="causal", K0=0.8, K1=1,
                             reference=1,
                             seed=101,
                             prior="jeffreys",
                             burnin=1000,
                             bbetween=NULL,
                             methodvar=NULL,
                             delta=c(1,2,3,4), dlag=c(.1,.1,.1,.05)
)

# view imputed values
impnooutcomes2mar %>% filter(id==idc) 
impnooutcomes2j2r %>% filter(id==idc) 
impnooutcomes2cir %>% filter(id==idc) 
impnooutcomes2lmcf %>% filter(id==idc) 
impnooutcomes2causal %>% filter(id==idc) 

