# analyse antidepressants data
# causal model with K initially 1 and multiplying by 0.5 each visit
# to compare between R and SAS (drop Catcov=PoolInv)
# R code - antidep.R
# SAS code - antidep.sas
# IW 19aug2020

library(mimix)
library(mice)
antidepressant<-read.table("http://raw.githubusercontent.com/UCL/mimix/master/inst/extdata/SASantidepmethodvar.csv",header=TRUE,sep=",",fileEncoding = "UTF-8-BOM")
table(antidepressant$VISIT.NUMBER,antidepressant$TREATMENT.NAME)

#recode covariates to numeric
antidepressant$PATIENT.SEX <- as.numeric(antidepressant$PATIENT.SEX)
antidepressant$TREATMENT.NAME <- as.numeric(antidepressant$TREATMENT.NAME)
# 1= drug, 2= placebo

table(antidepressant$VISIT.NUMBER,antidepressant$TREATMENT.NAME)

K0=1
K1=0.5
# ref = placebo
ref2 <- mimix("antidepressant",c("basval"),"change",
              "TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",
              100,2,"Causal",101,c("jeffreys"),1000,
              NULL,NULL,NULL,,,K0,K1)
fit2<-with(data= as.mids(ref2), expr = lm(change.7~TREATMENT.NAME+basval))
summary(pool(fit2))

# ref = drug
ref1 <- mimix("antidepressant",c("basval"),"change",
              "TREATMENT.NAME","PATIENT.NUMBER","VISIT.NUMBER",
              100,1,"Causal",101,c("jeffreys"),1000,
              NULL,NULL,NULL,,,K0,K1)
fit1<-with(data= as.mids(ref1), expr = lm(change.7~TREATMENT.NAME+basval))
summary(pool(fit1))

# > summary(pool(fit2))
            # term   estimate std.error statistic       df     p.value
# 1    (Intercept) -4.2881817 2.7754327 -1.545050 135.9273 0.124659251
# 2 TREATMENT.NAME  2.3260549 1.1391975  2.041836 139.5817 0.043051131
# 3         basval -0.2881528 0.1040768 -2.768654 137.6372 0.006405028

# > summary(pool(fit1))
            # term   estimate  std.error statistic       df      p.value
# 1    (Intercept) -3.1567852 2.59866343 -1.214773 149.3354 0.2263701001
# 2 TREATMENT.NAME  2.1051106 1.09305565  1.925895 143.9810 0.0560882527
# 3         basval -0.3654207 0.09807523 -3.725922 148.7030 0.0002758032
