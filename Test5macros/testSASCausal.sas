************************ test  5macros  causal model ***************
*   testng 0.5**(12-max(timevar))as k1power
*    against the R program parameters k0=1 K1=0.5   
*    includes R code at bottom
*    results copied into spreadsheet Comparison of reference-based imputation ...  
********************************************************************
* Hamdvisit0.csv can be downloaded from https://github.com/ucl/mimix;  
Proc import out = work.hamdcsv DATAFILE="\\ad.ucl.ac.uk\home0\rmjlkm0\Documents\GitHub\Five_Macros20171010\Hamdvisit0.csv"
   DBMS= csv REPLACE; 
RUN;
Proc print data=work.hamdcsv(obs=10) noobs; run;

* From suupplemental material in Royes paper p46;
* vist 7 last visit;
*  need find max where means max non missing ;
data Hamdtest;
   set Hamdcsv;
   if ~missing(HAMA_TOTAL) then 
      visit= visit_number ;
run;

proc sql;
create table HAMD2 as
   select *, 0.5**(7-max(Visit)) as k1power /*k1power=0.5**(tmax-da)*/
   from Hamdtest
   group by Patient_Number;
quit;
* reorder visit_number in ascending order;
proc sort data = HAMD2 out=Hamd3;
 by patient_number visit_Number;
run;

* but we want HAMD17.TOTAL as response so need basval for visit0;
data HAMD4;
  set HAMD3;
   if VISIT_NUMBER=0 then HAMD17_TOTAL = basval;
run; 

* J2R ref DRUG  No  Covbytime= basval;
*%part1A(Jobname=noBasCaus,Data=hamdno3618csvC,Subject=Patient_number,Response=change,Time=Visit_number, Treat=treatment_name,covgroup=treatment_name,id=k1power);


%part1A(Jobname=noBasCaus,Data=hamd4,Subject=Patient_number,Response=HAMD17_TOTAL,Time=Visit_number, Treat=treatment_name,covgroup=treatment_name,Id=k1power);
%part1B(Jobname=noBasCaus,Ndraws=1000,thin=100,seed=12345);
* Build the means based on sampled values of posterior for the linear model parameters.;
%part2A(Jobname=noBasCaus_D,inname=noBasCaus,method=Causal,causalk=k1power,ref=DRUG,vcmethod=REF);
* Build the imputed values, sampling from the predicted distribution;
%part2b(Jobname=noBasCaus_D,seed=54321);
* Lsmeans differences are comparisons to 0 (not all possible combinations);
%part3(Jobname=noBasCaus_D,anref=DRUG,ancovgroup=1,Label=Causal with sigma from reference);


proc print data=work.noBasCaus_D_Out noobs;
run;
******* investigate individual patient ids;
data work.noBasCaus_D_kmfinalout;
   set  noBasCaus_D_datafinal;
     if (visit_number=0) then HAMD17_TOTAL0=HAMD17_TOTAL;
      if (visit_number=4) then HAMD17_TOTAL4=HAMD17_TOTAL;
       if (visit_number=5) then HAMD17_TOTAL5=HAMD17_TOTAL;
	    if (visit_number=6) then HAMD17_TOTAL6=HAMD17_TOTAL;
		 if (visit_number=7) then HAMD17_TOTAL7=HAMD17_TOTAL;
run;

* interim 3618;
%macro analyse( id);
data data&id;
   set  work.noBasCaus_D_kmfinalout;
  if Patient_number = &id;
   Title "Causal ref Drug   Patient_number =  &id  ";
Proc means data=work.data&id;
    var  HAMD17_TOTAL0 HAMD17_TOTAL4-HAMD17_TOTAL7;
	run;
%mend  analyse;
%analyse(3618);
%analyse(1513);
%analyse(1517);

******* then compare with test_mimix_ref.R program output

******* try asthma dataset

* to read asthma 
* \ad.ucl.ac.uk/home0/rmjlkm0/Documents/GitHub/mimix/mimixR";
Proc import out = work.asthmacsv DATAFILE="\\ad.ucl.ac.uk\home0\rmjlkm0\Documents\GitHub\mimix\mimix\inst\extdata\asthma.csv"
   DBMS= csv REPLACE; 
RUN;

* need to make basval response value at time 0;
data asthmatime0;
   set asthmacsv;
    
* to apply proc mianalyse merge on basval;
data asthmatime2;
   set asthmacsv;
    fev=base;
   	if time =2;
 run;
data asthmatime0;
   set asthmatime2;
     time =0;
 run;  

Proc append base= asthmacsv data= asthmatime0;
run;

Proc sort data =  asthmacsv;
   by id time;
   run;



* From suupplemental material in Royes paper p46;
* vist 12 last visit;
*  need find max where means max non missing ;
data asthmatest;
   set asthmacsv;
   if ~missing(fev) then 
      timevar = time ;
run;


proc sql;
create table asthmaK as
   select *, 0.5**(12-max(timevar))as k1power /*k1power=0.5**(tmax-da)*/
   from asthmatest
   group by id;
quit;
* reorder visit_number in ascending order;
proc sort data = asthmaK out=asthmaKS;
 by id time;
run;

%part1A(Jobname=noBasCausAsthma,Data=asthmaKS,Subject=id,Response=fev,Time=time, Treat=treat,covgroup=treat,Id=k1power);
%part1B(Jobname=noBasCausAsthma,Ndraws=1000,thin=100,seed=12345);
* Build the means based on sampled values of posterior for the linear model parameters.;
%part2A(Jobname=noBasCaus_DAsthma,inname=noBasCausAsthma,method=Causal,causalk=k1power,ref=2,vcmethod=REF);
* Build the imputed values, sampling from the predicted distribution;
%part2b(Jobname=noBasCaus_DAsthma,seed=54321);
* Lsmeans differences are comparisons to 0 (not all possible combinations);
%part3(Jobname=noBasCaus_DAsthma,anref=2,ancovgroup=1,Label=Causal with sigma from reference);

proc print data=work.noBasCaus_DAsthma_Out noobs;
run;
******* investigate individual patient ids;
data work.noBasCaus_DAsthma_kmfinalout;
   set  noBasCaus_DAsthma_datafinal;
      if (time=2) then fev2=fev;
       if (time=4) then fev4=fev;
	    if (time=8) then fev8=fev;
		 if (time=12) then fev12=fev;
run;

* interim 3618;
%macro analyseasthma( id);
data data&id;
   set  work.noBasCaus_DAsthma_kmfinalout;
  if id = &id;
   Title "Causal ref Drug   Patient_number =  &id  ";
Proc means data=work.data&id;
    var  fev2 fev4 fev8 fev12;
	run;
%mend  analyseasthma;
%analyseasthma(5017);
%analyseasthma(5051);
%analyseasthma(5115);
%analyseasthma(5059);
%analyseasthma(5333);
%analyseasthma(5456);
%analyseasthma(5454);


*** below R code 
# want to analyse with  no base but using base as resp at  time 0
testasthma<-asthma 
testasthma2<-(subset(testasthma,testasthma$time==2))
testasthma2$fev<- testasthma2$base
testasthma2$time <- 0
testasthma3<-(rbind(testasthma2,testasthma))
head(testasthma3[order(testasthma3[,"id"],testasthma3[,"time"]),])
# need to order  sts4Dpatt[order(sts4Dpatt[,treatvar],sts4Dpatt[,"patt"]),]
asthmadepVisit0<- testasthma3[order(testasthma3[,"id"],testasthma3[,"time"]),]
impasthmaCausal_NoDeltanocov <- mimix("asthmadepVisit0",NULL,"fev","treat","id","time",1000,1,"Causal",101,c("jeffreys"),1000,NULL,NULL,NULL,NULL,NULL,K0=1,K1=0.5,mle=0 )

#impdatasetJ2Rref1Causal54321<-mimix("asthma",c("base"),"fev","treat","id","time",10,1,"Causal",54321,"jeffreys",1000,NULL,NULL,NULL,NULL,NULL,K0=1,K1=0.5,mle=0 )

#cat('
fit<-with(data= as.mids(impasthmaCausal_NoDeltanocov), expr = lm(fev.12~treat))
summary(pool(fit))
varlist <- c("fev.2","fev.4","fev.8","fev.12")
analyselist(5051,impasthmaCausal_NoDeltanocov,varlist)
analyselist(5115,impasthmaCausal_NoDeltanocov,varlist)
analyselist(5333,impasthmaCausal_NoDeltanocov,varlist)
analyselist(5059,impasthmaCausal_NoDeltanocov,varlist)
analyselist(5017,impasthmaCausal_NoDeltanocov,varlist)
analyselist(5456,impasthmaCausal_NoDeltanocov,varlist)
analyselist(5454,impasthmaCausal_NoDeltanocov,varlist)

pttestf<- function(n1,n2,mn1,sd1,mn2,sd2) {
  pttest = stats::pt((((mn1 - mn2) -0) /sqrt(sd1^2/n1+sd2^2/n2)),(n1+n2-2))
  return(pttest)
}
