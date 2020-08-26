/*Causal model estimation using the SAS 5-macros*/
/*An example data set (Chapter15_example.sas7bdat) is available on the DIA page of www.missingdata.org.uk*/


*%let Path1 = H:\rj\datasets\HAMD17;
*%let Path1 = H:\royes\datasets\HAMD17 data;
%let Path1 = H:\missing\Royes Joseph\data;
libname library "&Path1";

*Load the 5-macros;
*%let Path = H:\royes\RBI methods\Five_Macros20161111\TheFiveMacros;
%let Path = H:\missing\Royes Joseph\macros\Five_Macros20161111\TheFiveMacros;
%let Path2 = H:\missing\Royes Joseph\macros\Five_Macros20161111\Causal;

%PUT &Path;

%include "&Path\part1a_32.sas";
%include "&Path\part1b_46.sas";
%include "&Path2\part2A_38_causal.sas"; /*Here is the modified part 2A macro*/
%include "&Path\part2b_30.sas";
%include "&Path\part3_54.sas";

\\ad.ucl.ac.uk\home0\rmjlkm0\Documents\Five_Macros20161111\Five_Macros20161111
* km path%let Path= N:\Documents\GitHub\Five_Macros20171010\Five_Macros20171010\TheFiveMacros\part3_54.sas.
*%let Path= N:\Documents\GitHub\Five_Macros20171010\Five_Macros20171010\TheFiveMacros\part3_54.sas.
%let Pathkm= N:\Documents\Five_Macros20161111\Five_Macros20161111\TheFiveMacros;
%let Path2km= N:\Documents\Five_Macros20161111\Five_Macros20161111\Causal;

* from most uptodate macros;
%let Pathkm=N:\Documents\GitHub\Five_Macros20171010\Five_Macros20171010\TheFiveMacros;
%include "&Pathkm\part1a_33.sas";
%include "&Pathkm\part1b_47.sas";
%include "&Pathkm\part2a_40.sas";
%include "&Pathkm\part2b_31.sas";
%include "&Pathkm\part3_55.sas";
%include "&Pathkm\plotter_7.sas";


libname library "&Pathkm";
%PUT &Pathkm;
%PUT &PathkmCausal;

%include "&Pathkm\part1a_32.sas";
%include "&Pathkm\part1b_46.sas";
%include "&Path2km\part2A_38_causal.sas"; /*Here is the modified part 2A macro*/
%include "&Pathkm\part2b_30.sas";
%include "&Pathkm\part3_54.sas";



options mprint nomlogic;

* assumption K0=1, K1=0.5;

* put in macro to cycle thru values K1;

proc sql;
create table HAMD as
   select *, 0.8**(7-max(visit))as k1power /*k1power=0.5**(tmax-da)*/
   from library.Chapter15_example
   group by patient;
quit;




* try ;
* this from causal_eample run.sas;
data hamd;
length patient change visit basval 8;
set library.Chapter15_example;
k1power=0.5;/*here k is specified by a variable name*/
run;

proc sql;
create table HAMD as
   select *, 0.5**(7-max(visit))as k1power /*k1power=0.5**(tmax-da)*/
   from library.Chapter15_example
   group by patient;
quit;


%macro cycleK1_NocatcovD(K1);
proc sql;
create table HAMD as
   select *, &K1**(7-max(visit))as k1power /*k1power=0.5**(tmax-da)*/
   from library.Chapter15_example
   group by patient;
quit;
*Call macros;
* HAMDLT17 not closer than change; 
*;
*HAMATOTL not closer at all;

* remove Catcov  not much differne covgroup=Therapy;
* try wout Covbytime;
%part1A(Jobname=example,Data=hamd,Subject=Patient,Response=change,Time=Visit, Treat=Therapy,covgroup=Therapy,id=k1power);
*part1A(Jobname=example,Data=hamd,Subject=Patient,Response=change,Time=Visit, Treat=Therapy,Covbytime= Basval,covgroup=Therapy,id=k1power);

*%part1A(Jobname=example,Data=hamd,Subject=Patient,Response=change,Time=Visit, Treat=Therapy,Catcov=PoolInv,Covbytime= Basval,covgroup=Therapy,id=k1power);

%part1B(Jobname=example,Ndraws=1000,thin=100,seed=12345);

* Build the means based on sampled values of posterior for the linear model parameters.;
* try J2R;
*%part2A(Jobname=example_causal,inname=example,method=J2R ,causalk=k1power,ref=DRUG,vcmethod=REF);
*PLACEO/DRUG;
%part2A(Jobname=example_causal,inname=example,method=j2r,causalk=k1power,ref=PLACEBO,vcmethod=REF);
*%part2A(Jobname=example_causal,inname=example,method=CAUSAL,causalk=k1power,ref=DRUG,vcmethod=REF);
* Build the imputed values, sampling from the predicted distribution;
%part2b(Jobname=example_causal,seed=54321);
* Lsmeans differences are comparisons to 0 (not all possible combinations);
%part3(Jobname=example_causal,anref=PLACEBO,Label=j2r with sigma from reference);


************************* Provide summary table  ***********************;
title1 "Summary Table_ Causal results with Sigma from REF arm no catcov  0.5**(7-max(visit))as k1power";
*title1 "Summary Table_ Causal results with Sigma from REF arm  1 as k1power";

title1 "Summary table J2R Ndraws=1000,thin=100, ref arm PLACEBO No covbytime";
proc print data=work.example_causal_Out noobs;
run;


*extracting data from work.example_causal_Out;
*extracting data from work.example_causal_Out;

%let sufx=%sysevalf(&K1*10);
data testnameDIAnocatcov&sufx; 
 set work.example_causal_Out;
   K1=&K1;
   keep K1 diff LCL_Diff UCL_Diff ;
   if (row=8 ) then output; 
run;
%mend  cycleK1_NocatcovD;

%cycleK1_NocatcovD(0);
%cycleK1_NocatcovD(0.1);
%cycleK1_NocatcovD(0.2);
%cycleK1_NocatcovD(0.3);
%cycleK1_NocatcovD(0.4);
%cycleK1_NocatcovD(0.5);
%cycleK1_NocatcovD(0.6);
%cycleK1_NocatcovD(0.7);
%cycleK1_NocatcovD(0.8);
%cycleK1_NocatcovD(0.9);
%cycleK1_NocatcovD(1);

data testnameDIAnocatcovallDrug; 
   set 
   testnameDIAnocatcov0 
   testnameDIAnocatcov1 
   testnameDIAnocatcov2 
   testnameDIAnocatcov3 
   testnameDIAnocatcov4 
   testnameDIAnocatcov5 
   testnameDIAnocatcov6 
   testnameDIAnocatcov7 
   testnameDIAnocatcov8 
   testnameDIAnocatcov9 
   testnameDIAnocatcov10;
run ;



data test;
 set work.example_causal_Out;
  keep diff LCL_Diff UCL_Diff ;
 if (row=8 ) then output; 
run;
proc print data =test; run;

%macro Mixeddraw(draw,visit);
  
title "proc mixed example_causal_datafull(where=(draw=&draw and visit=&visit)";
Proc Mixed data= work.example_causal_datafinal (where=(draw=&draw and visit=&visit))  Method=REML;
  Class patient  Therapy visit Poolinv Therapy  ;
 Model change =  Therapy Poolinv basval /  DDFM=KenwardRoger;
 *Random intercept/ subject=patient Group=therapy;
 Repeated intercept / subject=patient Group=therapy ;
 LSMeans therapy / Diff;
run; 
%mend;
%Mixeddraw(1,5);

* test for J2R
*Call macros;
%part1A(Jobname=example,Data=hamd,Subject=Patient,Response=change,Time=Visit, Treat=Therapy,Catcov=PoolInv,Covbytime= Basval,covgroup=Therapy,id=k1power);

* try w/out covbytime to compare with stata same result ;
%part1A(Jobname=example,Data=hamd,Subject=Patient,Response=change,Time=Visit, Treat=Therapy,Catcov=PoolInv,Cov= Basval,covgroup=Therapy,id=k1power);
%part1B(Jobname=example,Ndraws=100,thin=100,seed=12345);

* Build the means based on sampled values of posterior for the linear model parameters.;
%part2A(Jobname=example_J2Rdrug,inname=example,method=J2R ,ref=DRUG,vcmethod=REF);
%part2b(Jobname=example_J2Rdrug,seed=54321);

%part2A(Jobname=example_J2Rplacebo,inname=example,method=J2R ,ref=PLACEBO,vcmethod=REF);
* Build the imputed values, sampling from the predicted distribution;
%part2b(Jobname=example_J2Rplacebo,seed=54321)
* Lsmeans differences are comparisons to 0 (not all possible combinations);
%part3(Jobname=example_J2Rplacebo,anref=PLACEBO,Label=J2R with sigma from reference);

title1 "Summary Table_ J2R results with Sigma from REF arm";
proc print data=work.example_J2Rplacebo_Out noobs;
run;
data work.example_J2R_kmfinalout;
   set  example_J2Rplacebo_datafull;
      if (visit=4) then change4=change;
       if (visit=5) then change5=change;
	    if (visit=6) then change6=change;
		 if (visit=7) then change7=change;
run;

data work.example_J2R_kmfinalout;
   set  example_J2Rdrug_datafull;
      if (visit=4) then change4=change;
       if (visit=5) then change5=change;
	    if (visit=6) then change6=change;
		 if (visit=7) then change7=change;
run;
%macro analyse( id);
   data  work.example_J2R_kmfinalout&id;
     set work.example_J2R_kmfinalout;
   if Patient = &id;
 Title "J2R Drug &id  ";
Proc means data=work.example_J2R_kmfinalout&id;
   var  change4-change7;
  run;
%mend  analyse;



%analyse(3618);
%analyse(3712);

Proc Mixed data= work.example_J2Rplacebo_datafull(where=(draw=7 and visit=7))  Method=REML;
  Class  Therapy visit Poolinv patient ;
 Model change =Therapy Poolinv basval /DDFM=KenwardRoger;
 *Random / subject=patient;
 Repeated / subject=patient Group=therapy;
 LSMeans therapy / Diff;
run; 

*magic sequence !!;
*); */; /*’*/ /*”*/; %mend;

* from Dia__run3;



* model 5;

%macro cycleDiarunNObasval(K1);
data hamd;
length patient change visit basval 8;
set library.Chapter15_example;
k1=&K1;/*here k is specified by a variable name*/
run;



*Call macros;
%part1A(Jobname=example
,Data=hamd
,Subject=Patient
,Response=change
,Time=Visit
,Treat=Therapy
,Catcov=PoolInv
,Covbytime= Basval
,covgroup=Therapy
,id=k1
);

%part1B(Jobname=example
,Ndraws=10
,thin=100
,seed=123456
);
*Part 1A and 1B need not to be repeated for different values of k;


* Build the means based on sampled values of posterior for the linear model parameters.;
%part2A(Jobname=example_causal
,inname=example
,method=CAUSAL /*if specify the method as CAUSAL, k should also be specified*/
,causalk=k1 /*By default CausalK takes value 1 and result CIR estimate*/
,ref=PLACEBO
,vcmethod=REF
);

* Build the imputed values, sampling from the predicted distribution;
%part2b(Jobname=example_causal
,seed=654321
);

* Run the analysis section, adding a Label to the results, for when we replay them all together at the end.;
* Lsmeans differences are comparisons to 0 (not all possible combinations);
%part3(Jobname=example_causal
,anref=PLACEBO
,Label=Causal with sigma from reference & k=&K1
);

options linesize=200;
proc print data=work.example_causal_Out noobs; run;

*extracting data from work.example_causal_Out;
* for negative suffixes;
%let sufx=%sysevalf(&K1*-10);
data testnameDIArun3&sufx; 
 set work.example_causal_Out;
   K1=&K1;
   keep K1 diff LCL_Diff UCL_Diff ;
   if (row=8 ) then output; 
run;
%mend cycleDiarun3;

%cycleDiarun3(-0.5);
%cycleDiarun3(-0.4);
%cycleDiarun3(-0.3);
%cycleDiarun3(-0.2);
%cycleDiarun3(-0.1);
* doesnt work for negative !;


%cycleDiarun3(0);
%cycleDiarun3(0.1);
%cycleDiarun3(0.2);
%cycleDiarun3(0.3);
%cycleDiarun3(0.4);
%cycleDiarun3(0.5);
%cycleDiarun3(0.6);
%cycleDiarun3(0.7);
%cycleDiarun3(0.8);
%cycleDiarun3(0.9);
%cycleDiarun3(1);
%cycleDiarun3(1.1);
%cycleDiarun3(1.2);
%cycleDiarun3(1.3);
%cycleDiarun3(1.4);
%cycleDiarun3(1.5);
%cycleDiarun3(1.6);
%cycleDiarun3(1.7);
%cycleDiarun3(1.8);
%cycleDiarun3(1.9);
%cycleDiarun3(2);
%cycleDiarun3(2.1);
%cycleDiarun3(2.2);
%cycleDiarun3(2.3);
%cycleDiarun3(2.4);
%cycleDiarun3(2.5);

data testDIArun25all;
  set 
  testnameDIArun30
  testnameDIArun31
  testnameDIArun32
  testnameDIArun33
  testnameDIArun34
  testnameDIArun35
  testnameDIArun36
  testnameDIArun37
  testnameDIArun38
  testnameDIArun39
  testnameDIArun310
  testnameDIArun312
  testnameDIArun313
  testnameDIArun314
  testnameDIArun315
  testnameDIArun316
  testnameDIArun317
  testnameDIArun318
  testnameDIArun319
  testnameDIArun320
  testnameDIArun321
  testnameDIArun322
  testnameDIArun323
  testnameDIArun324
  testnameDIArun325
  ;
run;

* then add the negative ones ;
data testDIArun25allneg;
 set
  testnameDIArun35
  testnameDIArun34
  testnameDIArun33
  testnameDIArun32
  testnameDIArun31
  testDIArun25all;
run;

data testDIArun3all;
  set 
  testnameDIArun30
  testnameDIArun31
  testnameDIArun32
  testnameDIArun33
  testnameDIArun34
  testnameDIArun35
  testnameDIArun36
  testnameDIArun37
  testnameDIArun38
  testnameDIArun39
  testnameDIArun310;
run;

* now plot in r ;

%macro kloop;

%cycleK1(0);
%cycleK1(0.1);
%cycleK1(0.2);
%cycleK1(0.3);
%cycleK1(0.4);
%cycleK1(0.5);
%cycleK1(0.6);
%cycleK1(0.7);
%cycleK1(0.8);
%cycleK1(0.9);
%cycleK1(1);
************************* Provide summary table  ***********************;
* Rank makes the order of records follow the order specified on the set statement, if you want;
data work.causal_ref;
set  example_causal_Out;
rank=_n_;
run;
title1 "Summary Table_ Causal results with Sigma from REF arm";
proc print data=work.causal_ref(drop=row rank) noobs;
run;



* goll standard w/out poolinv;
*Call macros;


proc sql;
create table HAMD as
   select *, 0.5**(7-max(visit))as k1power /*k1power=0.5**(tmax-da)*/
   from library.Chapter15_example
   group by patient;
quit;

*18/08;
%macro cycleDiarunNObasval(K1);
proc sql;
create table HAMD as
   select *, &K1**(7-max(visit))as k1power /*k1power=0.5**(tmax-da)*/
   from library.Chapter15_example
   group by patient;
  quit;
data Hamdx;
  set hamd;
   if gender= "F" then sex=1;
  if gender ="M" then sex=2;
 run;
*Call macros
* take Catcov Poolinv out to establish gold standard causal to compare stata/r;
* just above didnt produce same as R so try taking out cpvbytime=basval as well ;

%part1A(Jobname=example,Data=HAMDx,Subject=Patient,Response=change,Time=Visit, cov=sex,Treat=Therapy,covgroup=Therapy,id=k1power);
%part1B(Jobname=example,Ndraws=1000,thin=100,seed=12345);
* Build the means based on sampled values of posterior for the linear model parameters.;
%part2A(Jobname=example_causal,inname=example,method=CAUSAL ,causalk=k1power,ref=PLACEBO,vcmethod=REF);
* Build the imputed values, sampling from the predicted distribution;
%part2b(Jobname=example_causal,seed=54321);
* Lsmeans differences are comparisons to 0 (not all possible combinations);
%part3(Jobname=example_causal,anref=PLACEBO,Label=Causal with sigma from reference);


************************* Provide summary table  ***********************;
title1 "Summary Table_ Causal results with Sigma from REF arm PLACEBO, cov=sex, no Catcov=PoolInv, no covbytime=basval, K1power ";

options linesize=200;
proc print data=work.example_causal_Out noobs; run;

*extracting data from work.example_causal_Out;
* for negative suffixes;
* need mult by 10 to remove dec pt in suffix; 
%let sufx=%sysevalf(&K1*10);
data testnameDIArun3&sufx; 
 set work.example_causal_Out;
   K1=&K1;
   keep K1 diff LCL_Diff UCL_Diff ;
   if (row=8 ) then output; 
run;
%mend cycleDiarunNObasval;
%cycleDiarunNObasval(0);
%cycleDiarunNObasval(0.1);
%cycleDiarunNObasval(0.2);
%cycleDiarunNObasval(0.3);
%cycleDiarunNObasval(0.4);
%cycleDiarunNObasval(0.5);
%cycleDiarunNObasval(0.6);
%cycleDiarunNObasval(0.7);
%cycleDiarunNObasval(0.8);
%cycleDiarunNObasval(0.9);
%cycleDiarunNObasval(1);

data testDIArunNobasvalPLACEBO;
  set 
  testnameDIArun30
  testnameDIArun31
  testnameDIArun32
  testnameDIArun33
  testnameDIArun34
  testnameDIArun35
  testnameDIArun36
  testnameDIArun37
  testnameDIArun38
  testnameDIArun39
  testnameDIArun310;
run;



* try DRUG;
%part2A(Jobname=example_causal,inname=example,method=CAUSAL ,causalk=k1power,ref=DRUG,vcmethod=REF);
* Build the imputed values, sampling from the predicted distribution;
%part2b(Jobname=example_causal,seed=54321);
* Lsmeans differences are comparisons to 0 (not all possible combinations);
%part3(Jobname=example_causal,anref=PLACEBO,Label=Causal with sigma from reference)
title1 "Summary Table_ Causal results with Sigma from REF (DRUG)arm no Catcov=PoolInv, K1power ";
proc print data=work.example_causal_Out noobs;
run;
options linesize=200;
proc print data=work.example_causal_Out noobs; run;




%part1A(Jobname=example
,Data=hamd
,Subject=Patient
,Response=change
,Time=Visit
,Treat=Therapy
,Covbytime= Basval
,covgroup=Therapy
,id=k1
);

%part1B(Jobname=example
,Ndraws=10
,thin=100
,seed=123456
);
*Part 1A and 1B need not to be repeated for different values of k;


* Build the means based on sampled values of posterior for the linear model parameters.;
%part2A(Jobname=example_causal
,inname=example
,method=CAUSAL /*if specify the method as CAUSAL, k should also be specified*/
,causalk=k1 /*By default CausalK takes value 1 and result CIR estimate*/
,ref=PLACEBO
,vcmethod=REF
);

* Build the imputed values, sampling from the predicted distribution;
%part2b(Jobname=example_causal
,seed=654321
);

* Run the analysis section, adding a Label to the results, for when we replay them all together at the end.;
* Lsmeans differences are comparisons to 0 (not all possible combinations);
%part3(Jobname=example_causal
,anref=PLACEBO
,Label=Causal with sigma from reference & k=&K1
);


************************* Provide summary table  ***********************;
title1 "Summary Table_ Causal results with Sigma from REF arm no poolinv, K1power ";
proc print data=work.example_causal_Out noobs;
run;
options linesize=200;
proc print data=work.example_causal_Out noobs; run;
