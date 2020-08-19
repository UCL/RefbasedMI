/*
analyse antidepressants data
causal model with K initially 1 and multiplying by 0.5 each visit
to compare between R and SAS (drop Catcov=PoolInv)
R code - antidep.R
SAS code - antidep.sas
IW 19aug2020
*/

* corresponding SAS code;

*Load the 5-macros;
%let Path = N:\Home\missing\offtrt\Royes Joseph\macros\Five_Macros20161111\TheFiveMacros;
%let Pathc = N:\Home\missing\offtrt\Royes Joseph\macros\Five_Macros20161111\Causal;
%include "&Path\part1a_32.sas";
%include "&Path\part1b_46.sas";
%include "&Pathc\part2A_38_causal.sas"; /*Here is the modified part 2A macro*/
%include "&Path\part2b_30.sas";
%include "&Path\part3_54.sas";

libname library "N:\Home\missing\offtrt\Royes Joseph\data";

proc contents data=library.Chapter15_example;
run;

proc print data=library.Chapter15_example (obs=20);
run;

proc sql;
create table HAMD as
select *, 0.5**(7-max(visit)) as k1power /*k1power=0.5**(tmax-da)*/
  from library.Chapter15_example
group by patient;
quit;

%part1A(Jobname=example, Data=hamd, Subject=Patient, Response=change, Time=Visit, 
        Treat=Therapy, Covbytime=Basval, covgroup=Therapy, id=k1power);
%part1B(Jobname=example,Ndraws=100,thin=100,seed=12345);
* Build the means based on sampled values of posterior for the linear model parameters.;
%part2A(Jobname=example_causal,inname=example,method=CAUSAL,causalk=k1power,ref=PLACEBO,vcmethod=REF);
* Build the imputed values, sampling from the predicted distribution;
%part2b(Jobname=example_causal,seed=54321);
* Lsmeans differences are comparisons to 0 (not all possible combinations);
%part3(Jobname=example_causal,anref=PLACEBO,Label=Causal with sigma from reference);

* analysis not supplied;
