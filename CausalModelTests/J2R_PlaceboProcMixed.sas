/*
This example imputes under the causal model by specifying option Method = CAUSAL. 
The value of CausalK is specified by CausalK= or defaults to 1 giving the CIR estimate. 
CausalK can be specified using a scalar (a constant $CausalK$ for all participants; e.g.CausalK=0.5) or 
a variable that specifies $CausalK$ for each participant (e.g. CausalK=variable_name}). 
If CausalK is specified using a variable, that variable must be identified in macro 1A using option id=variable_name.
*/


%let path = N:\Documents\GitHub\Five_Macros20171010\Five_Macros20171010\TheFiveMacros;

libname library "&Path";

*%let path=/folders/myfolders/GSKCourse2017/TheFiveMacros/;
%include "&Path\part1a_33.sas";
%include "&Path\part1b_47.sas";
%include "&Path\part2a_40.sas";
%include "&Path\part2b_31.sas";
%include "&Path\part3_55.sas";
%include "&Path\plotter_7.sas";


* Make data set local;
data hamd;
length patient change visit basval 8;
set library.Chapter15_example;
*k1=0.5;/*here k is specified by a variable name*/
* K1=0 equal to J2R;
K1=0;
run;

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
,Ndraws=100
,thin=100
,seed=123456
);


*Part 1A and 1B need not to be repeated for different values of k;



* no need to rerun 1a 1b;

* Build the means based on sampled values of posterior for the linear model parameters.;
*********** Now use the J2R method;
%part2A(Jobname=DIA_J2R_Placebo
,inname=example
,method=J2R
,ref=PLACEBO
)

* Build the imputed values, sampling from the predicted distribution;
%part2b(Jobname=DIA_J2R_Placebo
  ,seed=123456
);


* Lsmeans differences are comparisons to 0 (not all possible combinations);
%part3(Jobname=DIA_J2R_Placebo,anref=PLACEBO,Label=Causal with sigma from reference);

data work.DIA_J2R_Placebo_kmout;
   set  DIA_J2R_Placebo_datafull;
      if (visit=4) then change4=change;
       if (visit=5) then change5=change;
	    if (visit=6) then change6=change;
		 if (visit=7) then change7=change;
run;
%macro analyse( id);
   data  work.DIA_J2R_Placebo_kmout&id;
     set work.DIA_J2R_Placebo_kmout;
   if Patient = &id;

   Title "Causal Placebo &id  ";
Proc means data=work.example_causal_kmout&id;
   var  change4-change7;
  run;
%mend  analyse;


%analyse(3618);



* try proc mixed %mianalyze(Data=Myoutput%mianalyze(Data=Myoutput(where=(visit=7))));

* this is syntax to obtain output as Royes example;
* just tried on one instance of a draw so depends on which draw!  obviously better to use 
* Rubins rules but not clear how to !; 
title "try proc mixed example_causal_kmout(where=(draw=1 and visit=7)";
Proc Mixed data= work.DIA_J2R_Placebo_kmout (where=(draw=52 and visit=7))  Method=REML;
  Class  Therapy visit Poolinv patient ;
 Model change =Therapy Poolinv basval /DDFM=KenwardRoger;
 *Random / subject=patient;
 Repeated / subject=patient Group=therapy;
 LSMeans therapy / Diff;
run; 

dm 'odsresults; clear';
