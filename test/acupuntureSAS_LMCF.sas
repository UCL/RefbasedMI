* This program for testing  between SAS 5mscros and Stata and R  mimix
* acupunture dasta and LMCF  when not using the covbytime option
* In the second model basval is head_base treated as the response change at time 0      

*  12/10/2020

* from most uptodate macros;
%let Pathkm=N:\Documents\GitHub\Five_Macros20171010\Five_Macros20171010\TheFiveMacros;
%include "&Pathkm\part1a_33.sas";
%include "&Pathkm\part1b_47.sas";
%include "&Pathkm\part2a_40.sas";
%include "&Pathkm\part2b_31.sas";
%include "&Pathkm\part3_55.sas";
%include "&Pathkm\plotter_7.sas";


libname library "&Pathkm";

*try acupuncture data for lmcf test?;
Proc import out = work.acupuncture DATAFILE="\\ad.ucl.ac.uk\home0\rmjlkm0\Downloads\acupuncture.csv"
   DBMS= csv REPLACE; 
RUN;
%part1A(Jobname=AcupLMCF,Data=acupuncture,Subject=id,Response=head,Time=time, Treat=treat,covgroup=treat);
%part1B(Jobname=AcupLMCF,Ndraws=1000,thin=100,seed=12345);
* Build the means based on sampled values of posterior for the linear model parameters.;
%part2A(Jobname=AcupLMCF_D,inname=AcupLMCF,method=ALMCF,vcmethod=this);
* Build the imputed values, sampling from the predicted distribution;
%part2b(Jobname=AcupLMCF_D,seed=54321);
* Lsmeans differences are comparisons to 0 (not all possible combinations);
%part3(Jobname=AcupLMCF_D,ancovgroup=1,Label=ALMCF with sigma );


ODS HTML  path="N:\Documents" file="acupunctureLMCFSAS.xls" RS=none style=MINIMAL; 
title1 "Summary Table_ LMCF 10 records in dataset with  basval as covariate ";
proc print data=work.AcupLMCF_D_Out noobs; 
%analyseLMCF(104);
%analyseLMCF(333);
%analyseLMCF(386);
%analyseLMCF(787);
%analyseLMCF(435);
%analyseLMCF(697);
%analyseLMCF(100);
%analyseLMCF(101);
run;
ODS HTML CLOSE;
run;
run;
dm 'odsresults; clear';

data work.AcupLMCF_D_data;
   set  AcupLMCF_D_datafinal;
      if (time=0) then head0=head;
       if (time=3) then head3=head;
	     if (time=12) then head12=head;
run;
%macro analyseLMCF( id);
   data  work.example_kmfinalout&id;
     set work.AcupLMCF_D_data;
   if id = &id;
 Title "C Drug Patient_number  &id  ";
Proc means data=work.example_kmfinalout&id;
   var  head0 head3 head12;
  run;
%mend  analyseLMCF;
%analyseLMCF(104);
%analyseLMCF(333);
%analyseLMCF(386);
%analyseLMCF(787);
%analyseLMCF(435);
%analyseLMCF(697);
%analyseLMCF(100);
%analyseLMCF(101);

ODS HTML CLOSE;
run;
