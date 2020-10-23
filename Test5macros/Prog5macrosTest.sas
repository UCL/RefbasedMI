
* This program displays the code which highlights discrepancy between SAS 5mscros and Stata mimix
* in the J2R DRUG reference arm when using the covbytime option
* In the first model basval is treated as a covariate with the response change variable having values at visit4 5,6 visit7 
* In the second model basval is treated as the response change at visit 0      

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

* Hamdvisit0.csv can be downloaded from https://github.com/ucl/mimix;  
Proc import out = work.hamdcsv DATAFILE="\\ad.ucl.ac.uk\home0\rmjlkm0\Documents\GitHub\Five_Macros20171010\Hamdvisit0.csv"
   DBMS= csv REPLACE; 
RUN;
/* Alternatively Fetch the file directly from the GitHub web site */
filename hamdcsvG  temp;
proc http
 url="https://raw.githubusercontent.com/ucl/mimix/master/data/Hamdvisit0.csv"
 method="GET"
 out=probly;
run;
Proc import 
     file=probly
	 out= hamdcsvG replace
	 dbms=csv;
  run;
Proc print data=work.hamdcsvG(obs=10) noobs; run;


* analyse data with no interims (i.e. 3618) and with basval as covariate  (ie take visit 0 out) ;  
data Hamdno3618csv;
  set hamdcsvG;
  if not (patient_number = 3618);
  if not(visit_number=0);
run; 

* to clear ersults viewer run below;
dm 'odsresults; clear'; run;

title1 "1st 10 records in dataset with basval as covariate ";
proc print data=Hamdno3618csv(obs=10); run; 

******************* corresponds to lines 6-13 in spdsht ***********************;  
* J2R ref DRUG Covbytime= basval;
%part1A(Jobname=inBasval,Data=hamdno3618csv,Subject=Patient_number,Response=change,Time=Visit_number, Treat=treatment_name,Covbytime= Basval,covgroup=treatment_name);
%part1B(Jobname=inBasval,Ndraws=1000,thin=100,seed=12345);
* Build the means based on sampled values of posterior for the linear model parameters.;
%part2A(Jobname=inBasval_J2RD,inname=inBasval,method=j2r,ref=DRUG,vcmethod=REF);
* Build the imputed values, sampling from the predicted distribution;
******************* corresponds line 6 spdsht ********************************;
%part2b(Jobname=inBasval_J2RD,seed=54321);
* Lsmeans differences are comparisons to 0 (not all possible combinations);
%part3(Jobname=inBasval_j2rD,anref=DRUG,ancovgroup=1,Label=j2r with sigma from reference);

******************* corresponds line 6 spdsht  (tment effect = 2.1118)  ********************************;
title1 "Summary table J2R Ndraws=1000,thin=100, ref arm drug, basval by covbytime";
proc print data=work.inBasval_j2rd_out(obs=10);
run;
* 6/10/20;
* try alternatve method using reg and mianalyze ;

* need to recode draw to __Imputation_ to use in proc mianalyze;
data inBasval_j2rD_datafinal_imp; 
  set inBasval_j2rD_datafinal ;
    _Imputation_ = draw;
	if visit_number =7;
 run;
%ODSOff
proc reg  data = inBasval_j2rD_datafinal_imp ;
    model change =  basval I_treatment_name;
	by _Imputation_;
	ods output parameterestimates=regparms_inBasval_j2rD;
run;
%ODSOn
******************* corresponds line 7 spdsht ********************************;
title1 "Proc mianalyse after regressing  change=basval + treatment by imputation";  
Proc mianalyze parms=regparms_inBasval_j2rD;
   modeleffects  I_treatment_name ;
   ;
run;

************* should be same if we analyse data as                                   ; 
************* now no basval  but use basval as visit 0 instesad *********************;
data HamdVisit0;
  set hamdcsvG;
  if not (patient_number = 3618);
run;  
title1 "1st 10 records in dataset with basval as value for response variable change at visit 0; ";
proc print data=HamdVisit0(obs=10); run; 

* J2R ref   basval is value for change at visit 0  so NO covbytime;
*********** corresponds to  lines 16-23 in spdsheet ************************;                             ;
%part1A(Jobname=nobasval,Data=HamdVisit0,Subject=Patient_number,Response=change,Time=Visit_number, Treat=treatment_name,covgroup=treatment_name);
%part1B(Jobname=nobasval,Ndraws=1000,thin=100,seed=12345);
* Build the means based on sampled values of posterior for the linear model parameters.;
********* corresponds to line 16 spdshet *************************************;
%part2A(Jobname=nobasval_J2RD,inname=nobasval,method=j2r,ref=DRUG,vcmethod=REF);
* Build the imputed values, sampling from the predicted distribution;
%part2b(Jobname=NObasval_J2RD,seed=54321);
* Lsmeans differences are comparisons to 0 (not all possible combinations);
%part3(Jobname=nobasval_j2rD,anref=PLACEBO,ancovgroup=1,Label=j2r with sigma from reference);
title1 "Summary table J2R Ndraws=1000,thin=100, ref arm Drug  visit 0 and  not covbytime";
proc print data=work.nobasval_j2rD_out noobs;
run;

* Now try fitting basval  nobasval_j2rp_finaldata needs basval ;
* need to recode draw to __Imputation_ to use in proc mianalyze;
title1 " 1st 10 obs datafinal imputed dataset";
proc print data=nobasval_j2rD_datafinal(obs=10); run; 
data nobasval_j2rD_datafinal_vis7; 
  set nobasval_j2rD_datafinal ;
    _Imputation_ = draw;
	if visit_number =7;
 run;
Proc sort Data= nobasval_j2rD_datafinal_vis7 ;
  BY patient_number;
run;
* obtain basval, make unique, in fact just need first ;
data hamdbasval; 
  set hamdno3618csv;
  keep  patient_number basval;
  by patient_number;
  if first.Patient_number then do;
   output;
 end;
run;
DATA nobasval_j2rD_vis7_fin;
  merge nobasval_j2rD_datafinal_vis7  hamdbasval;
  by patient_number;
run;
Proc sort data=nobasval_j2rD_vis7_fin ;BY _imputation_ patient_number ; run;
title1 " 1st 10 obs datafinal imputed dataset with basval";
proc print data=nobasval_j2rD_vis7_fin(obs=10); run; 
%ODSOff
proc reg  data = nobasval_j2rD_vis7_fin ;
    model change =  basval I_treatment_name;
	by _Imputation_;
	ods output parameterestimates=regparms_nobasval_j2rD;
run;
%ODSOn

********* corresponds to line 17 spdshet (tment effect = 1.94) *************************************;
title1 "Proc mianalyse after regressing  change=basval + treatment by imputation";   
Proc mianalyze parms=regparms_nobasval_j2rD;
   modeleffects  I_treatment_name ;
   ;
run;

******* investigate individual patient ids;
data work.nobasval_J2RD_kmfinalout;
   set  nobasval_j2rD_datafinal;
      if (visit_number=4) then change4=change;
       if (visit_number=5) then change5=change;
	    if (visit_number=6) then change6=change;
		 if (visit_number=7) then change7=change;
run;
%macro analyse( id);
   data  work.example_J2R_kmfinalout&id;
     set work.nobasval_J2RD_kmfinalout;
   if Patient_number = &id;
 Title "J2R Drug Patient_number  &id  ";
Proc means data=work.example_J2R_kmfinalout&id;
   var  change4-change7;
  run;
%mend  analyse;
%analyse(3712);


 
%macro ODSoff(); /* calc prior to by processing; */
ods graphics off;
ods exclude all;
ods noresults;
%mend;

%macro ODSon(); /*calc after by processing; */
ods graphics on;
ods exclude none;
ods results;
%mend;



**************** try causal models *********************;
data Hamdno3618csvC;
  set hamdcsvG;
  if not (patient_number = 3618);
  if not(visit_number=0);
k1power=0.5;/*here k is specified by a variable name*/
run; 

* J2R ref DRUG Covbytime= basval;
%part1A(Jobname=inBasCaus,Data=hamdno3618csvC,Subject=Patient_number,Response=change,Time=Visit_number, Treat=treatment_name,Covbytime= Basval,covgroup=treatment_name,id=k1power);
%part1B(Jobname=inBasCaus,Ndraws=1000,thin=100,seed=12345);
* Build the means based on sampled values of posterior for the linear model parameters.;
%part2A(Jobname=inBasCaus_J2RD,inname=inBasCaus,method=Causal,causalk=k1power,ref=DRUG,vcmethod=REF);
* Build the imputed values, sampling from the predicted distribution;
%part2b(Jobname=inBasCaus_J2RD,seed=54321);
* Lsmeans differences are comparisons to 0 (not all possible combinations);
%part3(Jobname=inBasCaus_j2rD,anref=DRUG,ancovgroup=1,Label=j2r with sigma from reference);


title1 "Summary Table_ Causal results with Sigma from REF arm basval k1power=0.5" ;

proc print data=work.inBasCaus_j2rD_Out noobs;
run;


proc print data=hamdno3618csvC(obs=10) noobs ; run;

* J2R ref DRUG  No  Covbytime= basval;
%part1A(Jobname=noBasCaus,Data=hamdno3618csvC,Subject=Patient_number,Response=change,Time=Visit_number, Treat=treatment_name,covgroup=treatment_name,id=k1power);
%part1B(Jobname=noBasCaus,Ndraws=1000,thin=100,seed=12345);
* Build the means based on sampled values of posterior for the linear model parameters.;
%part2A(Jobname=noBasCaus_D,inname=noBasCaus,method=Causal,causalk=k1power,ref=DRUG,vcmethod=REF);
* Build the imputed values, sampling from the predicted distribution;
%part2b(Jobname=noBasCaus_D,seed=54321);
* Lsmeans differences are comparisons to 0 (not all possible combinations);
%part3(Jobname=noBasCaus_D,anref=DRUG,ancovgroup=1,Label=j2r with sigma from reference);


title1 "Summary Table_ Causal results with Sigma from REF arm no basval k1power=0.5" ;

proc print data=work.noBasCaus_D_Out noobs;
run;

* run causal model on visit  data;

data Hamdno3618csvV0;
  set hamdcsvG;
  if not (patient_number = 3618);
  *if not(visit_number=0);
k1power=0.5;/*here k is specified by a variable name*/
run; 

* J2R ref DRUG  No  Covbytime= basval;
%part1A(Jobname=noBasCausV,Data=hamdno3618csvV0,Subject=Patient_number,Response=change,Time=Visit_number, Treat=treatment_name,covgroup=treatment_name,id=k1power);
%part1B(Jobname=noBasCausV,Ndraws=1000,thin=100,seed=12345);
* Build the means based on sampled values of posterior for the linear model parameters.;
%part2A(Jobname=noBasCausV_D,inname=noBasCausV,method=Causal,causalk=k1power,ref=DRUG,vcmethod=REF);
* Build the imputed values, sampling from the predicted distribution;
%part2b(Jobname=noBasCausV_D,seed=54321);
* Lsmeans differences are comparisons to 0 (not all possible combinations);
%part3(Jobname=noBasCausV_D,anref=DRUG,ancovgroup=1,Label=j2r with sigma from reference);

title1 "Summary Table_ Causal results with Sigma from REF arm  basval at visit0 k1power=0.5" ;



Proc print data= noBasCausV_D_datafinal(obs=10) noobs;

* to apply proc mianalyse merge on basval;
data noBasCausV_D_datafinal_vis7; 
  set noBasCausV_D_datafinal ;
    _Imputation_ = draw;
	if visit_number =7;
 run;
Proc sort Data= noBasCausV_D_datafinal_vis7 ;
  BY patient_number;
run;
DATA noBasCausV_D_datafinal_D_mg;
  merge noBasCausV_D_datafinal_vis7  hamdbasval;
  by patient_number;
run;
Proc sort data=noBasCausV_D_datafinal_D_mg ;BY _imputation_ patient_number ; run;
title1 " 1st 10 obs datafinal imputed dataset with basval";
proc print data=noBasCausV_D_datafinal_D_mg(obs=10); run; 
%ODSOff
proc reg  data = noBasCausV_D_datafinal_D_mg ;
    model change =  basval I_treatment_name;
	by _Imputation_;
	ods output parameterestimates=regparms_nobasval_VD;
run;
%ODSOn
title1 "Proc mianalyse after regressing  change=basval + treatment by imputation";   
Proc mianalyze parms=regparms_nobasval_VD;
   modeleffects  I_treatment_name ;
   ;
run;

proc reg  data = nobasval_j2rD_vis7_fin ;
    model change =  basval I_treatment_name;
	by _Imputation_;
	ods output parameterestimates=regparms_nobasval_j2rD;
run;
%ODSOn
title1 "Proc mianalyse after regressing  change=basval + treatment by imputation";   
Proc mianalyze parms=regparms_nobasval_j2rD;
   modeleffects  I_treatment_name ;
   ;
run;
