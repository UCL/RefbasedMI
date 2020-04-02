/*
This macro calculates the modified mean for each subject based on the drawn parameter values.

INPUT PARAMETERS
Jobname=         Identification text that links all the data sets together. [Default is the name of the data set].
                      [It is this jobname that allows the different macros to share information between themselves.]
INName=          Jobname used for input. [This allows several subsequent analyses with different jobnames.]
                    A copy of the file <INname>_Master is made as <jobname>_Master containing pointers.
Method=          Method used for imputation for every record [Default=MAR]
MethodV=         Variable holding the imputation method used for each individual record (Character only
                      If MethodV exist then Method= is ignored.
                      [Method must be one of the predefined methods or start with letter "U".]
Ref=             The treatment level used as the reference arm for every record
RefV=            Variable holding the treatment level used as the reference arm for each individual record.
                      If RefV= exists then Ref= is ignored. (Numeric or character. Must be same type as Treatment variable.)
VCmethod=        Method used for combining variance covariance matrix for the reference and actual arm.
VCmethodv=       Variable holding Vcmethod value for each subject. If VCmethodV exists then VCMethod= is ignored.
VCRef=           The reference arm used in building covariance matrix.
VCRefv=          Variable holding the reference used in building the covariance matrix for this subject.
CausalK=         (rj_change) CausalK can be specified using a scalar (if CausalK is constant across dropout patterns) 
                 or a variable that specifies CausalK for each participant (if CausalK is not constant across dropout patterns)(by default, CausalK=1)
Debug=0          Adds additional output to the log and stops intermediate data sets being deleted at the end.


INPUT DATA SETS
<....>_Master   Data set that holds information about the parameter estimation model and variable characteristics.
<....>_DataP    Data in parallel with pattern information.
<....>_PostLP   Sample from the posterior distribution for Linear Predictor parameters


OUTPUT
<....>_datapM   Means for each subject with missing data after modification based on drawn model parameter values. 
                     [Data with Predicted Mean for each subject with missing data based on drawn values.]

Actions:

Noted limitations:

UPDATES:
20 Aug 2011    JHR: Initial code completed
29 Nov 2011    JHR: Add error message when master file is open (Part2A_11)
03 DEC 2011    JHR: Handling of LMCF type issues by centering on avearge subject effect for covariates.
04 DEC 2011    JHR: Fix error in check of type when using REF as default for VCRef. Rename CDR as CIR (Part2A_14)
05 DEC 2011    JHR: Fix MethodV RefV VCMethodV and VCRefV as never worked. (Part2A_15)
06 DEC 2011    JHR: Add output of MAR predicted values for handling intermediate missing. (Part2A_16)
21 DEC 2011    JHR: Fix bug when no covariates at all. (Part2A_16)
30 Dec 2011    JHR: Add datetime stamp to master and remove Master records (2B etc.) for later Master entries.(Part2A_18)
01 Jan 2012    JHR: Change Pattern, Nmiss and Pindex to have Super_ prefix. (Part2A_19)
04 Jan 2012    JHR: Fix issue with case for VCMethod. (Part2A_20)
24 Jan 2012    JHR: Remove VCMethod SPLIT. (Part2A_24)
26 Jan 2012    JHR: Fix minor bug with VCRef= when there is only one VC matrix. (Part2A_25)
28 Feb 2012    JHR: Fix MethodV failing as method was defaulting to MAR (Part2A_27)
03 Mar 2012    JHR: Delete any Master material from Part2A onwards at start.  (Part2A_28)
11 May 2012    JHR: Better error message when dropping into Otherwise clause.  (Part2A_29)
18 May 2012    JHR: Error correctly when ther are no missing data at all to impute.  (Part2A_30)
10 Jul 2012    JHR: Add methods OFCMCF and AFCMCF. 
16 Jul 2012    JHR: Fix issue when use ALMCF or OLMCF where subject has no data. (Part2A_31)
17 Jul 2012    JHR: Make FCMP subroutines into function so they can reurn error messages. 
18 Jul 2012    JHR: Add Warning message handling, so any Warnings are flagged in the success box at the end. (Part2A_32)
31 Mar 2013    JHR: Fix issue with CatCovbyTime (Part2A_33)
14 APR 2013    JHR: Improve diagnostics when MethodV or RefV variables contain illegal values. [Part1A_34].
18 Apr 2013    JHR: Fix illegal text in diagnostic. [Part2A35]
01 May 2013    JHR: Add B. in SQL to handle case where MehtodV= variable has been included (unecessarily) in the ID= in Part1A. [Part2A36]
22 Apr 2014    JHR: Add code to CIR to stop error when subject has no data. [Part2A37]
29 Jul 2014    JHR: Fix out of bounds in OFCMCF and AFCMCF. [Part2A38]
Modifications by Royes Joseph, MRC BSU (search with key word "rj_change" for locating the changes)
30 Nov 2015     RJ: Defined function "CAUSAL" within proc fcmp [Here 'K_d' is defined as 'CausalK'], a parameter "CausalK" is specified in all function definition for the consistency; 
				    "CausalK" in CAUSAL function takes a default value 1, thus result CIR estimate by default.
					(Now the macro can run the causal method in addition to the RBI methods (MAR, CIR, CR, J2R, OLMCF, ALMCF,OFCMCF and AFCMCF); the default is MAR.) 
18 Dec 2015     RJ: CausalK can be specified using a scalar (a constant CausalK for all participants; e.g. CausalK=0.5) 
					or a variable that specifies CausalK for each participant (e.g. CausalK=variable_name)
04 Jan 2016    RJ: If CausalK is specified a variable name, will check whether the variable name exist in the dataset and it is numeric. 
21 Sep 2017    JHR: Change Jobname to Job2 in proc append of Newbits {part2A40}
                   
*/

*Flush anything remaining;
;run;quit;

options mprint nomlogic;


* These functions define the standard imputation methods;
* They return a null character string if they are successful;
proc fcmp outlib=work.funcs.trial;
* Missing at Random;
function MAR(CausalK,Last,This[*],Ref[*],MyCov[*],Out[*]) $ 128; /*rj_change: CausalK is included*/
	outargs This,ref,Out,Mycov;
	do i=1 to dim(out);
		out[i]=This[i]+Mycov[i];
	end;
	return("");
endsub;

* Copy Increment from Reference;
function CIR(CausalK,Last,This[*],Ref[*],MyCov[*],Out[*]) $ 128;  /*rj_change: CausalK is included*/
	outargs This,ref,Out,Mycov;
	if last=0 then do;
		do i=1 to dim(out);
			out[i]=Ref[i]+Mycov[i];
		end;
		return("Subject has no valid data. Increment for method CIR is set as zero for this subject");
	end;
	else do;
		do i=1 to last;
			out[i]=This[i]+Mycov[i];
		end;
		do i=last+1 to dim(out);
			out[i]=Ref[i]-ref[last]+this[last]+Mycov[i];
		end;
		return("");
	end;
endsub;

* Jump to Reference;
function J2R(CausalK,Last,This[*],Ref[*],MyCov[*],Out[*]) $ 128;  /*rj_change: CausalK is included*/
	outargs This,ref,Out,Mycov;
	do i=1 to last;
		out[i]=this[i]+Mycov[i];
	end;
	do i=last+1 to dim(out);
		out[i]=Ref[i]+Mycov[i];
	end;
	return("");
endsub;

* Copy Reference;
function CR(CausalK,Last,This[*],Ref[*],MyCov[*],Out[*])  $ 128;  /*rj_change: CausalK is included*/
	outargs This,ref,Out,Mycov;
	do i=1 to dim(out);
		out[i]=ref[i]+Mycov[i];
	end;
	return("");
endsub;

* Own Last Mean Carried Forward;
function OLMCF(CausalK,Last,This[*],Ref[*],MyCov[*],Out[*]) $ 128;  /*rj_change: CausalK is included*/
	outargs This,ref,Out,Mycov;
	do i=1 to last;
		out[i]=This[i]+Mycov[i];
	end;
	do i=last+1 to dim(out);
		* Max handles case where subject has no rsponse data at all;
		out[i]=this[max(1,last)]+Mycov[max(1,last)];
	end;
	if last>0 then return("");
	else return("Subject has no valid data. OLMCF imputed using marginal mean at first visit");
endsub;

*Average Last Mean carried forward;
function ALMCF(CausalK,Last,This[*],Ref[*],MyCov[*],Out[*]) $ 128;  /*rj_change: CausalK is included*/
	outargs This,ref,Out,Mycov;
	do i=1 to last;
		out[i]=This[i]+Mycov[i];
	end;
	do i=last+1 to dim(out);
		* Max handles case where subject has no rsponse data at all;
		out[i]=this[max(1,last)]+Mycov[i];
	end;
	if last>0 then return("");
	else return("Subject has no valid data. ALMCF imputed using marginal mean at first visit");
endsub;

* Own First Conditional Mean Carried Forward;
function OFCMCF(CausalK,Last,This[*],Ref[*],MyCov[*],Out[*])  $ 128;  /*rj_change: CausalK is included*/
	outargs This,ref,Out,Mycov;
	do i=1 to min(dim(out),last +1);
		out[i]=This[i]+Mycov[i];
	end;
	do i=last+2 to dim(out);
		out[i]=this[last+1]+Mycov[last+1];
	end;
	return("");
endsub;

*Average First Conditional Mean carried forward;
function AFCMCF(CausalK,Last,This[*],Ref[*],MyCov[*],Out[*]) $ 128;  /*rj_change: CausalK is included*/
	outargs This,ref,Out,Mycov;
	do i=1 to min(dim(out),last +1);
		out[i]=This[i]+Mycov[i];
	end;
	do i=last+2 to dim(out);
		out[i]=this[last+1]+Mycov[i];
	end;
	return("");
endsub;

 /*rj_change: New function "causal" is defined*/
* Causal model 1 with specific value of CausalK, CausalK independent on time since withrawal, CausalK_d=CausalK;
*But one can specify different CausalK for different pattern to account how long been on treatment or reason for discontinuation;
function CAUSAL(CausalK,Last,This[*],Ref[*],MyCov[*],Out[*]) $ 128;
	outargs This,ref,Out,Mycov;
	if last=0 then do;
		do i=1 to dim(out);
			out[i]=Ref[i]+Mycov[i];
		end;
		return("Subject has no valid data. Increment for method CAUSAL is set as zero for this subject");
	end;
	else do;
		do i=1 to last;
			out[i]=This[i]+Mycov[i];
		end;
		do i=last+1 to dim(out);
			out[i]=Ref[i]+ CausalK*(this[last]- ref[last])+Mycov[i];
		end;
		return("");
	end;
endsub;

quit;
options cmplib=work.funcs;


%macro part2A(jobname= 
,INName=
,method=
,methodv=
,ref=
,refv=
,VCmethod= 
,VCmethodv=
,VCref=
,VCrefv= 
,debug=0
,CausalK=1   /*rj_change: CausalK included and assigned default value*/
);

%**************** Approved list defined here *************************;
%* Note that any methods starting with the letter "U" (e.g.User1) are also allowed;
%*    ... but must have been compiled into work.funcs;
%local approved;
%let approved=MAR CIR CR J2R OLMCF ALMCF OFCMCF AFCMCF CAUSAL;  /*rj_change: 'CAUSAL' is included*/

* Defaults;
%if %length(&Debug) = 0 %then %let debug=0;
%if %length(&method)=0 and %length(&methodV)=0 %then %let method=MAR;
%let method=%upcase(&Method);

%* Local manipulation macro names;
%local macname Stopflag WarnFlag Subject PEWhere Treat Treattype Cov STD_Cov CovbyTime STD_CovbyTime CatCov I_CatCov
	CatcovbyTime I_CatcovbyTime Covgroup Covgptype;
%local Ntreat Ntimes sumcheck1 sumcheck2 sumcheck3 sumcheck4;
%local i txt decarr store calc keep loc1A loc1B job1 job2 Dsname errtext1 errtext2;



%* Macname       : Name of macro for printing in error messages etc;
%* Stopflag      : Flag indicate that at end of the Data step macro should end immediately;
%* Subject       : Name of subject variable;
%* Treat         : Name of treatment variable;
%* Treattype     : Type of the treatment variable;
%* Cov           : Cov= parameter from Part1A;
%* STD_Cov       : List of regression covariates prefixed by STD_;
%* CatCov        : CatCov= parameter from Part1A;
%* I_CatCov      : List of Catcov variables prefixed by I_;
%* CovbyTime     : CovbyTime= parameter from Part1A;
%* STD_CovbyTime : List of CovbyTime variates prefixed by STD_;
%* CatCovbyTime  : CatCovbyTime= parameter from Part1A;
%* I_CatCovbyTime: List of CatCovbyTime variabls prefixed by I_;
%* Covgroup      : Covgroup= parameter from Part1A;
%* CovgpType     : Type of the Covgrooupp variable from Part1A;
%* Ntreat        : Number of treatment levels;
%* Ntimes        : Number of levels for Time;

%let errtext1= ;
%let errtext2= ;

%let macname=Part2A;
%let stopflag=0;
%let warnflag=0;
%let Cov=;
%let CatCov=;
%let Covbytime=;
%let Catcovbytime=;
%let STD_Cov=;
%let I_CatCov=;
%let STD_Covbytime=;
%let I_Catcovbytime=;


%* Set Inname to jobname for standard case (continue in Jobname);
%if %length(&INName) =0 %then %do;
	%* Read from the jobname location;
	%let INName=&jobname;
%end;


%* Check if jobname has a valid master data set; 
%if not %sysfunc(exist(&Inname._master)) %then %do;
	%let errtext1=Master Data set does not exist with name <&Inname._master>; 
	%goto errtrap;
%end;

%* If this is a new tree branch then create new master data set;
%if &INName ^= &Jobname %then %do;
	* Copy the master data set; 
	data &jobname._master;
	set &inname._Master end=Myend;
	run;
	%if &syserr %then %do;
		%let errtext1=Error creating <&jobname._master> from <&inname._Master>: Have you got it open for editing?; 
		%goto errtrap;
	%end;
%end;

%* Loc1A and Loc1B are the locations for input and for output data sets. This allows one to use the Innmae facility;
%let loc1A=;
%let loc1B=;
proc sql noprint;
* Set status of Part1B as failed in case of error;
delete from &jobname._master where Role="Master2A";
insert into &jobname._master
set Vlabel="Failed", Role="Master2A";
quit;
%if &syserr %then %do;
	%if &syserr = 1012 %then %do;
		%let errtext1=Error accessing <&jobname._master>: Have you got it open for editing?; 
		%goto errtrap;
	%end;
	%else %goto myabort;
%end; 

proc sql noprint;
* Remove any debris from existing Master data set;
delete from &jobname._master where Role in(
	"Method","MethodV","Ref","RefV","VCRef","VCRefV","VCMeth","VCMethV",
	"Master2B",
	"Master3","ANWhere","ANTimes","ANRef","ANTreat","ANCov","ANCatCov","ANCovgp","Delta","DeltaV","DLag","Dgroups","DgroupsV");
delete from &jobname._master where Role="Seed" and Vlabel="Imputation";


* Check if previous sections (PART1A and PART1B) ran correctly;
select Vlabel into :loc1A from &inname._master(where=(Role="Master1A"));
%if %length(&loc1A) = 0 %then %do;
	%let errtext1=Part 1A has not been run for <&Inname.>; 
	%goto errtrap;
%end;
%if &loc1A = Failed %then %do;
	%let errtext1=Part 1A has been run but failed for <&Inname.>; 
	%goto errtrap;
%end;
%* Strip off datetime stamp;
%let Loc1A=%scan(&Loc1A,1,@);

select Vlabel into :loc1B from &inname._master(where=(Role="Master1B"));
%if %length(&loc1B) = 0 %then %do;
	%let errtext1=Part 1B has not been run for <&Inname.>; 
	%goto errtrap;
%end;
%if &loc1B = Failed %then %do;
	%let errtext1=Part 1B has been run but failed for <&Inname.>; 
	%goto errtrap;
%end;
%* Strip of datetime stamp;
%let Loc1B=%scan(&Loc1B,1,@);

* Reload the macro variables;
select vlabel into :DSName from &inname._master(where=(Role="DSName")); 
select vname into :Subject from &inname._master(where=(Role="Subject")); 
select vname into :PEWhere from &inname._master(where=(Role="PEWhere")); 
select vtype, vname into :Treattype, :Treat from &inname._master(where=(Role="Treat" and index=1)); 
select vname, "STD_"||Vname into :Cov separated by " " , :STD_cov separated by " " 
	from &inname._master(where=(Role="Cov")); 
select vname, "STD_"||Vname into :CovbyTime separated by " " , :STD_CovbyTime separated by " "
	from &inname._master(where=(Role="Covbytim")); 
select vname, "I_"||Vname into :Catcov separated by " " , :I_Catcov separated by " " 
	from &inname._master(where=(Role="CatCov" and index=1)); 
select vname, "I_"||Vname into :Catcovbytime separated by " " , :I_Catcovbytime separated by " " 
	from &inname._master(where=(Role="CatCovby" and index=1)); 
select vtype, vname into :Covgptype, :Covgroup from &inname._master(where=(Role="Covgp")); 
select max(index) into :Ntimes from &inname._master(where=(Role="Time")); 
select max(index) into :NTreat from &inname._master(where=(Role="Treat")); 
select max(draw) into :Ndraws from &loc1B._postlp; 
quit;

%* Left justify and trim;
%let Dsname=&DsName;
%let Treat=&Treat;
%let Ntreat=&Ntreat;
%let Ntimes=&Ntimes;
%let Subject=&Subject;
%let Ndraws=&Ndraws;

%* Put out recovered values;
%let Ntimes=&Ntimes;
%let DSname=&Dsname;
%let Subject=&Subject;

%if &debug %then %do;
	%put Treat=&Treat;
	%put Treattype=&Treattype;
	%put Cov=&Cov;
	%put STD_Cov=&STD_Cov;
	%put Covbytime=&Covbytime;
	%put STD_Covbytime=&STD_Covbytime;
	%put CatCov=&CatCov;
	%put I_CatCov=&I_CatCov;
	%put CatCovbytime=&CatCovbytime;
	%put I_CatCovbytime=&I_CatCovbytime;
	%put Covgroup=&Covgroup;
	%put Covgptype=&Covgptype;
	%put NTreat=&NTreat;
	%put NDraws=&NDraws;
%end;

%let job1=%scan(&jobname.,1,%str(.));
%let job2=%scan(&jobname.,2,%str(.));
%let txt=%scan(&jobname.,3,%str(.));
%if %length(&txt) %then %do;
	%let errtext1=The parameter Jobname <&Jobname.> contains more than one dot (.) so has more than two parts.; 
	%goto errtrap;
%end;

%let job1=&job1;
%let job2=&job2;
%if &debug %then %put job1=&job1;
%if &debug %then %put job2=&job2;

%if %length(&job2) %then %do;
	%if %sysfunc(libref(&job1)) %then %do;
		%let errtext1=The parameter Jobname <&Jobname.> contains a nonexistent libref; 
		%goto errtrap;
	%end;
%end;
%else %do;
	%let job2=&job1;
	%let job1=WORK;
%end;

* Check if jobname is a valid sas name; 
data _null_;
&job2=0;
run;
%if &syserr %then %do;
	%let errtext1=The parameter Jobname <&jobname> is not a legal SAS name; 
	%goto errtrap;
%end;


proc datasets library=&job1 nolist;
delete &job2._temp21 &job2._temp22 &job2._temp39 &job2._temp40 &job2._temp41 &job2._temp42 &job2._temp51
	&job2._datapM &job2._newbits;
%if &syserr %then %goto myabort; 
quit;
%if &syserr %then %goto myabort; 



* Check required variables are in the original data set;
data &job2._newbits;
keep Vname vlabel vnum vtype Vlen Vformat Role Clevel Nlevel Index;
length txt $256 Vname $32 vlabel $256 vnum 8 vtype $1 Role $8 Clevel $120 vformat $200;
drop txt dsid;
Clevel=" ";
Nlevel=.;
index=.;
dsid=open("&DSName","i");
if dsid<=0 then do;
	call symput("Stopflag","1");
	call symput("Errtext1","The original Data set <&DSName.> no longer exists");
	stop;
end;
%if %length(&Methodv) %then %do;
	vnum=varnum(dsid,"&MethodV");
	if vnum<=0 then do;
		call symput("Stopflag","1");
		call symput("Errtext1","MethodV= parameter: Variable <&MethodV> not found in Data set &DSName.");
		stop;
	end;
	vtype=vartype(dsid,vnum);
	if vtype ^= "C" then do;
		call symput("Stopflag","1");
		call symput("Errtext1","MethodV= parameter: Variable <&MethodV> in Data set <&DSName.> must be of type Character.");
		stop;
	end;
	Role="MethodV";
	Vname=varname(dsid,Vnum);
	vlen=varlen(dsid,Vnum);
	vlabel=varlabel(dsid,vnum);
	vformat=varfmt(dsid,vnum);
	output;	
	%if %length(&Method) %then %do;
		call symput("Stopflag","1");
		call symput("Errtext1","Method= parameter: A Methodv= has been declared. Method= should not be used. ");
		stop;
	%end;
%end;
%else %do;
	%if %length(&Method) = 0 %then %do;
		%let Method=MAR;
	%end;
	Role="Method";
	Vname=" ";
	vnum=.;
	Vtype=" ";
	vlen=.;
	vlabel="&Method";
	vformat=" ";
	output;	
%end;

%if %length(&Refv) %then %do;
	%if &Refv =&Treat %then %do;
		%let errtext1= Refv setting <&refv.> is the treatment variable. This is both silly and not allowed.; 
		%let errtext2= Change the setting for RefV; 
		%goto errtrap;
	%end;
	vnum=varnum(dsid,"&RefV");
	if vnum<=0 then do;
		call symput("Stopflag","1");
		call symput("Errtext1","RefV= parameter: Variable <&RefV> not found in Data set <&DSName.>");
		stop;
	end;
	vtype=vartype(dsid,vnum);
	if vtype ^= "&Treattype"  then do;
		call symput("Stopflag","1");
		call symput("Errtext1","RefV= parameter: Variable <&RefV.> in Data set <&DSName.> must be the same type as treatment variable <&Treat.>");
		stop;
	end;
	Role="RefV";
	Vname=varname(dsid,Vnum);
	vlen=varlen(dsid,Vnum);
	vlabel=varlabel(dsid,vnum);
	vformat=varfmt(dsid,vnum);
	output;	
	%if %length(&Ref) %then %do;
		call symput("Stopflag","1");
		call symput("Errtext1","Ref= parameter: A Refv= has been declared. Ref= should not be used. ");
		stop;
	%end;
%end;
%else %do;
	%if %length(&Ref) %then %do;
		Role="Ref";
		Vname=" ";
		vnum=.;
		Vtype=" ";
		vlen=.;
		vlabel="&Ref";
		vformat=" ";
		output;	
	%end;
%end;


%if %length(&VCMethodv) %then %do;
	vnum=varnum(dsid,"&VCMethodV");
	if vnum<=0 then do;
		call symput("Stopflag","1");
		call symput("Errtext1","MethodV= parameter: Variable <&VCMethodV> not found in Data set <&DSName.>");
		stop;
	end;
	vtype=vartype(dsid,vnum);
	if vtype ^= "C" then do;
		call symput("Stopflag","1");
		call symput("Errtext1","VCMethodV= parameter: Variable <&VCMethodV> in Data set <&DSName.> must be of type Character.");
		stop;
	end;
	Role="VCMethV";
	Vname=varname(dsid,Vnum);
	vlen=varlen(dsid,Vnum);
	vlabel=varlabel(dsid,vnum);
	vformat=varfmt(dsid,vnum);
	output;	
	%if %length(&VCMethod) %then %do;
		call symput("Stopflag","1");
		call symput("Errtext1","VCMethod= parameter: A VCMethodV= has been declared. VCMethod= should not be used. ");
		stop;
	%end;
%end;
%else %do;
	%if %length(&VCMethod) %then %do;
	if Upcase("&VCMethod") notin("THIS", "REF","ZERO"  ) then do;
		call symput("Stopflag","1");
		call symput("Errtext1","VCMethod= parameter: Must be one of THIS, REF, ZERO");
		stop;
	end;
	Role="VCMeth";
	Vname=" ";
	vnum=.;
	Vtype=" ";
	vlen=.;
	vlabel="&VCMethod";
	vformat=" ";
	output;	
	%end;
%end;

%if %length(&VCRefv) %then %do;
	vnum=varnum(dsid,"&VCRefV");
	if vnum<=0 then do;
		call symput("Stopflag","1");
		call symput("Errtext1","VCRefV= parameter: Variable <&VCRefV> not found in Data set &DSName.");
		stop;
	end;
	vtype=vartype(dsid,vnum);
	if vtype ^= "&Covgptype"  then do;
		call symput("Stopflag","1");
		call symput("Errtext1","VCRefV= parameter: Variable <&VCRefV.> in Data set &DSName. must be the same type as Covariance Group variable <&Covgroup.>");
		stop;
	end;
	Role="VCRefV";
	Vname=varname(dsid,Vnum);
	vlen=varlen(dsid,Vnum);
	vlabel=varlabel(dsid,vnum);
	vformat=varfmt(dsid,vnum);
	output;	
	%if %length(&VCRef) %then %do;
		call symput("Stopflag","1");
		call symput("Errtext1","VCRef= parameter: A VCRefV= has been declared. VCRef= should not be used. ");
		stop;
	%end;
%end;
%else %do;
	%if %length(VCRef) %then %do;
	Role="VCRef";
	Vname=" ";
	vnum=.;
	Vtype=" ";
	vlen=.;
	vlabel="&VCRef";
	vformat=" ";
	output;	
	%end;
%end;
/*rj_change : check if CausalK is specified by a numeric variable*/
%if %length(&Method) %then %do; 
if upcase("&method") = "CAUSAL" and anyalpha("&CausalK")^=0 then do;
	vnum=varnum(dsid,"&CausalK");
	if vnum<=0 then do;
		call symput("Stopflag","1");
		call symput("Errtext1","CausalK= parameter: Variable <&CausalK> not found in Data set <&DSName.>");
		stop;
	end;
	vtype=vartype(dsid,vnum);
	if vtype = "C" then do;
		call symput("Stopflag","1");
		call symput("Errtext1","CausalK= parameter: Variable <&CausalK> in Data set <&DSName.> must be of type Numeric.");
		stop;
	end;
	Role="CausalK";
	Vname=varname(dsid,Vnum);
	vlen=varlen(dsid,Vnum);
	vlabel=varlabel(dsid,vnum);
	vformat=varfmt(dsid,vnum);
	output;
end;	
%end;
rc=close(dsid);
run;
%if &stopflag %then %goto errtrap; 
%if &syserr %then %goto myabort; 


* Extract a standard set of linear predictor parameter values;
proc sort data=&loc1B._postlp out=&jobname._temp51(where= (draw=1));
by vname icat itime;
run;
%if &syserr %then %goto myabort; 

* Build text to 1) declare arrays, 2) store parameter values, 3) calculate linear predictor, 4) keep in data set;
data _Null_;
set &jobname._temp51 end=myend;
by vname icat itime;
length txt1 $ 1024;
length txt2 $ 1024;
length txt3 $ 1024;
length txt4 $ 1024;
retain txt1 txt2 txt3 txt4;
if last.vname then do;
	if upcase(Vname) ^= upcase("I_&Treat") then do;
		if icat^=. and itime^=. then do;
			*CATCOVbyTime;
			txt1=trim(txt1) || "array P_" || trim(Vname) || " [&Ndraws," || trim(left(put(icat,6.0))) ||"," ||
				trim(left(put(itime,6.0))) || "] _temporary_;" ;
			txt2=trim(txt2) || "if Vname='" || trim(Vname) || "' then P_" || trim(Vname) || "[draw,icat,itime]=Value;" ;
			txt3=trim(txt3) || "+P_" || trim(Vname) || "[draw," || trim(Vname) ||",itime]" ;
			txt4=trim(txt4) || ",B." || trim(Vname);
		end;
		else do;
			if icat=. and itime^=. then do;
				*CovbyTime;
				txt1=trim(txt1) || "array P_" || trim(Vname) || " [&Ndraws," || trim(left(put(itime,6.0))) || "] _temporary_;" ;
				txt2=trim(txt2) || "if Vname='" || trim(Vname) || "' then P_" ||trim(Vname) || "[draw,itime]=Value;" ;
				txt3=trim(txt3) || "+P_" || trim(Vname) || "[draw,itime]*" || trim(Vname) ;
				txt4=trim(txt4) || ",B." || trim(Vname);
			end;
			else do;
				if icat^=. and itime=. then do;
					*CatCov;
					txt1=trim(txt1) || "array P_" || trim(Vname) || " [&Ndraws," || trim(left(put(icat,6.0))) || "] _temporary_;" ;
					txt2=trim(txt2) || "if Vname='" || trim(Vname) || "' then P_" ||trim(Vname) || "[draw,icat]=Value;" ;
					txt3=trim(txt3) || "+P_" || trim(Vname) || "[draw," || trim(Vname) ||"]" ;
					txt4=trim(txt4) || ",B." || trim(Vname);
				end;
				else do;
					*Cov;
					txt1=trim(txt1) || "array P_" || trim(Vname) || " [&Ndraws] _temporary_;" ;
					txt2=trim(txt2) || "if Vname='" || trim(Vname) || "' then P_" ||trim(Vname) || "[draw]=Value;" ;
					txt3=trim(txt3) || "+P_" || trim(Vname) || "[draw]*" || trim(Vname);
					txt4=trim(txt4) || ",B." || trim(Vname);
				end;
			end;
		end;
	end;
end;
if myend then do;
	call symput("decarr",txt1);
	call symput("store",txt2);
	call symput("calc",txt3);
	call symput("keep",txt4);
end;
run;
%if &syserr %then %goto myabort; 

%if &debug %then %do;
	%put decarr=&decarr;
	%put store=&store;
	%put calc=&calc;
	%put keep=&keep;
%end;

* Check that the MethodV=, RefV=, VCMethodv= and VCRefV= variables are unique with subject;
* And extract from original data set;
%if %length(&MethodV.&RefV.&VCMethodV.&VCRefV.) %then %do;
	%let sumcheck1=0;
	%let sumcheck2=0;
	%let sumcheck3=0;
	%let sumcheck4=0;

	proc sql noprint;
	create table &jobname._temp39 as
	select &Subject 
		%if %length(&MethodV) %then %do; ,Max(&MethodV) as &MethodV, Min(&MethodV) as &MethodV._check  %end;
		%if %length(&RefV) %then %do; ,Max(&RefV) as &RefV, Min(&RefV) as &RefV._check %end;
		%if %length(&VCMethodV) %then %do; ,Max(&VCMethodV) as &VCMethodV, Min(&VCMethodV) as &VCMethodV._check  %end;
		%if %length(&VCRefV) %then %do; ,Max(&VCRefV) as VCRefV, Min(&VCRefV) as &VCRefV._check %end;
	from &dsname
	group by &Subject;
	%if &sqlrc %then %goto myabort;

	select 0
		%if %length(&MethodV) %then %do; ,sum(&MethodV ^= &MethodV._check) %end;
		%if %length(&RefV) %then %do; ,sum(&RefV ^= &RefV._check) %end;
		%if %length(&VCMethodV) %then %do; ,sum(&VCMethodV ^= &VCMethodV._check %end;
		%if %length(&VCRefV) %then %do; ,sum(&VCRefV ^= &VCRefV._check %end;
		into :dummy
		%if %length(&MethodV) %then %do; ,:sumcheck1 %end;
		%if %length(&RefV) %then %do; ,:sumcheck2 %end;
		%if %length(&VCMethodV) %then %do; ,:sumcheck3 %end;
		%if %length(&VCRefV) %then %do; ,:sumcheck4 %end;
	from &jobname._temp39
	%if &sqlrc %then %goto myabort;
	quit;
	%if &syserr %then %goto myabort; 

	%if %eval(&sumcheck1 + &sumcheck2 + &sumcheck3 + &sumcheck4) %then %do;
		%let errtext1= The following variables need to be unique with subject: ; 
		%if &sumcheck1 %then %let errtext1= &errtext1. MethodV=&Methodv. ;
		%if &sumcheck2 %then %let errtext1= &errtext1. RefV=&Refv. ;
		%if &sumcheck3 %then %let errtext1= &errtext1. VCMethodV=&VCMethodv. ;
		%if &sumcheck4 %then %let errtext1= &errtext1. VCRefV=&VCRefv. ;
		%goto errtrap;
	%end;

	* Add the required variables to the data;
	proc sql;
	create table &jobname._temp40 as
	select A.*
		%if %length(&MethodV) %then %do; ,B.&MethodV  %end;
		%if %length(&RefV) %then %do; ,B.&RefV %end;
		%if %length(&VCMethodV) %then %do; ,B.&VCMethodV %end;
		%if %length(&VCRefV) %then %do; ,B.&VCRefV %end;
	from &loc1B._datap(where=(Super_Nmiss>0)) A left join &jobname._temp39 B
	on A.&subject=B.&subject;
	quit;

	* Build required control variables;
	data &jobname._temp41;
	set &jobname._temp40;
%end;
%else %do;
* Build required control variables;
	data &jobname._temp41;
	set &loc1B._datap(where=(Super_Nmiss>0));
%end;

	keep &subject Subject_index I_&Treat  &I_Catcov &I_CatcovbyTime &STD_Cov &STD_CovbyTime
	%if &Covgroup ^= 1 %then %do; I_&Covgroup  %end;
	%if %length(&Methodv) %then %do; &Methodv  %end;
	%if %length(&Refv) %then %do; &Refv %end;
	%if %length(&VCMethodv) %then %do; &VCMethodv %end;
	Super_Response1-Super_Response&Ntimes
	Super_Nmiss Super_Pattern Super_Last Super_k /*rj_change : super_k is added and then declared below */
	Super_method Super_Ref  Super_VCMethod Super_VCRef;
 
* Declare Super_ref and super_VCRef with correct type;
length Super_method $8 Super_VCMethod $8 Super_k 8 ;

%* If methodv is there then copy it, and if method exist read it into new variable (default is MAR);
* All methods translated to upper case for internal use;
%if %length(&Methodv) %then Super_Method= UPCASE(&Methodv);
%if %length(&Methodv) = 0 and %length(&Method) %then Super_Method=UPCASE("&Method");
%if %length(&Methodv) = 0 and %length(&Method) = 0 %then Super_Method="MAR";;

%if %length(&CausalK) %then Super_k=&CausalK;; /*rj_change*/

* If reference variable exist copy it in, and if not then read it in from Ref. (Default is a missing value);
%if %length(&Refv) %then Super_Ref=&Refv; 
%if %length(&Refv) = 0 and %length(&Ref) %then %do;
	%if &Treattype = C %then %do;
		length Super_ref $%if %length(&Ref) <= 1 %then 2; %else %length(&Ref);;
		Super_Ref="&Ref";
	%end;
	%else %do;
		Super_Ref=&Ref;
	%end;
%end;
%if %length(&Refv) = 0 and %length(&Ref)=0 %then %do;
	%if &Treattype = C %then %do;
		Super_Ref=" ";
	%end;
	%else %do;
		Super_Ref=.;
	%end;
%end; ;

%if %length(&VCMethodv)%then Super_VCMethod=upcase(&VCMethodV);
%if %length(&VCMethodv)=0 and %length(&VCMethod) %then Super_VCMethod=upcase("&VCMethod");
%if %length(&VCMethodv)=0 and %length(&VCMethod)=0 %then Super_VCMethod=" ";;

%if %length(&VCRefv) %then Super_VCRef=&VCRefv; 
%if %length(&VCRefv) = 0 and %length(&VCRef) %then %do;
	%if &Covgptype = C %then %do;
		length Super_VCref $%length(&VCRef);
		Super_VCRef="&VCRef";
	%end;
	%else %do;
		Super_VCRef=&VCRef;
	%end;
%end;
%* If VCRef not set then default it to REF as long as Type of variables are consistent;
%*    but no ned if only useing one covariance group;
%if %length(&VCRefv) = 0 and %length(&VCRef)=0 %then %do;
	%if &Covgroup = I_1 %then %do;
		Super_VCRef=I_1;
	%end;
	%else %do;
		%if &Covgptype = &Treattype %then %do;
			Super_VCRef=Super_Ref;
		%end;
		%else %do;
			%let errtext1=Cannot use REF as default for VCRef when Covgroup= and Treat= variables are of dfferent type; 
			%goto errtrap;
		%end;
	%end;
%end;

run;
%if &syserr %then %goto myabort; 


* Now look up methods used;
* But first check we found something;
proc sql noprint;
select count(*) into :checkn from &jobname._temp41;
%if &checkn <= 0 %then %do;
	%let errtext1=There are no missing data found at all.; 
	%let errtext2=The macros will not proceed.; 
	%goto errtrap;
%end;
select distinct upcase(super_method) into :methods separated by " " 
from &jobname._temp41;
%if &sqlrc %then %goto myabort;
%if &debug %then %put Methods=&methods;

* and index the reference values;
create table &jobname._temp42 as
select A.*, index as I_ref
from &jobname._temp41 A left join &inname._master(where=(role="Treat")) B
%if &Treattype = N %then on A.super_ref = B.Nlevel;
%else on trim(left(A.super_ref)) = trim(left(B.Clevel)); ;
%if &sqlrc %then %goto myabort;
quit;
%if &syserr %then %goto myabort; 

%let i=1;
%let txt=%scan(&methods,&i,%str( ));
%do %while(%length(&txt));
	%let errtext1=Method <&txt> not defined. [New codes should be added to Approved list in Part1A];
	%if %upcase(%substr(&txt,1,1)) = U %then %do;
		%let errtext1=; 
	%end;
	%else %do;
		%let j=1;
		%let txt2=%scan(&approved,&j,%str( ));
		%do %while(%length(&txt2));
			%if &txt = &txt2 %then %do;
				%let errtext1=; 
			%end;
			%let j=%eval(&j +1);
			%let txt2=%scan(&approved,&j,%str( ));
		%end;
	%end;
	%if %length(&errtext1)%then %do;
		%goto errtrap;
	%end;
	%let i=%eval(&i +1);
	%let txt=%scan(&methods,&i,%str( ));
%end;

data &jobname._datapm;
set &loc1B._postlp(in=inpostlp)  &loc1B._datap(
	%if %length(&PEWhere) %then %do; where=(&PEWhere) %end;
	in=pass1) &jobname._temp42(in=pass2) ;
array values[&ndraws,&ntreat,&Ntimes] _temporary_;
array avgmeanp[&Ndraws,&Ntimes]  _temporary_ (%eval(&Ndraws*&Ntimes)*0);
retain counter 0;
%* declare arrays;
&decarr;

array super_mean[&Ntimes] super_mean1-super_mean&Ntimes;
array super_mar[&Ntimes] super_mar1-super_mar&Ntimes;
array A[&Ntimes] _temporary_;
array B[&Ntimes] _temporary_;
array C[&Ntimes] _temporary_;
array MyCov[&Ntimes] _temporary_;
keep &subject draw super_mean1-super_mean&Ntimes super_mar1-super_mar&Ntimes
	subject_index Super_Response1-Super_Response&Ntimes Super_Pattern Super_Nmiss Super_Last 
		Super_VCMethod Super_VCRef Super_k /*rj_change*/
		%if &Covgroup ^= 1 %then %do; I_&Covgroup  %end;;
length warntxt $128;
retain warncount 0;

if inpostlp then do;
	* Copy parameter values into the values array;
	if upcase(Vname) = upcase("I_&Treat") then do;
		values[draw,icat,itime]=value;
	end;
	else do;
		%* Store away the rest;
		&store
	end;
end;

if pass1 then do;
	counter=counter+1;
	do draw=1 to &Ndraws;
		* &calc contains code to calculate the rest of the linear predictor.;
		do itime=1 to &Ntimes;
			tempval=avgmeanp[draw,itime];
			avgmeanp[draw,itime]= tempval + ((0 &calc) - tempval)/counter ;
		end;
	end;
end;

if pass2 then do;
	* Set up defaults for VCMethod;
	if Super_VCMethod=" " then do;
		if super_method in("MAR","OLMCF","ALMCF","OFCMCF","AFCMCF") then Super_VCMethod="THIS";
		else if super_method in("CIR","J2R", "CAUSAL") then super_VCMethod="REF"; /*rj_change : CAUSAL is included*/
		else if super_method="CR" then super_VCMethod="REF";
		* If method not found then use "REF" if there is an Reference is set;
		*  or if no reference set then use "This";
		else if I_Ref=. then super_VCMethod="THIS";
		else super_VCMethod="REF";
	end;
	* Build in checks that REF= is set for standard methods that need it;
	if (super_method in("CIR","J2R","CR", "CAUSAL")) and (I_Ref < 1 or I_Ref > &Ntreat) then do; /*rj_change : CAUSAL is included*/
		call symput("Stopflag","1");
		call symput("Errtext1","Illegal value <Index= " || I_Ref || "> for Ref= or for content of Refv= variable. Method is "||
			left(trim(super_method)) || " Subject ID= <" || &subject || ">.");
		stop;
	end;

	* This makes user only one warning for each record in data set; 
	Warn=0;
	do draw=1 to &Ndraws;
		do itime=1 to &Ntimes;
			* We hand in the Treat*Time predicted mean for the average subject;
			A[itime]= values[draw,I_&treat,itime] +  AvgMeanp[draw,itime];
			if I_ref>.z then B[itime]= values[draw,I_Ref,itime] +  AvgMeanp[draw,itime];
			else B[itime]= .;
			* Mycov is the difference between this subjects predicted mean and that for an average usbject;
			MyCov[itime]= (0 &calc) - AvgMeanp[draw,itime];
		end;
		select(Super_Method);
			%let i=1;
			%let txt=%scan(&methods,&i);
			%do %while(%length(&txt));
				when("&txt") warntxt= &txt(Super_k,Super_Last,A,B,Mycov,C);/*rj_change :added super_K*/
				%let i=%eval(&i +1);
				%let txt=%scan(&methods,&i);
			%end;
			otherwise do;
				put "ERROR: In otherwise clause";
				put "ERROR: Trying to use Method <" Super_method ">";
				call symput("Stopflag","1");
				call symput("Errtext1","There is some problem with definition of Method");
				call symput("Errtext2","Please check this aspect of your macro call. See above for details.");
				stop;
			end;
		end;
		if warn=0 and (warntxt ^= "") then do;
			if warncount<10 then do;
				put "WARNING: For subject, &Subject=" &Subject;
				put "WARNING: " warntxt;
				call symput("Warnflag","1");
			end;
			warncount +1;
			warn +1;
		end;
		do itime=1 to &Ntimes;
			super_mean[itime]=C[itime];
			super_mar[itime]=A[itime]+MyCov[itime];
		end;
		output;
	end;
end;
run;
%if &stopflag %then %goto errtrap; 
%if &syserr %then %do;
		%let errtext1= Problems in calculating the imputation stage.;
		%if %length(&methodv) %then %let errtext2= Most likley problem is illegal Method in &Methodv variable; 
		%if %length(&refv) %then %let errtext2= &errtext2 or illegal Ref in RefV variable.; 
		%goto errtrap;
%end;


* Append this new information into the master data set;
proc append base=&jobname._master  data=&job2._newbits;
run;
%if &syserr %then %goto errtrap;

%if &debug=0 %then %do;
proc datasets library=&job1 nolist;
	delete &job2._temp21 &job2._temp22  &job2._temp39 &job2._temp40 &job2._temp41 &job2._temp42 &job2._temp51 
	&job2._newbits;
	%if &syserr %then %goto myabort; 
	quit;
	%if &syserr %then %goto myabort; 
%end;

proc sql;
update &jobname._master
set Vlabel="&jobname @"||put(datetime(),datetime.) where Role="Master2A"; 
delete from &jobname._master where Role in("Master2B","Master3");
%if &sqlrc %then %goto myabort;
quit;
%if &syserr %then %goto myabort; 


%put ################################;
%put Macro &macname ended succesfully;
%if &warnflag %then %do;
	%put WARNING: ... but see above for WARNINGS.;
%end;
%put ################################;
%goto myend;

%myabort:
%let errtext1=Untrapped error [Report to James Roger];

%errtrap:
	* Additional quit and run to flush any existing data step or proc;
	quit;
	run;
	* Now put error message and stop macro;	
	%PUT %str(ERR)OR: ; 
	%PUT %str(ERR)OR: In macro &macname with Jobname= &jobname; 
	%PUT %str(ERR)OR: &errtext1; 
	%PUT %str(ERR)OR: &errtext2; 
	%PUT %str(ERR)OR: ;
%myend:
%mend Part2A;

************************ End of macro Part2A *********************************************;
