/**********************************************************************************************
|
| Macro Name:  MIWithd
|
| Macro Version: 1
|
| SAS Version: 9.03
|
| Created By: James Roger
|
| Date: 06 October 2008
|
| Macro Purpose: Multiple imputation for Repeated Measures with specific
|		strategies for Missing data
|
| Macro Design:
|	Proc MI is used to generate multiple imputed means and VarCovar matrices
|		(from Posterior). This is done separately for each (treatment) group.
|	IML is used to simulate data to replace missing values based on selected
|		strategy. These are conditional MV Normal distributions conditioning
|		on the observed data and the covariates.
|	Proc MIANALZE is used to put these together and summarise result.
|	Srategies:
|		MAR  Use means and covariances for this group
|				[Equivalant to MMRM analysis] 
|		CC	 Use means and covariances for the control group instead (throughout) and its covariance.
|		CDC  Use actual group means up to withdrawal and then mean at last visit
|               plus differences in repeated means taken from control group and
|                covariancefrom control group. [If no data at all Start at control group mean.]
|		J2C  (Jump to Control) Use actual group means up to withdrawal and then jump to means
|				from control group after withdrawal. Use covariance from control group.
|		LMCF Use last previous observed group mean as mean, using same group for
|				both means and covariances.
|				[For LMCF, first mean is used throughout for those with no observations.]
|		Different strategies can be used by individual subjects (determined by a Variable in the data set)
|
| Input Parameters:
|
| NAME			DESCRIPTION																DEFAULT
| Data			Input data set name 
| Response		Name of the Response Variable
| Subject		Name of the variable identifying the Subject							Subject
| Time			Name of the variable identify time										Time
| Group			Name of the variable identifying the (treatment) Group
|				   This may be either numeric of Character variable.
| Covariates	List of any Covariates
|				 (which must be quantitative, that is numeric)
|				 (also covariates should be constant within a subject
|				   If not the maximum value across times will be used throughout)
|				   For classifcatory covariate try generating dummy variables.
| Nimpute		Number of imputations required											5
|				 (at least 50 is recommended)
| MISeed		Seed for random number generator in Proc MI								0
| IMLSeed		Seed for random number generator in IML									0
|				 (0 uses date-time)
| Method		Name of character Variable in data set determining the Method 
|				   for this subject (see above),
|					Possible values for variable are MAR, LMCF, CDC, CC					Method
| Control		Variable in the data set determining the level chosen for the control
|					This is only used when the Method variable has the value CC or CDC	[None]
| Out			Name of data set for estimated contrasts								MIwithd_pe
| OutData		Name of data set of imputed values
|				(When this parameter is used no analysis is done)	
| NBIter		The no. of iterations in the MCMC Burn in								1000
|				 (The default is 1000 which is much more than the Proc MI default of 200)
|				 Also we use the INITIAL=EM option to improve burnin.
| NIter			The no. of iterations between pulls from the posterior in MCMC			500				
|				 (The default is 500 which is much more than the Proc MI default of 100)
|				
| Output:
|
|
| Global macro variables created: None
|
| Macros called: None
|
| Example:
|
|	%MIWithd(data=psimv, response=madrs, subject=subject, time=week, Group=Trt,
|		Covariates=Baseline, Nimpute=50 );
|
|**********************************************************************************************
| Change Log
| 2008 Jan 04	James Roger	: Add code to remove records with missing covariate information.	
| 2008 Jan 04	James Roger	: Add check and warning for character covariates (not allowed).
| 2008 Jan 04	James Roger	: Remove bug using Control= with MAR. 
| 2008 Jan 24	James Roger : Add parameter NBIter and NIter.
| 2008 Jan 30   James Roger : Fix covariate means rather than sample
| 2008 Jan 30   James Roger : Separate Seed parameter into MISeed and IMLSeed
| 2008 Jan 31   James Roger : Allow control= parameter to be a list
| 2008 Feb 04   james Roger : Method and Control become variables in the data set. 
| 2008 may 01   James Roger : Fix bug when using OUTDATA
| 2008 Jun 05   James Roger : Add checks that variables exist in the data set (stops funny errors)
| 2008 Jun 06   James Roger : Change IML readng of data to make efficient (single scan)
| 2008 Sep 22   James Roger : Fix CC (IML_mc appeared as iml_)
| 2008 Oct 06   james Roger : Add J2C method and tidy diagnostics.
| 2009 Apr 26   James Roger : Fix bug when there are no covariates. (Version 23)
| 2009 Aug 10   james Roger : Change LOCF to LMCF. (Version 24).
| 2009 Aug 10   james Roger : New covariance matrix for for CDC and J2C (Version 24).
| 2009 Aug 17   james Roger : Add ERF to Proc analyze based ond DF from MIXED (Version 25) 
| 2009 Aug 19   James Roger : Change means for covariate to random rather than fixed.
| 2009 sep 02   James Roger : Add error message if all records in data are discarded.
| 2009 sep 02   James Roger : Add error message if IML fails (e.g. covariates are linearly related).
| 2010 aug 11   james Roger : Fix error trap when all data missing and no covariates exist (lnonmiss=0)
| 2010 oct 12   James Roger : Remove some redundent code. (MIWITHD30)
| 2010 feb 03   james Roger : Add formatted value to output label
| 2011 feb 03   james Roger : Add lsmeans as well as differences as default output. (MIWITHD31)
| 2011 nov 14   james Roger : Add df t value to default output. (MIWITHD31)
| 2013 feb 26   james roger : Fix error in selection of Sigma matrix for methods J2C and CDC 
                                 (used Active when should use reference).
| 2013 may 24   james Roger : Change ods listing close; to ods select none; etc.
**********************************************************************************************/

options nocenter mprint nomlogic;
%* options mprint mlogic symbolgen;
%* options nomprint nomlogic nosymbolgen;


%macro MIWithd(Data=, Response=, Subject=Subject, Time=Time ,Group=Group ,Covariates=, Nimpute=5,
	MISeed=0, IMLSeed=0, Method=Method, Control= ,Out=MIwithd_pe ,Outdata= ,NBIter=1000 ,NITER=500,
	imljrbug=0, macjrbug=0);

	%local i nlevels ntimes txt txt1 txt2 qcov ncov sqlcov sqlcov2 Clabel Clevel flag covflag flag dsid rc maxdf
			nrec1 nrec2  maxlsmdf;
	%let covflag=0;

	%* imljrbug macjrbug are Flags that turn on intermediate printing for diagnosing problems;

	%* nlevels is the number of levels for the group variable;
	%* ntimes is the number of levels for the time variable;
	%* ncov is the number of covariate variables;

	%* Build three lists of covariates to use in code later;
	%* qcov is list separated by spaces with variable names in quotes for reading variables in iml;
	%* sqlcov is list separated by commas with max(...) to pull covariates to higher stratum level
		using max( ) as the summarizing function;
	%* sqlcov2 is list separated by commas with B. added to pull covariates from second data set in
		the join (alias is B);
	%let i=1;
	%let qcov=;
	%let sqlcov=;
	%let sqlcov2=;
	%let cleancov= 1 ;
	%do %while(%length(%scan(&Covariates,&i)));
		%let txt=%scan(&Covariates,&i);
		%let qcov= &qcov "&txt";
		%let sqlcov= &sqlcov ,max( &txt ) as &txt;
		%let sqlcov2= &sqlcov2 ,B.&txt;
		%let cleancov= &cleancov & &txt > .z;
	%let i=%eval(&i +1);
	%end;
	%let ncov=%eval(&i -1);
	%if &macjrbug %then %put NCov= &Ncov;
	%if &macjrbug %then  %put sqlcov= &sqlcov;

	* Deleting working data sets before we start;
	proc datasets lib=work nolist;
	delete MIWithd_temp0 MIWithd_temp1 MIWithd_temp2 MIWithd_temp3 MIWithd_temp4 MIWithd_temp5
		MIWithd_temp6 MIWithd_temp7 MIWithd_temp8 MIWithd_temp9 MIWithd_temp10
		MIWithd_levels MIWithd_times
		MIWithd_Par1 MIWithd_par2 MIWithd_par3 MIWithd_par4
		MIWithd_outest1
		MIWithd_simdata MIWithd_pe MIWithd_finaldata MIWithd_diffs MIWithd_diffs2 MIWithd_lsm;
	quit;
 
	%* Check that data set exsts;
	%let dsid = %sysfunc(open(&data));
	%if &dsid %then %do;
		%* Check that required variables exist;
		%let flag=0;
		%if (%sysfunc(varnum(&dsid,&response)) = 0) OR
			(%sysfunc(varnum(&dsid,&subject)) = 0) OR
			(%sysfunc(varnum(&dsid,&time)) = 0) OR
			(%sysfunc(varnum(&dsid,&group)) = 0) OR
			(%sysfunc(varnum(&dsid,&method)) = 0)
		%then %let flag=1;
		%if %length(&Control) %then %if %sysfunc(varnum(&dsid,&control)) = 0 %then %let flag=1;
		%let i=1;
		%let txt=%scan(&Covariates,&i);
		%do %while(%length(&txt));
			%if %sysfunc(varnum(&dsid,&txt)) = 0 %then %let flag=1;
			%let i=%eval(&i +1);
			%let txt=%scan(&Covariates,&i);
		%end;
		%* Close it quickly;
		%let rc=%sysfunc(close(&dsid));
		%if &flag %then %do;
			%put %str(ERR)OR:;
			%put %str(ERR)OR: ==============================================================================;;
			%put %str(ERR)OR: ===   MIWithd macro (Sensitivity analysis using MI for NMAR assumptions)   ===;
			%put %str(ERR)OR: The input data set [ &Data ] does not contain all of these variables;
			%put %str(ERR)OR:    Response: &Response;
			%put %str(ERR)OR:     Subject: &Subject;
			%put %str(ERR)OR:        Time: &Time;
			%put %str(ERR)OR:       Group: &Group;
			%put %str(ERR)OR:      Method: &Method;
			%if %length(&Control) %then %put %str(ERR)OR:     Control: &Control;
			%let i=1;
			%let txt=%scan(&Covariates,&i);
			%do %while(%length(&txt));
				%put %str(ERR)OR: Covariate &i.: &txt;
				%let i=%eval(&i +1);
				%let txt=%scan(&Covariates,&i);
			%end;
			%put %str(ERR)OR:;
			%let flag=1;
			%goto errlabel;
		%end;
	%end;
	%else %do;
		%put %str(ERR)OR:;
		%put %str(ERR)OR: ==============================================================================;;
		%put %str(ERR)OR: ===   MIWithd macro (Sensitivity analysis using MI for NMAR assumptions)   ===;
		%put %str(ERR)OR: The input data set [ &Data ] does not exist;
		%put %str(ERR)OR:;
		%let flag=1;
		%goto errlabel;
	 %end;

	* Copy into local data set and sort data;
	proc sort data=&data out=MIWithd_temp0;
	where &cleancov;
	by &subject &time;
	run;
	%if &syserr = 3000 %then %do;
		%put %str(ERR)OR:;
		%put %str(ERR)OR: ==============================================================================;;
		%put %str(ERR)OR: ===   MIWithd macro (Sensitivity analysis using MI for NMAR assumptions)   ===;
		%put %str(ERR)OR: The previous error concerning "Where clause operator requires compatible variables.";
		%put %str(ERR)OR:    was generated because you have specified a Character variable in the Covariate list;
		%put %str(ERR)OR: All covariates must be numeric. You need to generate dummy variables for factor covariates;
		%put %str(ERR)OR:;
		%let flag=1;
		%goto errlabel;
	%end;

	* Add Warning if Data records were removed because of missing covariate values;
	data _NULL_;
	set &data NOBS=Nrec1;
	call symput("Nrec1",left(put(Nrec1,8.0)));
	stop;
	run;
	data _NULL_;
	set MIWithd_temp0 NOBS=Nrec2;
	call symput("Nrec2",left(put(Nrec2,8.0)));
	stop;
	run;
	%if %length(&nrec2) = 0 %then %do;
		%put %str(ERR)OR:;
		%put %str(ERR)OR:: There are no records selected from the data set;
		%put %str(ERR)OR:: Macro MIWITHD is aborted;
		%put %str(ERR)OR::;
		%let flag=1;
		%goto errlabel;
	%end;

	%if %eval(&Nrec1 - &Nrec2) %then %do;
		%put %str(WARN)ING:;
		%put %str(WARN)ING: %trim(&covflag) records were removed because of missing covariate information;
		%put %str(WARN)ING:;
	%end;

	* NODUPKEY removes any dupliate times for a subject;
	* Allows us to test for duplicate records in data;
	proc sort data=MIWithd_temp0 out=MIWithd_temp1 NODUPKEY;
	by &subject &time;
	run;

	* Stop if there were duplicates within a subject;
	data _NULL_;
	set MIWithd_temp0 NOBS=Nrec1;
	set MIWithd_temp1 NOBS=Nrec2;
	call symput("flag",put(Nrec1 ^= Nrec2,1.0));
	stop;
	run;
	%if &flag %then %do;
		%put %str(ERR)OR:;
		%put %str(ERR)OR:: There are duplicate times for the same level of subject;
		%put %str(ERR)OR:: Macro MIWITHD is aborted;
		%put %str(ERR)OR::;
		%let flag=1;
		%goto errlabel;
	%end;

	* Check if method is legal;
	%let flag=;
	data _NULL_;
	set MIWithd_temp1;
		%if %length(&control) %then %do; 
			if &method notin("MAR","CC","J2C","CDC","LMCF") then do;
			put "ERROR:";
			put "ERROR:: There are illegal values in the Method variable (&Method).";
			put "ERROR:: Value is " &method;  
			put "ERROR::    should be one of MAR, CC, J2C, CDC or LMCF";  
		%end;
		%* if control parameter is not set then only MAR and LMCF are allowed;
		%else %do;
			if &method notin("MAR", "LMCF") then do;
			put "ERROR:";
			put "ERROR:: There are illegal values in the Method variable (&Method).";
			put "ERROR:: Value is " &method " and no Control variable is specified";  
			put "ERROR::    should be one of MAR, or LMCF";  
		%end;

		put "ERROR:: Macro MIWITHD is aborted";
		put "ERROR::";
		call symput("flag",method);
		stop;
	end;
	run;
	%if %length(&flag) %then %do;
		%let flag=1;
		%goto errlabel;
	%end;

	proc sql noprint;
	* get distinct levels for groups;
	create table MIWithd_levels as
	select distinct(&group) as &group
	from MIWithd_temp1
	order by &group;

	* get distinct levels for time;
	create table MIWithd_times as
	select distinct(&time) as &time
	from MIWithd_temp1
	order by &time;

	select count(*) into :nlevels
	from MIWithd_levels;

	select count(*) into :ntimes
	from MIWithd_times;
	quit;

	* Trim spaces off the values;
	%let nlevels=&nlevels;
	%let ntimes=&ntimes;
 
	%if &macjrbug %then  %put &nlevels;
	%if &macjrbug %then  %put &ntimes;

	* Generate names for the response variables when we go parallel;
	data MIWithd_temp2;
	set MIWithd_times;
	My_time="Time_"||left(trim(put(_N_,3.0)));
	run;

	* ... and index for the group levels;
	%let Clevel= ;

	data MIWithd_temp3;
	set MIWithd_levels;
	My_group=_N_;
	run;

	proc sql noprint;
	* Now merge these back into the original data;
	create table MIWithd_temp4 as
	select MIWithd_temp1.*, My_time
	from MIWithd_temp1 A left join MIWithd_temp2 B
	on A.&time=B.&time;

	create table MIWithd_temp5 as
	select MIWithd_temp4.*, My_group 
	from MIWithd_temp4 A left join MIWithd_temp3 B
	on A.&Group=B.&Group;

	%if %length(&Control) %then %do;
		create table MIWithd_temp6 as
		select A.*, B.My_group as My_control 
		from MIWithd_temp5 A left join MIWithd_temp3 B
		on A.&Control = B.&Group
		order by &subject, My_time ;
	%end;
	%else %do;
		create table MIWithd_temp6 as
		select *, 0 as My_control 
		from MIWithd_temp5 A
		order by &subject, My_time ;
	%end;

	select min(My_control) into :flag
	from MIWithd_temp6;
	quit;

	%* trim flag;
	%let flag=&flag;
	%if .&flag = .. or .&flag = . %then %do;
		%put %str(ERR)OR:;
		%put %str(ERR)OR:: Control variable [&Control] takes a value that is not a Group value;
		%put %str(ERR)OR:: Macro MIWITHD is aborted;
		%put %str(ERR)OR::;
		%let flag=1;
		%goto errlabel;
	%end;

	* Make the data horizontal using the variable names from My_Time;
	proc transpose data=MIWithd_temp6 out= MIWithd_par1;
	var &response;
	id My_time ;
	by &Subject;
	run;

	%if &Ncov %then %do;
		* Use Max( ) to summarise data for covariates to higher stratum level,
		   Note this takes largest value from any record for a subject;
		proc sql;
		create table MIWithd_temp7 as
		select &subject, My_group, &method, my_control, &group &sqlcov
		from MIWithd_temp6
		group by &Subject, My_group, &method, my_control, &group;
		quit;
	%end;
	%else %do;
	%* As there are no covariayes use max(Mygroup) to force the GROUP BY to work;
		proc sql;
		create table MIWithd_temp7 as
		select &subject, max(My_group) as My_group, &method, my_control, &group
		from MIWithd_temp6
		group by &Subject, My_group, &method, my_control, &group;
		quit;
	%end;
	* Merge back in the rest of the variables;
	data MIWithd_par2;
	merge MIWithd_temp7 MIWithd_par1;
	by &subject;
	drop i k ;
	array columns[&ntimes] Time_1-Time_&Ntimes;
	*Encode pattern of missing data;
	k=1;
	pattern=0;
	do i=1 to &ntimes;
		if columns[i]<=.z then pattern=pattern+k;
		k=k*2;
	end;
	run;

	* Sort ready to do MI by group;
	proc sort data=MIWithd_par2 OUT=MIWithd_PAR4;
	by my_group &method my_control pattern &subject;
	run;

	* Now run the MI step and extract the simulated means and Var-Covar matrices;
	proc mi data= MIWithd_par4 nimpute=&nimpute %if &miseed ^= 0 %then %do; seed=&miseed %end; ;
	var Time_1-Time_&ntimes &covariates;
	mcmc OUTEST=MIWithd_outest1 NBITER=&NBITER NITER=&NITER INITIAL=EM;
	by my_group;


	* Build a list of the patterns for each group;
	* Count is used to control reading of data set;
	proc sql;
	create table MIWithd_par3 as
	select my_group, &method, my_control, pattern as idpat, Count(*) as Count
	from MIWithd_par2
	group by my_group, &method, my_control, pattern;
	run;

	*This is the main IML section where we build the imputed data;
	proc iml symsize=2097024;
		*Set control variables;
		ntimes=&Ntimes;
		nlevels=&nlevels;
		ncov=&Ncov;
		seed=&imlSeed;
		nimpute=&nimpute;
		* NCT is Length of the data including covariates;
		nct=ntimes+ncov;
* ===>;	%if &imljrbug %then print nct ntimes ncov;;
		Cols = concat("Time_",left(putN((1:ntimes),"2.0")));

		* C holds a list of the columns in the data that we need to read;
		%if %length(&qcov) %then %do;
			c=Cols||{&qcov};
		%end;
		%else %do;
			c=Cols;
		%end;

		* Set up control of output data set;
		ncols=c||"&Subject"||"Imputation"||"My_group";
		newy=j(1,nct+3,0);
		create MIWithd_simdata from newy [colname=nCols];

		* Read imputed means and VarCovars int IML matrices;
		VC=SHAPE(CONCAT('VC',LEFT(CHAR(1:(&Nimpute#&Nlevels)))),&Nimpute,&Nlevels);
		mean=SHAPE(CONCAT('Mean',LEFT(CHAR(1:(&Nimpute#&Nlevels)))),&Nimpute,&Nlevels);

		* Read in the imputed means and covariance matrices;
		USE MIWithd_outest1;
		do i=1 to NImpute;
			do j=1 to nlevels;
				read all var C into means where (_Imputation_=i & my_group=j & _Type_="PARM");
* ===>;			%if &imljrbug %then print means ;;

				newmeans=means[,1:nct];
* ===>;			%if &imljrbug %then print newmeans;;

				call valset(mean[i,j],newmeans);
				read all var c into Sigma where (_Imputation_=i & my_group=j & _Type_="COV" );
* ===>;			%if &imljrbug  %then print Sigma;;
				call valset(VC[i,j],sigma);
			end;
		end;
		%if &imljrbug  %then print "Closing Outest1" ;;
		CLOSE MIWithd_outest1;

		USE MIWithd_par3;
		read all var _num_ into nrdgroup;
		read all var _char_ into crdgroup;
* ===>;	%if &imljrbug  %then print nrdgroup crdgroup ;;
		CLOSE MIWithd_par3;
		* Ready to Read the raw data values;
		USE MIWithd_par4; 
		* Loop through the groups;
		rstart=1;
		do iml_rec=1 to nrow(nrdgroup);
			iml_mg=nrdgroup[iml_rec,1];
			iml_meth=crdgroup[iml_rec,1];
			if ncol(nrdgroup)=4 then do;
				iml_mc=nrdgroup[iml_rec,2];
				pat=nrdgroup[iml_rec,3];
				count=nrdgroup[iml_rec,4];
			end;
			else do;
				iml_mc=crdgroup[iml_rec,2];
				pat=nrdgroup[iml_rec,2];
				count=nrdgroup[iml_rec,3];
			end;
* ===>;		%if &imljrbug %then print iml_mg iml_meth iml_mc pat count;;

			* Decode pattern;
			flag=j(ntimes,1,0);
			k=2**ntimes;
			mypat=pat;
			do j= ntimes to 1 by -1;
				k=k/2;	
				if mypat >= k then do;
					mypat=mypat-k;
					flag[j]=1;
				end;
			end;
	
            %* Miss identifies columns in data which are missing;
            %* NonMiss identies columns in data which are not missing or are covariates;
            miss=loc(flag);
            if ncov>0 then nonmiss=loc(1-flag)||(ntimes+(1:ncov));
			else nonmiss=loc(1-flag);
			lnonmiss=nleng(nonmiss);

			* Read all var cols into Y where (pattern=pat & my_group=iml_mg & &method=iml_meth);
			* Note 1:Ntimes is raw data and remaining Ncov are covariates;
			CC=c||{"&Subject" "My_Control"};
			rend=rstart+count-1;
			read point (rstart:rend) var CC into Y;
			rstart=rend+1;
			* Make a matrix ready to poke results into;
			newy=j(nrow(Y),3+nct,.);

			* Now for each imputation generate imputed data; 
			do imputation=1 to nimpute;
				* print imputation;
				if pat=0 then do;
					*No missing data, so copy data into matrix ready for output;
					newy[,1:nct]=Y[,1:nct];
				end;
				else do;
					*Load values for conditional distribution;;
					if iml_meth = "MAR" then do;
						means=repeat(value(mean[Imputation,iml_mg]),nrow(Y),1);
						Sigma=value(VC[Imputation,iml_mg]);
						if lnonmiss then do;
							S11=sigma[nonmiss,nonmiss];
							S12=sigma[nonmiss,miss];
							S22=sigma[miss,miss];
						end;
					end;
					else if iml_meth = "CC" then do;
						means=repeat(value(mean[Imputation,iml_mc]),nrow(Y),1);
						Sigma=value(VC[Imputation,iml_mc]);
						if lnonmiss then do;
							S11=sigma[nonmiss,nonmiss];
							S12=sigma[nonmiss,miss];
							S22=sigma[miss,miss];
						end;
					end;
					else if iml_meth = "J2C" then do;
*						means=repeat(value(mean[Imputation,iml_mc]),nrow(Y),1);
						means=repeat(value(mean[Imputation,iml_mg]),nrow(Y),1);
						meansC=repeat(value(mean[Imputation,iml_mc]),nrow(Y),1);
* ===>;					%if &imljrbug %then print means meansc;;
						cc=miss`;
						call sort(cc,1);
						cc=cc`;
						do i=1 to ncol(cc);
							means[,cc[i]]=meansC[,cc[i]];
						end;
						Sigma=value(VC[Imputation,iml_mc]);
						if lnonmiss then do;
							S11=sigma[nonmiss,nonmiss];
							S12=sigma[nonmiss,miss];
							S22=sigma[miss,miss];
						end;
					end;
					else if iml_meth = "CDC" then do;
						means=repeat(value(mean[Imputation,iml_mg]),nrow(Y),1);
						meansC=repeat(value(mean[Imputation,iml_mc]),nrow(Y),1);
* ===>;					%if &imljrbug %then print means meansc;;
						cc=miss`;
						call sort(cc,1);
						cc=cc`;
						* Use do loop do make sure the values carry forward properly;
						do i=1 to ncol(cc);
							*Condition is that if all are missing and this is start then use start for control;
							*otherwise just step forward;
							if cc[i] = 1 then means[,1]=meansC[,1];
							else means[,cc[i]]=means[,cc[i]-1]+meansC[,cc[i]]-meansC[,cc[i]-1];
						end;
						Sigma=value(VC[Imputation,iml_mc]);
						if lnonmiss then do;
							S11=sigma[nonmiss,nonmiss];
							S12=sigma[nonmiss,miss];
							S22=sigma[miss,miss];
						end;
					end;
					else if iml_meth = "LMCF" then do;
						means=repeat(value(mean[Imputation,iml_mg]),nrow(Y),1);
						cc=miss`;
						call sort(cc,1);
						cc=cc`;
						* Use do loop do make sure the values carry forward properly;
						* If first time is missing then use first mean for this arm;
						do i=1 to ncol(cc);
							if cc[i] > 1 then means[,cc[i]]=means[,cc[i]-1];
						end;
						Sigma=value(VC[Imputation,iml_mg]);
						if lnonmiss then do;
							S11=sigma[nonmiss,nonmiss];
							S12=sigma[nonmiss,miss];
							S22=sigma[miss,miss];
						end;
					end;
* ===>;				%if &imljrbug %then print means Sigma;;
					if lnonmiss=0 then do;
						*Need to simulate complete data;
						U=root(sigma);
					    newy[,1:ntimes]=means[,1:ntimes]+rannor(repeat(seed,nrow(Y),ntimes))*U;
					end;
					else do;
					*Simulate missing data only;
						m1=means[,nonmiss];
						m2=means[,miss];
						Raw1=Y[,nonmiss];
* ===>;					%if &imljrbug %then print Y nonmiss raw1;;
* ===>;					%if &imljrbug %then print m1 m2 s11 s12 s22;;
						* Conditional distribution;
						t=solve(s11,s12);
* ===>;					%if &imljrbug %then print t;;
						Conds=s22-(s12`)*t;
						*print conds;
						meanval=m2 + (raw1 - m1)*t;
* ===>;					%if &imljrbug %then print meanval;;
						U=root(Conds);
					    Y1=meanval+rannor(repeat(seed,nrow(Y),ncol(miss)))*U;
* ===>;					%if &imljrbug %then print y1;;
						*Copy into matrix ready for output;
						newy[,nonmiss]=raw1;
						newy[,miss]=y1;
					end;
				end;

				* Fill out rest of matrix;
				newy[,1+nct]=Y[,nct+1];
				newy[,2+nct]=repeat(imputation,nrow(Y),1);
				newy[,3+nct]=repeat(iml_mg,nrow(Y),1);
* ===>;			%if &imljrbug %then print y;;
* ===>;			%if &imljrbug %then print newy;;
				append from newy ;
			end;
		end;
		CLOSE MIWithd_par4;
	quit;
	%if &syserr = 1012 %then %do;
		%put %str(ERR)OR:;
		%put %str(ERR)OR:: There is a problem in the IML code;
		%put %str(ERR)OR:: Matrix being inverted should be singular;
		%put %str(ERR)OR:: This is most likely because the covariates are linearly related;
		%put %str(ERR)OR:: Have you accidentally included dummy variables for all possible levels?;
		%put %str(ERR)OR:;
		%put %str(ERR)OR:: Macro MIWITHD is aborted;
		%put %str(ERR)OR::;
		%let flag=1;
		%goto errlabel;
	%end;
	%if &syserr %then %do;
		%put %str(ERR)OR:;
		%put %str(ERR)OR:: There is an unkown problem in the IML code;
		%put %str(ERR)OR:: Please Report to James Roger;
		%put %str(ERR)OR:;
		%put %str(ERR)OR:: Macro MIWITHD is aborted;
		%put %str(ERR)OR::;
		%let flag=1;
		%goto errlabel;
	%end;

	* Merge imputed data back onto original data set;
	* Note we could calculate any function of data here for univariate analysis;
	* We use last time point to define RESPONSE;
	%if %length(&Outdata) %then %do;

		proc sort data=MIWithd_simdata;
		by imputation &subject ;
		run;

		proc transpose data=MIWithd_simdata out=  MIWithd_temp8 prefix=Response;
		var time_1-time_&ntimes;
		by imputation &subject;
		run;

		proc sql;
		* Select data from last observation on subject as there covariate data;
		create table MIWithd_temp9 as
		select A.*
		from (select &subject , max(my_time) as my_time
			from MIWithd_temp4 group by &Subject) B
			left join MIWithd_temp4(drop=&time) A
		on A.my_time=B.my_time and A.&subject = B.&subject;

		*Get actual time data for new records;
		create table MIWithd_temp10 as
		select A.*, B.&time
		from MIWithd_temp8 A left join MIWithd_temp2 B
		on A._name_ = B.My_Time; 

		* Merge it all together;
		create table &Outdata as
		select A.Imputation, A.Response1 as &Response, A.&time, B.* 
		from  MIWithd_temp10 A left join miwithd_temp9(drop=my_time &Response) B
		on A.&subject = B.&subject
		order by imputation, &subject, &time;
		quit;

		%put NOTE:;
		%put NOTE: Macro MIWITHD. Imputed data created in variable Response in data set: &Outdata.;
		%put NOTE: No analysis carried out.;
		%put NOTE:;

	%end;
	%else %do;
		proc sql;
			create table MIWithd_finaldata as
			select A.&Subject, A.Imputation, A.My_group, A.Time_&ntimes as Response, B.&group &sqlcov2
			from MIWithd_simdata A left join MIWithd_temp7 B
			on A.&subject=B.&subject
			order by imputation, my_group, &subject;
		quit;


		ods select none;
		ods results off;

		* Run proc mixed on imputed data;
		proc mixed data=MIWithd_finaldata;
		class  &subject &group;
		model response = &group &covariates /ddfm=KR;
		lsmeans &group /diffs;
		repeated intercept /subject=&subject group=&group;
		by imputation;
		ods output diffs=miwithd_diffs;
		ods output lsmeans=miwithd_lsm;
		run;
		ods output clear;
		ods select all;
		ods results on;
	
		proc sql noprint;
		select max(df) into :maxdf  
		from miwithd_diffs;
		select max(df) into :maxlsmdf
		from miwithd_lsm;
		quit;

		* Generate label for proc Mianalyze;
		data miwithd_diffs2;
		set miwithd_lsm(in=inlsm) miwithd_diffs ;
		Length Response $ 16 Method $ 16 Label $ 32;
		Response="&Response";
		Method="&Method";
		if inlsm then do;
			if Vtype(&group)="N" and length(vformatn(&group))>0 then do;
				label=Effect||":: "||trim(left(putn(&group,vformatn(&group))));
			end;
			else do;
				label=Effect||":: "||trim(left(&group));
			end;
		end;
		else do;
			if Vtype(&group)="N" and Vtype(_&group)="N" and length(vformatn(&group))>0 
				and length(vformatn(_&group))>0 then do;
				label=Effect||": "||trim(left(putn(&group,vformatn(&group))))||" - "||
					trim(left(putn(_&Group,vformatn(_&Group))));
			end;
			else do;
				label=Effect||": "||trim(left(&group))||" - "||trim(left(_&Group));
			end;
		end;
		run;

		proc sort data=miwithd_diffs2;
		by label;
		run;

		* Finally run MIANALYZE for differences;
		title1  "Principled Multiple Imputation: Data set= &Data,  Response Var.= &Response, Method Var.= &Method"
			%if %length(&control) %then ", Control Var.=&Control"; ;
		title2 "Covariates= &Covariates";
		proc mianalyze data=miwithd_diffs2 edf=&maxdf;
		modeleffects estimate;
		stderr stderr;
		by label response method;
		ods output parameterestimates=&out ;
		run;
		ods output clear;

		proc print data=&out noobs;
		var Method Response Label Estimate Stderr LCLMean UCLMean df tvalue Probt;
		run;

		title1 "";
	%end;

	%let flag=0;
%errlabel:


	%* Deleting working data sets after we finish for routine running;
	%if (not &imljrbug) and (not &macjrbug) %then %do;
		proc datasets lib=work nolist;
		delete MIWithd_temp0 MIWithd_temp1 MIWithd_temp2 MIWithd_temp3 MIWithd_temp4 MIWithd_temp5
			MIWithd_temp6 MIWithd_temp7 MIWithd_temp8  MIWithd_temp9 MIWithd_temp10
			MIWithd_levels MIWithd_times
			MIWithd_Par1 MIWithd_par2 MIWithd_par3 MIWithd_par4
			MIWithd_outest1
			MIWithd_simdata MIWithd_pe MIWithd_finaldata MIWithd_diffs MIWithd_diffs2 MIWithd_lsm;
		quit;
	%end;

	%if &covflag %then %do;
		%put %str(WARN)ING: ==============================================================================;
		%put %str(WARN)ING: ===   MIWithd macro (Sensitivity analysis using MI for NMAR assumptions)   ===;
		%put %str(WARN)ING: Note that %trim(&covflag) records have been removed because of missing covariate information;
		%put %str(WARN)ING:;
	%end;

	%if &flag %then %do;
		%put %str(ERR)OR:;
		%put %str(ERR)OR: ==============================================================================;;
		%put %str(ERR)OR: ===   MIWithd macro (Sensitivity analysis using MI for NMAR assumptions)   ===;
		%put %str(ERR)OR: Macro MIWITHD has aborted. See above for details.;
		%put %str(ERR)OR:;
	%end;

%mend MIWithd;
;
