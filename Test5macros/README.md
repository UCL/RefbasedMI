Discrepancy found in treatment effects when running 5macros v mimix   

SAS
%part1A(Jobname=inBasval,Data=hamdno3618csv,Subject=Patient_number,Response=change,Time=Visit_number, Treat=treatment_name,Covbytime= Basval,covgroup=treatment_name);

Stata
mimix change tmentno,id(patient_number) time(visit_number) covar(basval) method(J2R) ref(1) m(1000) seed(301)
