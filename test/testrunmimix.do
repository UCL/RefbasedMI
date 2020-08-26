use "N:\Documents\GitHub\mimix\mimixR\asthma.dta", clear

use "N:\Documents\asthma.dta", clear
mimixOldv2 fev treat, id(id) time(time) method(cir) refgroup(2) covariates(base) clear m(2)  seed(101)

mimix fev treat, id(id) time(time) method(cir) refgroup(3) covariates(base) clear m(100)  seed(101)
mimixOldv2 fev treat, id(id) time(time) method(cir) refgroup(2) covariates(base) clear m(2)  seed(101)


preserve
keep if id==5099 & _mi_m >0
univar fev2 fev4 fev8 fev12 , dec(3)
restore
mi estimate , mcerror: regress fev treat  base  if time ==12  



mimix head treat,id(id) time(time) covariates(sex age head_base) method(j2r) refgroup(1) m(1000) seed(101) clear
* try covars diff order

mimix head treat,id(id) time(time) covariates(head_base  age sex) method(j2r) refgroup(1) m(1000) seed(101) clear


mi estimate ,mcerror: regress head sex age head_base if time==12



*keep if id==5333
use "\\ad.ucl.ac.uk\home0\rmjlkm0\Downloads\acupuncture.dta", clear
mimixOld  head treat, id(id) time(time) covar(head_base sex age) method(MAR)  m(2) seed(201) clear 

mimix  head treat, id(id) time(time) covar(head_base sex age) method(J2R) refgroup(0)  m(5) seed(201) clear 
mimix  head treat, id(id) time(time)  method(MAR) covar(head_base)  m(5) seed(201) clear  
 
 * try covars diff ord
 mimix  head treat, id(id) time(time) covar(age sex  head_base) method(J2R) refgroup(0)  m(5) seed(201) clear 
 
* create head3 head12

gen head3 = head if time ==3
gen head12 = head if time ==12

preserve
keep if id==100 & _mi_m>0
univar head3 head12 sex age head_base, dec(3)
 
restore
preserve
keep if id==101 & _mi_m>0
univar head3 head12 sex age head_base, dec(3)
restore
 
*

*replace fev2=.
gen fev2=.
gen fev4 =.
gen fev8=.
gen fev12=.
replace fev2=fev if time==2 
replace fev4=fev if time==4
replace fev8=fev if time==8 
replace fev12=fev if time==12   

preserve
keep if id==5099 & _mi_m >0
univar fev2 fev4 fev8 fev12 , dec(3)
restore
preserve
keep if id==5456 & _mi_m >0
univar fev2 fev4 fev8 fev12 , dec(3)
restore
preserve
keep if id==5127 & _mi_m >0
univar fev2 fev4 fev8 fev12 , dec(3)
restore
preserve
keep if id==5051 & _mi_m >0
univar fev2 fev4 fev8 fev12 , dec(3)
restore
preserve
keep if id==5074 & _mi_m >0
univar fev2 fev4 fev8 fev12 , dec(3)
restore
preserve
keep if id==5094 & _mi_m >0
univar fev2 fev4 fev8 fev12 , dec(3)
restore
preserve
keep if id==5333 & _mi_m >0
univar fev2 fev4 fev8 fev12 , dec(3)
restore





preserve
keep if id==100 & _mi_m>0
univar head3 head12 sex age head_base, dec(3)
 
 
 


univar head3 head12 sex age head_base, dec(3)

* test program to compare SAS 5macros
log using J2rCIRlog_refDRUG, text

import delimited "\\ad.ucl.ac.uk\home0\rmjlkm0\Documents\SASantidepmethodvar.csv", clear 
*ren treatmentname tmentno
* best recode drug/treatment otherwise drug ges recoded as 1 because 1st?
*drop tmentno
gen tmentno=.
replace tmentno =2 if (treatmentname=="PLACEBO") 
replace tmentno =1 if (treatmentname=="DRUG")
label define tmentlabel 2 "PLACEBO" 1 "DRUG"
label values tmentno tmentlabel

* testing for comparing SAS oututs Drug arm differs to R
*  M=2000  Causes error
mimix change tmentno,id(patientnumber) time(visitnumber)  method(J2R) ref(2) m(1000) seed(301)  regress clear
mimix change treatmentname,id(patientnumber) time(visitnumber) covar(basval  ) method(CIR) ref(DRUG) m(1000) seed(101)  regress clear
mimix change tmentno,id(patientnumber) time(visitnumber) covar(basval  ) method(J2R) ref(1) m(1000) seed(301)  regress clear

mi estimate , mcerror: regress change tmentno  if visitnumber ==7
mi estimate , mcerror: regress change treatmentname basval  if visitnumber ==7

mimix change tmentno,id(patientnumber) time(visitnumber) covar(basval  ) method(CIR) ref(2) m(500) seed(101)  clear
mi estimate , mcerror: regress change tmentno basval  if visitnumber ==7


mimix change tmentno,id(patientnumber) time(visitnumber) covar(basval ) methodvar(methodvar) refgroupvar(referencevar) m(2) seed(101) regress clear

* without methodvar
* create dummies from factor
tabulate pooledinvestigator, generate(pooledinvestigator)

mimix change tmentno,id(patientnumber) time(visitnumber) covar(basval pooledinvestigator1 pooledinvestigator2 pooledinvestigator3 pooledinvestigator4 pooledinvestigator5  ///
             pooledinvestigator6  ///
             pooledinvestigator7 pooledinvestigator8 pooledinvestigator9  pooledinvestigator10 pooledinvestigator11 pooledinvestigator12   pooledinvestigator13 ///
			 pooledinvestigator14               pooledinvestigator15 pooledinvestigator16  ) method(J2R) ref(2)  m(100) seed(101) regress clear
*mimixOld

gen sex =.
replace sex=1 if (patientsex == "F")
replace sex=2  if (patientsex == "M") 

log using J2rlog_DRUG, text
mimix change tmentno,id(patientnumber) time(visitnumber) covar(basval ) method(J2R) ref(1) m(200) seed(101) clear


mimix hamd17total tmentno,id(patientnumber) time(visitnumber) covar(basval  sex  ) method(MAR) ref(1) m(100) seed(101) clear
mi estimate, mcerror: regress hamd17total tmentno basval  sex  if visitnumber==7
mimix hamd17total tmentno,id(patientnumber) time(visitnumber) covar(basval  sex  ) method(LMCF) ref(1) m(100) seed(101) clear
mi estimate, mcerror: regress hamd17total tmentno basval  sex  if visitnumber==7
mimix hamd17total tmentno,id(patientnumber) time(visitnumber) covar(basval  sex  ) method(CIR) ref(1) m(100) seed(101) clear
mi estimate, mcerror: regress hamd17total tmentno basval  sex  if visitnumber==7

import delimited "\\ad.ucl.ac.uk\home0\rmjlkm0\Documents\SASantidep.csv", clear 

mimix change treatmentname,id(patientnumber) time(visitnumber) covar(basval  ) method(J2R) ref(DRUG) m(1000) seed(101)  clear
ren treatmentname tment
mimix change tment,id(patientnumber) time(visitnumber) covar(basval  pooledinvestigator) method(J2R) ref(DRUG) m(2) seed(101) mixed clear

#mimix change tment,id(patientnumber) time(visitnumber) covar(basval  pooledinvestigator) method(J2R) ref(DRUG) m(2) seed(101) mixed clear


mimix change treatmentname,id(patientnumber) time(visitnumber) covar(basval) method(J2R) ref(DRUG) m(2) seed(101) regress clear

gen change4=.
gen change5=.
gen change6=.
gen change7=.
replace change4 =change if visitnumber ==4
replace change5 =change if visitnumber ==5
replace change6 =change if visitnumber ==6
replace change7 =change if visitnumber ==7
preserve
keep if patientnumber == 3618 & _mi_m>0
univar change4 change5 change6 change7
restore
preserve
keep if patientnumber == 2721 & _mi_m>0
univar change4 change5 change6 change7
restore
preserve
keep if patientnumber == 2104 & _mi_m>0
univar change4 change5 change6 change7
restore
preserve
keep if patientnumber == 2230 & _mi_m>0
univar change4 change5 change6 change7
restore
preserve
keep if patientnumber == 1513 & _mi_m>0
univar change4 change5 change6 change7
restore

log off

gen hamd17total4=.
gen hamd17total5=.
gen hamd17total6=.
gen hamd17total7=.

replace hamd17total4 =hamd17total if visitnumber ==4
replace hamd17total5 =hamd17total if visitnumber ==5
replace hamd17total6 =hamd17total if visitnumber ==6
replace hamd17total7 =hamd17total if visitnumber ==7
preserve
keep if patientnumber == 2721 & _mi_m>0
univar hamd17total4 hamd17total5 hamd17total6 hamd17total7

preserve
keep if patientnumber == 2721 & _mi_m>0
preserve
keep if patientnumber == 2104 & _mi_m>0
univar hamd17total4 hamd17total5 hamd17total6 hamd17total7
restore
preserve
keep if patientnumber == 2230 & _mi_m>0
univar hamd17total4 hamd17total5 hamd17total6 hamd17total7
restore
preserve
keep if patientnumber == 1513 & _mi_m>0
univar hamd17total4 hamd17total5 hamd17total6 hamd17total7
restore




Figure 42 Figure 41 mimix hamd17total treatmentname, id(patientnumber) time(visitnumber) method(j2r) refgroup(DRUG) covariates(basval ) clear m(50) regress seed(101)