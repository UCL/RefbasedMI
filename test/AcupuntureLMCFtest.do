use "\\ad.ucl.ac.uk\home0\rmjlkm0\Downloads\acupuncture.dta", clear

mimix  head treat, id(id) time(time) covar(head_base) method(LMCF) refgroup(1)   m(1000) seed(54321) clear  



* error caused by interim option
mimix  head treat, id(id) time(time) covar(head_base) method(CIR) refgroup(1) interim(MAR)  m(10) seed(54321) clear 
mimix  head treat, id(id) time(time) covar(head_base) method(LMCF) refgroup(1)   m(1000) seed(54321) clear 
* create head3 head12

gen head3 = head if time ==3
gen head12 = head if time ==12

replace head3=head if time==3 
replace head12=head if time==12

* no head_base !!
mi estimate, mcerror: regress head12 treat  if time==12
mi estimate, mcerror: regress head12 treat head_base if time==12
* interims
preserve
keep if id==104 & _mi_m >0
univar head3 head12 , dec(3)
restore
preserve
keep if id==333 & _mi_m >0
univar head3 head12 , dec(3)
restore
preserve
keep if id==386 & _mi_m >0
univar head3 head12 , dec(3)
restore
preserve
keep if id==787 & _mi_m >0
univar head3 head12 , dec(3)
restore
preserve
keep if id==435 & _mi_m >0
univar head3 head12 , dec(3)
restore
preserve
keep if id==697 & _mi_m >0
univar head3 head12 , dec(3)
restore
preserve
keep if id==100 & _mi_m >0
univar head3 head12 , dec(3)
restore
preserve
keep if id==101 & _mi_m >0
univar head3 head12 , dec(3)
restore
log close


import delimited "C:\GitHub\mimix\mimixR\acupunturetime0.csv", numericcols(11) clear  
mimix  head treat, id(id) time(time)  method(LMCF) refgroup(1)   m(1000) seed(54321) clear 
mi estimate, mcerror: regress head12 treat  if time==12
mi estimate, mcerror: regress head12 treat head_base if time==12
* interims
preserve
keep if id==104 & _mi_m >0
univar head3 head12 , dec(3)
restore
preserve
keep if id==333 & _mi_m >0
univar head3 head12 , dec(3)
restore
preserve
keep if id==386 & _mi_m >0
univar head3 head12 , dec(3)
restore
preserve
keep if id==787 & _mi_m >0
univar head3 head12 , dec(3)
restore
preserve
keep if id==435 & _mi_m >0
univar head3 head12 , dec(3)
restore
preserve
keep if id==697 & _mi_m >0
univar head3 head12 , dec(3)
restore
preserve
keep if id==100 & _mi_m >0
univar head3 head12 , dec(3)
restore
preserve
keep if id==101 & _mi_m >0
univar head3 head12 , dec(3)
restore
log close


