/*
check_treatvar_changetype.do

Changing treatvar from numeric to string changes results in Stata
but only if it changes the sort order

Below, string trtcharPA and numeric treat have the same sort order 
and give the same imputations

IW 22mar2021
*/


use "C:\ado\ian\test_Rmimix\asthma.dta", clear
tab treat
decode treat, gen(trtcharAP)
gen trtcharPA=string(treat)+". "+trtcharAP
label val treat
tab1 tr*
replace fev=fev*1000

mimix fev treat, id(id) time(time) method(J2R) refgroup(2) covariates(base) ///
	saving(z23,replace) m(2) seed(101)
mimix fev trtcharAP, id(id) time(time) method(J2R) refgroup("Placebo") covariates(base) ///
	saving(zAP,replace) m(2) seed(101)
mimix fev trtcharPA, id(id) time(time) method(J2R) refgroup("2. Placebo") covariates(base) ///
	saving(zPA,replace) m(2) seed(101)

use z23, clear
l if id==5017, sepby(_mi_m)
use zAP, clear
l if id==5017, sepby(_mi_m)
use zPA, clear
l if id==5017, sepby(_mi_m)
