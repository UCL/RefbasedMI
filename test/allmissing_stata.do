/*
What does mimix do when an observation is wholly missing?
When ref = that arm: all methods are the same, except LMCF
When ref = other arm: 
	MAR=same as when ref=that arm
	LMCF=MAR at first time only
	J2R=CIR at all times
	everything else completely different
My understanding:
	J2R = CR always and = MAR when ref=that arm
	CIR and LMCF are formally undefined
*/
foreach ref in 2 3 {
	di as input "Ref = `ref'"
	foreach method in mar j2r cr cir lmcf {
		di as input "Method = `method'" _c
		use c:\stata\asthma, clear
		qui replace fev=. if id==5001
		qui mimix fev treat, id(id) time(time) method(`method') ref(`ref') covariates(base) clear m(1) seed(101)
		l if id==5001 & _mi_m==1
	}
}
