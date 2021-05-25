/*
What does mimix do when an observation is wholly missing?
When ref = that arm: 
	LMCF = MAR at time 1 only
	J2R, CIR and CR all = MAR at all times
When ref = other arm: 
	MAR and LMCF are same as for "ref = that arm"
	J2R = CIR at all times but not equal to MAR
	CR different from everything else 
We can use this to test the R code
IW updated 29mar2021
Takes individual 5001 in placebo arm and sets to wholly missing
*/
label list treat1
foreach ref in 2 3 {
	di as input "Ref = `ref'"
	foreach method in mar lmcf j2r  cir cr {
		di as input "Method = `method'" _c
		use c:\stata\asthma, clear
		qui replace fev=. if id==5001
		qui mimix fev treat, id(id) time(time) method(`method') ref(`ref') covariates(base) clear m(1) seed(101)
		qui mi reshape wide fev, i(id) j(time)
		l if id==5001
	}
}
