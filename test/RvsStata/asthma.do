/* 
Ian to compare with R.do
mimix results for comparison with R 
Ian 8jul2020
updated 23nov2022: read R results in; add CC and MAR
*/

cd C:\ado\ian\RefBasedMI

use "data\asthma.dta", clear
replace treat=treat-1
replace treat=3 in 1/200
label val treat
tab treat
gen base2=base^2
replace fev=fev*1000
save c:\temp\asthmatemp, replace

reg fev i.treat base* if time==12

* compare with R
local bStata=_b[base2]
local seStata=_se[base2]
use test/RvsStata/R_CC, clear
di "Method CC (point est): Stata = " %9.4f `bStata' ", R = " %9.4f Estimate[5] 
di "Method CC (std error): Stata = " %9.4f `seStata' ", R = " %9.4f StdError[5] 
assert reldif(`=Estimate[5]',`bStata') < 1E-6
assert reldif(`=StdError[5]',`seStata')< 1E-6 

forvalues ref=1/2 {
	foreach method in MAR J2R CR CIR {
		di as input _new(3) "Method `method', ref `ref'"
		set seed 101
		use c:\temp\asthmatemp, clear
		mimix fev treat, id(id) time(time) method(`method') refgroup(`ref') covariates(base base2) ///
			clear m(100) seed(101)
		mi estimate, mcerror: reg fev i.treat base base2 if time==12
		mat b = e(b_mi)
		local bStata = b[1,"base2"]
		mat V = e(V_mi)
		local seStata = sqrt(V["base2","base2"])
		mat B = e(B_mi)
		local bStata_mcse = sqrt(B["base2","base2"]/e(M_mi))

		di as input "Compare with R"
		use test/RvsStata/R_`method'_ref`ref', clear
		l
		decode term, gen(termchar)
		assert termchar[5]=="base2"
		local bR = estimate[5]
		local diffz = (`bR' - `bStata') / `bStata_mcse' / sqrt(2) 
			// sqrt(2) assumes MC errors same in R and Stata
		di "Method `method', ref `ref': Stata = " %9.4f `bStata' ", R = " %9.4f `bR' _c
		di ", each MCSE = " %9.4f `bStata_mcse' ", z = " %9.2f `diffz'
	}
}
