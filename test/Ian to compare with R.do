/* 
Ian to compare with R.do
mimix results for comparison with R 
Ian 8jul2020
*/
forvalues ref=2/3 {
	foreach method in J2R CR CIR {
		di as input _new(3) "Method `method', ref `ref'"
		set seed 101
		use "C:\ado\ian\test_Rmimix\asthma.dta", clear
		replace treat=4 in 1/200
		tab treat, nol
		gen base2=base^2
		replace fev=fev*1000
		mimix fev treat, id(id) time(time) method(`method') refgroup(`ref') covariates(base base2) ///
			clear m(100) seed(101)
		mi estimate, mcerror: reg fev i.treat base base2 if time==12
	}
}
