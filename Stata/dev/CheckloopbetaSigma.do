* km use mimix_all_stdata1
mata: mimix_all= st_data(., .)

local m=2
local nct =5
di "m= " `m', "nct = " `nct'
	forvalues k=1/`m' {						
		mata: mean_group`i'_imp`k' = mimix_all[`k',1..`nct']
		}
		
		
		mata: mata_VAR_group`i'_imp`k'=J(`nct',`nct',0)
		local step = `nct'+ 1
		forvalues r = 1/`nct' {				
			forvalues j = 1/`nct'{ 
				if `j' <= `r' {
					mata: mata_VAR_group`i'_imp`k'[`r', `j'] = mimix_all[`k', `step']
					local step = `step' + 1
				}
			}
		}
		mata: mata_VAR_group`i'_imp`k' = makesymmetric(mata_VAR_group`i'_imp`k')
	}	 
	/* 
	
/* allnew4 is equilv to fev2  , analyse id 5456  , allnew8 is id

use "C:\Users\rmjlkm0\Documents\mimix\mata_all_newimp2.dta"  ,clear

summary((allnew8)

*create subset 5456
keep if allnew8 == 5456
