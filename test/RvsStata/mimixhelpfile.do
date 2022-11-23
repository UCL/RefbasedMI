/* 
mimixhelpfile.do
Run the mimix help file examples subject to modifications:
	removed the saving option
	increased M from 50 to 500
	added CC analysis
IW 23/11/2022
*/

cd C:\ado\ian\RefBasedMI

* Start with CC to check correct dataset loading
use "data\asthma.dta", clear
reg fev i.treat base if time==12

*    Multiple imputation assuming the response variable fev is MAR
use "data\asthma.dta", clear
mimix fev treat, id(id) time(time) method(mar) covariates(base) clear m(500) seed(101)
mi estimate, mcerror: reg fev i.treat base if time==12

*    Multiple imputation and regression analysis assuming last mean carried forward for the response variable fev
use "data\asthma.dta", clear
mimix fev treat, id(id) time(time) method(lmcf) covariates(base) clear m(500) regress seed(101)
mi estimate, mcerror: reg fev i.treat base if time==12

*    Multiple imputation and regression analysis assuming jump to reference for the response variable fev, with placebo=2 as the reference
use "data\asthma.dta", clear
mimix fev treat, id(id) time(time) method(j2r) refgroup(2) covariates(base) clear m(500) regress seed(101)
mi estimate, mcerror: reg fev i.treat base if time==12

*    Saving the imputed dataset with filename mimix_example, assuming copy increments in reference for the response variable fev, with placebo=2 as the reference
use "data\asthma.dta", clear
mimix fev treat, id(id) time(time) method(cir) refgroup(2) covariates(base) clear m(500) seed(101)
mi estimate, mcerror: reg fev i.treat base if time==12
