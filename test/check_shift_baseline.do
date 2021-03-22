/* 
check_shift_baseline.do
verify that no method depends on actual value of baseline
IW 22mar2021
*/

local id 5099
pda
cd c:\temp
use c:\stata\asthma, clear
replace fev=. if id==`id'
local ref 2
foreach method in mar j2r cr cir lmcf {
	mimix fev treat, id(id) time(time) method(`method') ref(`ref') cov(base) saving(z`method',replace) m(2) seed(101)
}
replace base = base+10
foreach method in mar j2r cr cir lmcf {
	mimix fev treat, id(id) time(time) method(`method') ref(`ref') cov(base) saving(z`method'p,replace) m(2) seed(101)
}
foreach method in mar j2r cr cir lmcf {
	di "Method = `method'"
	use z`method', clear
	l if id==`id', sepby(_mi_m)
	use z`method'p, clear
	l if id==`id', sepby(_mi_m)
}
