/* 
twelvetimes_make.do
Stata file to create data with 12 times and 4 arms
to test RefBasedMI in R
IW 22nov2021
*/
set seed 51650
clear
set obs 4
gen group = _n
expand 200
sort group
by group: gen gpid = _n
gen u = rnormal()
gen v = rnormal()
by group: gen d = 12*runiform() if runiform()<0.5
expand 12
sort group gpid
by group gpid: gen time = _n
gen yvar = 10 + group + u + v*time + rnormal() if time<=d
drop u v d
gen myid=1000*group + gpid
drop if mi(yvar)
summ

outsheet myid group time yvar if gpid<=20 using twelvetimes20.csv, replace comma
outsheet myid group time yvar using twelvetimes200.csv, replace comma
