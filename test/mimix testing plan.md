# Tests to do on mimix

IW drafted 18dec2020

Note: 

- I've written "similar" meaning "within MC error" for some comparisons

- I've written "exactly the same" for several comparisons of different options of the R package. But here "similar" is acceptable *if* we understand why.

- Stata = Stata mimix; SAS = SAS 5 macros

# Main methods - J2R, CIR, CR, MAR - after discontinuation

- show all methods give exactly the same imputed values in reference arm

- show similar results (imputed values and treatment effect) when used with a different seed 

- show similar results (imputed values and treatment effect) to Stata

- note: we think Stata may impute wrongly after discontinuation in individuals who also have interim missing values.


# Causal model

I suggest looking at just one case e.g. K0=1, K1=0.9

- show K0=0 corresponds to J2R

- show K0=K1=1 corresponds to CIR

- increment all observed outcomes at a particular time t in the non-reference arm, and show this causes all later imputed values in individuals who discontinue at t to increase by K0*K1^time, and does not affect other imputed values.

-- done by Ian

- show similar results (imputed values and treatment effect) to SAS 


# Interim missing, imputed as MAR

- show interim imputed values are exactly the same for all methods (J2R, CIR, CR, MAR, LMCF, causal) and all reference arms (active, placebo)

- show changing the interim method changes only the interim imputed values


# Delta method

- show that a simple specification of delta increases all post-discontinuation imputed values correctly and does not affect interim imputed values
