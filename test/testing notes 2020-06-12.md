# Rmimix: Ian's 2nd set of comments 
12jun2020

## output 

- drop the [1] and quotes be dropped in the outputs
- change X1 X1cum exid in "summary missing pattern" to something more meaningful as previously suggested
- at start of imputation, print "Starting imputation"

## mimix

- idvar shouldn't be changed to .id
- to reduce risk of errors, code should count how many missing values remain unimputed in data set e.g. using sum(is.na(subset(impJ2Rdelta,.imp>0))), & output with a warning if >0
- it's strange that method and refer are separate options (using fixed values) but methodvar specifies them both (using variables). I'd prefer to have a refervar option too. 
- change refer option name to ref or reference?
- after conversion to mids format, imputation methods are still pmm (we agreed to change to "mimix") and predictor matrix is still present (we agreed to change to N/A).

## delta

- works correctly for interim missings
- works correctly for simple treatment discontinuation 
- fails if delta is specified and dlag isn't. dlag should default to 1's.
- help should state the required length of delta and dlag vectors (=#time points)
- if delta or dlag is mis-specified, imputed values can be missing. OUCH!! Build in a check that delta method has not de-imputed a value?
- suggested @details help text for delta method: "Specifying delta and dlag allows imputations to differ sytematically from RBI methods. They provide an increment which is added on to all values imputed after treatment discontinuation, but not to interim missing values. Values of delta are cumulated after treatment discontinuation. For example, for an individual who  discontinued treatment at the 2nd time point, we take the vector of delta's starting at the 3rd time point and add their cumulative sums to the imputed values. Specifying dlag modifies this behaviour, so that the vector of delta's starting at the 3rd time point is multipled elementwise by the vector dlag. The formula for the increment at time k for an individual who discontinued after time p is b_1*a_{p+1} + b_2*a_{p+2} + ... + b_{k-p}*a_k where delta=(a_1,a_2,...) and dlag=(b_1,b_2,...). A common increment of 3 at all time points after treatment discontinuation is achieved by delta=c(3,3,3,...) and dlag=c(1,0,0,...)."
- the algorithm is wrong for individuals with both interim missing and discontinuation e.g. id=5051 pattern=XOXX: delta seems not to be added after discontinuation. OUCH!! 

## causal - initial comments

- i ran it with K0=1, K1=1. I reckon that this, compared with J2R, should increase all imputed values by an amount that depends only on when discontinuation occurred. But id=5017 (pattern OXXX) gets exactly the same imputations as J2R but should be different, while 5115 (pattern OOXO i.e. only interim missing) gets different imputations but should be same as J2R. OUCH!! 
- this will need more testing in future.
