# Ian's mimix testing notes, 8feb2021

## version tested

Runmimix.R committed to github 8feb2021

## Help file
- Description is weak: prefer "Multiply impute missing data using reference-based imputation".
- Details needs some clarification too
- Value, what does "impdatset" mean? sholdn't this say "an object of type xxx containing yyy"?
- Examples: why "Not run"?
- need to say "mle" option is not encouraged: "this option produces data that are NOT suitable for combination using Rubin's rules"
- dlag are "delta lags" not "delta values"; typo "Roger's", doc not paper, and give location
- harmonise capitalisation
- is there a vignette?

## Output

- typo "begining"
- suppress "interim at id=#" (except as debugging option)
- is "test pass2 in runmimx" correct?

# Main methods - J2R, CIR, CR, MAR - after discontinuation

- show all methods give exactly the same imputed values in reference arm
YES

- show similar results (imputed values and treatment effect) when used with a different seed 
YES VERY SIMILAR, AND ABOUT THE RIGHT DEGREE OF SIMILARITY

- show similar results (imputed values and treatment effect) to Stata
YES SIMILAR WITHIN MC ERROR

# Causal model

I suggest looking at just one case e.g. K0=1, K1=0.9

- show K0=0 corresponds to J2R

- show K0=K1=1 corresponds to CIR

- increment all observed outcomes at a particular time t in the non-reference arm, and show this causes all later imputed values in individuals who discontinue at t to increase by K0*K1^time, and does not affect other imputed values.

-- done by Ian (Ian_test_causal.r)

- show similar results (imputed values and treatment effect) to SAS 

??

# Interim missing, imputed as MAR

- show interim imputed values are exactly the same for all methods (J2R, CIR, CR, MAR, LMCF, causal) and all reference arms (active, placebo)

CHECKED FOR J2R/ref 1 = CIR/ref 2

- show changing the interim method changes only the interim imputed values

NOT NOW AVAILABLE - ONLY DOING MAR



# Delta method

- show that a simple specification of delta increases all post-discontinuation imputed values correctly and does not affect interim imputed values

CHECKED - CORRECT!
