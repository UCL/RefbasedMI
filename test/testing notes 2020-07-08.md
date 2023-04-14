# Ian's mimix testing notes, 8jul2020

updated to v0.0.7
- had to delete package ellipsis at N:\R\win-library\3.6.then reinstall it, in order for install_github() to work.

## help file

- see word doc

## front page

## causal

- *** need a timescale option: is k1 per unit time or per visit?

- i compared k0=k1=0 with J2R: perfect agreement

- *** i compared k0=k1=1 with CIR: disagreement after discontinuation when there were intermittent missing data before discontinuation.

- *** I was unclear how to verify the results more broadly. Do you have code for the examples in the paper?

## delta

- *** setting dlag of the wrong length leads to missing imputed values and "WARNING!!! unimputed data values, possibly due to mis-specified delta". Should catch this error and fail, as for delta of wrong length.

- tested for arbitrary values of delta and dlag and verified that the increases in imputed values are as expected after discontinuation and are zero before discontinuation

## main package
- 3 arms: changed first 200 obs in asthma data to a 3rd arm. Verified agreement with Stata package to within MC error with M=100. 

- Also noted R package is MUCH faster than Stata.

- added second covariate (=base^2). Verified agreement with Stata package to within MC error with M=100. 

- multiplied outcome by 1000. Results are exactly multipled by 1000.

- after conversion to mids format, attributes such as pmm are wrong. However this is NOT A PROBLEM as it is happening outside our function.

## to discuss

- allow incomplete covariates? or add a recommendation to use mean imputation or missing indicator?