# mimix
 *0.0.14*

<h2>Multiple imputation sensitivity analysis for longitudinal trials using R</h2> 

We aim to port the functionality of the Stata program **mimix**  into R 

The purpose of mimix is as described in the paper

Reference-based sensitivity analysis via multiple imputation for longitudinal trials with protocol deviation
by Suzie Cro, Tim P. Morris, Michael G. Kenward, and James R. Carpenter
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5796638/

The 5 methods available for sensitivity analysis are
 
|  Method         | option cmd             | reference group required |
| --------------- | --------------- | --------------------  |
| Randomized-arm                | MAR |  n |
| Jump to reference	            | J2R |  y |
| Copy increments in reference	| CIR |  y |
| Copy reference	              | CR  |  y |
| Last mean carried forward	    | LMCF|  n |


The R program does not provide the interim option available in Stata (where the individual has data observed later ) or the methodvar option (where  different imputation methods are required for different individuals). 

# installation

To install package from Github

library(githubinstall)
githubinstall("rmimixpackage")
