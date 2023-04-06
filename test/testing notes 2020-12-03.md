# Notes on mimix causal

## Successes in monotone imputation
I tested the impact of increasing fev by x units at time 4 in the treatment arm.
This increases the treatment effect at time 4 by 1 without impact on the treatment effects at other times.
I verified that this increases imputed values for patterns 0100 and 1100 only:
  by x*K0*K1^4 at time 8 
  by x*K0*K1^8 at time 12

## Limitations
I haven't tested non-monotone imputation

## Problems
Could allow K1 not specified if K0=0

Mistake at runmimix.R line 53: running mimix with options
  method="causal",K0=8,K1=0.5
wrongly fails with error 
  K1 Causal constant not in range 0..1 

Can idvar be >1 variable? not at present!

How do I carry another variable into the imputed data set? 

