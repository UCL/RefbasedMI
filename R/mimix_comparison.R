#' mimix: Comparisons with Stata and SAS   
#' 
#' @description 
#' mimix is available in Stata whils th Five_macros suite does similar in SAS  
#' 
#' @section  Comparison with Stata:
#' This mimix is based on the Stata version and has similar functionality while adding
#' the causal method and delta adjustment. 
#' As with the Stata version the input data requires the longitudinal input data in
#' long format with one record per individual at each timepoint.
#' The program differs in how interim missing cases - those cases which 
#' have a missing measurement at a timepoint previous to a later observed measurement - are treated.        
#' Under Stata by default, the interim missing are treated the same as for the post-discontinuation
#' missing unless the interim option is explicitly used. 
#' Here the interims are treated as under MAR, the post-discontinuations then imputed under
#' the specified method. There is no interim option as there is in Stata.
#' Unlike Stata an option is supplied whereby the prior used in the MCMC draws can be changed from the 
#' default jeffreys (as in Stata) to either the ridge or uniform  
#' 
#' 
#' @section Comparison with SAS:
#' 
#' Whilst this program is based on the Stata program, the latter is an adaptation of the SAS macro miwithd,
#' written by james Roger, subsequently updated to the Five_Macros suite of macros
#' This program uses the same approach for the delta adjustment as described in the Five_macros, 
#' in comparing outputs from our program with the Five_macros it is to be noted that interaction between treatment and covariates
#' is not allowed in the SAS macros, and comparisons are only valid for example in testing the Causal model by specifically not
#' not using the covbytime and catcovbytime options in the Five_macros 
#' Not using these options also means that the LMCF method can be compared with either ALMCF or OLMCF in the Five_macros.
#' When there is no observed data (common in the acupuncture data) the first mean is used in Stata,
#' a warning is given in the Five_macros     
#'
#'  
#'  
#' @references 
#'     Cro s, Morris T, Kenward G,Carpenter Joshttps://www.ncbi.nlm.nih.gov/pmc/articles/PMC5796638/
#'     White I, Royes J, Best N, https://arxiv.org/abs/1705.04506
#'     Roger J,
#'     URL: https://www.lshtm.ac.uk/research/centres-projects-groups/missing-data#sensitivity-analysis, 
#' @docType package
#' @name mimix_Comparison 
NULL



# if package then use this  "_PACKAGE"  
#> "_PACKAGE"

# #' @keywords internal 
# @docType package
# @name mimix_package
#NULL
# #> NULL
# #' mimix: A package porting the Stata mimix command  
# #' The mimix package provides the functionality of the Stata package plus
# #' delta and causal methods