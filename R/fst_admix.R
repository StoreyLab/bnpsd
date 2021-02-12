#' Calculate FST for the admixed individuals
#'
#' This function returns the generalized FST of the admixed individuals given their admixture proportion matrix, the coancestry matrix of intermediate subpopulations (or its special cases, see `coanc_subpops` parameter below), and optional weights for individuals.
#' This FST equals the weighted mean of the diagonal of the coancestry matrix (see `\link{coanc_admix}`).
#' Below there are `n` individuals and `k` intermediate subpopulations.
#' 
#' @param admix_proportions The `n`-by-`k` admixture proportion matrix
#' @param coanc_subpops Either the `k`-by-`k` intermediate subpopulation coancestry matrix (for the complete admixture model), or the length-`k` vector of intermediate subpopulation FST values (for the BN-PSD model; assumes zero coancestries between subpopulations), or a scalar FST value shared by all intermediate subpopulations (also assumes zero coancestry between subpopulations).
#' @param weights Optional length-`n` vector of weights for individuals that define their generalized FST (default uniform weights)
#'
#' @return The generalized FST of the admixed individuals
#'
#' @examples
#' # set desired parameters
#' # number of individuals
#' n_ind <- 1000
#' # number of intermediate subpopulations
#' k_subpops <- 10
#' 
#' # differentiation of intermediate subpopulations
#' coanc_subpops <- ( 1 : k_subpops ) / k_subpops
#'
#' # construct admixture proportions
#' admix_proportions <- admix_prop_1d_linear(n_ind, k_subpops, sigma = 1)
#'
#' # lastly, calculate Fst!!! (uniform weights in this case)
#' fst_admix(admix_proportions, coanc_subpops)
#'
#' @export
fst_admix <- function(admix_proportions, coanc_subpops, weights = NULL) {
    # die informatively...
    if (missing( admix_proportions ))
        stop('`admix_proportions` is required!')
    if (missing( coanc_subpops ))
        stop('`coanc_subpops` is required!')
    
    # calculate necessary intermediates
    # could probably be more efficient skipping un-needed matrix products, but such efficiency is not needed here because this function is usually called once only and n is big but not huge

    # coancestry matrix
    coancestry <- coanc_admix( admix_proportions, coanc_subpops )
    # vector of inbreeding coefficients of individuals
    inbreeding <- diag( coancestry )
    
    # return weighted mean inbreeding (Fst)
    if ( is.null(weights) ) {
        return( mean( inbreeding ) )
    } else {
        return( drop( inbreeding %*% weights ) )
    }
}
