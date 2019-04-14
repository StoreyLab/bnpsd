#' Calculate FST for the admixed individuals
#'
#' This function returns the \eqn{F_{ST}}{FST} of the admixed individuals given the admixture proportion matrix for \eqn{n} individuals and \eqn{k} intermediate subpopulations, the coancestry matrix of intermediate subpopulations (or its special cases, see \code{coanc_subpops} parameter below), and optional weights for individuals.
#' This \eqn{F_{ST}}{FST} equals the weighted mean of the diagonal of the coancestry matrix (see \code{\link{coanc_admix}}).
#' 
#' @param admix_proportions The \eqn{n \times k}{n-by-k} admixture proportion matrix
#' @param coanc_subpops Either the \eqn{k \times k}{k-by-k} intermediate subpopulation coancestry matrix (for the complete admixture model), or the length-\eqn{k} vector of intermediate subpopulation \eqn{F_{ST}}{FST} values (for the BN-PSD model), or a scalar \eqn{F_{ST}}{FST} value shared by all intermediate subpopulations.
#' @param weights The length-\eqn{n} vector of weights for individuals that define \eqn{F_{ST}}{FST} (default uniform weights)
#'
#' @return The \eqn{F_{ST}}{FST} of the admixed individuals
#'
#' @examples
#' # set desired parameters
#' n <- 1000 # number of individuals
#' k <- 10 # number of intermediate subpopulations
#' s <- 0.5 # desired bias coefficient
#' sigma <- 1 # for 1D admixture model
#' # differentiation of intermediate subpopulations
#' coanc_subpops <- (1:k)/k
#' # construct admixture proportions
#' admix_proportions <- admix_prop_1d_linear(n, k, sigma)
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
