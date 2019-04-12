#' Calculate FST for the admixed individuals
#'
#' Given the admixture proportion matrix \eqn{Q} for \eqn{n} individuals and \eqn{k} intermediate subpopulations, the vector of intermediate inbreeding coefficients \eqn{F} (per-subpopulation \eqn{F_{ST}}{FST}'s), and weights for individuals, this function returns the \eqn{F_{ST}}{FST} of the admixed individuals.
#' This \eqn{F_{ST}}{FST} equals the weighted mean of the diagonal of the coancestry matrix (see \code{\link{coanc_admix}}).
#' 
#' @param Q The \eqn{n \times k}{n-by-k} admixture proportion matrix
#' @param F The length-\eqn{k} vector of subpopulation inbreeding coefficients
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
#' F <- (1:k)/k
#' # construct final admixture proportions
#' Q <- q1d(n = n, k = k, sigma = sigma)
#' # lastly, calculate Fst!!! (uniform weights in this case)
#' F <- fst(Q, F)
#'
#' @export
fst <- function(Q, F, weights = NULL) {
    # die informatively...
    if (missing(Q))
        stop('Q is missing!')
    if (missing(F))
        stop('F is missing!')
    
    # calculate necessary intermediates
    # could probably be more efficient skipping un-needed matrix products, but such efficiency is not needed here because this function is usually called once only and n is big but not huge

    # coancestry matrix
    coancestry <- coanc_admix(Q, F)
    # vector of inbreeding coefficients of individuals
    inbreeding <- diag(coancestry)
    
    # return weighted mean inbreeding (Fst)
    if (is.null(weights)) {
        return( mean( inbreeding ) )
    } else {
        return( drop( inbreeding %*% weights ) )
    }
}

