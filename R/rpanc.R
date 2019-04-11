#' Draw uniform ancestral allele frequencies
#'
#' This is simply a wrapper around \code{\link[stats]{runif}} with different defaults and additional validations.
#'
#' @param m Number of loci
#' @param min Minimum allele frequency to draw
#' @param max Maximum allele frequency to draw
#'
#' @return A length-\eqn{m} vector of ancestral allele frequencies
#'
#' @examples
#' pAnc <- rpanc(m = 10)
#' 
#' @export
rpanc <- function(m, min = 0.01, max = 0.5) {
    # some restrictions for allele frequencies
    if (min < 0)
        stop('minimum allele frequency cannot be negative: ', min)
    if (max > 1)
        stop('maximum allele frequency cannot be greater than 1: ', max)
    if (min > max)
        stop('minimum allele frequency (', min, ') cannot be larger than the maximum allele frequency (', max, ')')
    
    # now draw them!
    stats::runif(m, min = min, max = max)
}
