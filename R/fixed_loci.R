#' Identify fixed loci
#'
#' A locus is "fixed" if the non-missing sub-vector contains all 0's or all 2's (the locus is completely homozygous for one allele or completely homozygous for the other allele).
#' This function tests each locus, returning a vector that is `TRUE` for each fixed locus, `FALSE` otherwise.
#' Loci with only missing elements (`NA`) are treated as fixed.
#' Below `m` is the number of loci, and `n` is the number of individuals.
#'
#' @param X The `m`-by-`n` genotype matrix
#'
#' @return A length-`m` boolean vector where the `i` element is `TRUE` if locus `i` is fixed or completely missing, `FALSE` otherwise.
#'
#' @examples
#' # here's a toy genotype matrix
#' X <- matrix(
#'        data = c(
#'               2, 2, NA,  # fixed locus (with one missing element)
#'               0, NA, 0,  # another fixed locus, for opposite allele
#'               1, 1, 1,   # NOT fixed (heterozygotes are not considered fixed)
#'               0, 1, 2,   # a completely variable locus
#'               NA, NA, NA # completely missing locus (will be treated as fixed)
#'              ),
#'        ncol = 3, byrow = TRUE)
#' 
#' # test that we get the desired values
#' stopifnot(
#'   fixed_loci(X) == c(TRUE, TRUE, FALSE, FALSE, TRUE)
#' )
#'
#' @export
fixed_loci <- function(X) {
    # stop if required parameters are missing
    if (missing(X))
        stop('`X` is required!')
    # ensure that things that should be matrices are so
    if (!is.matrix(X))
        stop('`X` must be a matrix!')
    
    # is this too slow? (what if we did it using Rcpp?)
    # main step is calculating allele frequencies per row, which automatically handles missingness
    # this returns NaN for completely missing rows
    p_anc_hat <- rowMeans(X, na.rm = TRUE) / 2
    
    # this is the return value we want
    is.na(p_anc_hat) | p_anc_hat == 0 | p_anc_hat == 1
}
