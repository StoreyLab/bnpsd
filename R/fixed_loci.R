#' Identify fixed loci
#'
#' A locus is said to be fixed if the non-missing sub-vector contains all 0's or all 2's (the locus is completely homozygous for one allele or completely homozygous for the other allele).
#' This function tests each locus, returning a vector that is TRUE for each fixed locus, FALSE otherwise.
#' A locus with only missing elements (NA) will also be marked as fixed (TRUE).
#' Below \eqn{m} is the number of loci, and \eqn{n} is the number of individuals.
#'
#' @param X The \eqn{m \times n}{m-by-n} genotype matrix
#'
#' @return A length-\eqn{m} boolean vector where the \eqn{i}th element is TRUE if locus \eqn{i} is fixed or completely missing, FALSE otherwise.
#'
#' @examples
#' # here's a toy matrix
#' X <- matrix(
#'        data = c(
#'               2, 2, NA, # fixed locus (with one missing element)
#'               0, NA, 0, # another fixed locus, for opposite allele
#'               1, 1, 1, # NOT fixed (heterozygotes are not considered fixed)
#'               0, 1, 2, # a completely variable locus
#'               NA, NA, NA # completely missing locus (will be treated as fixed)
#'              ),
#'        ncol = 3, byrow = TRUE)
#' # test that we get the desired values
#' stopifnot(
#'   fixed_loci(X) == c(TRUE, TRUE, FALSE, FALSE, TRUE)
#' )
#'
#' @export
fixed_loci <- function(X) {
    # is this too slow? (what if we did it using Rcpp?)
    # main step is calculating allele frequencies per row, which automatically handles missingness
    # this returns NaN for completely missing rows
    pHat <- rowMeans(X, na.rm = TRUE) / 2
    # this is the return value we want
    is.na(pHat) | pHat == 0 | pHat == 1
}
