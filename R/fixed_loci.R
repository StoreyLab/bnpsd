#' Identify fixed loci
#'
#' A locus is "fixed" if the non-missing sub-vector contains all 0's or all 2's (the locus is completely homozygous for one allele or completely homozygous for the other allele).
#' This function tests each locus, returning a vector that is `TRUE` for each fixed locus, `FALSE` otherwise.
#' Loci with only missing elements (`NA`) are treated as fixed.
#' The parameter `maf_min` extends the "fixed" definition to loci whose minor allele frequency is smaller or equal than this value.
#' Below `m` is the number of loci, and `n` is the number of individuals.
#'
#' @param X The `m`-by-`n` genotype matrix
#' @param maf_min The minimum minor allele frequency (default zero), to extend the working definition of "fixed" to include rare variants.
#' Loci with minor allele frequencies less than or *equal* to this value are marked as fixed.
#' Must be a scalar between 0 and 0.5.
#'
#' @return A length-`m` boolean vector where the `i` element is `TRUE` if locus `i` is fixed or completely missing, `FALSE` otherwise.
#' If `X` had row names, they are copied to the names of this output vector.
#'
#' @examples
#' # here's a toy genotype matrix
#' X <- matrix(
#'        data = c(
#'               2, 2, NA,  # fixed locus (with one missing element)
#'               0, NA, 0,  # another fixed locus, for opposite allele
#'               1, 1, 1,   # NOT fixed (heterozygotes are not considered fixed)
#'               0, 1, 2,   # a completely variable locus
#'               0, 0, 1,   # a somewhat "rare" variant
#'               NA, NA, NA # completely missing locus (will be treated as fixed)
#'              ),
#'        ncol = 3, byrow = TRUE)
#' 
#' # test that we get the desired values
#' stopifnot(
#'   fixed_loci(X) == c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE)
#' )
#'
#' # the "rare" variant gets marked as "fixed" if we set `maf_min` to its frequency
#' stopifnot(
#'   fixed_loci(X, maf_min = 1/6) == c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE)
#' )
#'
#' @export
fixed_loci <- function( X, maf_min = 0 ) {
    # stop if required parameters are missing
    if (missing(X))
        stop('`X` is required!')
    # ensure that things that should be matrices are so
    if (!is.matrix(X))
        stop('`X` must be a matrix!')

    # check `maf_min` range, etc
    if ( !is.numeric( maf_min ) )
        stop( '`maf_min` must be numeric!' )
    if ( length( maf_min ) != 1 )
        stop( '`maf_min` must be a scalar (have length 1)!' )
    if ( is.na( maf_min ) )
        stop( '`maf_min` cannot be `NA`!' )
    if ( maf_min < 0 )
        stop( '`maf_min` must be non-negative!' )
    if ( maf_min > 0.5 )
        stop( '`maf_min` must not exceed 0.5!' )
    
    # is this too slow? (what if we did it using Rcpp?)
    # main step is calculating allele frequencies per row, which automatically handles missingness
    # this returns NaN for completely missing rows
    # NOTE: if X had rownames, they become the names of p_anc_hat!
    p_anc_hat <- rowMeans(X, na.rm = TRUE) / 2
    
    # this is the return value we want
    # NOTE: rownames are again inherited into the output automatically!
    is.na(p_anc_hat) | p_anc_hat <= maf_min | p_anc_hat >= 1 - maf_min
}
