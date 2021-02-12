#' Draw random Uniform or Beta ancestral allele frequencies
#'
#' This is simply a wrapper around `\link[stats]{runif}` or `\link[stats]{rbeta}` (depending on parameters) with different defaults and additional validations.
#'
#' @param m_loci Number of loci to draw.
#' @param p_min Minimum allele frequency to draw (Uniform case only).
#' @param p_max Maximum allele frequency to draw (Uniform case only).
#' @param beta Shape parameter for a symmetric Beta.
#' If `NA` (default), Uniform(`p_min`, `p_max`) is used.
#' Otherwise, a Symmetric Beta is used and the user-specified range is ignored (values in \[0, 1\] will be returned).
#'
#' @return A length-`m` vector of random ancestral allele frequencies
#'
#' @examples
#' # Default is uniform with range between 0.01 and 0.5
#' p_anc <- draw_p_anc(m_loci = 10)
#'
#' # Use of `beta` triggers a symmetric Beta distribution.
#' # This parameter has increased density for rare minor allele frequencies,
#' # resembling the 1000 Genomes allele frequency distribution
#' p_anc <- draw_p_anc(m_loci = 10, beta = 0.03)
#' 
#' @export
draw_p_anc <- function(m_loci, p_min = 0.01, p_max = 0.5, beta = NA) {
    # make sure the first argument isn't missing
    if (missing(m_loci))
        stop('`m_loci` is required!')

    if ( is.na(beta) ) {
        # original "uniform" case, with potentially custom limits
        
        # some restrictions for allele frequencies
        if (p_min < 0)
            stop('minimum allele frequency cannot be negative: ', p_min)
        if (p_max > 1)
            stop('maximum allele frequency cannot be greater than 1: ', p_max)
        if (p_min > p_max)
            stop('minimum allele frequency (', p_min, ') cannot be larger than the maximum allele frequency (', p_max, ')')
        
        # now draw them!
        stats::runif(m_loci, min = p_min, max = p_max)
    } else {
        # new symmetric Beta distribution
        # in this case we don't restrict to lower half or impose any other range limits, we just return whatever rbeta returns
        stats::rbeta(m_loci, shape1 = beta, shape2 = beta)
    }
}
