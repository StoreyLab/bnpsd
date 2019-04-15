#' Draw random uniform ancestral allele frequencies
#'
#' This is simply a wrapper around \code{\link[stats]{runif}} with different defaults and additional validations.
#'
#' @param m_loci Number of loci to draw
#' @param p_min Minimum allele frequency to draw
#' @param p_max Maximum allele frequency to draw
#'
#' @return A length-\eqn{m} vector of random ancestral allele frequencies
#'
#' @examples
#' p_anc <- draw_p_anc(m_loci = 10)
#' 
#' @export
draw_p_anc <- function(m_loci, p_min = 0.01, p_max = 0.5) {
    # make sure the first argument isn't missing
    if (missing(m_loci))
        stop('`m_loci` is required!')
    
    # some restrictions for allele frequencies
    if (p_min < 0)
        stop('minimum allele frequency cannot be negative: ', p_min)
    if (p_max > 1)
        stop('maximum allele frequency cannot be greater than 1: ', p_max)
    if (p_min > p_max)
        stop('minimum allele frequency (', p_min, ') cannot be larger than the maximum allele frequency (', p_max, ')')
    
    # now draw them!
    stats::runif(m_loci, min = p_min, max = p_max)
}

# stick deprecated function name here

#' @title Draw uniform ancestral allele frequencies
#' @description Draw uniform ancestral allele frequencies
#' @param m Number of loci
#' @param min Minimum allele frequency to draw
#' @param max Maximum allele frequency to draw
#' @return A vector of ancestral allele frequencies
#'
#' @name rpanc-deprecated
#' @usage rpanc(m, min = 0.01, max = 0.5)
#' @seealso \code{\link{bnpsd-deprecated}}
#' @keywords internal
NULL

#' @rdname bnpsd-deprecated
#' @section \code{rpanc}:
#' For \code{rpanc}, use \code{\link{draw_p_anc}}.
#'
#' @export
rpanc <- function(m, min = 0.01, max = 0.5) {
    # mark as deprecated
    .Deprecated('draw_p_anc')
    # return as usual, to not break things just yet
    draw_p_anc(m, min, max)
}
