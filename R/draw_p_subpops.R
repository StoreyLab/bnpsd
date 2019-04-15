#' Draw allele frequencies for independent subpopulations
#'
#' Allele frequencies \eqn{p_i^{S_u}} for independent subpopulations \eqn{S_u} at locus \eqn{i} are drawn from the Balding-Nichols distribution with ancestral allele frequency \eqn{p_i^T} and \eqn{F_{ST}}{FST} parameter \eqn{f^T_{S_u}} as
#' \deqn{p_i^{S_u} \sim \mbox{Beta}(\nu_u p_i^T, \nu_u (1-p_i^T)),}
#' where \eqn{\nu_u = 1/f^T_{S_u} - 1}.
#' Below \eqn{m} is the number of loci and \eqn{k} is the number of subpopulations.
#'
#' @param p_anc The length-\eqn{m} vector of ancestral allele frequencies per locus
#' @param inbr_subpops The length-\eqn{k} vector of subpopulation \eqn{F_{ST}}{FST} values
#'
#' @return The \eqn{m \times k}{m-by-k} matrix of independent subpopulation allele frequencies
#'
#' @examples
#' m_loci <- 10 # number of loci
#' p_anc <- draw_p_anc(m_loci) # random vector of ancestral allele frequencies
#' inbr_subpops <- c(0.1, 0.3) # FST values for two subpops
#' p_subpops <- draw_p_subpops(p_anc, inbr_subpops) # matrix of intermediate subpop allele freqs
#'
#' @export
draw_p_subpops <- function(p_anc, inbr_subpops) {
    # basic param checking
    if (missing(p_anc))
        stop('ancestral allele frequencies `p_anc` are required!')
    if (missing(inbr_subpops))
        stop('`inbr_subpops` (FST) scalar or vector are required!')
    
    # number of loci and subpopulations to simulate
    m_loci <- length(p_anc)
    k_subpops <- length(inbr_subpops)
    
    # let's translate parameters for Balding-Nichols case
    nu <- 1 / inbr_subpops - 1 # nu is a vector or a scalar, same as inbr_subpops (whatever that is)
    p_anc_alt <- 1 - p_anc # precompute vector of "alternative" p_anc's, shared by all subpopulations below
    # vectorization makes a lot of sense for each subpopulation... (doing all SNPs together)
    p_subpops <- matrix(nrow = m_loci, ncol = k_subpops) # matrix of intermediate allele frequencies we want...
    for ( j in 1 : k_subpops ) {
        nuj <- nu[j]
        if (is.infinite(nuj)) {
            # there is no drift from p_anc in this special case (inbr_subpops == 0)
            # (coded separately because rbeta incorrectly returns 0.5 instead)
            p_subpops[, j] <- p_anc
        } else {
            p_subpops[, j] <- stats::rbeta(m_loci, nuj * p_anc, nuj * p_anc_alt) # draw all SNPs for this population, store immediately
        }
    }
    
    p_subpops # the only thing we want out of this
}

# stick deprecated function name here

#' @title Draw intermediate subpopulation allele frequencies
#' @description Draw intermediate subpopulation allele frequencies
#' @param pAnc The vector of ancestral allele frequencies per locus
#' @param F The vector of subpopulation \eqn{F_{ST}}{FST} values
#' @return The matrix of intermediate subpopulation allele frequencies
#'
#' @name rpint-deprecated
#' @usage rpint(pAnc, F)
#' @seealso \code{\link{bnpsd-deprecated}}
#' @keywords internal
NULL

#' @rdname bnpsd-deprecated
#' @section \code{rpint}:
#' For \code{rpint}, use \code{\link{draw_p_subpops}}.
#'
#' @export
rpint <- function(pAnc, F) {
    # mark as deprecated
    .Deprecated('draw_p_subpops')
    # return as usual, to not break things just yet
    draw_p_subpops(pAnc, F)
}
