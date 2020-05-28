#' Draw allele frequencies for independent subpopulations
#'
#' Allele frequencies \eqn{p_i^{S_u}} for independent subpopulations \eqn{S_u} at locus \eqn{i} are drawn from the Balding-Nichols distribution with ancestral allele frequency \eqn{p_i^T} and \eqn{F_{ST}}{FST} parameter \eqn{f^T_{S_u}} as
#' \deqn{p_i^{S_u} \sim \mbox{Beta}(\nu_u p_i^T, \nu_u (1-p_i^T)),}
#' where \eqn{\nu_u = 1/f^T_{S_u} - 1}.
#' Below \eqn{m} is the number of loci and \eqn{k} is the number of subpopulations.
#'
#' @param p_anc The scalar or length-\eqn{m} vector of ancestral allele frequencies per locus.
#' @param inbr_subpops The length-\eqn{k} vector of subpopulation \eqn{F_{ST}}{FST} values.
#' @param m_loci Optional.
#' The desired number of loci \eqn{m}, to be used if \code{p_anc} is a scalar.
#' Stops if both \code{length(p_anc) > 1} and \code{m_loci} are set and they disagree.
#' @param k_subpops Optional.
#' The desired number of subpopulations \eqn{k}, to be used if \code{inbr_subpops} is a scalar.
#' Stops if both \code{length(inbr_subpops) > 1} and \code{k_subpops} are set and they disagree.
#'
#' @return The \eqn{m \times k}{m-by-k} matrix of independent subpopulation allele frequencies
#'
#' @examples
#' # a typical, non-trivial example
#' # number of loci
#' m_loci <- 10
#' # random vector of ancestral allele frequencies
#' p_anc <- draw_p_anc(m_loci)
#' # FST values for two subpops
#' inbr_subpops <- c(0.1, 0.3)
#' # matrix of intermediate subpop allele freqs
#' p_subpops <- draw_p_subpops(p_anc, inbr_subpops)
#'
#' # special case of scalar p_anc
#' p_subpops <- draw_p_subpops(p_anc = 0.5, inbr_subpops, m_loci = m_loci)
#' stopifnot ( nrow( p_subpops ) == m_loci )
#'
#' # special case of scalar inbr_subpops
#' k_subpops <- 2
#' p_subpops <- draw_p_subpops(p_anc, inbr_subpops = 0.2, k_subpops = k_subpops)
#' stopifnot ( ncol( p_subpops ) == k_subpops )
#' 
#' # both main parameters scalars but return value still matrix
#' p_subpops <- draw_p_subpops(p_anc = 0.5, inbr_subpops = 0.2, m_loci = m_loci, k_subpops = k_subpops)
#' stopifnot ( nrow( p_subpops ) == m_loci )
#' stopifnot ( ncol( p_subpops ) == k_subpops )
#'
#' # passing scalar parameters without setting dimensions separately results in a 1x1 matrix
#' p_subpops <- draw_p_subpops(p_anc = 0.5, inbr_subpops = 0.2)
#' stopifnot ( nrow( p_subpops ) == 1 )
#' stopifnot ( ncol( p_subpops ) == 1 )
#'
#' @export
draw_p_subpops <- function(p_anc, inbr_subpops, m_loci = NA, k_subpops = NA) {
    # basic param checking
    if (missing(p_anc))
        stop('ancestral allele frequencies `p_anc` are required!')
    if (missing(inbr_subpops))
        stop('`inbr_subpops` (FST) scalar or vector are required!')
    
    # number of loci to simulate
    # allow p_anc to be a scalar, m_loci can be passed separately
    # actual length
    m_loci_p <- length(p_anc)
    # make sure both things were not set and contradict each other
    if (
        m_loci_p > 1 &&    # length(p_anc) > 1
        !is.na(m_loci) &&  # and m_loci was also set
        m_loci != m_loci_p # and they disagree
    )
        stop('length of `p_anc` (', m_loci_p, ') disagrees with passed parameter `m_loci` (', m_loci, ')')
    # if this is missing, always set to actual value (even if it is 1)
    if ( is.na( m_loci ) )
        m_loci <- m_loci_p
    
    # number of subpopulations to simulate
    # allow for scalar inbr_subpops here, k_subpops can be passed separately
    # actual length
    k_subpops_inbr <- length(inbr_subpops)
    # stop if both things were set and don't agree
    if (
        k_subpops_inbr > 1 &&       # length(inbr_subpops) > 1
        !is.na(k_subpops) &&        # and k_subpops was also set
        k_subpops != k_subpops_inbr # and they disagree
    )
        stop('length of `inbr_subpops` (', k_subpops_inbr, ') disagrees with passed parameter `k_subpops` (', k_subpops, ')')
    # if this is missing, always set to actual value (even if it is 1)
    if ( is.na( k_subpops ) )
        k_subpops <- k_subpops_inbr
    
    # let's translate parameters for Balding-Nichols case
    nu <- 1 / inbr_subpops - 1 # nu is a vector or a scalar, same as inbr_subpops (whatever that is)
    p_anc_alt <- 1 - p_anc # precompute vector of "alternative" p_anc's, shared by all subpopulations below
    # vectorization makes a lot of sense for each subpopulation... (doing all SNPs together)
    p_subpops <- matrix(nrow = m_loci, ncol = k_subpops) # matrix of intermediate allele frequencies we want...
    for ( j in 1 : k_subpops ) {
        # handle scalar inbr_subpops here
        nuj <- if (k_subpops_inbr == 1) nu else nu[j] 
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
