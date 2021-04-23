#' Draw allele frequencies for independent subpopulations
#'
#' The allele frequency matrix `P` for `m_loci` loci (rows) and `k_subpops` independent subpopulations (columns) are drawn from the Balding-Nichols distribution with ancestral allele frequencies `p_anc` and FST parameters `inbr_subpops` equivalent to
#' `P[ i, j ] <- rbeta( 1, nu_j * p_anc[i], nu_j * ( 1 - p_anc[i] ) )`,
#' where `nu_j <- 1 / inbr_subpops[j] - 1`.
#' The actual function is more efficient than the above code.
#'
#' @param p_anc The scalar or length-`m_loci` vector of ancestral allele frequencies per locus.
#' @param inbr_subpops The scalar or length-`k_subpops` vector of subpopulation FST values.
#' @param m_loci If `p_anc` is scalar, optionally provide the desired number of loci (lest only one locus be simulated).
#' Stops if both `length(p_anc) > 1` and `m_loci` is not `NA` and they disagree.
#' @param k_subpops If `inbr_subpops` is a scalar, optionally provide the desired number of subpopulations (lest a single subpopulation be simulated).
#' Stops if both `length(inbr_subpops) > 1` and `k_subpops` is not `NA` and they disagree.
#'
#' @return The `m_loci`-by-`k_subpops` matrix of independent subpopulation allele frequencies.
#' If `p_anc` is length-`m_loci` with names, these are copied to the row names of this output matrix.
#' If `inbr_subpops` is length-`k_subpops` with names, these are copied to the column names of this output matrix.
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
#' @seealso
#' [draw_p_subpops_tree()] for version for subpopulations related by a tree, which can therefore be non-independent.
#' 
#' @export
draw_p_subpops <- function(p_anc, inbr_subpops, m_loci = NA, k_subpops = NA) {
    # basic param checking
    if (missing(p_anc))
        stop('ancestral allele frequencies `p_anc` are required!')
    if (missing(inbr_subpops))
        stop('`inbr_subpops` (FST) scalar or vector are required!')
    
    # look at data ranges
    # all of these are probabilities so problems happen when they're out of range
    if ( any( p_anc < 0 ) )
        stop( '`p_anc` cannot be negative!' )
    if ( any( p_anc > 1 ) )
        stop( '`p_anc` cannot exceed 1!' )
    if ( any( inbr_subpops < 0 ) )
        stop( '`inbr_subpops` cannot be negative!' )
    if ( any( inbr_subpops > 1 ) )
        stop( '`inbr_subpops` cannot exceed 1!' )
    
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

    # add names now if it makes sense
    # p_anc can be scalar, so this only makes sense if the length matches (and names are defined)
    if ( !is.null( names( p_anc ) ) && length( p_anc ) == m_loci )
        rownames( p_subpops ) <- names( p_anc )
    # similarly for inbr_subpops
    if ( !is.null( names( inbr_subpops ) ) && length( inbr_subpops ) == k_subpops )
        colnames( p_subpops ) <- names( inbr_subpops )
    
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
    
    return( p_subpops ) # the only thing we want out of this
}
