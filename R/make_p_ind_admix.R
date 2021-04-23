#' Construct individual-specific allele frequency matrix under the PSD admixture model
#'
#' Here `m` is the number of loci, `n` the number of individuals, and `k` the number of intermediate subpopulations.
#' The `m`-by-`n` individual-specific allele frequency matrix `p_ind` is constructed from the `m`-by-`k` intermediate subpopulation allele frequency matrix `p_subpops` and the `n`-by-`k` admixture proportion matrix `admix_proportions` equivalent to
#' `p_ind <- p_subpops %*% t( admix_proportions )`.
#' This function is a wrapper around [tcrossprod()], but also ensures the output allele frequencies are in \[0, 1\] (otherwise not guaranteed due to limited machine precision).
#' 
#' If both `admix_proportions` and `p_subpops` have column names, and if they disagree, the function stops as a precaution, as this suggests the data is misaligned or inconsistent in some way.
#' 
#' @param p_subpops The `m`-by-`k` matrix of intermediate subpopulation allele frequencies.
#' @param admix_proportions The `n`-by-`k` matrix of admixture proportions.
#'
#' @return The `m`-by-`n` matrix of individual-specific allele frequencies `p_ind`.
#' Row names equal those from `p_subpops`, and column names equal rownames from `admix_proportions`, if present.
#'
#' @examples
#' # data dimensions
#' # number of loci
#' m_loci <- 10
#' # number of individuals
#' n_ind <- 5
#' # number of intermediate subpops
#' k_subpops <- 2
#'
#' # FST values for k = 2 subpops
#' inbr_subpops <- c(0.1, 0.3)
#' 
#' # non-trivial admixture proportions
#' admix_proportions <- admix_prop_1d_linear(n_ind, k_subpops, sigma = 1)
#' 
#' # random vector of ancestral allele frequencies
#' p_anc <- draw_p_anc(m_loci)
#' 
#' # matrix of intermediate subpop allele freqs
#' p_subpops <- draw_p_subpops(p_anc, inbr_subpops)
#' 
#' # matrix of individual-specific allele frequencies
#' p_ind <- make_p_ind_admix(p_subpops, admix_proportions)
#'
#' @export
make_p_ind_admix <- function(p_subpops, admix_proportions) {
    # make sure nothing is missing
    if (missing(p_subpops))
        stop('`p_subpops` is required!')
    if (missing(admix_proportions))
        stop('`admix_proportions` is required!')
    
    # ensure that things that should be matrices are so
    if (!is.matrix(p_subpops))
        stop('`p_subpops` must be a matrix!')
    if (!is.matrix(admix_proportions))
        stop('`admix_proportions` must be a matrix!')
    
    # validate data dimensions
    if (ncol(p_subpops) != ncol(admix_proportions))
        stop('`p_subpops` and `admix_proportions` are not compatible: ncol(p_subpops) == ', ncol(p_subpops), ' != ', ncol(admix_proportions), ' == ncol(admix_proportions)')

    # if both inputs have names, let's stop if they don't agree, explain if it appears that they are not aligned in particular
    compare_names(
        colnames( p_subpops ),
        colnames( admix_proportions ),
        'p_subpops',
        'admix_proportions'
    )
    
    # this is the main multiplication
    # NOTE: row and column names for `p_ind` are automatically inherited as desired!
    p_ind <- tcrossprod(p_subpops, admix_proportions)
    
    # sometimes p_ind has values slighly outside of [0,1], simply due to machine precision errors
    # fix that here!
    p_ind[p_ind < 0] <- 0
    p_ind[p_ind > 1] <- 1
    
    return(p_ind)
}
