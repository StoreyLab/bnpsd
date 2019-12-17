#' Simulate random allele frequencies and genotypes from the BN-PSD admixture model
#'
#' This function returns simulated ancestral, intermediate, and individual-specific allele frequencies and genotypes given the admixture structure, as determined by the admixture proportions and the vector of intermediate subpopulation \eqn{F_{ST}}{FST} values.
#' The function is a wrapper around \code{\link{draw_p_anc}}, \code{\link{draw_p_subpops}}, \code{\link{make_p_ind_admix}}, and \code{\link{draw_genotypes_admix}} with additional features such as requiring polymorphic loci.
#' Importantly, by default fixed loci are re-drawn from the start (starting from the ancestral allele frequencies) so no fixed loci are in the output and no biases are introduced by re-drawing genotypes conditional on any of the previous allele frequencies (ancestral, intermediate, or individual-specific).
#' Below \eqn{m} is the number of loci, \eqn{n} is the number of individuals, and \eqn{k} is the number of intermediate subpopulations.
#'
#' @param admix_proportions The \eqn{n \times k}{n-by-k} matrix of admixture proportions.
#' @param inbr_subpops The length-\eqn{k} vector (or scalar) of intermediate subpopulation \eqn{F_{ST}}{FST} values.
#' @param m_loci The number of loci to draw.
#' @param want_genotypes If `TRUE` (default), includes the matrix of random genotypes in the return list.
#' @param want_p_ind If `TRUE` (NOT default), includes the matrix of individual-specific allele frequencies in the return list.
#' Note that by default `p_ind` is not constructed in full at all, instead a fast low-memory algorithm constructs it in parts as needed only; beware that setting `want_p_ind = TRUE` increases memory usage in comparison.
#' @param want_p_subpops If `TRUE` (NOT default), includes the matrix of random intermediate subpopulation allele frequencies in the return list.
#' @param want_p_anc If `TRUE` (default), includes the matrix of random ancestral allele frequencies in the return list.
#' @param verbose If `TRUE`, prints messages for every stage in the algorithm.
#' @param require_polymorphic_loci If TRUE (default), returned genotype matrix will not include any fixed loci (loci that happened to be fixed are drawn again, starting from their ancestral allele frequencies, and checked iteratively until no fixed loci remain, so that the final number of polymorphic loci is exactly \eqn{m_loci}).
#' @param beta Shape parameter for a symmetric Beta for ancestral allele frequencies `p_anc`.
#' If `NA` (default), `p_anc` is uniform with range in \[0.01, 0.5\].
#' Otherwise, `p_anc` has a symmetric Beta distribution with range in \[0, 1\].
#'
#' @return A named list that includes the following randomly-generated data in this order:
#' \describe{
#'   \item{X:}{
#'     An \eqn{m \times n}{m-by-n} matrix of genotypes.
#'     Included if \code{want_genotypes = TRUE}.
#'   }
#'   \item{p_anc:}{
#'     A length-\eqn{m} vector of ancestral allele frequencies.
#'     Included if \code{want_p_anc = TRUE}.
#'   }
#'   \item{p_subpops:}{
#'     An \eqn{m \times k}{m-by-k} matrix of intermediate subpopulation allele frequencies
#'     Included if \code{want_p_subpops = TRUE}.
#'   }
#'   \item{p_ind:}{
#'     An \eqn{m \times n}{m-by-n} matrix of individual-specific allele frequencies.
#'     Included if \code{want_p_ind = TRUE}.
#'   }
#' }
#'
#' @examples
#' # dimensions
#' # number of loci
#' m_loci <- 10
#' # number of individuals
#' n_ind <- 5
#' # number of intermediate subpops
#' k_subpops <- 2
#'
#' # define population structure
#' # FST values for k = 2 subpopulations
#' inbr_subpops <- c(0.1, 0.3)
#' # admixture proportions from 1D geography
#' admix_proportions <- admix_prop_1d_linear(n_ind, k_subpops, sigma = 1)
#'
#' # draw all random allele freqs and genotypes
#' out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci)
#'
#' # return value is a list with these items:
#' 
#' # genotypes
#' X <- out$X
#' 
#' # ancestral AFs
#' p_anc <- out$p_anc
#' 
#' # # these are excluded by default, but would be included if ...
#' # # ... `want_p_subpops == TRUE`
#' # # intermediate subpopulation AFs
#' # p_subpops <- out$p_subpops
#' # 
#' # # ... `want_p_ind == TRUE`
#' # # individual-specific AFs
#' # p_ind <- out$p_ind
#' 
#' @export
draw_all_admix <- function(
                           admix_proportions,
                           inbr_subpops,
                           m_loci,
                           want_genotypes = TRUE,
                           want_p_ind = FALSE,
                           want_p_subpops = FALSE,
                           want_p_anc = TRUE,
                           verbose = FALSE,
                           require_polymorphic_loci = TRUE,
                           beta = NA
                           ) {
    # stop if required parameters are missing
    if (missing(admix_proportions))
        stop('`admix_proportions` is required!')
    if (missing(inbr_subpops))
        stop('`inbr_subpops` is required!')
    if (missing(m_loci))
        stop('`m_loci` is required!')

    # ensure that things that should be matrices are so
    if (!is.matrix(admix_proportions))
        stop('`admix_proportions` must be a matrix!')
    
    # get dimensions, test coherence
#    n_ind <- nrow(admix_proportions) # actually not used anywhere!?!
    k_subpops <- ncol(admix_proportions)
    k_subpops_inbr <- length(inbr_subpops)
    if (
        k_subpops_inbr > 1 &&       # it's ok if inbr_subpops is a scalar in this case
        k_subpops != k_subpops_inbr # but if it's not scalar, it must agree with admix_proportions
    )
        stop('`admix_proportions` and `inbr_subpops` are not compatible: ncol(admix_proportions) == ', k_subpops, ' != ', k_subpops_inbr, ' == length(inbr_subpops)')
    
    # generate the random ancestral allele frequencies, in usual range and with minimum threshold for simplicity
    if (verbose)
        message('drawing p_anc')
    p_anc <- draw_p_anc(m_loci, beta = beta)
    
    # draw intermediate allele frequencies from Balding-Nichols
    if (verbose)
        message('drawing p_subpops')
    # pass m_loci for consistency check (it already ought to match length(p_anc))
    # pass k_subpops in case inbr_subpops was a scalar (otherwise provides a redundant check)
    p_subpops <- draw_p_subpops(p_anc, inbr_subpops, m_loci = m_loci, k_subpops = k_subpops)
    
    if (want_p_ind) {
        # this triggers higher-memory path
        if (verbose)
            message('drawing p_ind')
        p_ind <- make_p_ind_admix(p_subpops, admix_proportions)
        
        # draw genotypes from p_ind
        if (want_genotypes) {
            if (verbose)
                message('drawing X')
            X <- draw_genotypes_admix(p_ind)
        }
    } else {
        # draw genotypes skipping p_ind
        # (always low memory now!)
        if (want_genotypes) {
            if (verbose)
                message('drawing X')
            X <- draw_genotypes_admix(p_subpops, admix_proportions)
        }
    }
    
    if (require_polymorphic_loci && want_genotypes) {
        # check for fixed loci, draw them again if needed
        # Note this only applies to want_genotypes==TRUE, since p_ind and p_subpops are continuous and therefore practically never truly fixed
        fixed_loci_indexes <- fixed_loci(X) # boolean vector identifies fixed loci
        m_loci_fixed <- sum(fixed_loci_indexes) # number of cases
        if (m_loci_fixed > 0) {
            # call self with desired number of loci, all the same parameters otherwise
            # note that since this is also called asking for no fixed loci, it will work recursively within itself and return when all loci desired are not fixed!
            if (verbose)
                message('re-drawing fixed loci')
            obj <- draw_all_admix(
                admix_proportions = admix_proportions,
                inbr_subpops = inbr_subpops,
                m_loci = m_loci_fixed,
                want_genotypes = want_genotypes,
                want_p_ind = want_p_ind,
                want_p_subpops = want_p_subpops,
                want_p_anc = want_p_anc,
                verbose = FALSE, # don't show more messages for additional iterations
                require_polymorphic_loci = require_polymorphic_loci,
                beta = beta
            )
            # overwrite fixed loci with redrawn polymorphic data
            X[fixed_loci_indexes, ] <- obj$X # guaranteed to be there
            if (want_p_ind)
                p_ind[fixed_loci_indexes, ] <- obj$p_ind
            if (want_p_subpops)
                p_subpops[fixed_loci_indexes, ] <- obj$p_subpops
            if (want_p_anc)
                p_anc[fixed_loci_indexes] <- obj$p_anc
        }
    }
    
    # now prepare output list
    out <- list()
    if (want_genotypes)
        out$X <- X
    if (want_p_anc)
        out$p_anc <- p_anc
    if (want_p_subpops)
        out$p_subpops <- p_subpops
    if (want_p_ind)
        out$p_ind <- p_ind
    return(out)
}
