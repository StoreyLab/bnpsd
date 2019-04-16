#' Construct individual-specific allele frequency matrix under the PSD admixture model
#'
#' Here \eqn{m} is the number of loci, \eqn{n} the number of individuals, and \eqn{k} the number of intermediate subpopulations.
#' The \eqn{m \times n}{m-by-n} Individual-specific Allele Frequency (IAF) matrix \eqn{P} is constructed from the \eqn{m \times k}{m-by-k} intermediate subpopulation allele frequency matrix \eqn{B} and the \eqn{n \times k}{n-by-k} admixture proportion matrix \eqn{Q} using
#' \deqn{P = B Q^T.}{P = B * Q^T.}
#' This function is a wrapper around \code{\link{tcrossprod}}, but also ensures the output allele frequencies are in [0, 1], as this is not guaranteed by \code{\link{tcrossprod}} due to limited machine precision.
#' 
#' @param p_subpops The \eqn{m \times k}{m-by-k} matrix of intermediate subpopulation allele frequencies.
#' @param admix_proportions The \eqn{n \times k}{n-by-k} matrix of admixture proportions.
#'
#' @return The \eqn{m \times n}{m-by-n} matrix of individual-specific allele frequencies.
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

    # this is the main multiplication
    p_ind <- tcrossprod(p_subpops, admix_proportions)
    
    # sometimes p_ind has values slighly outside of [0,1], simply due to machine precision errors
    # fix that here!
    p_ind[p_ind < 0] <- 0
    p_ind[p_ind > 1] <- 1
    
    return(p_ind)
}

# stick deprecated function name here

#' @title Construct individual-specific allele frequency matrix under the PSD admixture model
#' @description Construct individual-specific allele frequency matrix under the PSD admixture model
#' @param B The matrix of intermediate subpopulation allele frequencies.
#' @param Q The matrix of admixture proportions.
#' @return The matrix of individual-specific allele frequencies.
#'
#' @name rpiaf-deprecated
#' @usage rpiaf(B, Q)
#' @seealso \code{\link{bnpsd-deprecated}}
#' @keywords internal
NULL

#' @rdname bnpsd-deprecated
#' @section \code{rpiaf}:
#' For \code{rpiaf}, use \code{\link{make_p_ind_admix}}.
#'
#' @export
rpiaf <- function(B, Q) {
    # mark as deprecated
    .Deprecated('make_p_ind_admix')
    # return as usual, to not break things just yet
    make_p_ind_admix(B, Q)
}
