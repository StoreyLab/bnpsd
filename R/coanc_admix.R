#' Construct the coancestry matrix of an admixture model
#'
#' In the most general case, the \eqn{n \times n}{n-by-n} coancestry matrix \eqn{\Theta} of admixed individuals is determined by the \eqn{n \times k}{n-by-k} admixture proportion matrix \eqn{Q} and the \eqn{k \times k}{k-by-k} intermediate subpopulation coancestry matrix \eqn{\Psi}, given by
#' \deqn{\Theta = Q \Psi Q^T}{\Theta = Q * \Psi * Q^T}
#' In the BN-PSD model \eqn{\Psi} is a diagonal matrix (with \eqn{F_{ST}}{FST} values for the intermediate subpopulations along the diagonal, zero values off-diagonal).
#'
#' @param admix_proportions The \eqn{n \times k}{n-by-k} admixture proportion matrix
#' @param coanc_subpops Either the \eqn{k \times k}{k-by-k} intermediate subpopulation coancestry matrix (for the complete admixture model), or the length-\eqn{k} vector of intermediate subpopulation \eqn{F_{ST}}{FST} values (for the BN-PSD model), or a scalar \eqn{F_{ST}}{FST} value shared by all intermediate subpopulations.
#'
#' @return The \eqn{n \times n}{n-by-n} coancestry matrix.
#'
#' @examples
#' # a trivial case: unadmixed individuals from independent subpopulations
#' n <- 5 # number of individuals and subpops
#' admix_proportions <- diag(rep.int(1, n)) # unadmixed individuals
#' coanc_subpops <- 0.2 # equal Fst for all subpops
#' coancestry <- coanc(admix_proportions, coanc_subpops) # diagonal coancestry matryx
#'
#' # a more complicated admixture model
#' n <- 5 # number of individuals
#' k <- 2 # number of intermediate subpops
#' sigma <- 1 # dispersion parameter of intermediate subpops
#' admix_proportions <- q1d(n, k, sigma) # non-trivial admixture proportions
#' coanc_subpops <- c(0.1, 0.3) # different Fst for each of the k subpops
#' coancestry <- coanc(admix_proportions, coanc_subpops) # non-trivial coancestry matrix
#' 
#' @export 
coanc_admix <- function(admix_proportions, coanc_subpops) {
    # die if things are missing outright
    if (missing( admix_proportions ))
        stop('`admix_proportions` is required!')
    if (missing( coanc_subpops ))
        stop('`coanc_subpops` is required!')

    # behavior depends on dimensions of coanc_subpops:
    k <- ncol(admix_proportions) # dimension that matters the most
    if (is.matrix(coanc_subpops)) {
        
        # case 1 - coanc_subpops square matrix
        
        # check dimensions
        if (nrow(coanc_subpops) != k)
            stop('`admix_proportions` and `coanc_subpops` are not compatible: nrow(coanc_subpops) == ', nrow(coanc_subpops), ' != ', k, ' == ncol(admix_proportions)')
        if (ncol(coanc_subpops) != k)
            stop('`admix_proportions` and `coanc_subpops` are not compatible: ncol(coanc_subpops) == ', ncol(coanc_subpops), ' != ', k, ' == ncol(admix_proportions)')

        # compute in most general form!
        tcrossprod(admix_proportions %*% coanc_subpops, admix_proportions)
        
    } else if (length(coanc_subpops) == k) {
        
        # case 2 - coanc_subpops vector
        tcrossprod(admix_proportions %*% diag(coanc_subpops), admix_proportions)
        
    } else if (length(coanc_subpops) == 1) {
        
        # case 3 - coanc_subpops scalar
        tcrossprod(admix_proportions) * coanc_subpops
        
    } else
        stop('`admix_proportions` and `coanc_subpops` are not compatible: length(coanc_subpops) = ', length(coanc_subpops), ' != ', k, ' == ncol(admix_proportions)')
}

# stick deprecated function name here

#' @title Construct the coancestry matrix of an admixture model
#' @description Construct the coancestry matrix of an admixture model
#' @param Q The admixture proportion matrix
#' @param F Either the intermediate subpopulation coancestry matrix, or the vector of intermediate subpopulation \eqn{F_{ST}}{FST} values, or a scalar \eqn{F_{ST}}{FST} value shared by all intermediate subpopulations.
#' @return The coancestry matrix
#'
#' @name coanc-deprecated
#' @usage coanc(Q, F)
#' @seealso \code{\link{bnpsd-deprecated}}
#' @keywords internal
NULL

#' @rdname bnpsd-deprecated
#' @section \code{coanc}:
#' For \code{coanc}, use \code{\link{coanc_admix}}.
#'
#' @export
coanc <- function(Q, F) {
    # mark as deprecated
    .Deprecated('coanc_admix')
    # return as usual, to not break things just yet
    coanc_admix(Q, F)
}
