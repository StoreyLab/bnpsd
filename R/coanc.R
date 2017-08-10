#' Construct the coancestry matrix of an admixture model
#'
#' In the most general case, the \eqn{n \times n}{n-by-n} coancestry matrix \eqn{\Theta} of admixed individuals is determined by the \eqn{n \times k}{n-by-k} admixture proportion matrix \eqn{Q} and the \eqn{k \times k}{k-by-k} intermediate subpopulation coancestry matrix \eqn{\Psi}, given by
#' \deqn{\Theta = Q \Psi Q^T}{\Theta = Q * \Psi * Q^T}
#' In the BN-PSD model \eqn{\Psi} is a diagonal matrix (with \eqn{F_{ST}}{FST} values for the intermediate subpopulations along the diagonal, zero values off-diagonal).
#'
#' @param Q The \eqn{n \times k}{n-by-k} admixture proportion matrix
#' @param F Either the \eqn{k \times k}{k-by-k} intermediate subpopulation coancestry matrix (for the complete admixture model), or the length-\eqn{k} vector of intermediate subpopulation \eqn{F_{ST}}{FST} values (for the BN-PSD model), or a scalar \eqn{F_{ST}}{FST} value shared by all intermediate subpopulations.
#'
#' @return The \eqn{n \times n}{n-by-n} coancestry matrix \eqn{\Theta}
#'
#' @examples
#' # a trivial case: unadmixed individuals from independent subpopulations
#' n <- 5 # number of individuals/subpops
#' Q <- diag(rep.int(1, n)) # unadmixed individuals
#' F <- 0.2 # equal Fst for all subpops
#' Theta <- coanc(Q, F) # diagonal coancestry matryx
#'
#' # a more complicated admixture model
#' n <- 5 # number of individuals
#' k <- 2 # number of intermediate subpops
#' sigma <- 1 # dispersion parameter of intermediate subpops
#' Q <- q1d(n, k, sigma) # non-trivial admixture proportions
#' F <- c(0.1, 0.3) # different Fst for each of the k subpops
#' Theta <- coanc(Q, F) # non-trivial coancestry matrix
#' 
#' @export 
coanc <- function(Q, F) {
    k <- ncol(Q) # dimension that matters the most
    if (is.matrix(F)) {
        if (nrow(F) != k) stop('Fatal: Q and F are not compatible: nrow(F) == ', nrow(F), ' != ', k, ' == ncol(Q)')
        if (ncol(F) != k) stop('Fatal: Q and F are not compatible: ncol(F) == ', ncol(F), ' != ', k, ' == ncol(Q)')
        tcrossprod(Q %*% F, Q) # most general version, F is an arbitrary matrix of coancestries!
    } else if (length(F) == k) {
        tcrossprod(Q %*% diag(F), Q) # version for vector F of non-equal values (intermediate generality)
    } else if (length(F) == 1) {
        tcrossprod(Q) * F # this is the matrix we want
    } else stop('Fatal: Q and F are not compatible: length(F) = ', length(F), ' != ', k, ' == ncol(Q)')
}

