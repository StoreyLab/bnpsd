#' Transform coancestry matrix to kinship matrix
#'
#' Let
#' \eqn{\Theta = (\theta_{jk})} be the coancestry matrix and
#' \eqn{\Phi = (\varphi_{jk})} be the kinship matrix.
#' These matrices agree off-diagonal, but the diagonal gets transformed as
#' \deqn{\phi_{jj} = \fracc{1 + \theta_{jj}}{2}.}
#' Below \eqn{n} is the number of individuals.
#'
#' This function starts by copying the input matrix, so it preserves column and row names.
#'
#' Note that this function is the inverse of popkin::inbrDiag
#' 
#' @param Theta The \eqn{n \times n}{n-by-n} coancestry matrix
#'
#' @return The \eqn{n \times n}{n-by-n} kinship matrix
#'
#' @examples
#' # a trivial case: unadmixed individuals from independent subpopulations
#' n <- 5 # number of individuals/subpops
#' Q <- diag(rep.int(1, n)) # unadmixed individuals
#' F <- 0.2 # equal Fst for all subpops
#' Theta <- coanc(Q, F) # diagonal coancestry matryx
#' Phi <- coanc_to_kinship(Theta)
#'
#' # a more complicated admixture model
#' n <- 5 # number of individuals
#' k <- 2 # number of intermediate subpops
#' sigma <- 1 # dispersion parameter of intermediate subpops
#' Q <- q1d(n, k, sigma) # non-trivial admixture proportions
#' F <- c(0.1, 0.3) # different Fst for each of the k subpops
#' Theta <- coanc(Q, F) # non-trivial coancestry matrix
#' Phi <- coanc_to_kinship(Theta)
#'
#' @export
coanc_to_kinship <- function(Theta) {
    # first copy matrix
    Phi <- Theta
    # transform diagonal (only difference)
    diag(Phi) <- (diag(Phi)+1)/2
    # return
    return( Phi )
}
