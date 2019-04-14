#' Transform coancestry matrix to kinship matrix
#'
#' Let
#' \eqn{\Theta = (\theta_{jk})} be the coancestry matrix and
#' \eqn{\Phi = (\varphi_{jk})} be the kinship matrix.
#' These matrices agree off-diagonal, but the diagonal gets transformed as
#' \deqn{\phi_{jj} = \frac{1 + \theta_{jj}}{2}.}
#' Below \eqn{n} is the number of individuals.
#'
#' @param coancestry The \eqn{n \times n}{n-by-n} coancestry matrix
#'
#' @return The \eqn{n \times n}{n-by-n} kinship matrix, preserving column and row names.
#'
#' @examples
#' # a trivial case: unadmixed individuals from independent subpopulations
#' n <- 5 # number of individuals/subpops
#' admix_proportions <- diag(rep.int(1, n)) # unadmixed individuals
#' inbr_subpops <- 0.2 # equal Fst for all subpops
#' coancestry <- coanc_admix(admix_proportions, inbr_subpops) # diagonal coancestry matryx
#' kinship <- coanc_to_kinship(coancestry)
#'
#' # a more complicated admixture model
#' n <- 5 # number of individuals
#' k <- 2 # number of intermediate subpops
#' sigma <- 1 # dispersion parameter of intermediate subpops
#' admix_proportions <- admix_prop_1d_linear(n, k, sigma) # non-trivial admixture proportions
#' inbr_subpops <- c(0.1, 0.3) # different Fst for each of the k subpops
#' coancestry <- coanc_admix(admix_proportions, inbr_subpops) # non-trivial coancestry matrix
#' kinship <- coanc_to_kinship( coancestry )
#'
#' @seealso
#' The inverse function is given by \code{\link[popkin]{inbr_diag}}.
#' 
#' @export
coanc_to_kinship <- function(coancestry) {
    # first copy matrix
    kinship <- coancestry
    # transform diagonal (only difference)
    diag(kinship) <- ( diag(kinship) + 1 ) / 2
    # return
    return( kinship )
}
