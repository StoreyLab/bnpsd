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
#' # number of individuals/subpops
#' n_ind <- 5
#' # unadmixed individuals
#' admix_proportions <- diag(rep.int(1, n_ind))
#' # equal Fst for all subpops
#' inbr_subpops <- 0.2
#' # diagonal coancestry matryx
#' coancestry <- coanc_admix(admix_proportions, inbr_subpops)
#' kinship <- coanc_to_kinship(coancestry)
#'
#' # a more complicated admixture model
#' # number of individuals
#' n_ind <- 5
#' # number of intermediate subpops
#' k_subpops <- 2
#' # non-trivial admixture proportions
#' admix_proportions <- admix_prop_1d_linear(n_ind, k_subpops, sigma = 1)
#' # different Fst for each of the k subpops
#' inbr_subpops <- c(0.1, 0.3)
#' # non-trivial coancestry matrix
#' coancestry <- coanc_admix(admix_proportions, inbr_subpops)
#' kinship <- coanc_to_kinship( coancestry )
#'
#' @seealso
#' The inverse function is given by \code{\link[popkin]{inbr_diag}}.
#' 
#' @export
coanc_to_kinship <- function(coancestry) {
    # die if this is missing
    if (missing(coancestry))
        stop('`coancestry` is required!')
    # ensure that things that should be matrices are so
    if (!is.matrix(coancestry))
        stop('`coancestry` must be a matrix!')
    # make sure it is numeric
    if (!is.numeric(coancestry))
        stop('`coancestry` must be numeric!')
    # check dimensions
    m <- nrow(coancestry)
    n <- ncol(coancestry)
    if (n != m)
        stop('`coancestry` must be a square matrix!  (nrow ', m, ' != ncol ', n, ')')
    
    # first copy matrix
    kinship <- coancestry
    # transform diagonal (only difference)
    diag(kinship) <- ( diag(kinship) + 1 ) / 2
    # return
    return( kinship )
}
