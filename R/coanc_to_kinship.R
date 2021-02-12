#' Transform coancestry matrix to kinship matrix
#'
#' If `Theta` is the coancestry matrix and `Phi` is the kinship matrix (both are `n`-by-`n` symmetric), then these matrices agree off-diagonal, but the diagonal gets transformed as
#' `diag( Phi ) = ( 1 + diag( Theta ) ) / 2`.
#'
#' @param coancestry The `n`-by-`n` coancestry matrix
#'
#' @return The `n`-by-`n` kinship matrix, preserving column and row names.
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
#' The inverse function is given by `\link[popkin]{inbr_diag}`.
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
