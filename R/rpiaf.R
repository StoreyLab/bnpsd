## PUBLIC?: used only by ~/storey/fst/software/ParseTeraStructure.R

#' Construct individual-specific allele frequency matrix
#'
#' Here \eqn{m} is the number of loci, \eqn{n} the number of individuals, and \eqn{k} the number of intermediate populations.
#' The \eqn{m \times n}{m-by-n} Individual-specific Allele Frequency (IAF) matrix \eqn{P} is constructed from the \eqn{m \times k}{m-by-k} intermediate population allele frequency matrix \eqn{B} and the \eqn{n \times k}{n-by-k} admixture proportion matrix \eqn{Q} using
#' \deqn{P = B Q^T.}{P = B * Q^T.}
#' 
#' @param B The \eqn{m \times k}{m-by-k} intermediate population allele frequency matrix
#' @param Q The \eqn{n \times k}{n-by-k} admixture proportion matrix
#'
#' @return The \eqn{m \times n}{m-by-n} IAF matrix \eqn{P}
#'
#' @examples
#' m <- 10 # number of loci
#' n <- 5 # number of individuals
#' k <- 2 # number of intermediate subpops
#' pAnc <- rpanc(m) # random vector of ancestral allele frequencies
#' F <- c(0.1, 0.3) # FST values for k=2 subpopulations
#' B <- rpint(pAnc, F) # matrix of intermediate subpop allele freqs
#' sigma <- 1 # dispersion parameter of intermediate subpops
#' Q <- q1d(n, k, sigma) # non-trivial admixture proportions
#' P <- rpiaf(B,Q)
#'
#' @export
rpiaf <- function(B,Q) {
    ## constructs the individual allele frequency matrix implied by PSD (given ancestral allele frequencies and admixture coefficients, both of which may be estimates)
    ## it's a simple matrix product, but I always forget the matrix orientations, here I just get it right once and remember forever...
    ## Q is (n,k), B is (m,k), P is (m,n)
    ## so want P <- B %*% t(Q), below does this more efficiently

    ## validate data dimensions, to have more reasonable error messages
    ## these values must match or we can't multiply B and Q
    if (ncol(B) != ncol(Q)) stop('Fatal: B and Q are not compatible: ncol(B) == ', ncol(B), ' != ', ncol(Q), ' == ncol(Q)')

    ## this is the main multiplication
    P <- tcrossprod(B, Q)
    
    ## unfortunately, in very extreme cases, P may have values slighly outside of [0,1] simply due to machine precision errors, fix that here!
    P[P<0] <- 0
    P[P>1] <- 1
    P # return!
}
