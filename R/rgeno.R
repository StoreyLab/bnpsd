#' Draw genotypes from the admixture model
#'
#' Given the Individual-specific Allele Frequency (IAF) \eqn{\pi_{ij}} for locus \eqn{i} and individual \eqn{j}, genotypes are drawn binomially:
#' \deqn{x_{ij}|\pi_{ij} \sim \mbox{Binomial}(2, \pi_{ij}).}
#' Below \eqn{m} is the number of loci, \eqn{n} the number of individuals, and \eqn{k} the number of intermediate subpopulations.
#' If an admixture proportion matrix \eqn{Q} is provided as the second argument, the first argument \eqn{P} is treated as the intermediate subpopulation allele frequency matrix and the IAF matrix is given by
#' \deqn{P Q^T.}{P \%*\% t(Q).}
#' If \eqn{Q} is missing, then \eqn{P} is treated as the IAF matrix.
#' 
#' To reduce memory, set \code{lowMem=TRUE} to draw genotypes one locus at the time from \eqn{P} and \eqn{Q} (both must be present).
#' This low-memory algorithm prevents the construction of the entire IAF matrix, but is considerably slower than the standard algorithm.
#'
#' @param P The \eqn{m \times n}{m-by-n} IAF matrix (if \eqn{Q} is missing) or the \eqn{m \times k}{m-by-k} intermediate subpopulation allele frequency matrix (if \eqn{Q} is present)
#' @param Q The optional \eqn{n \times k}{n-by-k} admixture proportion matrix
#' @param lowMem If \code{TRUE}, the low-memory algorithm is used (Q must be present)
#'
#' @return The \eqn{m \times n}{m-by-n} genotype matrix
#'
#' @examples
#' # dimensions
#' m <- 10 # number of loci
#' n <- 5 # number of individuals
#' k <- 2 # number of intermediate subpops
#'
#' # define population structure
#' F <- c(0.1, 0.3) # FST values for k=2 subpopulations
#' sigma <- 1 # dispersion parameter of intermediate subpops
#' Q <- q1d(n, k, sigma) # non-trivial admixture proportions
#'
#' # draw alelle frequencies
#' pAnc <- rpanc(m) # random vector of ancestral allele frequencies
#' B <- rpint(pAnc, F) # matrix of intermediate subpop allele freqs
#' P <- rpiaf(B,Q) # matrix of individual-specific allele frequencies
#'
#' # draw genotypes from intermediate populations (one individual each)
#' Xb <- rgeno(B)
#' # and genotypes for admixed individuals
#' Xp <- rgeno(P)
#'
#' @export
rgeno <- function(P, Q=NULL, lowMem=FALSE) {
    
    ## let's determine the data dimensions, validate inputs
    m <- nrow(P) # number of SNPs in all cases!
    if (is.null(Q)) {
        n <- ncol(P) # easy case, it's just P and nothing is really validated
    } else {
        n <- nrow(Q)
        ## these values must match or we can't multiply P and Q
        if (ncol(P) != ncol(Q)) stop('Fatal: P and Q are not compatible: ncol(P) == ', ncol(P), ' != ', ncol(Q), ' == ncol(Q)')
    }
    
    ## calculations
    if (lowMem) {
        if (is.null(Q)) stop('Fatal: Q matrix must be provided to rgeno if lowMem=TRUE')
        
        ## 2016-06-13:
        ## This version is designed to handle larger numbers of individuals by avoiding having P*Q in memory the whole time
        ## In fact, we don't save P*Q either!
        ## Since P and Q are the percursors and may be fairly low-memory (for small k), we can draw genotypes more slowly, making the IAF vector one SNP at the time.
        ## For some simulations this is acceptable, if we don't care about the IAFs P*Q (we can still construct it outside given B and Q!)
        X <- matrix(0, nrow=m, ncol=n) # initialize the big matrix!
        ## navigate each SNP...
        for (i in 1:m) {
            Pi <- tcrossprod(P[i,], Q) # create IAF vector at this SNP only
            X[i,] <- stats::rbinom(n, 2, Pi) # draw genotypes and store!
        }
    } else {
        if (!is.null(Q)) {
            ## in this case we multiply P and Q to get the IAF matrix we want
            P <- rpiaf(P, Q) # overwrite P in this case
        }
        ## now P is the IAF matrix (either way), draw binomially from it!
        X <- matrix(stats::rbinom(P, 2, P), nrow=m, ncol=n)
    }

    ## return for all cases
    X
}
