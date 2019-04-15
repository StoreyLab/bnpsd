#' Draw genotypes from the admixture model
#'
#' Given the Individual-specific Allele Frequency (IAF) \eqn{\pi_{ij}} for locus \eqn{i} and individual \eqn{j}, genotypes are drawn binomially:
#' \deqn{x_{ij}|\pi_{ij} \sim \mbox{Binomial}(2, \pi_{ij}).}
#' Below \eqn{m} is the number of loci, \eqn{n} the number of individuals, and \eqn{k} the number of intermediate subpopulations.
#' If an admixture proportion matrix \eqn{Q} is provided as the second argument, the first argument \eqn{P} is treated as the intermediate subpopulation allele frequency matrix and the IAF matrix is given by
#' \deqn{P Q^T.}{P \%*\% t(Q).}
#' If \eqn{Q} is missing, then \eqn{P} is treated as the IAF matrix.
#' 
#' To reduce memory, set \code{lowMem = TRUE} to draw genotypes one locus at the time from \eqn{P} and \eqn{Q} (both must be present).
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
#' m_loci <- 10 # number of loci
#' n_ind <- 5 # number of individuals
#' k_subpops <- 2 # number of intermediate subpops
#'
#' # define population structure
#' # FST values for k = 2 subpops
#' inbr_subpops <- c(0.1, 0.3)
#' # non-trivial admixture proportions
#' admix_proportions <- admix_prop_1d_linear(n_ind, k_subpops, sigma = 1)
#'
#' # draw alelle frequencies
#' # vector of ancestral allele frequencies
#' p_anc <- draw_p_anc(m_loci)
#' # matrix of intermediate subpop allele freqs
#' p_subpops <- draw_p_subpops(p_anc, inbr_subpops)
#' # matrix of individual-specific allele frequencies
#' p_ind <- make_p_ind_admix(p_subpops, admix_proportions)
#'
#' # draw genotypes from intermediate subpops (one individual each)
#' X_subpops <- rgeno(p_subpops)
#' # and genotypes for admixed individuals
#' X_ind <- rgeno(p_ind)
#'
#' @export
rgeno <- function(P, Q = NULL, lowMem = FALSE) {
    
    # let's determine the data dimensions, validate inputs
    m <- nrow(P) # number of SNPs in all cases!
    if (is.null(Q)) {
        n <- ncol(P) # easy case, it's just P and nothing is really validated
    } else {
        n <- nrow(Q)
        # these values must match or we can't multiply P and Q
        if (ncol(P) != ncol(Q))
            stop('P and Q are not compatible: ncol(P) == ', ncol(P), ' != ', ncol(Q), ' == ncol(Q)')
    }
    
    # calculations
    if (lowMem) {
        if (is.null(Q))
            stop('Q matrix must be provided to rgeno if `lowMem = TRUE`')
        
        # This version handles larger numbers of individuals by avoiding having P*Q in memory the whole time
        # In fact, we don't save P*Q either!
        # Since P and Q are the percursors and may be fairly low-memory (for small k), we can draw genotypes more slowly, making the IAF vector one SNP at the time.
        # For some simulations this is acceptable, if we don't care about the IAFs P*Q (we can still construct it outside given B and Q!)
        X <- matrix(0, nrow = m, ncol = n) # initialize the big matrix!
        # navigate each SNP...
        for (i in 1:m) {
            Pi <- tcrossprod(P[i,], Q) # create IAF vector at this SNP only
            X[i,] <- stats::rbinom(n, 2, Pi) # draw genotypes and store!
        }
    } else {
        if (!is.null(Q)) {
            # in this case we multiply P and Q to get the IAF matrix we want
            P <- make_p_ind_admix(P, Q) # overwrite P in this case
        }
        # now P is the IAF matrix (either way), draw binomially from it!
        # it appears that rbinom reads P by column, and output is turned back into matrix also by column, therefore retaining consistency!
        X <- matrix(stats::rbinom(P, 2, P), nrow = m, ncol = n)
    }

    # return for all cases
    X
}
