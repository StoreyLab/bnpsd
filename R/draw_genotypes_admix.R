#' Draw genotypes from the admixture model
#'
#' Given the Individual-specific Allele Frequency (IAF) \eqn{\pi_{ij}} for locus \eqn{i} and individual \eqn{j}, genotypes are drawn binomially:
#' \deqn{x_{ij}|\pi_{ij} \sim \mbox{Binomial}(2, \pi_{ij}).}
#' Below \eqn{m} is the number of loci, \eqn{n} the number of individuals, and \eqn{k} the number of intermediate subpopulations.
#' If an admixture proportion matrix \eqn{Q} is provided as the second argument, the first argument \eqn{P} is treated as the intermediate subpopulation allele frequency matrix and the IAF matrix is given by
#' \deqn{P Q^T.}{P \%*\% t(Q).}
#' If \eqn{Q} is missing, then \eqn{P} is treated as the IAF matrix.
#' 
#' To reduce memory, set \code{low_mem = TRUE} to draw genotypes one locus at the time from \eqn{P} and \eqn{Q} (both must be present).
#' This low-memory algorithm prevents the construction of the entire IAF matrix, but is considerably slower than the standard algorithm.
#'
#' @param p_ind The \eqn{m \times n}{m-by-n} IAF matrix (if \code{admix_proportions} is missing) or the \eqn{m \times k}{m-by-k} intermediate subpopulation allele frequency matrix (if \code{admix_proportions} is present)
#' @param admix_proportions The optional \eqn{n \times k}{n-by-k} admixture proportion matrix
#' @param low_mem If \code{TRUE}, the low-memory algorithm is used (\code{admix_proportions} must be present)
#'
#' @return The \eqn{m \times n}{m-by-n} genotype matrix
#'
#' @examples
#' # dimensions
#' # number of loci
#' m_loci <- 10
#' # number of individuals
#' n_ind <- 5
#' # number of intermediate subpops
#' k_subpops <- 2
#'
#' # define population structure
#' # FST values for k = 2 subpops
#' inbr_subpops <- c(0.1, 0.3)
#' # non-trivial admixture proportions
#' admix_proportions <- admix_prop_1d_linear(n_ind, k_subpops, sigma = 1)
#'
#' # draw allele frequencies
#' # vector of ancestral allele frequencies
#' p_anc <- draw_p_anc(m_loci)
#' 
#' # matrix of intermediate subpop allele freqs
#' p_subpops <- draw_p_subpops(p_anc, inbr_subpops)
#' 
#' # matrix of individual-specific allele frequencies
#' p_ind <- make_p_ind_admix(p_subpops, admix_proportions)
#'
#' # draw genotypes from intermediate subpops (one individual each)
#' X_subpops <- draw_genotypes_admix(p_subpops)
#' 
#' # and genotypes for admixed individuals
#' X_ind <- draw_genotypes_admix(p_ind)
#' 
#' # draw genotypes for admixed individuals without p_ind intermediate
#' # (p_ind is computer internally and discarded when done)
#' X_ind <- draw_genotypes_admix(p_subpops, admix_proportions)
#' 
#' # use low-memory version (p_ind is computed by row, never fully in memory)
#' X_ind <- draw_genotypes_admix(p_subpops, admix_proportions, low_mem = TRUE)
#'
#' @export
draw_genotypes_admix <- function(p_ind, admix_proportions = NULL, low_mem = FALSE) {
    # die if this is missing
    if (missing(p_ind))
        stop('`p_ind` is required!')

    # ensure that things that should be matrices are so
    if (!is.matrix(p_ind))
        stop('`p_ind` must be a matrix!')

    # let's determine the data dimensions, validate inputs
    m <- nrow(p_ind) # number of SNPs in all cases!
    if (is.null(admix_proportions)) {
        n <- ncol(p_ind) # easy case, it's just p_ind and nothing is really validated
    } else {
        # ensure that things that should be matrices are so
        if (!is.matrix(admix_proportions))
            stop('`admix_proportions` must be a matrix!')
        
        n <- nrow(admix_proportions)
        # these values must match or we can't multiply p_ind and admix_proportions
        if (ncol(p_ind) != ncol(admix_proportions))
            stop('`p_ind` and `admix_proportions` are not compatible: ncol(p_ind) == ', ncol(p_ind), ' != ', ncol(admix_proportions), ' == ncol(admix_proportions)')
    }
    
    # calculations
    if (low_mem) {
        if (is.null(admix_proportions))
            stop('`admix_proportions` must be provided if `low_mem = TRUE`')
        
        # This version handles larger numbers of individuals by avoiding having p_ind*admix_proportions in memory the whole time
        # In fact, we don't save p_ind*admix_proportions either!
        # Since p_ind and admix_proportions are the percursors and may be fairly low-memory (for small k), we can draw genotypes more slowly, making the IAF vector one SNP at the time.
        # For some simulations this is acceptable, if we don't care about the IAFs p_ind*admix_proportions (we can still construct it outside given B and admix_proportions!)
        X <- matrix(0, nrow = m, ncol = n) # initialize the big matrix!
        # navigate each SNP...
        for (i in 1:m) {
            p_indi <- tcrossprod(p_ind[i,], admix_proportions) # create IAF vector at this SNP only
            X[i,] <- stats::rbinom(n, 2, p_indi) # draw genotypes and store!
        }
    } else {
        if (!is.null(admix_proportions)) {
            # in this case we multiply p_ind and admix_proportions to get the IAF matrix we want
            p_ind <- make_p_ind_admix(p_ind, admix_proportions) # overwrite p_ind in this case
        }
        # now p_ind is the IAF matrix (either way), draw binomially from it!
        # it appears that rbinom reads p_ind by column, and output is turned back into matrix also by column, therefore retaining consistency!
        X <- matrix(stats::rbinom(p_ind, 2, p_ind), nrow = m, ncol = n)
    }

    # return for all cases
    X
}

# stick deprecated function name here

#' @title Draw genotypes from the admixture model
#' @description Draw genotypes from the admixture model
#' @param P The IAF matrix (if Q is missing) or the intermediate subpopulation allele frequency matrix (if Q is present)
#' @param Q The optional admixture proportion matrix
#' @param lowMem If \code{TRUE}, the low-memory algorithm is used (Q must be present)
#' @return The genotype matrix
#'
#' @name rgeno-deprecated
#' @usage rgeno(P, Q = NULL, lowMem = FALSE)
#' @seealso \code{\link{bnpsd-deprecated}}
#' @keywords internal
NULL

#' @rdname bnpsd-deprecated
#' @section \code{rgeno}:
#' For \code{rgeno}, use \code{\link{draw_genotypes_admix}}.
#'
#' @export
rgeno <- function(P, Q = NULL, lowMem = FALSE) {
    # mark as deprecated
    .Deprecated('draw_genotypes_admix')
    # return as usual, to not break things just yet
    draw_genotypes_admix(P, Q, lowMem)
}
