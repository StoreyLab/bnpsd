#' Draw genotypes from the admixture model
#'
#' Given the Individual-specific Allele Frequency (IAF) \eqn{\pi_{ij}} for locus \eqn{i} and individual \eqn{j}, genotypes are drawn binomially:
#' \deqn{x_{ij}|\pi_{ij} \sim \mbox{Binomial}(2, \pi_{ij}).}
#' Below \eqn{m} is the number of loci, \eqn{n} the number of individuals, and \eqn{k} the number of intermediate subpopulations.
#' If an admixture proportion matrix \eqn{Q} is provided as the second argument, the first argument \eqn{P} is treated as the intermediate subpopulation allele frequency matrix and the IAF matrix is given by
#' \deqn{P Q^T.}{P \%*\% t(Q).}
#' However, in this case the IAF matrix is computed in parts only, never stored in full, greatly reducing memory usage.
#' If \eqn{Q} is missing, then \eqn{P} is treated as the IAF matrix.
#' 
#' @param p_ind The \eqn{m \times n}{m-by-n} IAF matrix (if \code{admix_proportions} is missing) or the \eqn{m \times k}{m-by-k} intermediate subpopulation allele frequency matrix (if \code{admix_proportions} is present)
#' @param admix_proportions The optional \eqn{n \times k}{n-by-k} admixture proportion matrix
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
#' # (p_ind is computed internally in parts, never stored in full,
#' # reducing memory use substantially)
#' X_ind <- draw_genotypes_admix(p_subpops, admix_proportions)
#' 
#' @export
draw_genotypes_admix <- function(p_ind, admix_proportions = NULL) {
    # die if this is missing
    if ( missing( p_ind ) )
        stop('`p_ind` is required!')

    # ensure that things that should be matrices are so
    if ( !is.matrix( p_ind ) )
        stop('`p_ind` must be a matrix!')

    # number of loci in all cases!
    m_loci <- nrow( p_ind )
    
    # calculations and additional validations
    if ( is.null( admix_proportions ) ) {
        
        # easy case, it's just p_ind and nothing is really validated
        n_ind <- ncol( p_ind )
        
        # p_ind is the IAF matrix, draw binomially from it!
        # it appears that rbinom reads p_ind by column, and output is turned back into matrix also by column, therefore retaining consistency!
        X <- matrix(
            stats::rbinom(p_ind, 2, p_ind),
            nrow = m_loci,
            ncol = n_ind
        )

    } else {
        
        # ensure that things that should be matrices are so
        if ( !is.matrix( admix_proportions ) )
            stop('`admix_proportions` must be a matrix!')
        
        n_ind <- nrow( admix_proportions )
        # these values must match or we can't multiply p_ind and admix_proportions
        if ( ncol( p_ind ) != ncol( admix_proportions ) )
            stop('`p_ind` and `admix_proportions` are not compatible: ncol( p_ind ) == ', ncol( p_ind ), ' != ', ncol( admix_proportions ), ' == ncol( admix_proportions )')
        
        # This version handles larger numbers of individuals by avoiding having p_ind * admix_proportions in memory the whole time
        # In fact, we don't save p_ind * admix_proportions either!
        # Since p_ind and admix_proportions are the percursors and may be fairly low-memory (for small k), we can draw genotypes more slowly, making the IAF vector one SNP at the time.
        # For some simulations this is acceptable, if we don't care about the IAFs p_ind*admix_proportions (we can still construct it outside given B and admix_proportions!)
        X <- matrix(0, nrow = m_loci, ncol = n_ind) # initialize the big matrix!
        # fill by whichever is the smallest dimension (vectorizes the most)
        if ( n_ind <= m_loci ) {
            # navigate each individual
            for ( j in 1 : n_ind ) {
                # create IAF vector for this individual only
                # admix_proportions[ j, ] is length (k)
                # p_ind is (m, k)
                # p_ind_j is length (m)
                p_ind_j <- drop( p_ind %*% admix_proportions[ j, ] )
                # draw genotypes and store!
                X[ , j ] <- stats::rbinom( m_loci, 2, p_ind_j )
            }
        } else {
            # navigate each locus
            for ( i in 1 : m_loci ) {
                # create IAF vector at this locus only
                # p_ind[ i, ] is length (k)
                # admix_proportions is (n, k)
                # p_ind_i is length (n)
                p_ind_i <- drop( admix_proportions %*% p_ind[ i, ] )
                # draw genotypes and store!
                X[ i, ] <- stats::rbinom( n_ind, 2, p_ind_i )
            }
        }
        
    }

    # return for all cases
    X
}

# stick deprecated function name here

#' @title Draw genotypes from the admixture model
#' @description Draw genotypes from the admixture model
#' @param P The IAF matrix (if Q is missing) or the intermediate subpopulation allele frequency matrix (if Q is present)
#' @param Q The optional admixture proportion matrix
#' @param lowMem OBSOLETE (the default algorithm is now the old lowMem=TRUE)
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
    draw_genotypes_admix(P, Q)
}
