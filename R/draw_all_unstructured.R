#' Simulate random allele frequencies and genotypes from the unstructured model
#'
#' This function returns simulated ancestral allele frequencies and genotypes without structure, meaning individuals draw their genotypes independently and identically from the Binomial distribution with the same ancestral allele frequency per locus.
#' The function is a wrapper around [draw_p_anc()] with additional features such as requiring polymorphic loci, mimicking [draw_all_admix()] in options as applicable.
#' Importantly, by default fixed loci (where all individuals were homozygous for the same allele) are re-drawn from the start (starting from the ancestral allele frequencies) so no fixed loci are in the output.
#' Below `m_loci` (also `m`) is the number of loci and `n_ind` is the number of individuals.
#'
#' @param n_ind The number of individuals to draw (required).
#' @param m_loci The number of loci to draw.  Required except when `p_anc` below is provided and is a vector, in which case the number of loci equals the length of `p_anc` and the value of `m_loci` passed is ignored.
#' @param beta Shape parameter for a symmetric Beta for ancestral allele frequencies `p_anc`.
#' If `NA` (default), `p_anc` is uniform with range in \[0.01, 0.5\].
#' Otherwise, `p_anc` has a symmetric Beta distribution with range in \[0, 1\].
#' Has no effect if `p_anc` option is non-`NULL`.
#' @param p_anc If provided, it is used as the ancestral allele frequencies (instead of drawing random ones).  Must either be a scalar or a length-`m_loci` vector.
#' If scalar, `m_loci` is required, and the returned `p_anc` is the scalar value repeated `m_loci` times.
#' If `p_anc` is a vector, its length is used to define `m_loci` and the value of `m_loci` passed is ignored.
#' If a locus was fixed and has to be redrawn, the ancestral allele frequency in `p_anc` is retained and only genotypes are redrawn.
#' @param require_polymorphic_loci If `TRUE` (default), returned genotype matrix will not include any fixed loci (loci that happened to be fixed are drawn again, starting from their ancestral allele frequencies, and checked iteratively until no fixed loci remain, so that the final number of polymorphic loci is exactly `m_loci`).
#' @param maf_min The minimum minor allele frequency (default zero), to extend the working definition of "fixed" above to include rare variants.
#' This helps simulate a frequency-based locus ascertainment bias.
#' Loci with minor allele frequencies less than or equal to this value are treated as fixed (passed to [fixed_loci()]).
#' This parameter has no effect if `require_polymorphic_loci` is `FALSE`.
#' @param verbose If `TRUE`, prints messages for every stage in the algorithm.
#'
#' @return A named list with the following items (which may be missing depending on options):
#'
#' - `X`: An `m_loci`-by-`n_ind` matrix of genotypes.
#' - `p_anc`: A length-`m_loci` vector of ancestral allele frequencies.
#'
#' @examples
#' # dimensions
#' # number of loci
#' m_loci <- 10
#' # number of individuals
#' n_ind <- 5
#'
#' # draw all random allele freqs and genotypes
#' out <- draw_all_unstructured( n_ind, m_loci )
#'
#' # return value is a list with these items:
#' 
#' # genotypes
#' X <- out$X
#' 
#' # ancestral AFs
#' p_anc <- out$p_anc
#' 
#' @export
draw_all_unstructured <- function(
                                  n_ind,
                                  m_loci = NA,
                                  beta = NA,
                                  p_anc = NULL,
                                  require_polymorphic_loci = TRUE,
                                  maf_min = 0,
                                  verbose = TRUE
                                  ) {
    # only n_ind is absolutely required
    if ( missing( n_ind ) )
        stop( '`n_ind` is required!' )
    
    # validate or generate p_anc
    p_anc_in <- p_anc # remember its original value, for later
    if ( !is.null( p_anc ) ) {
        if ( length( p_anc ) == 1 ) {
            # here we require m_loci as well
            if ( is.na( m_loci ) )
                stop( '`m_loci` must be non-NA when `p_anc` has length 1!' )
            # expand scalar into vector, so it looks as desired
            p_anc <- rep.int( p_anc, m_loci )
        } else {
            # use given length as number of loci, potentially overwriting whatever was passed
            m_loci <- length( p_anc )
        }
        
        # validate range
        if ( any( p_anc < 0 ) )
            stop('Provided `p_anc` has negative values!')
        if ( any( p_anc > 1 ) )
            stop('Provided `p_anc` has values exceeding 1!')
    } else {
        if ( is.na( m_loci ) )
            stop( '`m_loci` must be non-NA when `p_anc` is NULL!' )
        # generate the random ancestral allele frequencies, in usual range and with minimum threshold for simplicity
        if (verbose)
            message('drawing p_anc')
        p_anc <- draw_p_anc( m_loci, beta = beta )
    }
    
    # draw genotypes
    if (verbose)
        message('drawing X')
    X <- matrix(
        stats::rbinom( n_ind * m_loci, 2, p_anc ),
        nrow = m_loci,
        ncol = n_ind
    )
    
    if ( require_polymorphic_loci ) {
        # check for fixed loci, draw them again if needed
        # Note this only applies to want_genotypes==TRUE, since p_ind and p_subpops are continuous and therefore practically never truly fixed
        fixed_loci_indexes <- fixed_loci( X, maf_min = maf_min ) # boolean vector identifies fixed loci
        m_loci_fixed <- sum(fixed_loci_indexes) # number of cases
        while ( m_loci_fixed > 0 ) {
            # p_anc is tricky here
            # this is what we'll pass below
            p_anc_redo <- p_anc_in
            # if null or length 1, nothing changes, so look at rest
            if ( !is.null( p_anc_in ) && length( p_anc_in ) > 1 ) {
                # only pass cases that are getting redone, must be a vector of the right length
                p_anc_redo <- p_anc_in[ fixed_loci_indexes ]
            }
            # call self with desired number of loci, all the same parameters otherwise
            if (verbose)
                message('re-drawing fixed loci')
            obj <- draw_all_unstructured(
                n_ind = n_ind,
                m_loci = m_loci_fixed,
                beta = beta,
                p_anc = p_anc_redo,
                require_polymorphic_loci = FALSE, # loop here, avoid recursion in function, which in extreme cases leads to "node stack overflow" errors!
                maf_min = maf_min,
                verbose = FALSE # don't show more messages for additional iterations
            )
            # overwrite fixed loci with redrawn polymorphic data
            X[ fixed_loci_indexes, ] <- obj$X
            p_anc[ fixed_loci_indexes ] <- obj$p_anc
            
            # for next iteration, look for cases that remain fixed
            fixed_loci_indexes <- fixed_loci( X, maf_min = maf_min ) # boolean vector identifies fixed loci
            m_loci_fixed <- sum(fixed_loci_indexes) # number of cases
        }
    }

    # inherit names of p_anc, if any
    rownames( X ) <- names( p_anc )
    # there are never column names, individuals are IID because they're unstructured!
    
    # return desired objects
    return(
        list(
            X = X,
            p_anc = p_anc
        )
    )
}
