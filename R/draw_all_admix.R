#' Simulate random allele frequencies and genotypes from the BN-PSD admixture model
#'
#' This function returns simulated ancestral, intermediate, and individual-specific allele frequencies and genotypes given the admixture structure, as determined by the admixture proportions and the vector or tree of intermediate subpopulation FST values.
#' The function is a wrapper around [draw_p_anc()], [draw_p_subpops()]/[draw_p_subpops_tree()], [make_p_ind_admix()], and [draw_genotypes_admix()] with additional features such as requiring polymorphic loci.
#' Importantly, by default fixed loci (where all individuals were homozygous for the same allele) are re-drawn from the start (starting from the ancestral allele frequencies) so no fixed loci are in the output and no biases are introduced by re-drawing genotypes conditional on any of the previous allele frequencies (ancestral, intermediate, or individual-specific).
#' Below `m_loci` (also `m`) is the number of loci, `n` is the number of individuals, and `k` is the number of intermediate subpopulations.
#'
#' As a precaution, function stops if both the column names of `admix_proportions` and the names in `inbr_subpops` or `tree_subpops` exist and disagree, which might be because these two data are not aligned or there is some other inconsistency.
#'
#' @param admix_proportions The `n`-by-`k` matrix of admixture proportions.
#' @param inbr_subpops The length-`k` vector (or scalar) of intermediate subpopulation FST values.
#' Either this or `tree_subpops` must be provided (but not both).
#' @param m_loci The number of loci to draw.
#' @param tree_subpops The coancestry tree relating the `k` intermediate subpopulations.
#' Must be a `phylo` object from the `ape` package (see [ape::read.tree()]).
#' Either this or `inbr_subpops` must be provided (but not both).
#' @param want_genotypes If `TRUE` (default), includes the matrix of random genotypes in the return list.
#' @param want_p_ind If `TRUE` (NOT default), includes the matrix of individual-specific allele frequencies in the return list.
#' Note that by default `p_ind` is not constructed in full at all, instead a fast low-memory algorithm constructs it in parts as needed only; beware that setting `want_p_ind = TRUE` increases memory usage in comparison.
#' @param want_p_subpops If `TRUE` (NOT default), includes the matrix of random intermediate subpopulation allele frequencies in the return list.
#' @param want_p_anc If `TRUE` (default), includes the vector of random ancestral allele frequencies in the return list.
#' @param verbose If `TRUE`, prints messages for every stage in the algorithm.
#' @param require_polymorphic_loci If `TRUE` (default), returned genotype matrix will not include any fixed loci (loci that happened to be fixed are drawn again, starting from their ancestral allele frequencies, and checked iteratively until no fixed loci remain, so that the final number of polymorphic loci is exactly `m_loci`).
#' @param maf_min The minimum minor allele frequency (default zero), to extend the working definition of "fixed" above to include rare variants.
#' This helps simulate a frequency-based locus ascertainment bias.
#' Loci with minor allele frequencies less than or *equal* to this value are treated as fixed (passed to [fixed_loci()]).
#' This parameter has no effect if `require_polymorphic_loci` is `FALSE`.
#' @param beta Shape parameter for a symmetric Beta for ancestral allele frequencies `p_anc`.
#' If `NA` (default), `p_anc` is uniform with range in \[0.01, 0.5\].
#' Otherwise, `p_anc` has a symmetric Beta distribution with range in \[0, 1\].
#' Has no effect if either `p_anc` or `p_anc_distr` options are non-`NULL`.
#' @param p_anc If provided, it is used as the ancestral allele frequencies (instead of drawing random ones).  Must either be a scalar or a length-`m_loci` vector.
#' If scalar and `want_p_anc = TRUE`, then the returned `p_anc` is the scalar value repeated `m_loci` times (it is always a vector).
#' If a locus was fixed and has to be redrawn, the ancestral allele frequency in `p_anc` is retained and only downstream allele frequencies and genotypes are redrawn (contrast to `p_anc_distr` below).
#' @param p_anc_distr If provided, ancestral allele frequencies are drawn with replacement from this vector (which may have any length) or function, instead of from [draw_p_anc()].
#' If a function, must accept a single parameter specifying the number of loci to draw.
#' If a locus was fixed and has to be redrawn, the ancestral allele frequency is redrawn from the distribution (contrast to `p_anc` above).
#'
#' @return A named list with the following items (which may be missing depending on options):
#'
#' - `X`: An `m`-by-`n` matrix of genotypes.
#'   Included if `want_genotypes = TRUE`.
#' - `p_anc`: A length-`m` vector of ancestral allele frequencies.
#'   Included if `want_p_anc = TRUE`.
#' - `p_subpops`: An `m`-by-`k` matrix of intermediate subpopulation allele frequencies
#'   Included if `want_p_subpops = TRUE`.
#' - `p_ind`: An `m`-by-`n` matrix of individual-specific allele frequencies.
#'   Included if `want_p_ind = TRUE`.
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
#' # FST values for k = 2 subpopulations
#' inbr_subpops <- c(0.1, 0.3)
#' # admixture proportions from 1D geography
#' admix_proportions <- admix_prop_1d_linear(n_ind, k_subpops, sigma = 1)
#'
#' # draw all random allele freqs and genotypes
#' out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci)
#'
#' # return value is a list with these items:
#' 
#' # genotypes
#' X <- out$X
#' 
#' # ancestral AFs
#' p_anc <- out$p_anc
#' 
#' # # these are excluded by default, but would be included if ...
#' # # ... `want_p_subpops == TRUE`
#' # # intermediate subpopulation AFs
#' # p_subpops <- out$p_subpops
#' # 
#' # # ... `want_p_ind == TRUE`
#' # # individual-specific AFs
#' # p_ind <- out$p_ind
#' 
#' @export
draw_all_admix <- function(
                           admix_proportions,
                           inbr_subpops = NULL,
                           m_loci,
                           tree_subpops = NULL,
                           want_genotypes = TRUE,
                           want_p_ind = FALSE,
                           want_p_subpops = FALSE,
                           want_p_anc = TRUE,
                           verbose = FALSE,
                           require_polymorphic_loci = TRUE,
                           maf_min = 0,
                           beta = NA,
                           p_anc = NULL,
                           p_anc_distr = NULL
                           ) {
    # stop if required parameters are missing
    if ( missing( admix_proportions ) )
        stop( '`admix_proportions` is required!' )
    if ( is.null( inbr_subpops ) && is.null( tree_subpops ) )
        stop( 'Either `inbr_subpops` or `tree_subpops` is required!' )
    if ( missing( m_loci ) )
        stop( '`m_loci` is required!' )

    # stop if both *_subpops structures were provided
    if ( !is.null( inbr_subpops ) && !is.null( tree_subpops ) )
        stop( '`inbr_subpops` and `tree_subpops` cannot both be provided!' )

    # stop if both `p_anc` and `p_anc_distr` were provided
    if ( !is.null( p_anc ) && !is.null( p_anc_distr ) )
        stop( '`p_anc` and `p_anc_distr` cannot both be provided!' )
    
    # ensure that things that should be matrices are so
    if ( !is.matrix( admix_proportions ) )
        stop( '`admix_proportions` must be a matrix!' )

    # get dimensions, test coherence
    #n_ind <- nrow(admix_proportions) # actually not used anywhere!?!
    k_subpops <- ncol(admix_proportions)
    if ( is.null( tree_subpops ) ) {
        k_subpops_inbr <- length(inbr_subpops)
        if (
            k_subpops_inbr > 1 &&       # it's ok if inbr_subpops is a scalar in this case
            k_subpops != k_subpops_inbr # but if it's not scalar, it must agree with admix_proportions
        )
            stop('`admix_proportions` and `inbr_subpops` are not compatible: ncol(admix_proportions) == ', k_subpops, ' != ', k_subpops_inbr, ' == length(inbr_subpops)')

        # if both inputs have names, let's stop if they don't agree, explain if it appears that they are not aligned in particular
        compare_names(
            names( inbr_subpops ),
            colnames( admix_proportions ),
            'inbr_subpops',
            'admix_proportions'
        )
    } else {
        # run overall tree validation
        validate_coanc_tree( tree_subpops )
        # now check for coherence with `admix_proportions`
        k_subpops_tree <- length( tree_subpops$tip.label )
        if ( k_subpops != k_subpops_tree )
            stop('`admix_proportions` and `tree_subpops` are not compatible: ncol(admix_proportions) == ', k_subpops, ' != ', k_subpops_tree, ' == length( tree_subpops$tip.label )')

        # if both inputs have names, let's stop if they don't agree, explain if it appears that they are not aligned in particular
        compare_names(
            names_coanc( tree_subpops ),
            colnames( admix_proportions ),
            'tree_subpops',
            'admix_proportions'
        )
    }

    # validate or generate p_anc
    p_anc_in <- p_anc # remember its original value, for later
    if ( !is.null( p_anc ) ) {
        # validate length
        m_loci_p_anc <- length( p_anc )
        if ( m_loci_p_anc == 1 ) {
            # expand scalar into vector, so it looks as desired
            p_anc <- rep.int( p_anc, m_loci )
        } else if ( m_loci_p_anc != m_loci )
            # otherwise must agree with m_loci
            stop('Provided `p_anc` has length (', m_loci_p_anc, ') neither 1 nor m_loci (', m_loci, ')')
        
        # validate range
        if ( any( p_anc < 0 ) )
            stop('Provided `p_anc` has negative values!')
        if ( any( p_anc > 1 ) )
            stop('Provided `p_anc` has values exceeding 1!')
    } else if ( !is.null( p_anc_distr ) ) {
        # draw allele frequencies from a custom distribution!
        if (verbose)
            message('drawing p_anc from p_anc_distr')
        if ( is.function( p_anc_distr ) ) {
            p_anc <- p_anc_distr( m_loci ) # ask for these many samples
        } else if ( is.numeric( p_anc_distr ) ) {
            # sample from vector, with replacement to trivially handle cases that may otherwise be too small (and to ensure sufficient randomness and exchangeability as we redraw fixed loci).
            p_anc <- sample( p_anc_distr, m_loci, replace = TRUE )
        }
    } else {
        # generate the random ancestral allele frequencies, in usual range and with minimum threshold for simplicity
        if (verbose)
            message('drawing p_anc')
        p_anc <- draw_p_anc(m_loci, beta = beta)
    }
    
    # draw intermediate allele frequencies from Balding-Nichols
    if (verbose)
        message('drawing p_subpops')
    # both cases: pass m_loci for consistency check (required if p_anc is a scalar, otherwise it already ought to match length(p_anc))
    if ( is.null( tree_subpops ) ) {
        # pass k_subpops in case inbr_subpops was a scalar (otherwise provides a redundant check)
        p_subpops <- draw_p_subpops(p_anc, inbr_subpops, m_loci = m_loci, k_subpops = k_subpops)
        # if inbr_subpops was a scalar, then p_subpops will not have column names
        # try to get them from admix_proportions in this case (ok if this ends up being NULL anyway)
        if ( length( inbr_subpops ) == 1 )
            colnames( p_subpops ) <- colnames( admix_proportions )
    } else {
        p_subpops <- draw_p_subpops_tree( p_anc, tree_subpops, m_loci = m_loci )
    }
    
    if (want_p_ind) {
        # this triggers higher-memory path
        if (verbose)
            message('drawing p_ind')
        p_ind <- make_p_ind_admix(p_subpops, admix_proportions)
        
        # draw genotypes from p_ind
        if (want_genotypes) {
            if (verbose)
                message('drawing X')
            X <- draw_genotypes_admix(p_ind)
        }
    } else {
        # draw genotypes skipping p_ind
        # (always low memory now!)
        if (want_genotypes) {
            if (verbose)
                message('drawing X')
            X <- draw_genotypes_admix(p_subpops, admix_proportions)
        }
    }
    
    if (require_polymorphic_loci && want_genotypes) {
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
            obj <- draw_all_admix(
                admix_proportions = admix_proportions,
                inbr_subpops = inbr_subpops,
                tree_subpops = tree_subpops,
                m_loci = m_loci_fixed,
                want_genotypes = want_genotypes,
                want_p_ind = want_p_ind,
                want_p_subpops = want_p_subpops,
                want_p_anc = want_p_anc,
                verbose = FALSE, # don't show more messages for additional iterations
                require_polymorphic_loci = FALSE, # loop here, avoid recursion in function, which in extreme cases leads to "node stack overflow" errors!
                maf_min = maf_min,
                beta = beta,
                p_anc = p_anc_redo,
                p_anc_distr = p_anc_distr
            )
            # overwrite fixed loci with redrawn polymorphic data
            X[fixed_loci_indexes, ] <- obj$X # guaranteed to be there
            if (want_p_ind)
                p_ind[fixed_loci_indexes, ] <- obj$p_ind
            if (want_p_subpops)
                p_subpops[fixed_loci_indexes, ] <- obj$p_subpops
            # this works when we didn't specify p_anc
            # otherwise p_anc doesn't change, it is what we wanted it to be
            if ( want_p_anc && is.null( p_anc_in ) )
                p_anc[fixed_loci_indexes] <- obj$p_anc
            
            # for next iteration, look for cases that remain fixed
            fixed_loci_indexes <- fixed_loci( X, maf_min = maf_min ) # boolean vector identifies fixed loci
            m_loci_fixed <- sum(fixed_loci_indexes) # number of cases
        }
    }
    
    # now prepare output list
    out <- list()
    if (want_genotypes)
        out$X <- X
    if (want_p_anc)
        out$p_anc <- p_anc
    if (want_p_subpops)
        out$p_subpops <- p_subpops
    if (want_p_ind)
        out$p_ind <- p_ind
    return(out)
}
