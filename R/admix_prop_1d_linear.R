#' Construct admixture proportion matrix for 1D geography
#'
#' Assumes `k_subpops` intermediate subpopulations placed along a line at locations `1 : k_subpops` spread by random walks, then `n_ind` individuals equally spaced in \[`coord_ind_first`,`coord_ind_last`\] draw their admixture proportions relative to the Normal density that models the random walks of each of these intermediate subpopulations.
#' The spread of the random walks (the standard deviation of the Normal densities) is `sigma`.
#' If `sigma` is missing, it can be set indirectly by providing three variables: (1) the desired bias coefficient `bias_coeff`, (2) the coancestry matrix of the intermediate subpopulations `coanc_subpops` (up to a scalar factor), and (3) the final `fst` of the admixed individuals (see details below).
#'
#' If `sigma` is `NA`, its value is determined from the desired `bias_coeff`, `coanc_subpops` up to a scalar factor, and `fst`.
#' Uniform weights for the final generalized FST are assumed.
#' The scale of `coanc_subpops` is irrelevant because it cancels out in `bias_coeff`; after `sigma` is found, `coanc_subpops` is rescaled to give the desired final FST.
#' However, the function stops if any rescaled `coanc_subpops` values are greater than 1, which are not allowed since they are IBD probabilities.
#'
#' @param n_ind Number of individuals.
#' @param k_subpops Number of intermediate subpopulations.
#' @param sigma Spread of intermediate subpopulations (standard deviation of normal densities).
#' The edge cases `sigma = 0` and `sigma = Inf` are handled appropriately!
#' @param coord_ind_first Location of first individual (default `0.5`).
#' @param coord_ind_last Location of last individual (default `k_subpops + 0.5`).
#'
#' OPTIONS FOR BIAS COEFFICIENT VERSION
#' 
#' @param bias_coeff If `sigma` is `NA`, this bias coefficient is required.
#' @param coanc_subpops If `sigma` is `NA`, this intermediate subpops coancestry is required.
#' It can be provided as a `k_subpops`-by-`k_subpops` matrix, a length-`k_subpops` population inbreeding vector (for independent subpopulations, where between-subpop coancestries are zero) or scalar (if population inbreeding values are all equal and coancestries are zero).
#' This `coanc_subpops` can be in the wrong scale (it cancels out in calculations), which is returned corrected, to result in the desired `fst` (next).
#' @param fst If `sigma` is `NA`, this FST of the admixed individuals is required.
#'
#' @return If `sigma` was provided, returns the `n_ind`-by-`k_subpops` admixture proportion matrix (`admix_proportions`).
#' If `sigma` is missing, returns a named list containing:
#' - `admix_proportions`: the `n_ind`-by-`k_subpops` admixture proportion matrix.
#'   If `coanc_subpops` had names, they are copied to the columns of this matrix.
#' - `coanc_subpops`: the input `coanc_subpops` rescaled.
#' - `sigma`: the fit value of the spread of intermediate subpopulations
#' - `coanc_factor`: multiplicative factor used to rescale `coanc_subpops`
#'
#' @examples
#' # admixture matrix for 1000 individuals drawing alleles from 10 subpops
#' # simple version: spread of 2 standard deviations along the 1D geography
#' # (just set sigma)
#' admix_proportions <- admix_prop_1d_linear(n_ind = 1000, k_subpops = 10, sigma = 2)
#'
#' # as sigma approaches zero, admix_proportions approaches the independent subpopulations matrix
#' admix_prop_1d_linear(n_ind = 10, k_subpops = 2, sigma = 0)
#'
#' # advanced version: a similar model but with a bias coefficient of exactly 1/2
#' # (must provide bias_coeff, coanc_subpops, and fst in lieu of sigma)
#' k_subpops <- 10
#' # FST vector for intermediate independent subpops, up to a factor (will be rescaled below)
#' coanc_subpops <- 1 : k_subpops
#' obj <- admix_prop_1d_linear(
#'     n_ind = 1000,
#'     k_subpops = k_subpops,
#'     bias_coeff = 0.5,
#'     coanc_subpops = coanc_subpops,
#'     fst = 0.1 # desired final FST of admixed individuals
#' )
#' 
#' # in this case return value is a named list with three items:
#' # admixture proportions
#' admix_proportions <- obj$admix_proportions
#' 
#' # rescaled coancestry data (matrix or vector) for intermediate subpops
#' coanc_subpops <- obj$coanc_subpops
#' 
#' # and the sigma that gives the desired bias_coeff and final FST
#' sigma <- obj$sigma
#'
#' @export
admix_prop_1d_linear <- function(
                                 n_ind,
                                 k_subpops,
                                 sigma = NA,
                                 coord_ind_first = 0.5,
                                 coord_ind_last = k_subpops + 0.5,
                                 bias_coeff = NA,
                                 coanc_subpops = NULL,
                                 fst = NA
                                 ) {
    # stop if these required parameters are missing
    if (missing(n_ind))
        stop('`n_ind` is required!')
    if (missing(k_subpops))
        stop('`k_subpops` is required!')

    # figure out if we need the more complicated algorithm...
    fit_bias_coeff <- is.na(sigma) # remember after it was set
    if (fit_bias_coeff) {
        # this triggers version that fits bias coefficient
        
        # check for more required parameters
        if (is.na(bias_coeff))
            stop('`bias_coeff` is required when sigma is missing!')
        if (is.null(coanc_subpops))
            stop('`coanc_subpops` is required when sigma is missing!')
        if (is.na(fst))
            stop('fst` is required when sigma is missing!')

        # fit sigma!
        sigma <- bias_coeff_admix_fit(
            bias_coeff = bias_coeff,
            coanc_subpops = coanc_subpops,
            n_ind = n_ind,
            k_subpops = k_subpops,
            func = admix_prop_1d_linear,
            coord_ind_first = coord_ind_first,
            coord_ind_last = coord_ind_last
        )
    } else {
        # validate input sigma here (bias_coeff_admix_fit ought to return valid numbers)
        if ( sigma < 0 )
            stop('sigma must be non-negative!')
    }
    sigma2 <- - 2 * sigma^2 # square once for loop below (plus other Normal constants)
    
    # the x-coordinates of the n individuals based on [coord_ind_first, coord_ind_last] limits
    xs <- coord_ind_first + ( 0 : ( n_ind - 1 ) ) / ( n_ind - 1 ) * (coord_ind_last - coord_ind_first)
    # and subpopulations (same deal, except fixed [1, k] range)
    mus <- 1 : k_subpops
    
    # construct the coefficients of each person now!
    admix_proportions <- matrix( 0, nrow = n_ind, ncol = k_subpops )
    for ( i in 1 : n_ind ) {
        # this is the "normal" way (pun intended)
        if ( sigma > 0 ) {
            # collect the density values for each intermediate subpopulation at individual i's position
            # line implements super fast Normal without constant factors (which only involve constant sigma)
            # NOTE: sigma2 has negative sign built into it!
            # NOTE: sigma = Inf is correctly handled here (gives all admix_proportions == 1 before normalizing)
            admix_proportions[i,] <- exp( (xs[i] - mus)^2 / sigma2 )
        }
        
        # if sigma was zero, the entire row is still at zero
        # if the normal way failed, the row will also be all zeroes
        # so test for that and apply a discrete clustering either way
        if ( sigma == 0 || sum( admix_proportions[ i, ] ) == 0 ) {
            # let's handle this special case, the limit of which is the island model
            # compute distances to the subpopulations
            distances <- (xs[i] - mus)^2
            # find the minimum distance
            min_distance <- min( distances )
            # ok to set to booleans (turns numeric upon saving to matrix, or if sigma==0, upon normalization)
            # the minima will be TRUE (1), the rest FALSE (0)
            # this ensures ties get admix_proportions split evenly (after normalizing at the end)
            admix_proportions[i,] <- distances == min_distance
        }
    } 
    # normalize to have rows/coefficients sum to 1!
    admix_proportions_row_sums <- rowSums( admix_proportions )
    # in some extreme examples a whole row is zero, resulting in NAs
    # we should have handled that already, but check anyway
    if ( any( admix_proportions_row_sums == 0 ) )
        stop( 'Resulting `admix_proportions` had rows with zero sums, cannot normalize!  sigma = ', sigma )
    admix_proportions <- admix_proportions / admix_proportions_row_sums
    
    # check for issues so far
    if ( anyNA( admix_proportions ) )
        stop( '`admix_proportions` had NAs!  sigma = ', sigma )

    if (fit_bias_coeff) {
        # this triggers version that fits bias coefficient
        # let's rescale coanc_subpops now!
        obj <- rescale_coanc_subpops(admix_proportions, coanc_subpops, fst) 
        
        # get names from `coanc_subpops`, inherit into `admix_proportions` if non-NULL
        names_coanc_subpops <- names_coanc( coanc_subpops )
        if ( !is.null( names_coanc_subpops ) )
            colnames( admix_proportions ) <- names_coanc_subpops
        
        # return all this additional data!
        return(
            list(
                admix_proportions = admix_proportions,
                coanc_subpops = obj$coanc_subpops,
                sigma = sigma,
                coanc_factor = obj$factor
            )
        )
    } else {
        # in direct case, always return admix_proportions
        # colnames are always NULL in this case
        return(admix_proportions)
    }
}
