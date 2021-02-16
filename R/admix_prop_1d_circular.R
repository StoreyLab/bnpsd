#' Construct admixture proportion matrix for circular 1D geography
#'
#' Assumes `k_subpops` intermediate subpopulations placed along a circumference (the \[`0`, `2 * pi`\] line that wraps around) with even spacing spread by random walks (see details below), then `n_ind` individuals sampled equally spaced in \[`coord_ind_first`,`coord_ind_last`\] (default \[`0`, `2 * pi`\] with a small gap so first and last individual do not overlap) draw their admixture proportions relative to the Von Mises density that models the random walks of each of these intermediate subpopulations.
#' The spread of the random walks is `sigma = 1 / sqrt(kappa)` of the Von Mises density.
#' If `sigma` is missing, it can be set indirectly by providing three variables: (1) the desired bias coefficient `bias_coeff`, (2) the coancestry matrix of the intermediate subpopulations `coanc_subpops` (up to a scalar factor), and (3) the final `fst` of the admixed individuals (see details below).
#'
#' Assuming the full range of \[`0`, `2 * pi`\] is considered, and the first and last individuals do not overlap, the gap between individuals is `delta = 2 * pi / n`.
#' To not have any individuals on the edge, we place the first individual at `delta / 2` and the last at `2 * pi - delta / 2`.
#' The location of subpopulation `j` is `delta / 2 + ( j - 1/2 ) / k * (2 * pi - delta)`, chosen to agree with the default correspondence between individuals and subpopulations of the linear 1D geography admixture scenario ([admix_prop_1d_linear()]).
#'
#' If `sigma` is `NA`, its value is determined from the desired `bias_coeff`, `coanc_subpops` up to a scalar factor, and `fst`.
#' Uniform weights for the final generalized FST are assumed.
#' The scale of `coanc_subpops` is irrelevant because it cancels out in `bias_coeff`; after `sigma` is found, `coanc_subpops` is rescaled to give the desired final FST.
#' However, the function stops if any rescaled `coanc_subpops` values are greater than 1, which are not allowed since they are IBD probabilities.
#'
#' @param n_ind Number of individuals
#' @param k_subpops Number of intermediate subpopulations
#' @param sigma Spread of intermediate subpopulations (approximate standard deviation of Von Mises densities, see above)
#' The edge cases `sigma = 0` and `sigma = Inf` are handled appropriately!
#' @param coord_ind_first Location of first individual
#' @param coord_ind_last Location of last individual
#'
#' OPTIONS FOR BIAS COEFFICIENT VERSION
#' 
#' @param bias_coeff If `sigma` is `NA`, this bias coefficient is required.
#' @param coanc_subpops If `sigma` is `NA`, this intermediate subpops coancestry is required.
#' It can be provided as a `k_subpops`-by-`k_subpops` matrix, a length-`k_subpops` population inbreeding vector (for independent subpopulations, where between-subpop coancestries are zero) or scalar (if population inbreeding values are all equal and coancestries are zero).
#' This `coanc_subpops` can be in the wrong scale (it cancels out in calculations), which is returned corrected, to result in the desired `fst` (next).
#' @param fst If `sigma` is `NA`, this FST of the admixed individuals is required.
#'
#' @return If `sigma` was provided, the `n_ind`-by-`k_subpops` admixture proportion matrix (`admix_proportions`).
#' If `sigma` is missing, a named list is returned containing `admix_proportions`, the rescaled `coanc_subpops`, and the `sigma` (which together give the desired `bias_coeff` and `fst`).
#'
#' @examples
#' # admixture matrix for 1000 individuals drawing alleles from 10 subpops
#' # simple version: spread of about 2 standard deviations along the circular 1D geography
#' # (just set sigma)
#' admix_proportions <- admix_prop_1d_circular(n_ind = 1000, k_subpops = 10, sigma = 2)
#'
#' # advanced version: a similar model but with a bias coefficient of exactly 1/2
#' # (must provide bias_coeff, coanc_subpops, and fst in lieu of sigma)
#' k_subpops <- 10
#' # FST vector for intermediate independent subpops, up to a factor (will be rescaled below)
#' coanc_subpops <- 1 : k_subpops
#' obj <- admix_prop_1d_circular(
#'     n_ind = 1000,
#'     k_subpops = k_subpops,
#'     bias_coeff = 0.5,
#'     coanc_subpops = coanc_subpops,
#'     fst = 0.1 # desired final FST of admixed individuals
#' )
#' 
#' # in this case return value is a named list with three items:
#' admix_proportions <- obj$admix_proportions
#' 
#' # rescaled coancestry data (matrix or vector) for intermediate subpops
#' coanc_subpops <- obj$coanc_subpops
#' 
#' # and the sigma that gives the desired bias_coeff and final FST
#' sigma <- obj$sigma
#'
#' @export
admix_prop_1d_circular <- function(
                                   n_ind,
                                   k_subpops,
                                   sigma = NA,
                                   coord_ind_first = 2 * pi / (2 * n_ind),
                                   coord_ind_last = 2 * pi * (1 - 1 / (2 * n_ind) ),
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
            stop('`fst` is required when sigma is missing!')
        
        # fit sigma!
        sigma <- bias_coeff_admix_fit(
            bias_coeff = bias_coeff,
            coanc_subpops = coanc_subpops,
            n_ind = n_ind,
            k_subpops = k_subpops,
            func = admix_prop_1d_circular,
            coord_ind_first = coord_ind_first,
            coord_ind_last = coord_ind_last
        )
    } else {
        # validate input sigma here (bias_coeff_admix_fit ought to return valid numbers)
        if ( sigma < 0 )
            stop('sigma must be non-negative!')
    }
    sigma2 <- sigma^2 # square once for loop below

    # the x-coordinates of the n individuals based on [coord_ind_first, coord_ind_last] limits
    xs <- coord_ind_first + ( 0 : ( n_ind - 1 ) ) / ( n_ind - 1 ) * (coord_ind_last - coord_ind_first)
    # and subpopulations (same deal, except fixed [0, 2*pi] range)
    mus <- ( ( (1 : k_subpops) - 0.5 ) / k_subpops * (1 - 1 / n_ind) + 1 / (2 * n_ind) ) * 2 * pi
    
    # construct the coefficients of each person now!
    admix_proportions <- matrix(nrow = n_ind, ncol = k_subpops)
    for (i in 1 : n_ind) {
        if (sigma == 0) {
            # let's handle this special case, the limit of which is the island model
            # compute distances to the subpopulations
            distances <- 1 - cos(xs[i] - mus)
            # find the minimum distance
            min_distance <- min( distances )
            # ok to set to booleans (normalization will turn numeric)
            # the minima will be TRUE (1), the rest FALSE (0)
            # this ensures ties get admix_proportions split evenly (after normalizing at the end)
            admix_proportions[i,] <- distances == min_distance
        } else {
            # collect the density values for each intermediate subpopulation at individual i's position
            # line implements super fast Von Mises without constant factors (which only involve constant sigma)
            # NOTE: sigma = Inf is correctly handled here (gives all admix_proportions == 1 before normalizing)
            admix_proportions[i,] <- exp( cos(xs[i] - mus) / sigma2 )
        }
    }
    # normalize to have rows/coefficients sum to 1!
    admix_proportions <- admix_proportions / rowSums(admix_proportions)

    if (fit_bias_coeff) {
        # this triggers version that fits bias coefficient
        coanc_subpops <- rescale_coanc_subpops(admix_proportions, coanc_subpops, fst) # let's rescale coanc_subpops now!
        return( list(admix_proportions = admix_proportions, coanc_subpops = coanc_subpops, sigma = sigma) ) # return all this additional data!
    } else {
        return(admix_proportions) # in direct case, always return admix_proportions
    }
}
