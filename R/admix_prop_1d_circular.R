# 1/kappa = sigma^2

#' Construct admixture proportion matrix for circular 1D geography
#'
#' Assumes \eqn{k} intermediate subpopulations placed along a circumference (the \eqn{[0, 2\pi]} line that wraps around) with even spacing spread by random walks (see details below), then \eqn{n} individuals sampled equally spaced in \eqn{[a,b]} (default \eqn{[0, 2\pi]} with a small gap so first and last individual do not overlap) draw their admixture proportions relative to the Von Mises density that models the random walks of each of these intermediate subpopulations.
#' The spread of the random walks (the \eqn{\sigma=1/\sqrt{\kappa}} of the Von Mises densities) is set to \code{sigma} if not missing, otherwise \eqn{\sigma} is found numerically to give the desired bias coefficient \code{bias_coeff}, the coancestry matrix of the intermediate subpopulations \code{coanc_subpops} (up to a scalar factor), and the final \eqn{F_{ST}}{FST} of the admixed individuals (see details below).
#'
#' Assuming the full range of \eqn{[0, 2\pi]} is considered, and the first and last individuals do not overlap, the gap between individuals is \eqn{\Delta = 2 \pi / n}.
#' To not have any individuals on the edge, we place the first individual at \eqn{\Delta / 2} and the last at \eqn{2 \pi - \Delta / 2}.
#' The location of subpopulation \eqn{j} is
#' \deqn{\Delta / 2 + (j-1/2)/k (2 \pi - \Delta),}
#' chosen to agree with the default correspondence between individuals and subpopulations of the linear 1D geography admixture scenario (\code{\link{admix_prop_1d_linear}}).
#'
#' When \code{sigma} is missing, the function determines its value using the desired \code{bias_coeff}, \code{coanc_subpops} up to a scalar factor, and \code{fst}.
#' Uniform weights for the final generalized \eqn{F_{ST}}{FST} are assumed.
#' The scaling factor of the input \code{coanc_subpops} is irrelevant because it cancels out in \code{bias_coeff}; after \code{sigma} is found, \code{coanc_subpops} is rescaled to give the desired final \eqn{F_{ST}}{FST}.
#' However, the function stops with a fatal error if the rescaled \code{coanc_subpops} takes on any values greater than 1, which are not allowed since \code{coanc_subpops} are IBD probabilities.
#'
#' @param n_ind Number of individuals
#' @param k_subpops Number of intermediate subpopulations
#' @param sigma Spread of intermediate subpopulations (approximate standard deviation of Von Mises densities, see above)
#' The edge cases \code{sigma = 0} and \code{sigma = Inf} are handled appropriately!
#' @param coord_ind_first Location of first individual
#' @param coord_ind_last Location of last individual
#'
#' OPTIONS FOR BIAS COEFFICIENT VERSION
#' 
#' @param bias_coeff The desired bias coefficient, which specifies \eqn{\sigma} indirectly.
#' Required if \code{sigma} is missing.
#' @param coanc_subpops The length-\eqn{k} vector of inbreeding coefficients (or \eqn{F_{ST}}{FST}'s) of the intermediate subpopulations, up to a scaling factor (which cancels out in calculations).
#' Required if \code{sigma} is missing.
#' @param fst The desired final \eqn{F_{ST}}{FST} of the admixed individuals.
#' Required if \code{sigma} is missing.
#'
#' @return If \code{sigma} was provided, the \eqn{n \times k}{n-by-k} admixture proportion matrix.
#' If \code{sigma} is missing, a named list is returned containing \code{admix_proportions}, the rescaled \code{coanc_subpops}, and the \code{sigma} that together give the desired \eqn{bias_coeff} and final \eqn{F_{ST}}{FST} of the admixed individuals.
#'
#' @examples
#' # admixture matrix for 1000 individuals drawing alleles from 10 subpops
#' # and a spread of about 2 standard deviations along the circular 1D geography
#' admix_proportions <- admix_prop_1d_circular(n_ind = 1000, k_subpops = 10, sigma = 2)
#'
#' # a similar model but with a bias coefficient of exactly 1/2
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
                                   bias_coeff,
                                   coanc_subpops,
                                   fst
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
        if (missing(bias_coeff))
            stop('`bias_coeff` is required when sigma is missing!')
        if (missing(coanc_subpops))
            stop('`coanc_subpops` is required when sigma is missing!')
        if (missing(fst))
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
# stick deprecated function name here

#' @title Construct admixture proportion matrix for 1D geography
#' @description Construct admixture proportion matrix for 1D geography
#'
#' @param n Number of individuals
#' @param k Number of intermediate subpopulations
#' @param sigma Spread of intermediate subpopulations (approximate standard deviation of Von Mises densities, see above)
#' @param a Location of first individual
#' @param b Location of last individual
#'
#' OPTIONS FOR BIAS COEFFICIENT VERSION
#' 
#' @param s The desired bias coefficient, which specifies \eqn{\sigma} indirectly.  Required if \code{sigma} is missing
#' @param F The vector of inbreeding coefficients of the intermediate subpopulations, up to a scaling factor (which cancels out in calculations).  Required if \code{sigma} is missing
#' @param Fst The desired final \eqn{F_{ST}}{FST} of the admixed individuals.  Required if \code{sigma} is missing
#' @return If \code{sigma} was provided, the \eqn{n \times k}{n-by-k} admixture proportion matrix \eqn{Q}.
#' If \code{sigma} is missing, a named list is returned containing \code{Q}, the rescaled \code{F}, and the \code{sigma} that together give the desired \eqn{s} and final \eqn{F_{ST}}{FST} of the admixed individuals.
#'
#' @name q1dc-deprecated
#' @usage q1dc(n, k, sigma, a = 0, b = 2 * pi, s, F, Fst)
#' @seealso \code{\link{bnpsd-deprecated}}
#' @keywords internal
NULL

#' @rdname bnpsd-deprecated
#' @section \code{q1d}:
#' For \code{q1d}, use \code{\link{admix_prop_1d_circular}}.
#'
#' @export
q1dc <- function(n, k, sigma, a = 0, b = 2 * pi, s, F, Fst) {
    # mark as deprecated
    .Deprecated('admix_prop_1d_circular')
    # return as usual, to not break things just yet
    admix_prop_1d_circular(
        n_ind = n,
        k_subpops = k,
        sigma = sigma,
        coord_ind_first = a,
        coord_ind_last = b,
        bias_coeff = s,
        coanc_subpops = F,
        fst = Fst
    )
}
