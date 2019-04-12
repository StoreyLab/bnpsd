# 1/kappa = sigma^2

#' Construct admixture proportion matrix for circular 1D geography
#'
#' Assumes \eqn{k} intermediate subpopulations placed along a circumference (the \eqn{[0, 2\pi]} line that wraps around) with even spacing spread by random walks, then \eqn{n} individuals sampled equally spaced in \eqn{[a,b]} (default \eqn{[0, 2\pi]}) draw their admixture proportions relative to the Von Mises density that models the random walks of each of these intermediate subpopulations.
#' The spread of the random walks (the \eqn{\sigma=1/\sqrt{\kappa}} of the Von Mises densities) is set to \code{sigma} if not missing, otherwise \eqn{\sigma} is found numerically to give the desired bias coefficient \code{s}, the vector \code{F} of \eqn{F_{ST}}{FST}'s for the intermediate subpopulations up to a scalar factor, and the final \eqn{F_{ST}}{FST} of the admixed individuals (see details below).
#'
#' When \code{sigma} is missing, the function determines its value using the desired \code{s}, \code{F} up to a scalar factor, and \code{Fst}.
#' Uniform weights for the final generalized \eqn{F_{ST}}{FST} are assumed.
#' The scaling factor of the input \code{F} is irrelevant because it cancels out in \code{s}; after \code{sigma} is found, \code{F} is rescaled to give the desired final \eqn{F_{ST}}{FST}.
#' However, the function stops with a fatal error if the rescaled \code{F} takes on any values greater than 1, which are not allowed since \code{F} are IBD probabilities.
#'
#' @param n Number of individuals
#' @param k Number of intermediate subpopulations
#' @param sigma Spread of intermediate subpopulations (approximate standard deviation of Von Mises densities, see above)
#' The edge cases \code{sigma = 0} and \code{sigma = Inf} are handled appropriately!
#' @param a Location of first individual
#' @param b Location of last individual
#'
#' OPTIONS FOR BIAS COEFFICIENT VERSION
#' 
#' @param s The desired bias coefficient, which specifies \eqn{\sigma} indirectly.  Required if \code{sigma} is missing
#' @param F The length-\eqn{k} vector of inbreeding coefficients (or \eqn{F_{ST}}{FST}'s) of the intermediate subpopulations, up to a scaling factor (which cancels out in calculations).  Required if \code{sigma} is missing
#' @param Fst The desired final \eqn{F_{ST}}{FST} of the admixed individuals.  Required if \code{sigma} is missing
#' @param interval Restrict the search space of \eqn{\sigma} to this interval
#' @param tol The numerical tolerance used to declare the solution found
#'
#' @return If \code{sigma} was provided, the \eqn{n \times k}{n-by-k} admixture proportion matrix \eqn{Q}.  If \code{sigma} is missing, a named list is returned containing \code{Q}, the rescaled \code{F}, and the \code{sigma} that together give the desired \eqn{s} and final \eqn{F_{ST}}{FST} of the admixed individuals.
#'
#' @examples
#' # admixture matrix for 1000 individuals drawing alleles from 10 subpops
#' # and a spread of about 2 standard deviations along the circular 1D geography
#' Q <- q1dc(n = 1000, k = 10, sigma = 2)
#'
#' # a similar model but with a bias coefficient "s" of exactly 1/2
#' k <- 10
#' F <- 1:k # Fst vector for intermediate subpops, up to a factor (will be rescaled below)
#' Fst <- 0.1 # desired final Fst of admixed individuals
#' obj <- q1dc(n = 1000, k = k, s = 0.5, F = F, Fst = Fst)
#' # in this case return value is a named list with three items:
#' Q <- obj$Q # admixture proportions
#' F <- obj$F # rescaled Fst vector for intermediate subpops
#' sigma <- obj$sigma # and the sigma that gives the desired s and final Fst
#'
#' @export
q1dc <- function(n, k, sigma, a = 0, b = 2 * pi, s, F, Fst, interval = c(0.1,10), tol = .Machine$double.eps) {
    
    # figure out if we need the more complicated algorithm...
    sigmaMissing <- missing(sigma) # remember after it was set
    if (sigmaMissing) { # this triggers s version
        if (missing(s))
            stop('s is required when sigma is missing!')
        if (missing(F))
            stop('F is required when sigma is missing!')
        if (missing(Fst))
            stop('Fst is required when sigma is missing!')
        sigma <- bias_coeff_admix_fit(s, F, n, q1dc, interval, tol)
    } else {
        # validate input sigma here (bias_coeff_admix_fit ought to return valid numbers)
        if ( sigma < 0 )
            stop('sigma must be non-negative!')
    }
    sigma2 <- sigma^2 # square once for loop below
    
    # the x-coordinates of the n individuals based on [a,b] limits
    xs <- a + (0:(n-1)) / n * (b - a)
    # and subpopulations (same deal, except fixed [0,2*pi] range)
    mus <- (0:(k-1)) / k * 2 * pi
    
    # construct the coefficients of each person now!
    Q <- matrix(nrow = n, ncol = k) # dimensions match that of makeQ
    for (i in 1:n) {
        if (sigma == 0) {
            # let's handle this special case, the limit of which is the island model
            # compute distances to the subpopulations
            distances <- cos(xs[i] - mus)
            # find the minimum distance
            min_distance <- min( distances )
            # ok to set to booleans (normalization will turn numeric)
            # the minima will be TRUE (1), the rest FALSE (0)
            # this ensures ties get admix_proportions split evenly (after normalizing at the end)
            Q[i,] <- distances == min_distance
        } else {
            # collect the density values for each intermediate subpopulation at individual i's position
            # line implements super fast Von Mises without constant factors (which only involve constant sigma)
            # NOTE: sigma = Inf is correctly handled here (gives all Q == 1 before normalizing)
            Q[i,] <- exp( cos(xs[i] - mus) / sigma2 )
        }
    }
    # normalize to have rows/coefficients sum to 1!
    Q <- Q / rowSums(Q)

    if (sigmaMissing) { # this triggers s version
        F <- rescaleF(Q, F, Fst) # let's rescale F now!
        return( list(Q = Q, F = F, sigma = sigma) ) # return all this additional data!
    } else {
        return(Q) # in direct case, always return Q
    }
}
