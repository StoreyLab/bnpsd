## PUBLIC: used by ~/storey/fst/software/simPSD02.*.R

## TODO: add other q* scenarios (fixed island, 1D circular, 2D, etc)

#' Construct admixture proportion matrix for 1D geography
#'
#' Assumes \eqn{k} ancestral populations placed along a line at locations \eqn{1:k} spread by random walks, then \eqn{n} individuals sampled equally spaced in \eqn{[a,b]} (default \eqn{[0.5, k+0.5]}) draw their admixture proportions relative to the Normal density that models the random walks of each of these ancestral populations.
#' The spread of the random walks (the \eqn{\sigma} of the Normal densities) is set to \code{sigma} if not missing, otherwise \eqn{\sigma} is found numerically to give the desired bias coefficient \code{s}, the vector \code{F} of \eqn{F_{ST}}{FST}'s for the intermediate subpopulations up to a scalar factor, and the final \eqn{F_{ST}}{FST} of the admixed individuals (see details below).
#'
#' When \code{sigma} is missing, the function determines its value using the desired \code{s}, \code{F} up to a scalar factor, and \code{Fst}.
#' Uniform weights for the final generalized \eqn{F_{ST}}{FST} are assumed.
#' The scaling factor of the input \code{F} is irrelevant because it cancels out in \code{s}; after \code{sigma} is found, \code{F} is rescaled to give the desired final \eqn{F_{ST}}{FST}.
#' However, the function stops with a fatal error if the rescaled \code{F} takes on any values greater than 1, which are not allowed since \code{F} are IBD probabilities.
#'
#' @param n Number of individuals
#' @param k Number of ancestral populations
#' @param sigma Spread of ancestral populations (standard deviation of normal densities)
#' @param a Location of first individual
#' @param b Location of last individual
#'
#' OPTIONS FOR BIAS COEFFICIENT VERSION
#' 
#' @param s The desired bias coefficient, which specifies \eqn{\sigma} indirectly.  Required if \code{sigma} is missing
#' @param F The length-\eqn{k} vector of inbreeding coefficients (or \eqn{F_{ST}}{FST}'s) of the ancestral populations, up to a scaling factor (which cancels out in calculations).  Required if \code{sigma} is missing
#' @param Fst The desired final \eqn{F_{ST}}{FST} of the admixed individuals.  Required if \code{sigma} is missing
#' @param interval Restrict the search space of \eqn{\sigma} to this interval
#' @param tol The numerical tolerance used to declare the solution found
#'
#' @return If \code{sigma} was provided, the \eqn{n \times k}{n-by-k} admixture proportion matrix \eqn{Q}.  If \code{sigma} is missing, a named list is returned containing \code{Q}, the rescaled \code{F}, and the \code{sigma} that together give the desired \eqn{s} and final \eqn{F_{ST}}{FST} of the admixed individuals.
#'
#' @examples
#' ## admixture matrix for 1000 individuals drawing alleles from 10 ancestral populations
#' ## and a spread of 2 standard deviations along the 1D geography
#' Q <- q1d(n=1000, k=10, sigma=2)
#'
#' ## a similar model but with a bias coefficient "s" of exactly 1/2
#' k <- 10
#' F <- 1:k # Fst vector for intermediate subpops, up to a factor (will be rescaled below)
#' Fst <- 0.1 # desired final Fst of admixed individuals
#' obj <- q1d(n=1000, k=k, s=0.5, F=F, Fst=Fst)
#' ## in this case return value is a named list with three items:
#' Q <- obj$Q # admixture proportions
#' F <- obj$F # rescaled Fst vector for intermediate subpopulations
#' sigma <- obj$sigma # and the sigma that gives the desired s and final Fst
#'
#' @export
q1d <- function(n, k, sigma, a=0.5, b=k+0.5, s, F, Fst, interval=c(0.1,10), tol=.Machine$double.eps) {
    ## figure out if we need the more complicated algorithm...
    sigmaMissing <- missing(sigma) # remember after it was set
    if (sigmaMissing) { ## this triggers s version
        if (missing(s)) stop('Fatal in q1d: s is required when sigma is missing!')
        if (missing(F)) stop('Fatal in q1d: F is required when sigma is missing!')
        if (missing(Fst)) stop('Fatal in q1d: Fst is required when sigma is missing!')
        sigma <- solveSigma(s, F, n, interval, tol)
    }
    
    ## the x-coordinates of the n individuals based on [a,b] limits
    xs <- a + (0:(n-1))/(n-1)*(b-a)
    
    ## construct the coefficients of each person now!
    Q <- matrix(nrow=n, ncol=k) # dimensions match that of makeQ
    for (i in 1:n) {
        ## collect the density values for each ancestral population at individual i's position
        Q[i,] <- stats::dnorm(xs[i], mean=1:k, sd=sigma)
    }
    ## normalize to have rows/coefficients sum to 1!
    Q <- Q/rowSums(Q) # return!

    if (sigmaMissing) { ## this triggers s version
        F <- rescaleFst(Q, F, Fst) # let's rescale F now!
        return( list(Q=Q, F=F, sigma=sigma) ) # return all this additional data!
    } else {
        return(Q) # in direct case, always return Q
    }
}

