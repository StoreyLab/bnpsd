## Find the sigma that gives the desired bias coefficient
##
## This function uses a numerical solver to find the value of \eqn{\sigma} (the spread of the ancestral populations, more specifically the standard deviation of the normal densities that model this spread) that give the desired bias coefficient \eqn{s} given the other PSD parameters.
##
## Implicitly assumes our 1D admixture scenario, where the major unknown is the \eqn{\sigma} parameter.
## Therefore the structure of the admixture coefficients is not explicitly stated in the inputs.
## Also assumes uniform weights in calculating s.
##
## @param s The desired bias coefficient
## @param F The length-\eqn{k} vector of inbreeding coefficients (or \eqn{F_{ST}}'s) of the ancestral populations, up to a scaling factor (which cancels out in calculations)
## @param n The number of individuals
## @param inverval Restrict the search space of \eqn{\sigma} to this interval
## @param tol The numerical tolerance used to declare the solution found.
##
## @return The desired value of \eqn{\sigma}
##
## @examples
## # number of ancestral populations
## k <- 10
## # differentiation of ancestral populations
## F <- (1:k)/k
## # set bias coeff. of 1/2, F and n, implicitly assumes our 1D scenario
## sigma <- solveSigma(s=0.5, F=F, n=1000) 
solveSigma <- function(s, F, n, interval=c(0.1,10), tol=.Machine$double.eps) {
    ## finds the value of sigma that give a desired bias coefficient "s"
    ## we use a very generic numeric method, since the function is complicated at best
    ## interval shouldn't matter since it gets extended automatically, though we probably want to stay away from sigma=0 because it leads to NaNs (in that limit we're supposed to approach the island model)
    ## psd2s is monotonically increasing with sigma (hence extendInt='upX' gives hint of where to extend, if needed)
    ## default tolerance of .Machine$double.eps^0.25 (~ 1e-4) was kinda crummy, reduced to ~2e-16
    sigmaRoot <- stats::uniroot(function(sigma) psd2s(sigma, F, n)-s, interval, extendInt='upX', tol=tol)
    sigmaRoot$root # return this only
}

