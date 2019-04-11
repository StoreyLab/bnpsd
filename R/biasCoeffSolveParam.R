# Find the param that gives the desired bias coefficient
#
# This function uses a numerical solver to find the value of \code{param} in \code{func} that give the desired bias coefficient \eqn{s} given the other PSD parameters.
#
# Implicitly assumes our 1D admixture scenario, where the major unknown is the \code{param} parameter.
# Therefore the structure of the admixture coefficients is not explicitly stated in the inputs.
# Also assumes uniform weights in calculating s.
#
# @param s The desired bias coefficient
# @param F The length-\eqn{k} vector of inbreeding coefficients (or \eqn{F_{ST}}'s) of the intermediate subpopulations, up to a scaling factor (which cancels out in calculations)
# @param n The number of individuals
# @param func A function that accepts \code{(n, k=length(F), param)} as inputs and returns the admixture matrix Q
# @param inverval Restrict the search space of \code{param} to this interval
# @param tol The numerical tolerance used to declare the solution found.
# @param extendInt hints root solver about monotonicity (default 'upX' assumes \code{func} is monotonically increasing with \code{param})
#
# @return The desired value of \code{param}
#
# @examples
# # number of intermediate subpops
# k <- 10
# # differentiation of subpops
# F <- (1:k)/k
# # set bias coeff. of 1/2, F and n, implicitly assumes our 1D scenario
# sigma <- biasCoeffSolveParam(s=0.5, F=F, n=1000, q1d) 
biasCoeffSolveParam <- function(s, F, n, func, interval = c(0.1,10), tol = .Machine$double.eps, extendInt = 'upX') {
    # finds the value of param that give a desired bias coefficient "s"
    # we use a very generic numeric method, since the function is complicated at best

    k <- length(F) # a constant in the optimization...
    
    # function whose zero we want!
    param2sDelta <- function(param) {
        # first, construct the admixture coefficients, which critically depend on param
        Q <- func(n, k, param)
        # get bias coefficient, return difference from desired value!
        bias_coeff(Q, F) - s
    }
    
    # interval shouldn't matter since it gets extended automatically, though we probably want to stay away from param=0 because it leads to NaNs (in that limit we're supposed to approach the island model)
    # default tolerance of .Machine$double.eps^0.25 (~ 1e-4) was kinda crummy, reduced to ~2e-16
    paramRoot <- stats::uniroot(param2sDelta, interval, extendInt = extendInt, tol = tol)
    paramRoot$root # return this only
}
