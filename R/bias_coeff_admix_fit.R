# Find the param that gives the desired bias coefficient
#
# This function uses a numerical solver to find the value of \code{param} in \code{func} that give the desired bias coefficient \eqn{s} given the other PSD parameters.
#
# Implicitly assumes our 1D admixture scenario, where the major unknown is the \code{param} parameter.
# Therefore the structure of the admixture coefficients is not explicitly stated in the inputs.
# Also assumes uniform weights in calculating s.
#
# @param bias_coeff The desired bias coefficient
# @param inbr_subpops The length-\eqn{k} vector of inbreeding coefficients (or \eqn{F_{ST}}'s) of the intermediate subpopulations, up to a scaling factor (which cancels out in calculations)
# @param n_ind The number of individuals
# @param func A function that accepts \code{(n, k = length(inbr_subpops), param)} as inputs and returns the admixture matrix admix_proportions
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
# inbr_subpops <- (1:k)/k
# # set bias coeff. of 1/2, inbr_subpops and n, implicitly assumes our 1D scenario
# sigma <- bias_coeff_admix_fit(bias_coeff = 0.5, inbr_subpops = inbr_subpops, n_ind = 1000, func = q1d) 
bias_coeff_admix_fit <- function(
                                 bias_coeff,
                                 inbr_subpops,
                                 n_ind,
                                 func,
                                 interval = c(0.1, 10),
                                 tol = .Machine$double.eps,
                                 extendInt = 'upX'
                                 ) {
    # finds the value of param that give a desired bias coefficient "s"
    # we use a very generic numeric method, since the function is complicated at best
    if (missing(bias_coeff))
        stop('Desired `bias_coeff` is required!')
    if (missing(inbr_subpops))
        stop('`inbr_subpops` is required!')
    if (missing(n_ind))
        stop('Number of individuals `n_ind` is required!')
    if (missing(func))
        stop('Admixture proportion function `func` is required!')
    
    k <- length(inbr_subpops) # a constant in the optimization...
    
    # function whose zero we want!
    # "sigma" is the only parameter to optimize
    # need it inside here as it depends on extra parameters
    bias_coeff_admix_objective <- function(sigma) {
        # first, construct the admixture coefficients, which critically depend on sigma
        admix_proportions <- func(n_ind, k, sigma)
        # get bias coefficient, return difference from desired value!
        bias_coeff_admix(admix_proportions, inbr_subpops) - bias_coeff
    }
    
    # interval shouldn't matter since it gets extended automatically, though we probably want to stay away from sigma = 0 because it leads to NaNs (in that limit we're supposed to approach the island model)
    # default tolerance of .Machine$double.eps^0.25 (~ 1e-4) was not good enough, reduced to ~ 2e-16
    obj <- stats::uniroot(bias_coeff_admix_objective, interval, extendInt = extendInt, tol = tol)
    obj$root # return this only
}
