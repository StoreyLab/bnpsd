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
                                 func
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

    # validate range
    if (bias_coeff > 1)
        stop('Desired `bias_coeff` must be <= 1, passed ', bias_coeff)
    if (bias_coeff < 0)
        stop('Desired `bias_coeff` must be >= 0, passed ', bias_coeff)
    
    k <- length(inbr_subpops) # a constant in the optimization...
    
    # actual minimum is greater than zero, but it seems that message should be different for non-negative values out of this range
    # set sigma = 0 here, in both cases that we have (q1d, q1dc) it's the minimum bias (independent subpopulations)!
    admix_prop_bias_coeff_min <- func(n_ind, k, sigma = 0)
    bias_coeff_min <- bias_coeff_admix(admix_prop_bias_coeff_min, inbr_subpops)
    if (bias_coeff < bias_coeff_min)
        stop('Desired `bias_coeff` must be greater than ', bias_coeff_min, ' (the minimum achievable with `sigma = 0`), passed ', bias_coeff)
    
    # function whose zero we want!
    # "sigma" is the only parameter to optimize
    # need this function inside here as it depends on extra parameters
    bias_coeff_admix_objective <- function( x ) {
        # this is a transformation that helps us explore up to `sigma = Inf` in a compact space instead:
        sigma <- x / (1 - x)
        # x == 0 => sigma == 0
        # x == 1 => sigma == Inf
        # first, construct the admixture coefficients, which critically depend on sigma
        admix_proportions <- func(n_ind, k, sigma)
        # get bias coefficient, return difference from desired value!
        bias_coeff_admix(admix_proportions, inbr_subpops) - bias_coeff
    }
    
    # default tolerance of .Machine$double.eps^0.25 (~ 1e-4) was not good enough, reduced to ~ 2e-16
    obj <- stats::uniroot(
                      bias_coeff_admix_objective,
                      interval = c(0, 1),
                      extendInt = "no",
                      tol = .Machine$double.eps
                  )
    # this is the value in terms of `x`
    x_root <- obj$root
    # transform back to sigma, return this
    x_root / (1 - x_root)
}
