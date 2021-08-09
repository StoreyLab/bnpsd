# Find the sigma that gives the desired bias coefficient
#
# This function uses a numerical solver to find the value of `sigma` in `func` that give the desired bias coefficient `s` given the other PSD parameters.
# The admixture scenarios guaranteed to work are the 1D linear and 1D circular geographies as provided in `func`.
# Therefore the structure of the admixture coefficients is not explicitly stated in the inputs.
# Assumes uniform weights in calculating the bias coefficient.
#
# @param bias_coeff The desired bias coefficient.
# @param coanc_subpops The coancestry matrix of the intermediate subpopulations (or equivalent vector or scalar), up to a scaling factor (which is irrelevant as it cancels out in the bias coefficient).
# @param n_ind The number of individuals.
# @param k_subpops The number of subpopulations.
# @param func A function that accepts `(n_ind, k_subpops, sigma, coord_ind_first = coord_ind_first, coord_ind_last = coord_ind_last)` as inputs and returns the admixture proportions matrix.
# @param coord_ind_first Location of first individual (to pass to func).
# @param coord_ind_last Location of last individual (to pass to func).
#
# @return The desired value of `sigma`
#
# @examples
# # number of intermediate subpops
# k <- 10
# # set bias coeff. of 1/2, coanc_subpops and n, implicitly assumes our 1D scenario
# sigma <- bias_coeff_admix_fit(
#     bias_coeff = 0.5,
#     coanc_subpops = (1:k)/k,
#     n_ind = 1000,
#     k_subpops = k,
#     func = admix_prop_1d_linear,
#     coord_ind_first = 0.5,
#     coord_ind_last = k + 0.5
# ) 
bias_coeff_admix_fit <- function(
                                 bias_coeff,
                                 coanc_subpops,
                                 n_ind,
                                 k_subpops,
                                 func,
                                 coord_ind_first,
                                 coord_ind_last
                                 ) {
    # finds the value of param that give a desired bias coefficient "s"
    # we use a very generic numeric method, since the function is complicated at best
    if (missing(bias_coeff))
        stop('Desired `bias_coeff` is required!')
    if (missing(coanc_subpops))
        stop('`coanc_subpops` is required!')
    if (missing(n_ind))
        stop('Number of individuals `n_ind` is required!')
    if (missing(k_subpops))
        stop('Number of subpopulations `k_subpops` is required!')
    if (missing(func))
        stop('Admixture proportion function `func` is required!')
    if (missing(coord_ind_first))
        stop('`coord_ind_first` is required!')
    if (missing(coord_ind_last))
        stop('`coord_ind_last` is required!')
    
    # validate range
    if (bias_coeff > 1)
        stop('Desired `bias_coeff` must be <= 1, passed ', bias_coeff)
    if (bias_coeff < 0)
        stop('Desired `bias_coeff` must be >= 0, passed ', bias_coeff)

    # tolerance should be precision-dependent, but when configured with --disable-long-double the value of .Machine$double.eps doesn't change!
    # so here we test for that config difference (by testing if .Machine$sizeof.longdouble is zero) and reduce the tolerance accordingly
    tol <- .Machine$double.eps
    if ( .Machine$sizeof.longdouble == 0 )
        tol <- sqrt( tol )
    
    # handle an edge case that is problematic in some systems (looking at you Apple M1)
    if ( bias_coeff == 1 ) {
        # in the existing cases the answer is:
        sigma <- Inf
        # before returning that value, let's check that we get the right answer
        # (to ensure that this will not make incorrect assumption if applied to other "func"s other than the two official ones
        admix_proportions <- func(
            n_ind,
            k_subpops,
            sigma,
            coord_ind_first,
            coord_ind_last
        )
        # get actual bias coefficient of this case
        bias_coeff2 <- bias_coeff_admix(admix_proportions, coanc_subpops)
        # equality is not always feasible due to machine precision, but if it's close enough let's do this
        if ( abs( bias_coeff2 - bias_coeff ) < tol )
            return( sigma )
    }
    
    # actual minimum is greater than zero, but it seems that message should be different for non-negative values out of this range
    # set sigma = 0 here, in both cases that we have (admix_prop_1d_linear, admix_prop_1d_circular) it's the minimum bias (independent subpopulations)!
    admix_prop_bias_coeff_min <- func(
        n_ind,
        k_subpops,
        sigma = 0,
        coord_ind_first = coord_ind_first,
        coord_ind_last = coord_ind_last
    )
    bias_coeff_min <- bias_coeff_admix(admix_prop_bias_coeff_min, coanc_subpops)
    if (bias_coeff < bias_coeff_min)
        stop('Desired `bias_coeff` must be greater than ', bias_coeff_min, ' (the minimum achievable with `sigma = 0`), passed ', bias_coeff, '.  Tip: This minimum depends most strongly on the input `coanc_subpops`, so the main alternative to increasing `bias_coeff` is to change the shape of `coanc_subpops` (its overall scale does not change the minimum value)')

    # function whose zero we want!
    # "sigma" is the only parameter to optimize
    # need this function inside here as it depends on extra parameters
    bias_coeff_admix_objective <- function( x ) {
        # this is a transformation that helps us explore up to `sigma = Inf` in a compact space instead:
        sigma <- x / (1 - x)
        # x == 0 => sigma == 0
        # x == 1 => sigma == Inf
        # first, construct the admixture coefficients, which critically depend on sigma
        admix_proportions <- func(
            n_ind,
            k_subpops,
            sigma,
            coord_ind_first = coord_ind_first,
            coord_ind_last = coord_ind_last
        )
        # get bias coefficient, return difference from desired value!
        delta <- bias_coeff_admix(admix_proportions, coanc_subpops) - bias_coeff
        
        ## # hack to prevent issues when attempting to fit extreme values (at boundaries)
        ## # in particular, without this (and in particular in lower precision settings, such as --disable-long-double), uniroot below complains because delta has the same sign on both extrema (due to machine error, not a real sign issue).  If the sign wasn't checked then the tolerance would find that the solution is at the boundary, so here we force the check earlier.
        ## if ( abs(delta) < tol )
        ##     delta <- 0
        
        # return delta now
        return( delta )
    }

    # to avoid uniroot dying, check opposite ends and make sure sings are opposite (which is what that function checks and dies of if it fails)
    # I believe the issue is handled earlier, but just in case, directly assess this problematic case and return hack answers if that happens
    v0 <- bias_coeff_admix_objective( 0 ) # sigma =   0, s = admix_prop_bias_coeff_min
    v1 <- bias_coeff_admix_objective( 1 ) # sigma = Inf, s = 1
    # this function is monotonically increasing, so v0 should be negative and v1 positive
    # handle cases when this is not so
    if ( v0 >= 0 )
        return( 0 )
    if ( v1 <= 0 )
        return( Inf )
    
    # default tolerance of .Machine$double.eps^0.25 (~ 1e-4) was not good enough, reduced to ~ 2e-16
    obj <- stats::uniroot(
                      bias_coeff_admix_objective,
                      interval = c(0, 1),
                      extendInt = "no",
                      tol = tol
                  )
    # this is the value in terms of `x`
    x_root <- obj$root
    # transform back to sigma, return this
    x_root / (1 - x_root)
}
