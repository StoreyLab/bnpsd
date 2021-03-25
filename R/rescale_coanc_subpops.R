# Rescale subpopulation inbreeding vector to give a desired FST for the admixed individuals
#
# Given a final admixture proportion matrix for `n` individuals and `k` intermediate subpopulations, the coancestry of intermediate subpopulations assumed to be incorrectly scaled, and the desired final FST of the admixed individuals, this function returns the correctly-scaled intermediate coancestry matrix (or its special cases).
# 
# The idea is that we have determined the admixture structure, and `coanc_subpops` in a relative scale, but now we want to get `coanc_subpops` in the correct scale so the final FST is reasonable.
# After rescaling, the function stops with an error if `coanc_subpops` has a maximum value beyond 1.
#
# @param admix_proportions The `n`-by-`k` admixture proportion matrix.
# @param coanc_subpops The length-`k` vector of intermediate subpopulation inbreeding coefficients, or `k`-by-`k` intermediate coancestry matrix, assumed to be in the wrong scale.
# @param fst The desired final FST of the admixed individuals.
# @param weights The length-`n` vector of weights for individuals that define FST (default uniform weights).
#
# @return A named list containing
# - `coanc_subpops`: The rescaled intermediate subpopulation inbreeding coefficient vector
# - `factor`: the multiplicative factor used to rescale `coanc_subpops`
#
# @examples
# # set desired parameters
# # number of individuals
# n_ind <- 1000
# # number of intermediate subpops
# k_subpops <- 10
# # desired FST
# Fst <- 0.1
# 
# # differentiation of subpops (relative (wrong) scale)
# coanc_subpops <- ( 1 : k_subpops ) / k_subpops
#
# # construct admixture proportions
# admix_proportions <- admix_prop_1d_linear(n_ind, k_subpops, sigma = 1)
#
# # lastly, rescale coanc_subpops to give desired FST!!!
# obj <- rescale_coanc_subpops(admix_proportions, coanc_subpops, Fst)
# # rescaled vector or matrix
# obj$coanc_subpops
# # and factor used
# obj$factor
#
rescale_coanc_subpops <- function(admix_proportions, coanc_subpops, fst, weights = NULL) {
    # die informatively...
    if (missing(admix_proportions))
        stop('`admix_proportions` is required!')
    if (missing(coanc_subpops))
        stop('`coanc_subpops` is required!')
    if (missing(fst))
        stop('`fst` is required!')

    # NOTE: we let fst_admix further validate inputs
    
    # calculate necessary intermediates
    # could probably be more efficient skipping un-needed matrix products, but such efficiency is not needed here because this function is usually called once only and n is big but not huge

    # get fst under this model (wrong scale, yields adjustment!)
    fst_0 <- fst_admix( admix_proportions, coanc_subpops, weights)
    # calculate multiplicative factor that corrects coanc_subpops (will return separately)
    factor <- fst / fst_0
    # this fixes scale
    coanc_subpops <- coanc_subpops * factor
    
    # unfortunately, some coanc_subpops may be rescaled to values greater than 1, which would be disallowed under the probabilistic inbreeding framework.
    # Check and stop if needed!
    if ( max( coanc_subpops ) > 1 )
        stop( 'Rescaling for `fst = ', fst, '` resulted in max(coanc_subpops) = ', max(coanc_subpops), ' > 1.  This can be fixed by changing input `coanc_subpops` or lowering `fst`.' )

    # return rescaled version if things were good!
    # factor too
    return(
        list(
            coanc_subpops = coanc_subpops,
            factor = factor
        )
    )
}

