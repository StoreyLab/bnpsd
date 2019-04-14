# Rescale subpopulation inbreeding vector to give a desired FST for the admixed individuals
#
# Given a final admixture proportion matrix for \eqn{n} individuals and \eqn{k} intermediate subpopulations, the coancestry of intermediate subpopulations assumed to be incorrectly scaled, and the desired final \eqn{F_{ST}}{FST} of the admixed individuals, this function returns the correctly-scaled intermediate coancestry matrix (or its special cases).
# 
# The idea is that we have determined the admixture structure, and \eqn{coanc_subpops} in a relative scale, but now we want to get \eqn{coanc_subpops} in the correct scale so the final \eqn{F_{ST}}{FST} is reasonable.
# After rescaling, the function stops with an error if \eqn{coanc_subpops} has a maximum value beyond 1.
#
# @param admix_proportions The \eqn{n \times k}{n-by-k} admixture proportion matrix
# @param coanc_subpops The length-\eqn{k} vector of intermediate subpopulation inbreeding coefficients, assumed to be in the wrong scale
# @param fst The desired final \eqn{F_{ST}}{FST} of the admixed individuals
# @param weights The length-\eqn{n} vector of weights for individuals that define \eqn{F_{ST}}{FST} (default uniform weights)
#
# @return The rescaled intermediate subpopulation inbreeding coefficient vector
#
# @examples
# # set desired parameters
# n <- 1000 # number of individuals
# k <- 10 # number of intermediate subpops
# sigma <- 1 # ...
# Fst <- 0.1 # desired FST
# # differentiation of subpops (relative (wrong) scale)
# coanc_subpops <- (1:k) / k
# # construct admixture proportions
# admix_proportions <- admix_prop_1d_linear(n, k, sigma)
# # lastly, rescale coanc_subpops to give desired FST!!!
# coanc_subpops <- rescale_coanc_subpops(admix_proportions, coanc_subpops, Fst)
rescale_coanc_subpops <- function(admix_proportions, coanc_subpops, fst, weights = NULL) {
    # die informatively...
    if (missing(admix_proportions))
        stop('`admix_proportions` is required!')
    if (missing(coanc_subpops))
        stop('`coanc_subpops` is required!')
    if (missing(fst))
        stop('`fst` is required!')
    
    # calculate necessary intermediates
    # could probably be more efficient skipping un-needed matrix products, but such efficiency is not needed here because this function is usually called once only and n is big but not huge
    fst_0 <- fst_admix( admix_proportions, coanc_subpops, weights) # get fst under this model (wrong scale, yields adjustment!)
    coanc_subpops <- coanc_subpops * fst / fst_0 # this fixes scale
    # unfortunately, some coanc_subpops may be rescaled to values greater than 1, which would be disallowed under the probabilistic inbreeding framework.
    # Check and stop if needed!
    if (max(coanc_subpops) > 1)
        stop('Rescaling for `fst = ', fst, '` resulted in max(coanc_subpops) = ', max(coanc_subpops), ' > 1')
    coanc_subpops # return rescaled version if things were good!
}

