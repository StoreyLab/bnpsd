## Rescale ancestral inbreeding to give desired final FST
##
## Given a final admixture proportion matrix \eqn{Q} for \eqn{n} individuals and \eqn{k} ancestral populations, the vector of ancestral inbreeding coefficients \eqn{F} (ancestral \eqn{F_{ST}}{FST}'s) that is assumed to be incorrectly scaled, and the desired final \eqn{F_{ST}}{FST} of the admixed individuals, this function returns the correctly-scaled \eqn{F} vector.
##
## The idea is that we have determined the admixture structure \eqn{Q}, and \eqn{F} in a relative scale, but now we want to get \eqn{F} in the correct scale so the final \eqn{F_{ST}}{FST} is reasonable.  Note a given \eqn{F} may be rescaled to have a maximum value beyond 1, the function stops with an error if this happens.
## The function assumes uniform weights for \eqn{F_{ST}}{FST}.
##
## @param Q The \eqn{n \times k}{n-by-k} admixture proportion matrix
## @param F The length-\eqn{k} vector of ancestral inbreeding coefficients assumed to be in the wrong scale
## @param Fst The desired final \eqn{F_{ST}}{FST} of the admixed individuals
##
## @return The rescaled ancestral inbreeding coefficients
##
## @examples
## # set desired parameters
## n <- 1000 # number of individuals
## k <- 10 # number of ancestral populations
## s <- 0.5 # desired bias coefficient
## Fst <- 0.1 # desired FST
## # differentiation of ancestral populations (relative (wrong) scale)
## F <- (1:k)/k
## # set bias coeff. of 1/2, F and n, implicitly assumes our 1D scenario
## sigma <- solveSigma(s=s, F=F, n=n)
## # construct final admixture proportions
## Q <- q1d(n=n, k=k, sigma=sigma)
## # lastly, rescale F to give desired Fst!!!
## F <- rescaleFst(Q, F, Fst)
rescaleFst <- function(Q, F, Fst) {
    ## given a PSD model with all but the scale wrong, and a desired Fst, returns rescaled ancestral Fsts!
    ## Q is final, but F is the weights with wrong scale (ancestral differentiations), Fst is desired Fst
    ## assumes uniform weights are desired (a common assumption for most of the code in this package)
    ## also assumes vector F
    fstWrong <- drop( F %*% colMeans(Q^2) ) # current Fst
    F <- F * Fst/fstWrong # this fixes scale
    ## unfortunately, some F (particularly very non-uniform vectors) may be rescaled to values greater than 1, which would be disallowed under the probabilistic inbreeding framework.  Check if this is so, with a fatal message if needed!
    if (max(F)>1) stop('Fatal: after rescaling for Fst=', Fst, ', some ancestral Fsts became greater than 1:', "\n", paste(F, collapse=' '))
    F # return rescaled version if things were good!
}

