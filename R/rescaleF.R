## Rescale subpopulation inbreeding vector to give a desired FST for the admixed individuals
##
## Given a final admixture proportion matrix \eqn{Q} for \eqn{n} individuals and \eqn{k} intermediate subpopulations, the vector of intermediate inbreeding coefficients \eqn{F} (per-subpopulation \eqn{F_{ST}}{FST}'s) assumed to be incorrectly scaled, and the desired final \eqn{F_{ST}}{FST} of the admixed individuals, this function returns the correctly-scaled \eqn{F} vector.
## 
## The idea is that we have determined the admixture structure \eqn{Q}, and \eqn{F} in a relative scale, but now we want to get \eqn{F} in the correct scale so the final \eqn{F_{ST}}{FST} is reasonable.
## After rescaling, the function stops with an error if \eqn{F} has a maximum value beyond 1.
##
## @param Q The \eqn{n \times k}{n-by-k} admixture proportion matrix
## @param F The length-\eqn{k} vector of ancestral inbreeding coefficients assumed to be in the wrong scale
## @param Fst The desired final \eqn{F_{ST}}{FST} of the admixed individuals
## @param w The length-\eqn{n} vector of weights for individuals that define \eqn{F_{ST}}{FST} (default uniform weights)
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
## F <- rescaleF(Q, F, Fst)
rescaleF <- function(Q, F, Fst, w) {
    ## die informatively...
    if (missing(Q)) stop('Fatal: Q is missing!')
    if (missing(F)) stop('Fatal: F is missing!')
    if (missing(Fst)) stop('Fatal: Fst is missing!')
    ## calculate necessary intermediates
    ## could probably be more efficient skipping un-needed matrix products, but such efficiency is not needed here because this function is usually called once only and n is big but not huge
    fstWrong <- fst(Q,F,w) # get Fst under this model (wrong scale, yields adjustment!)
    ## fstWrong <- drop( F %*% colMeans(Q^2) ) # current Fst (this only works for uniform weights)
    F <- F * Fst/fstWrong # this fixes scale
    ## unfortunately, some F (particularly very non-uniform vectors) may be rescaled to values greater than 1, which would be disallowed under the probabilistic inbreeding framework.  Check and stop if needed!
    if (max(F)>1) stop('Fatal: after rescaling for Fst=', Fst, ', some ancestral Fsts became greater than 1:', "\n", paste(F, collapse=' '))
    F # return rescaled version if things were good!
}

