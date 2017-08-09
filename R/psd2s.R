## internal function, returns the bias coefficient "s" given other q1d PSD params
## assumes uniform weights for Fst and thetaBar!!!
## F is vector of Fst's up to a scaling constant (it cancels out anyway)
psd2s <- function(sigma, F, n) {
    ## this function returns the bias coefficient "s" for my admixture model given sigma (which we hope to optimize later) and other parameters
    ## internally we use the construction that gives uniform weights, so they are assumed!

    ## derived parameter
    k <- length(F)

    ## first, construct the admixture coefficients, which critically depend on sigma
    Q <- q1d(n, k, sigma)

    ## COMMENTED OUT: clearer but slower code (due to matrix products), optimized below...
    ## construct coancestry matrix, needed for next steps
    ##    Theta <- coanc(Q,F)
    ## now construct numerator and denominator of "s"
    ##    thetaBar <- mean(Theta) # mean coancestry using uniform weights
    ##    fst <- mean(diag(Theta)) # Fst is mean inbreeding (coancestry diagonal) using uniform weights

    ## mean coancestry using uniform weights
    thetaBar <- drop( F %*% colMeans(Q)^2 )
    ## Fst is mean inbreeding (coancestry diagonal) using uniform weights
    fst <- drop( F %*% colMeans(Q^2) )
    thetaBar/fst # this is desired "s"
}

