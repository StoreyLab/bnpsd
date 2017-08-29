biasCoeff <- function(Q, F, w) {
    ## this function returns the bias coefficient "s" for arbitrary BN_PSD admixture models
    ## includes an optimized implementation for most common case
    if (missing(w) && !is.matrix(F)) {
        ## optimized version for uniform weights and vector F
        
        ## mean coancestry using uniform weights
        thetaBar <- drop( F %*% colMeans(Q)^2 )
        ## Fst is mean inbreeding (coancestry diagonal) using uniform weights
        fst <- drop( F %*% colMeans(Q^2) )
        
    } else {
        ## more general and cleaner algorithm
        ## but slower due to matrix products
        
        ## construct coancestry matrix, needed for next steps
        Theta <- coanc(Q,F)
        ## now construct numerator and denominator of "s"
        ## drop both to not have represented as matrices...
        thetaBar <- drop( w %*% Theta %*% w ) # mean coancestry in general
        fst <- drop( diag(Theta) %*% w ) # Fst is mean inbreeding (coancestry diagonal)
    }
    
    thetaBar/fst # this is desired "s"
}

