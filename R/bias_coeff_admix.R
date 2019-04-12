bias_coeff_admix <- function(admix_proportions, coanc_subpops, weights = NULL) {
    # die if things are missing
    if (missing(admix_proportions))
        stop('`admix_proportions` is required!')
    if (missing(coanc_subpops))
        stop('`coanc_subpops` is required!')

    # in the interest of speed, does not validate matrix dimensions
    # (this is an internal function, so validations ought to happen earlier)
    
    # this function returns the bias coefficient "s" for arbitrary BN_PSD admixture models
    # includes an optimized implementation for most common case
    if (is.null(weights) && !is.matrix(coanc_subpops)) {
        # optimized version assumes:
        # - uniform weights
        # - vector coanc_subpops
        
        # mean coancestry using uniform weights
        mean_coancestry <- drop( coanc_subpops %*% colMeans(admix_proportions)^2 )
        # Fst is mean inbreeding (coancestry diagonal) using uniform weights
        Fst <- drop( coanc_subpops %*% colMeans(admix_proportions^2) )
        
    } else {
        # more general and cleaner algorithm
        # but slower due to matrix products

        # construct coancestry matrix, needed for next steps
        coancestry <- coanc_admix(admix_proportions, coanc_subpops)

        # now construct numerator and denominator of "s"
        # still have to handle the case of missing weights
        if (is.null(weights)) {
            mean_coancestry <- mean( coancestry )
            Fst <- mean( diag(coancestry) ) # Fst is mean inbreeding (coancestry diagonal)
        } else {
            # drop both to not have represented as matrices...
            mean_coancestry <- drop( weights %*% coancestry %*% weights ) # mean coancestry in general
            Fst <- drop( diag(coancestry) %*% weights ) # Fst is mean inbreeding (coancestry diagonal)
        }
    }

    # this is desired bias coefficient "s"
    mean_coancestry / Fst
}

