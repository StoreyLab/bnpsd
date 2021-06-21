#' Undifferentiate an allele distribution
#'
#' This function takes a vector of allele frequencies and an FST value, and returns a new distribution of allele frequencies that is consistent with reversing differentiation with the given FST, in the sense that the new distribution is more concentrated around the middle (0.5) than the input/original by an amount predicted from theory.
#' The new distribution is created by weighing the input distribution with a random mixing distribution with a lower variance.
#' An automatic method is provided that always selects a Beta distribution with just the right concentration to work for the given data and FST.
#' Explicit methods are also provided for more control, but are more likely to result in errors when mixing variances are not small enough (see details below).
#'
#' Model: Suppose we started from an allele frequency `p0` with expectation 0.5 and variance `V0`.
#' Differentiation creates a new allele frequency `p1` without changing its mean (`E(p1|p0) = p0`) and with a conditional variance given by FST `F`: `Var(p1|p0) = p0*(1-p0)*F`.
#' The total variance of the new distribution (calculated using the law of total variance) equals
#' `V1 = Var(p1) = F/4 + (1-F)*V0`.
#' (Also `E(p1) = 0.5`).
#' So the new variance is greater for `F>0` (note `V0 <= 1/4` for any distribution bounded in (0,1)).
#' Thus, given `V1` and `F`, the goal is to construct a new distribution with the original, lower variance of `V0 = (V1-F/4)/(1-F)`.
#' An error is thrown if `V1 < F/4` in input data, which is inconsistent with this assumed model.
#' 
#' Construction of "undifferentiated" allele frequencies:
#' `p3 = w*p + (1-w)*p2`, where `p` is the input with sample variance `V1` and `p2` is a random draw from the mixing distribution `distr` with expectation 0.5 and known variance `V2`.
#' The output variance is `V3 = w^2*V1 + (1-w)^2*V2`, which we set to the desired `V0 = (V1-F/4)/(1-F)` and solve for `w` (the largest of the two quadratic roots is used).
#' An error is thrown if `V2 > V0` (the output variance must be larger than the mixing variance).
#' This error is avoided by manually adjusting choice of `distr` and `alpha` (for `distr = "beta"`), or preferably with `distr = "auto"` (default), which selects a Beta distribution with `alpha = (1/(4*V0)-1)/2` that is guaranteed to work for any valid `V0` (assuming `V1 < F/4`).
#' 
#'
#' @param p A vector of observed allele frequencies.
#' @param F The FST value of the differentiation to reverse.
#' @param distr Name of the mixing distribution to use.
#' - "auto" picks a symmetric Beta distribution with parameters that ensure a small enough variance to succeed.
#' - "beta" is a symmetric Beta distribution with parameter `alpha` as provided below.
#' - "uniform" is a uniform distribution (same as "beta" with `alpha = beta = 1`).
#' - "point" is a distribution fully concentrated/fixed at 0.5 (same as the limit of "beta" with `alpha = beta = Inf`, which has zero variance).
#' @param alpha Shape parameter for `distr = "beta"`, ignored otherwise.
#'
#' @return A list with two named elements:
#' - `p`: A new vector of allele frequencies with the same length as `p`, with the desired variance (see below) obtained by weighing the input `p` with new random data from distribution `distr`.
#' - `w`: The weight used for the input data (`1-w` for the mixing distribution).
#'
#' @examples
#' # create random uniform data for these many loci
#' m <- 100
#' p <- runif( m )
#' # differentiate the distribution using Balding-Nichols model
#' F <- 0.1
#' nu <- 1 / F - 1
#' p2 <- rbeta( m, p * nu, (1 - p) * nu )
#'
#' # now undifferentiate with this function
#' # default "automatic" distribution recommended
#' # (avoids possible errors for specific distributions)
#' p3 <- undiff_af( p2, F )$p
#' 
#' # note p3 does not equal p (original is unrecoverable)
#' # but variances (assuming expectation is 0.5 for all) should be close to each other,
#' # and both be lower than p2's variance:
#' V1 <- mean( ( p - 0.5 )^2 )
#' V2 <- mean( ( p2 - 0.5 )^2 )
#' V3 <- mean( ( p3 - 0.5 )^2 )
#' # so p3 is stochastically consistent with p as far as the variance is concerned
#' 
#' @export
undiff_af <- function( p, F, distr = c('auto', 'uniform', 'beta', 'point'), alpha = 1 ) {
    # check inputs
    if ( missing( p ) )
        stop( '`p` is required!' )
    if ( missing( F ) )
        stop( '`F` is required!' )
    # this one has some choices
    distr <- match.arg( distr )

    # continue validations
    if ( !is.numeric( p ) )
        stop( '`p` must be numeric!' )
    if ( anyNA( p ) )
        stop( '`p` cannot have missing values!' )
    if ( min( p ) < 0 )
        stop( '`p` must have non-negative values!  Observed min: ', min( p ) )
    if ( max( p ) > 1 )
        stop( '`p` cannot exceed 1!  Observed max: ', max( p ) )

    if ( length( F ) != 1 )
        stop( '`F` must be scalar!' )
    if ( !is.numeric( F ) )
        stop( '`F` must be numeric!' )
    if ( is.na( F ) )
        stop( '`F` cannot be NA!' )
    if ( F < 0 )
        stop( '`F` must be non-negative!  Observed: ', F )
    if ( F > 1 )
        stop( '`F` cannot exceed 1!  Observed: ', F )
    
    # handle some potential trivial cases
    # nothing to do if asked for no differentiation
    if ( F == 0 )
        return( list( p = p, w = 1 ) )
    
    # variance of observed MAF distribution (assuming symmetry)
    # satisfied by construction:
    # 0 <= Vp <= 1/4
    Vp <- mean( ( p - 0.5 )^2 )
    # for assumed model to be correct, `Vp >= F / 4` is required!
    # make sure a higher `F` value wasn't passed, that has no solution
    if ( F > 4 * Vp )
        stop( '`F` (', F, ') cannot be larger than 4 times the observed MAF variance (', 4 * Vp, ')!  This violates model assumptions and results in negative output variance.  Please select a lower value for `F`!' )
    # desired output variance:
    Vo <- ( Vp - F / 4 ) / ( 1 - F )
    # assuming F is correct, this now satisfies
    # 0 <= Vo <= Vp <= 1/4
    # F / 4 <= Vp

    # "auto" is a hack where a Beta that is guaranteed to work is chosen, avoiding errors!
    if ( distr == 'auto' ) {
        # this is beta with `alpha` such that V2 == Vo
        alpha <- ( 1 / ( 4 * Vo ) - 1 ) / 2
        distr <- 'beta' # treat as this case from now on
    }

    # define second distribution's known variance
    # needed for calculating weight
    if ( distr == 'point' ) {
        V2 <- 0
    } else if ( distr == 'uniform' ) {
        V2 <- 1/12
    } else if ( distr == 'beta' ) {
        # beta-specific validations
        if ( alpha < 0 )
            stop( '`alpha` must be non-negative!' )
        # now actual variance formula
        V2 <- 1 / ( 4 * ( 2 * alpha + 1 ) )
    } else
        stop( 'Unimplemented distribution: ', distr )

    # catch an early issue if `F` is so large that the mixing distribution doesn't work
    # in other words, output variance must be at least as large as mixing variance
    if ( Vo < V2 )
        stop( 'Output variance (', Vo, ') is smaller than mixing variance (', V2, ')!  This violates model assumptions and can result in imaginary weights.  This can be fixed by either reducing `F` (to increase output variance) or by picking a mixing distribution with a lower variance than the desired output variance.' )
    # now we're at
    # 0 <= V2 <= Vo <= Vp <= 1/4
    # F / 4 <= Vp
    
    # normalize all variances by Vp now:
    Vo <- Vo / Vp
    V2 <- V2 / Vp

    # after this normalization, the inequalities satisfied are
    # 0 <= V2 <= Vo <= 1 <= 1/4 / Vp <= 1 / F
    
    # want to find `w` that solves:
    # Vo = w^2 + ( 1 - w )^2 * V2
    # Vo = w^2 + ( 1 + w^2 - 2 *w ) * V2
    # w^2 * ( 1 + V2 ) - 2 * w * V2 + (V2 - Vo) = 0
    ## a <- 1 + V2
    ## b <- - 2 * V2
    ## c <- V2 - Vo
    ## # manual quadratic solution
    ## det <- b^2 - 4 * a * c
    ## if ( det < 0 )
    ##     stop( 'determinant is negative!' )
    ## w <- ( -b + sqrt( det ) ) / ( 2 * a )

    # version with more direct variance notation, to see signs and potential issues
    # NOTES:
    # - this root gives weights in [0,1], the other root is mostly negative
    # - (this det is actually det/4, meh)
    # we have already required Vo >= V2, so determinant is guaranteed to be non-negative!
    ## det <- V2^2 - (1 + V2) * (V2 - Vo)
    det <- Vo - V2 * ( 1 - Vo )
    #det <- Vo * ( 1 + V2 ) - V2
    w <- ( V2 + sqrt( det ) ) / ( 1 + V2 )
    # more inequalities:
    # Vo^2 <= det <= Vo
    # (min attained at `V2 = Vo`)
    # (max attained at `V2 = 0`)
    # plots suggest w(V2) is overall strictly decreasing (V2 decreases `det`, but increases the other bound otherwise, so it's not obvious)
    # if so, then upper bound is at V2 = 0, and lower bound is at V2 = Vo:
    # 0 <= 2 * Vo / ( 1 + Vo ) <= w <= sqrt( Vo ) <= 1
    # verified correctness empirically!
    # these are tightest bounds!
    
    # I'm a bit surprised when `V2 = Vo` we don't just have `w = 0`
    # saw other root has this behavior, but otherwise other roon is inadvisable because it is mostly negative, so better stick with this branch
    # so this branch/root favors weighing input data non-zero generally, which is best when desired variance exceeds the mixing variance (expected)
    
    # once weight is determined and no errors occur, now we can bother to draw `p2` (could be a large number of them, so do last)
    # number of loci
    m <- length( p )
    if ( distr == 'point' ) {
        p2 <- 0.5 # ok as scalar
    } else if ( distr == 'uniform' ) {
        p2 <- stats::runif( m )
    } else if ( distr == 'beta' ) {
        p2 <- stats::rbeta( m, alpha, alpha )
    } else
        stop( 'Unimplemented distribution: ', distr )
    
    # undifferentiate by mixing!
    q <- p * w + ( 1 - w ) * p2
    
    return(
        list(
            p = q,
            w = w
        )
    )
}

# for internal tests
# the forward differentiation function, in this case a simple BN step (could be much more complicated, meh)
diff_af <- function( p, F ) {
    # redifferentiate for test
    nu <- 1 / F - 1
    m <- length( p )
    # tests suggests that this vectorization (with alpha and beta vectors) works correctly!
    p2 <- stats::rbeta( m, p * nu, (1 - p) * nu )
    return( p2 )
}

