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
#' `p_out = w*p_in + (1-w)*p_mix`, where `p_in` is the input with sample variance `V_in` (`V1` in above model) and `p_mix` is a random draw from the mixing distribution `distr` with expectation 0.5 and known variance `V_mix`.
#' The output variance is `V_out = w^2*V_in + (1-w)^2*V_mix`, which we set to the desired `V_out = (V_in-F/4)/(1-F)` (`V0` in above model) and solve for `w` (the largest of the two quadratic roots is used).
#' An error is thrown if `V_mix > V_out` (the output variance must be larger than the mixing variance).
#' This error is avoided by manually adjusting choice of `distr` and `alpha` (for `distr = "beta"`), or preferably with `distr = "auto"` (default), which selects a Beta distribution with `alpha = (1/(4*V_out)-1)/2 + eps` that is guaranteed to work for any valid `V_out` (assuming `V_in < F/4`).
#' 
#'
#' @param p A vector of observed allele frequencies.
#' @param F The FST value of the differentiation to reverse.
#' @param distr Name of the mixing distribution to use.
#' - "auto" picks a symmetric Beta distribution with parameters that ensure a small enough variance to succeed.
#' - "beta" is a symmetric Beta distribution with parameter `alpha` as provided below.
#' - "uniform" is a uniform distribution (same as "beta" with `alpha = 1`).
#' - "point" is a distribution fully concentrated/fixed at 0.5 (same as the limit of "beta" with `alpha = Inf`, which has zero variance).
#' @param alpha Shape parameter for `distr = "beta"`, ignored otherwise.
#' @param eps If `distr = "auto"`, this small value is added to the calculated `alpha` to avoid roundoff errors and ensuring that the mixing variance is smaller than the maximum allowed.
#'
#' @return A list with the new distribution and several other informative statistics, which are named elements:
#' - `p`: A new vector of allele frequencies with the same length as input `p`, with the desired variance (see details) obtained by weighing the input `p` with new random data from distribution `distr`.
#' - `w`: The weight used for the input data (`1-w` for the mixing distribution).
#' - `F_max`: The maximum FST possible for undifferentiating this data (equals four times the input variance (see details), which results in zero output variance).
#' - `V_in`: sample variance of input `p`, assuming its expectation is 0.5.
#' - `V_out`: target variance of output `p`.
#' - `V_mix`: variance of mixing distribution.
#' - `alpha`: the value of `alpha` used for symmetric Beta mixing distribution, informative if `distr = "auto"`.
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
#' # note p3 does not equal p exactly (original is unrecoverable)
#' # but variances (assuming expectation is 0.5 for all) should be close to each other,
#' # and both be lower than p2's variance:
#' V1 <- mean( ( p - 0.5 )^2 )
#' V2 <- mean( ( p2 - 0.5 )^2 )
#' V3 <- mean( ( p3 - 0.5 )^2 )
#' # so p3 is stochastically consistent with p as far as the variance is concerned
#' 
#' @export
undiff_af <- function(
                      p,
                      F,
                      distr = c('auto', 'uniform', 'beta', 'point'),
                      alpha = 1,
                      eps = 10 * .Machine$double.eps
                      ) {
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
    # 0 <= V_in <= 1/4
    V_in <- mean( ( p - 0.5 )^2 )
    # for assumed model to be correct, `V_in >= F / 4` is required!
    # make sure a higher `F` value wasn't passed, that has no solution
    F_max <- 4 * V_in # return too, useful info
    if ( F > F_max )
        stop( '`F` (', F, ') cannot be larger than 4 times the observed MAF variance (', F_max, ')!  This violates model assumptions and results in negative output variance.  Please select a lower value for `F`!' )
    # desired output variance:
    V_out <- ( V_in - F / 4 ) / ( 1 - F )
    # assuming F is correct, this now satisfies
    # 0 <= V_out <= V_in <= 1/4
    # F / 4 <= V_in

    # "auto" is a hack where a Beta that is guaranteed to work is chosen, avoiding errors!
    if ( distr == 'auto' ) {
        # this is beta with `alpha` such that V_mix == V_out.
        # to avoid roundoff errors, it is best to choose a slightly larger alpha
        # (more concentrated around 0.5, i.e. lower V_mix; V_mix <= V_out is strictly required!)
        alpha <- ( 1 / ( 4 * V_out ) - 1 ) / 2 + eps
        distr <- 'beta' # treat as this case from now on
    }

    # define second distribution's known variance
    # needed for calculating weight
    if ( distr == 'point' ) {
        V_mix <- 0
        alpha <- Inf
    } else if ( distr == 'uniform' ) {
        V_mix <- 1/12
        alpha <- 1
    } else if ( distr == 'beta' ) {
        # beta-specific validations
        if ( alpha < 0 )
            stop( '`alpha` must be non-negative!' )
        # now actual variance formula
        V_mix <- 1 / ( 4 * ( 2 * alpha + 1 ) )
    } else
        stop( 'Unimplemented distribution: ', distr )

    # catch an early issue if `F` is so large that the mixing distribution doesn't work
    # in other words, output variance must be at least as large as mixing variance
    if ( V_out < V_mix )
        stop( 'Output variance (', V_out, ') is smaller than mixing variance (', V_mix, ')!  This violates model assumptions and can result in imaginary weights.  This can be fixed by either reducing `F` (to increase output variance) or by picking a mixing distribution with a lower variance than the desired output variance.' )
    # now we're at
    # 0 <= V_mix <= V_out <= V_in <= 1/4
    # F / 4 <= V_in
    
    # normalize all variances by V_in now:
    V_out <- V_out / V_in
    V_mix <- V_mix / V_in

    # after this normalization, the inequalities satisfied are
    # 0 <= V_mix <= V_out <= 1 <= 1/4 / V_in <= 1 / F
    
    # want to find `w` that solves:
    # V_out = w^2 + ( 1 - w )^2 * V_mix
    # V_out = w^2 + ( 1 + w^2 - 2 *w ) * V_mix
    # w^2 * ( 1 + V_mix ) - 2 * w * V_mix + (V_mix - V_out) = 0
    ## a <- 1 + V_mix
    ## b <- - 2 * V_mix
    ## c <- V_mix - V_out
    ## # manual quadratic solution
    ## det <- b^2 - 4 * a * c
    ## if ( det < 0 )
    ##     stop( 'determinant is negative!' )
    ## w <- ( -b + sqrt( det ) ) / ( 2 * a )

    # version with more direct variance notation, to see signs and potential issues
    # NOTES:
    # - this root gives weights in [0,1], the other root is mostly negative
    # - (this det is actually det/4, meh)
    # we have already required V_out >= V_mix, so determinant is guaranteed to be non-negative!
    ## det <- V_mix^2 - (1 + V_mix) * (V_mix - V_out)
    det <- V_out - V_mix * ( 1 - V_out )
    #det <- V_out * ( 1 + V_mix ) - V_mix
    w <- ( V_mix + sqrt( det ) ) / ( 1 + V_mix )
    # more inequalities:
    # V_out^2 <= det <= V_out
    # (min attained at `V_mix = V_out`)
    # (max attained at `V_mix = 0`)
    # plots suggest w(V_mix) is overall strictly decreasing (V_mix decreases `det`, but increases the other bound otherwise, so it's not obvious)
    # if so, then upper bound is at V_mix = 0, and lower bound is at V_mix = V_out:
    # 0 <= 2 * V_out / ( 1 + V_out ) <= w <= sqrt( V_out ) <= 1
    # verified correctness empirically!
    # these are tightest bounds!
    
    # I'm a bit surprised when `V_mix = V_out` we don't just have `w = 0`
    # saw other root has this behavior, but otherwise other roon is inadvisable because it is mostly negative, so better stick with this branch
    # so this branch/root favors weighing input data non-zero generally, which is best when desired variance exceeds the mixing variance (expected)
    
    # once weight is determined and no errors occur, now we can bother to draw `p_mix` (could be a large number of them, so do last)
    # number of loci
    m <- length( p )
    if ( distr == 'point' ) {
        p_mix <- 0.5 # ok as scalar
    } else if ( distr == 'uniform' ) {
        p_mix <- stats::runif( m )
    } else if ( distr == 'beta' ) {
        p_mix <- stats::rbeta( m, alpha, alpha )
    } else
        stop( 'Unimplemented distribution: ', distr )
    
    # undifferentiate by mixing!
    p_out <- p * w + ( 1 - w ) * p_mix
    
    return(
        list(
            p = p_out,
            w = w,
            F_max = F_max,
            V_in = V_in,
            V_out = V_out * V_in, # unnormalize
            V_mix = V_mix * V_in, # unnormalize
            alpha = alpha
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

