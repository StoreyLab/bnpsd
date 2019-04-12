context('bnpsd')

# start with lower-level/internal tests, more informative that higher-level function errors

test_that("fixed_loci works in toy cases", {
    # here's a toy matrix
    X <- matrix(
        data = c(
            2, 2, NA, # fixed locus (with one missing element)
            0, NA, 0, # another fixed locus, for opposite allele
            1, 1, 1, # NOT fixed (heterozygotes are not considered fixed)
            0, 1, 2, # a completely variable locus
            NA, NA, NA # completely missing locus (will be treated as fixed)
        ),
        ncol = 3,
        byrow = TRUE
    )
    # test that we get the desired values
    expect_equal(
        fixed_loci(X), c(TRUE, TRUE, FALSE, FALSE, TRUE)
    )
})

test_that("coanc_admix works in toy cases", {
    # set up vars
    Q <- diag(c(1, 1)) # an independent subpops model with two subpops
    F <- c(0.1, 0.3)
    
    # die when the important parameters are missing
    expect_error( coanc_admix() ) # all missing 
    expect_error( coanc_admix(Q) ) # coanc_subpops missing
    expect_error( coanc_admix(coanc_subpops = F) ) # admix_proportions missing
    
    ThetaExp <- diag(F) # the Theta we expect for this setup
    Theta <- coanc_admix(Q, F)
    expect_equal(Theta, ThetaExp)

    # same Q, scalar F
    F <- 0.2
    ThetaExp <- diag(c(F, F)) # the Theta we expect for this setup
    Theta <- coanc_admix(Q, F)
    expect_equal(Theta, ThetaExp)

    # same Q, matrix F
    Fv <- c(0.1, 0.4) # vector version
    F <- diag(Fv) # matrix version
    Theta <- coanc_admix(Q, F)
    expect_equal(Theta, F) # F is the theta we expect in this case

    # most complex case, just a general math check
    Q <- matrix(c(0.7, 0.3, 0.2, 0.8), nrow = 2, byrow = TRUE)
    F <- matrix(c(0.3, 0.1, 0.1, 0.3), nrow = 2, byrow = TRUE)
    ThetaExp <- Q %*% F %*% t(Q) # the Theta we expect for this setup (slower but more explicit version)
    Theta <- coanc_admix(Q, F)
    expect_equal(Theta, ThetaExp)
}) 

test_that("coanc_to_kinship works in toy cases", {
    Q <- diag(c(1, 1)) # an IS model with two subpops
    F <- c(0.1, 0.3)
    PhiExp <- diag( (1+F) / 2 ) # the Phi we expect for this setup
    Theta <- coanc_admix(Q, F)
    Phi <- coanc_to_kinship(Theta)
    expect_equal(Phi, PhiExp)

    # same Q, scalar F
    F <- 0.2
    K <- (1+F) / 2 # transform at this stage
    PhiExp <- diag(c(K, K)) # the Phi we expect for this setup
    Theta <- coanc_admix(Q, F)
    Phi <- coanc_to_kinship(Theta)
    expect_equal(Phi, PhiExp)

    # same Q, matrix F
    Fv <- c(0.1, 0.4) # vector version
    F <- diag(Fv) # matrix version
    PhiExp <- diag((1+Fv) / 2)
    Theta <- coanc_admix(Q, F)
    Phi <- coanc_to_kinship(Theta)
    expect_equal(Phi, PhiExp)

    # most complex case, just a general math check
    Q <- matrix(c(0.7, 0.3, 0.2, 0.8), nrow = 2, byrow = TRUE)
    F <- matrix(c(0.3, 0.1, 0.1, 0.3), nrow = 2, byrow = TRUE)
    PhiExp <- Q %*% F %*% t(Q) # the Theta we expect for this setup (slower but more explicit version)
    diag(PhiExp) <- (1 + diag(PhiExp))/2 # explicit transformation to kinship
    Theta <- coanc_admix(Q, F)
    Phi <- coanc_to_kinship(Theta)
    expect_equal(Phi, PhiExp)
}) 

test_that("fst works in toy cases", {
    Q <- diag(c(1, 1)) # an IS model with two subpops
    F <- c(0.1, 0.3)
    fst1 <- mean(F) # the Theta we expect for this setup
    fst2 <- fst(Q, F)
    expect_equal(fst1, fst2)

    # same Q, scalar F
    F <- 0.2
    fst2 <- fst(Q, F)
    expect_equal(F, fst2)

    # same Q, matrix F
    Fv <- c(0.1, 0.4) # vector version
    F <- diag(Fv) # matrix version
    fst1 <- mean(Fv)
    fst2 <- fst(Q, F)
    expect_equal(fst1, fst2) # F is the theta we expect in this case

    # most complex case, just a general math check
    Q <- matrix(c(0.7, 0.3, 0.2, 0.8), nrow = 2, byrow = TRUE)
    F <- matrix(c(0.3, 0.1, 0.1, 0.3), nrow = 2, byrow = TRUE)
    fst1 <- mean(diag(Q %*% F %*% t(Q))) # the Fst we expect for this setup (slower but more explicit version)
    fst2 <- fst(Q, F)
    expect_equal(fst1, fst2)
})

test_that("bias_coeff_admix agrees with explicitly calculated bias coeff s", {
    # set up some simulated data
    n <- 5
    k <- 2
    sigma <- 1
    Q <- q1d(n, k, sigma)
    F <- 1:k # scale doesn't matter right now...

    Theta <- coanc_admix(Q, F) # in wrong scale but meh
    sWant <- mean(Theta) / mean(diag(Theta)) # this is the correct bias coeff, with uniform weights
    s <- bias_coeff_admix(Q, F) # calculation to compare to
    expect_equal(s, sWant)
    expect_true(s > 0) # other obvious properties...
    expect_true(s <= 1)

    # repeat with matrix F, missing (uniform) weights
    F_mat <- diag(F)
    s <- bias_coeff_admix(Q, F_mat) # calculation to compare to
    expect_equal(s, sWant)
    expect_true(s > 0) # other obvious properties...
    expect_true(s <= 1)
    
    # repeat with non-uniform weights...
    weights <- runif(n) # random weights for given number of individuals
    weights <- weights / sum(weights) # normalize to add up to 1! # NOTE: should check sum(weights)!= 0, meh...
    sWant <- drop(weights %*% Theta %*% weights) / drop( diag(Theta) %*% weights ) # this is the correct bias coeff, with uniform weights
    s <- bias_coeff_admix(Q, F, weights) # calculation to compare to
    expect_equal(s, sWant)
    expect_true(s > 0) # other obvious properties...
    expect_true(s <= 1)

    # repeat with matrix F and non-uniform weights
    s <- bias_coeff_admix(Q, F_mat, weights) # calculation to compare to
    expect_equal(s, sWant)
    expect_true(s > 0) # other obvious properties...
    expect_true(s <= 1)
})

test_that("qis returns valid admixture coefficients", {
    labs <- c(1, 2, 2, 3, 3, 3, 4, 4, 4, 4)
    n <- length(labs)
    k <- length(unique(labs))
    Q <- qis(labs)
    # general tests for admixture matrices
    expect_equal(nrow(Q), n) # n rows
    expect_equal(ncol(Q), k) # k columns
    expect_true(all(Q >= 0)) # all are non-negative
    expect_true(all(Q <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(Q), rep.int(1, n)) # rows sum to 1, vector length n
    # specific tests for qis
    expect_true(all(Q %in% c(TRUE, FALSE)))
    expect_true(all(colnames(Q) == sort(unique(labs))))
    
    # test with provided subpops
    subpops <- 4:1
    Q <- qis(labs, subpops)
    # general tests for admixture matrices
    expect_equal(nrow(Q), n) # n rows
    expect_equal(ncol(Q), k) # k columns
    expect_true(all(Q >= 0)) # all are non-negative
    expect_true(all(Q <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(Q), rep.int(1,n)) # rows sum to 1, vector length n
    # specific tests for qis
    expect_true(all(Q %in% c(TRUE, FALSE)))
    expect_true(all(colnames(Q) == subpops))
    
    # test with provided subpops (additional labels)
    k <- 10
    subpops <- 1:k
    Q <- qis(labs, subpops)
    # general tests for admixture matrices
    expect_equal(nrow(Q), n) # n rows
    expect_equal(ncol(Q), k) # k columns
    expect_true(all(Q >= 0)) # all are non-negative
    expect_true(all(Q <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(Q), rep.int(1, n)) # rows sum to 1, vector length n
    # specific tests for qis
    expect_true(all(Q %in% c(TRUE, FALSE)))
    expect_true(all(colnames(Q) == subpops))

    # test with provided subpops (missing labels, must die!)
    subpops <- 1:3 # missing 4!
    expect_error( qis(labs, subpops) )
})

test_that("q1d returns valid admixture coefficients", {
    n <- 10
    k <- 2
    Q <- q1d(n, k, sigma = 1)
    expect_equal(nrow(Q), n) # n rows
    expect_equal(ncol(Q), k) # k columns
    expect_true(all(Q >= 0)) # all are non-negative
    expect_true(all(Q <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(Q), rep.int(1, n)) # rows sum to 1, vector length n

    # test with sigma == 0 (special case that makes usual formula break)
    Q <- q1d(n, k, sigma = 0)
    expect_equal(nrow(Q), n) # n rows
    expect_equal(ncol(Q), k) # k columns
    expect_true(all(Q >= 0)) # all are non-negative
    expect_true(all(Q <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(Q), rep.int(1, n)) # rows sum to 1, vector length n
    # in this case it should equal independent subpopulations
    labs <- c( rep.int(1, 5), rep.int(2, 5) ) # two subpops
    Q2 <- qis(labs)
    dimnames(Q2) <- NULL # before comparing, must toss column names
    expect_equal(Q, Q2)

    # test s version
    obj <- q1d(n, k, s = 0.5, F = 1:k, Fst = 0.1)
    Q <- obj$Q # returns many things in this case, get Q here
    expect_equal(nrow(Q), n) # n rows
    expect_equal(ncol(Q), k) # k columns
    expect_true(all(Q >= 0)) # all are non-negative
    expect_true(all(Q <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(Q), rep.int(1, n)) # rows sum to 1, vector length n
})

test_that("q1dc returns valid admixture coefficients", {
    n <- 10
    k <- 2
    Q <- q1dc(n, k, sigma = 1)
    expect_equal(nrow(Q), n) # n rows
    expect_equal(ncol(Q), k) # k columns
    expect_true(all(Q >= 0)) # all are non-negative
    expect_true(all(Q <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(Q), rep.int(1, n)) # rows sum to 1, vector length n

    # test with sigma == 0 (special case that makes usual formula break)
    Q <- q1dc(n, k, sigma = 0)
    expect_equal(nrow(Q), n) # n rows
    expect_equal(ncol(Q), k) # k columns
    expect_true(all(Q >= 0)) # all are non-negative
    expect_true(all(Q <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(Q), rep.int(1, n)) # rows sum to 1, vector length n
    # though the result is nearly island-like, there is an annoying shift I'd rather not try to figure out for this test...

    # test s version
    obj <- q1dc(n, k, s = 0.5, F = 1:k, Fst = 0.1)
    Q <- obj$Q # returns many things in this case, get Q here
    expect_equal(nrow(Q), n) # n rows
    expect_equal(ncol(Q), k) # k columns
    expect_true(all(Q >= 0)) # all are non-negative
    expect_true(all(Q <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(Q), rep.int(1, n)) # rows sum to 1, vector length n
})

test_that("bias_coeff_admix_fit agrees with reverse func", {
    n <- 1000
    inbr_subpops <- c(0.1, 0.2, 0.3)
    k <- length(inbr_subpops)
    s_want <- 0.5

    # test with q1d
    sigma <- bias_coeff_admix_fit(bias_coeff = s_want, inbr_subpops = inbr_subpops, n_ind = n, func = q1d) # get sigma
    # construct everything and verify s == s_want
    Q <- q1d(n, k, sigma) # now get Q from there
    Theta <- coanc_admix(Q, inbr_subpops)
    s <- mean(Theta) / mean(diag(Theta)) # this is the correct bias coeff, with uniform weights
    expect_equal(s, s_want)
    # since we set 0 < s_want < 1, nothing else to test
    
    # test with q1dc
    sigma <- bias_coeff_admix_fit(bias_coeff = s_want, inbr_subpops = inbr_subpops, n_ind = n, func = q1dc) # get sigma
    # construct everything and verify s == s_want
    Q <- q1dc(n, k, sigma) # now get Q from there
    Theta <- coanc_admix(Q, inbr_subpops)
    s <- mean(Theta) / mean(diag(Theta)) # this is the correct bias coeff, with uniform weights
    expect_equal(s, s_want)
    # since we set 0 < s_want < 1, nothing else to test

    # test extreme case s_want = 1
    s_want <- 1

    # test with q1d
    sigma <- bias_coeff_admix_fit(bias_coeff = s_want, inbr_subpops = inbr_subpops, n_ind = n, func = q1d) # get sigma
    expect_true( is.infinite(sigma) ) # only `sigma = Inf` should achieve the max
    # construct everything and verify s == s_want
    Q <- q1d(n, k, sigma) # now get Q from there
    Theta <- coanc_admix(Q, inbr_subpops)
    s <- mean(Theta) / mean(diag(Theta)) # this is the correct bias coeff, with uniform weights
    expect_equal(s, s_want)
    
    # test with q1dc
    sigma <- bias_coeff_admix_fit(bias_coeff = s_want, inbr_subpops = inbr_subpops, n_ind = n, func = q1dc) # get sigma
    expect_true( is.infinite(sigma) ) # only `sigma = Inf` should achieve the max
    # construct everything and verify s == s_want
    Q <- q1dc(n, k, sigma) # now get Q from there
    Theta <- coanc_admix(Q, inbr_subpops)
    s <- mean(Theta) / mean(diag(Theta)) # this is the correct bias coeff, with uniform weights
    expect_equal(s, s_want)

    # test extreme case s_want = minimum
    # test with q1d
    # construct directly
    admix_prop_bias_coeff_min <- q1d(n, k, sigma = 0)
    # this is the mminimum s_want
    s_want <- bias_coeff_admix(admix_prop_bias_coeff_min, inbr_subpops)
    sigma <- bias_coeff_admix_fit(bias_coeff = s_want, inbr_subpops = inbr_subpops, n_ind = n, func = q1d) # get sigma
    expect_equal( sigma, 0 ) # only `sigma = 0` should achieve the min
    # construct everything and verify s == s_want
    Q <- q1d(n, k, sigma) # now get Q from there
    Theta <- coanc_admix(Q, inbr_subpops)
    s <- mean(Theta) / mean(diag(Theta)) # this is the correct bias coeff, with uniform weights
    expect_equal(s, s_want)
    
    # test with q1dc
    # construct directly
    admix_prop_bias_coeff_min <- q1dc(n, k, sigma = 0)
    # this is the mminimum s_want
    s_want <- bias_coeff_admix(admix_prop_bias_coeff_min, inbr_subpops)
    sigma <- bias_coeff_admix_fit(bias_coeff = s_want, inbr_subpops = inbr_subpops, n_ind = n, func = q1dc) # get sigma
    expect_equal( sigma, 0 ) # only `sigma = 0` should achieve the min
    # construct everything and verify s == s_want
    Q <- q1dc(n, k, sigma) # now get Q from there
    Theta <- coanc_admix(Q, inbr_subpops)
    s <- mean(Theta) / mean(diag(Theta)) # this is the correct bias coeff, with uniform weights
    expect_equal(s, s_want)
})

test_that("rescaleF agrees with explicitly Fst calculation", {
    n <- 5
    k <- 2
    sigma <- 1
    Fst <- 0.1
    Q <- q1d(n, k, sigma)
    F <- 1:k # scale doesn't matter right now...
    F2 <- rescaleF(Q, F, Fst) # calculation to compare to
    Theta <- coanc_admix(Q, F2) # in wrong scale but meh
    Fst2 <- mean(diag(Theta)) # this is the actual Fst, with uniform weights
    expect_equal(Fst, Fst2)
    # since 0 < Fst=0.1 < 1, there's nothing else to test
})

test_that("rpanc is in range", {
    m <- 1000
    pAnc <- rpanc(m)
    expect_equal(length(pAnc), m)
    expect_true(all(pAnc >= 0)) # all are non-negative
    expect_true(all(pAnc <= 1)) # all are smaller or equal than 1
})

test_that("rpint is in range", {
    m <- 1000
    F <- c(0.1, 0.2, 0.3)
    k <- length(F)
    pAnc <- rpanc(m)
    B <- rpint(pAnc, F)
    expect_equal(nrow(B), m)
    expect_equal(ncol(B), k)
    expect_true(all(B >= 0)) # all are non-negative
    expect_true(all(B <= 1)) # all are smaller or equal than 1
})

test_that("rpiaf is in range", {
    m <- 1000
    F <- c(0.1, 0.2, 0.3)
    k <- length(F)
    Q <- diag(rep.int(1, k)) # island model for test...
    Q <- rbind(Q, Q, Q) # repeat so we have multiple people per island
    n <- nrow(Q) # number of individuals (3 * k)
    pAnc <- rpanc(m)
    B <- rpint(pAnc, F)
    P <- rpiaf(B, Q)
    expect_equal(nrow(P), m)
    expect_equal(ncol(P), n)
    expect_true(all(P >= 0)) # all are non-negative
    expect_true(all(P <= 1)) # all are smaller or equal than 1
})

test_that("rgeno is in range", {
    m <- 1000
    F <- c(0.1, 0.2, 0.3)
    k <- length(F)
    Q <- diag(rep.int(1, k)) # island model for test...
    Q <- rbind(Q, Q, Q) # repeat so we have multiple people per island
    n <- nrow(Q) # number of individuals (3 * k)
    pAnc <- rpanc(m)
    B <- rpint(pAnc, F)
    P <- rpiaf(B, Q)
    X <- rgeno(P) # direct test
    expect_equal(nrow(X), m)
    expect_equal(ncol(X), n)
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    X <- rgeno(B, Q) # indirect draw test
    expect_equal(nrow(X), m)
    expect_equal(ncol(X), n)
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    X <- rgeno(B, Q, lowMem = TRUE) # indirect draw test with low-mem algo
    expect_equal(nrow(X), m)
    expect_equal(ncol(X), n)
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
})

test_that("rbnpsd works well", {
    m <- 1000
    F <- c(0.1, 0.2, 0.3)
    k <- length(F)
    Q <- diag(rep.int(1, k)) # island model for test...
    Q <- rbind(Q, Q, Q) # repeat so we have multiple people per island
    n <- nrow(Q) # number of individuals (3*k)
    # run rbnpsd
    out <- rbnpsd(Q, F, m)
    X <- out$X # genotypes
    P <- out$P # IAFs
    B <- out$B # Intermediate AFs
    pAnc <- out$Pa # Ancestral AFs
    # test X
    expect_equal(nrow(X), m)
    expect_equal(ncol(X), n)
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    # test P
    expect_equal(nrow(P), m)
    expect_equal(ncol(P), n)
    expect_true(all(P >= 0)) # all are non-negative
    expect_true(all(P <= 1)) # all are smaller or equal than 1
    # test B
    expect_equal(nrow(B), m)
    expect_equal(ncol(B), k)
    expect_true(all(B >= 0)) # all are non-negative
    expect_true(all(B <= 1)) # all are smaller or equal than 1
    # test pAnc
    expect_equal(length(pAnc), m)
    expect_true(all(pAnc >= 0)) # all are non-negative
    expect_true(all(pAnc <= 1)) # all are smaller or equal than 1
})

test_that("rbnpsd `noFixed = TRUE` works well", {
    m <- 1000
    F <- c(0.1, 0.2, 0.3)
    k <- length(F)
    Q <- diag(rep.int(1, k)) # island model for test...
    Q <- rbind(Q, Q, Q) # repeat so we have multiple people per island
    n <- nrow(Q) # number of individuals (3*k)
    # run rbnpsd
    out <- rbnpsd(Q, F, m, noFixed = TRUE)
    X <- out$X # genotypes
    P <- out$P # IAFs
    B <- out$B # Intermediate AFs
    pAnc <- out$Pa # Ancestral AFs
    # test X
    expect_equal(nrow(X), m)
    expect_equal(ncol(X), n)
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    expect_true(!any(fixed_loci(X))) # we don't expect any loci to be fixed
    # test P
    expect_equal(nrow(P), m)
    expect_equal(ncol(P), n)
    expect_true(all(P >= 0)) # all are non-negative
    expect_true(all(P <= 1)) # all are smaller or equal than 1
    # test B
    expect_equal(nrow(B), m)
    expect_equal(ncol(B), k)
    expect_true(all(B >= 0)) # all are non-negative
    expect_true(all(B <= 1)) # all are smaller or equal than 1
    # test pAnc
    expect_equal(length(pAnc), m)
    expect_true(all(pAnc >= 0)) # all are non-negative
    expect_true(all(pAnc <= 1)) # all are smaller or equal than 1
})

test_that("rbnpsd `lowMem = TRUE` works well", {
    m <- 1000
    F <- c(0.1, 0.2, 0.3)
    k <- length(F)
    Q <- diag(rep.int(1, k)) # island model for test...
    Q <- rbind(Q, Q, Q) # repeat so we have multiple people per island
    n <- nrow(Q) # number of individuals (3*k)
    # run rbnpsd
    out <- rbnpsd(Q, F, m, lowMem = TRUE)
    X <- out$X # genotypes
    B <- out$B # Intermediate AFs
    pAnc <- out$Pa # Ancestral AFs
    # test X
    expect_equal(nrow(X), m)
    expect_equal(ncol(X), n)
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    # test B
    expect_equal(nrow(B), m)
    expect_equal(ncol(B), k)
    expect_true(all(B >= 0)) # all are non-negative
    expect_true(all(B <= 1)) # all are smaller or equal than 1
    # test pAnc
    expect_equal(length(pAnc), m)
    expect_true(all(pAnc >= 0)) # all are non-negative
    expect_true(all(pAnc <= 1)) # all are smaller or equal than 1
})
