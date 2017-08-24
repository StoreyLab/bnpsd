context('bnpsd')

## start with lower-level/internal tests, more informative that higher-level function errors

test_that("coanc works in toy cases", {
    Q <- diag(c(1,1)) # an IS model with two subpops
    F <- c(0.1, 0.3)
    ThetaExp <- diag(F) # the Theta we expect for this setup
    Theta <- coanc(Q, F)
    expect_equal(Theta, ThetaExp)

    ## same Q, scalar F
    F <- 0.2
    ThetaExp <- diag(c(F,F)) # the Theta we expect for this setup
    Theta <- coanc(Q, F)
    expect_equal(Theta, ThetaExp)

    ## same Q, matrix F
    Fv <- c(0.1,0.4) # vector version
    F <- diag(Fv) # matrix version
    Theta <- coanc(Q, F)
    expect_equal(Theta, F) # F is the theta we expect in this case

    ## most complex case, just a general math check
    Q <- matrix(c(0.7, 0.3, 0.2, 0.8), nrow=2, byrow=TRUE)
    F <- matrix(c(0.3, 0.1, 0.1, 0.3), nrow=2, byrow=TRUE)
    ThetaExp <- Q %*% F %*% t(Q) # the Theta we expect for this setup (slower but more explicit version)
    Theta <- coanc(Q, F)
    expect_equal(Theta, ThetaExp)
}) 

test_that("fst works in toy cases", {
    Q <- diag(c(1,1)) # an IS model with two subpops
    F <- c(0.1, 0.3)
    fst1 <- mean(F) # the Theta we expect for this setup
    fst2 <- fst(Q, F)
    expect_equal(fst1, fst2)

    ## same Q, scalar F
    F <- 0.2
    fst2 <- fst(Q, F)
    expect_equal(F, fst2)

    ## same Q, matrix F
    Fv <- c(0.1,0.4) # vector version
    F <- diag(Fv) # matrix version
    fst1 <- mean(Fv)
    fst2 <- fst(Q, F)
    expect_equal(fst1, fst2) # F is the theta we expect in this case

    ## most complex case, just a general math check
    Q <- matrix(c(0.7, 0.3, 0.2, 0.8), nrow=2, byrow=TRUE)
    F <- matrix(c(0.3, 0.1, 0.1, 0.3), nrow=2, byrow=TRUE)
    fst1 <- mean(diag(Q %*% F %*% t(Q))) # the Fst we expect for this setup (slower but more explicit version)
    fst2 <- fst(Q, F)
    expect_equal(fst1, fst2)
}) 

test_that("qis returns valid admixture coefficients", {
    labs <- c(1,2,2,3,3,3,4,4,4,4)
    n <- length(labs)
    k <- length(unique(labs))
    Q <- qis(labs)
    ## general tests for admixture matrices
    expect_equal(nrow(Q), n) # n rows
    expect_equal(ncol(Q), k) # k columns
    expect_true(all(Q >= 0)) # all are non-negative
    expect_true(all(Q <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(Q), rep.int(1,n)) # rows sum to 1, vector length n
    ## specific tests for qis
    expect_true(all(Q %in% c(TRUE,FALSE)))
    expect_true(all(colnames(Q) == sort(unique(labs))))
    
    ## test with provided subpops
    subpops <- 4:1
    Q <- qis(labs, subpops)
    ## general tests for admixture matrices
    expect_equal(nrow(Q), n) # n rows
    expect_equal(ncol(Q), k) # k columns
    expect_true(all(Q >= 0)) # all are non-negative
    expect_true(all(Q <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(Q), rep.int(1,n)) # rows sum to 1, vector length n
    ## specific tests for qis
    expect_true(all(Q %in% c(TRUE,FALSE)))
    expect_true(all(colnames(Q) == subpops))
    
    ## test with provided subpops (additional labels)
    k <- 10
    subpops <- 1:k
    Q <- qis(labs, subpops)
    ## general tests for admixture matrices
    expect_equal(nrow(Q), n) # n rows
    expect_equal(ncol(Q), k) # k columns
    expect_true(all(Q >= 0)) # all are non-negative
    expect_true(all(Q <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(Q), rep.int(1,n)) # rows sum to 1, vector length n
    ## specific tests for qis
    expect_true(all(Q %in% c(TRUE,FALSE)))
    expect_true(all(colnames(Q) == subpops))

    ## test with provided subpops (missing labels, must die!)
    subpops <- 1:3 # missing 4!
    expect_error( qis(labs, subpops) )
})

test_that("q1d returns valid admixture coefficients", {
    n <- 5
    k <- 2
    sigma <- 1
    Q <- q1d(n, k, sigma)
    expect_equal(nrow(Q), n) # n rows
    expect_equal(ncol(Q), k) # k columns
    expect_true(all(Q >= 0)) # all are non-negative
    expect_true(all(Q <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(Q), rep.int(1,n)) # rows sum to 1, vector length n
})

test_that("q1d_sigma2s agrees with explicitly calculated bias coeff s", {
    n <- 5
    k <- 2
    sigma <- 1
    Q <- q1d(n, k, sigma)
    F <- 1:k # scale doesn't matter right now...
    Theta <- coanc(Q, F) # in wrong scale but meh
    sWant <- mean(Theta)/mean(diag(Theta)) # this is the correct bias coeff, with uniform weights
    s <- q1d_sigma2s(sigma, F, n) # calculation to compare to
    expect_equal(s, sWant)
    expect_true(s > 0) # other obvious properties...
    expect_true(s <= 1)
})

test_that("q1d_s2sigma agrees with reverse func", {
    n <- 1000
    F <- c(0.1,0.2,0.3)
    k <- length(F)
    sWant <- 0.5
    sigma <- q1d_s2sigma(s=sWant, F=F, n=n) # get sigma
    ## construct everything and verify s == sWant
    Q <- q1d(n, k, sigma) # now get Q from there
    Theta <- coanc(Q, F)
    s <- mean(Theta)/mean(diag(Theta)) # this is the correct bias coeff, with uniform weights
    expect_equal(s, sWant)
    ## since we set 0<sWant<1, nothing else to test
})

test_that("rescaleF agrees with explicitly Fst calculation", {
    n <- 5
    k <- 2
    sigma <- 1
    Fst <- 0.1
    Q <- q1d(n, k, sigma)
    F <- 1:k # scale doesn't matter right now...
    F2 <- rescaleF(Q, F, Fst) # calculation to compare to
    Theta <- coanc(Q, F2) # in wrong scale but meh
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
    F <- c(0.1,0.2,0.3)
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
    F <- c(0.1,0.2,0.3)
    k <- length(F)
    Q <- diag(rep.int(1,k)) # island model for test...
    Q <- rbind(Q,Q,Q) # repeat so we have multiple people per island
    n <- nrow(Q) # number of individuals (3*k)
    pAnc <- rpanc(m)
    B <- rpint(pAnc, F)
    P <- rpiaf(B,Q)
    expect_equal(nrow(P), m)
    expect_equal(ncol(P), n)
    expect_true(all(P >= 0)) # all are non-negative
    expect_true(all(P <= 1)) # all are smaller or equal than 1
})

test_that("rgeno is in range", {
    m <- 1000
    F <- c(0.1,0.2,0.3)
    k <- length(F)
    Q <- diag(rep.int(1,k)) # island model for test...
    Q <- rbind(Q,Q,Q) # repeat so we have multiple people per island
    n <- nrow(Q) # number of individuals (3*k)
    pAnc <- rpanc(m)
    B <- rpint(pAnc, F)
    P <- rpiaf(B,Q)
    X <- rgeno(P) # direct test
    expect_equal(nrow(X), m)
    expect_equal(ncol(X), n)
    expect_true(all(X %in% c(0,1,2))) # only three values allowed!
    X <- rgeno(B,Q) # indirect draw test
    expect_equal(nrow(X), m)
    expect_equal(ncol(X), n)
    expect_true(all(X %in% c(0,1,2))) # only three values allowed!
    X <- rgeno(B,Q,lowMem=TRUE) # indirect draw test with low-mem algo
    expect_equal(nrow(X), m)
    expect_equal(ncol(X), n)
    expect_true(all(X %in% c(0,1,2))) # only three values allowed!
})

test_that("rbnpsd works well", {
    m <- 1000
    F <- c(0.1,0.2,0.3)
    k <- length(F)
    Q <- diag(rep.int(1,k)) # island model for test...
    Q <- rbind(Q,Q,Q) # repeat so we have multiple people per island
    n <- nrow(Q) # number of individuals (3*k)
    ## run rbnpsd
    out <- rbnpsd(Q, F, m)
    X <- out$X # genotypes
    P <- out$P # IAFs
    B <- out$B # Intermediate AFs
    pAnc <- out$Pa # Ancestral AFs
    ## test X
    expect_equal(nrow(X), m)
    expect_equal(ncol(X), n)
    expect_true(all(X %in% c(0,1,2))) # only three values allowed!
    ## test P
    expect_equal(nrow(P), m)
    expect_equal(ncol(P), n)
    expect_true(all(P >= 0)) # all are non-negative
    expect_true(all(P <= 1)) # all are smaller or equal than 1
    ## test B
    expect_equal(nrow(B), m)
    expect_equal(ncol(B), k)
    expect_true(all(B >= 0)) # all are non-negative
    expect_true(all(B <= 1)) # all are smaller or equal than 1
    ## test pAnc
    expect_equal(length(pAnc), m)
    expect_true(all(pAnc >= 0)) # all are non-negative
    expect_true(all(pAnc <= 1)) # all are smaller or equal than 1
})

test_that("rbnpsd lowMem=TRUE works well", {
    m <- 1000
    F <- c(0.1,0.2,0.3)
    k <- length(F)
    Q <- diag(rep.int(1,k)) # island model for test...
    Q <- rbind(Q,Q,Q) # repeat so we have multiple people per island
    n <- nrow(Q) # number of individuals (3*k)
    ## run rbnpsd
    out <- rbnpsd(Q, F, m, lowMem=TRUE)
    X <- out$X # genotypes
    B <- out$B # Intermediate AFs
    pAnc <- out$Pa # Ancestral AFs
    ## test X
    expect_equal(nrow(X), m)
    expect_equal(ncol(X), n)
    expect_true(all(X %in% c(0,1,2))) # only three values allowed!
    ## test B
    expect_equal(nrow(B), m)
    expect_equal(ncol(B), k)
    expect_true(all(B >= 0)) # all are non-negative
    expect_true(all(B <= 1)) # all are smaller or equal than 1
    ## test pAnc
    expect_equal(length(pAnc), m)
    expect_true(all(pAnc >= 0)) # all are non-negative
    expect_true(all(pAnc <= 1)) # all are smaller or equal than 1
})
