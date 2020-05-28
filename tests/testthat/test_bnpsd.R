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
    admix_proportions <- diag(c(1, 1)) # an independent subpops model with two subpops
    F <- c(0.1, 0.3)
    
    # die when the important parameters are missing
    expect_error( coanc_admix() ) # all missing 
    expect_error( coanc_admix(admix_proportions) ) # coanc_subpops missing
    expect_error( coanc_admix(coanc_subpops = F) ) # admix_proportions missing
    
    coancestry_expected <- diag(F) # the coancestry we expect for this setup
    coancestry <- coanc_admix(admix_proportions, F)
    expect_equal(coancestry, coancestry_expected)

    # same admix_proportions, scalar F
    F <- 0.2
    coancestry_expected <- diag(c(F, F)) # the coancestry we expect for this setup
    coancestry <- coanc_admix(admix_proportions, F)
    expect_equal(coancestry, coancestry_expected)

    # same admix_proportions, matrix F
    Fv <- c(0.1, 0.4) # vector version
    F <- diag(Fv) # matrix version
    coancestry <- coanc_admix(admix_proportions, F)
    expect_equal(coancestry, F) # F is the theta we expect in this case

    # most complex case, just a general math check
    admix_proportions <- matrix(c(0.7, 0.3, 0.2, 0.8), nrow = 2, byrow = TRUE)
    F <- matrix(c(0.3, 0.1, 0.1, 0.3), nrow = 2, byrow = TRUE)
    coancestry_expected <- admix_proportions %*% F %*% t(admix_proportions) # the coancestry we expect for this setup (slower but more explicit version)
    coancestry <- coanc_admix(admix_proportions, F)
    expect_equal(coancestry, coancestry_expected)
}) 

test_that("coanc_to_kinship works in toy cases", {
    admix_proportions <- diag(c(1, 1)) # an IS model with two subpops
    F <- c(0.1, 0.3)
    PhiExp <- diag( (1+F) / 2 ) # the Phi we expect for this setup
    coancestry <- coanc_admix(admix_proportions, F)
    Phi <- coanc_to_kinship(coancestry)
    expect_equal(Phi, PhiExp)

    # same admix_proportions, scalar F
    F <- 0.2
    K <- (1+F) / 2 # transform at this stage
    PhiExp <- diag(c(K, K)) # the Phi we expect for this setup
    coancestry <- coanc_admix(admix_proportions, F)
    Phi <- coanc_to_kinship(coancestry)
    expect_equal(Phi, PhiExp)

    # same admix_proportions, matrix F
    Fv <- c(0.1, 0.4) # vector version
    F <- diag(Fv) # matrix version
    PhiExp <- diag((1+Fv) / 2)
    coancestry <- coanc_admix(admix_proportions, F)
    Phi <- coanc_to_kinship(coancestry)
    expect_equal(Phi, PhiExp)

    # most complex case, just a general math check
    admix_proportions <- matrix(c(0.7, 0.3, 0.2, 0.8), nrow = 2, byrow = TRUE)
    F <- matrix(c(0.3, 0.1, 0.1, 0.3), nrow = 2, byrow = TRUE)
    PhiExp <- admix_proportions %*% F %*% t(admix_proportions) # the coancestry we expect for this setup (slower but more explicit version)
    diag(PhiExp) <- (1 + diag(PhiExp))/2 # explicit transformation to kinship
    coancestry <- coanc_admix(admix_proportions, F)
    Phi <- coanc_to_kinship(coancestry)
    expect_equal(Phi, PhiExp)
}) 

test_that("fst_admix works in toy cases", {
    admix_proportions <- diag(c(1, 1)) # an IS model with two subpops
    F <- c(0.1, 0.3)
    fst1 <- mean(F) # the coancestry we expect for this setup
    fst2 <- fst_admix(admix_proportions, F)
    expect_equal(fst1, fst2)

    # same admix_proportions, scalar F
    F <- 0.2
    fst2 <- fst_admix(admix_proportions, F)
    expect_equal(F, fst2)

    # same admix_proportions, matrix F
    Fv <- c(0.1, 0.4) # vector version
    F <- diag(Fv) # matrix version
    fst1 <- mean(Fv)
    fst2 <- fst_admix(admix_proportions, F)
    expect_equal(fst1, fst2) # F is the theta we expect in this case

    # most complex case, just a general math check
    admix_proportions <- matrix(c(0.7, 0.3, 0.2, 0.8), nrow = 2, byrow = TRUE)
    F <- matrix(c(0.3, 0.1, 0.1, 0.3), nrow = 2, byrow = TRUE)
    fst1 <- mean(diag(admix_proportions %*% F %*% t(admix_proportions))) # the FST we expect for this setup (slower but more explicit version)
    fst2 <- fst_admix(admix_proportions, F)
    expect_equal(fst1, fst2)
})

test_that("bias_coeff_admix agrees with explicitly calculated bias_coeff", {
    # set up some simulated data
    n <- 5
    k_subpops <- 2
    sigma <- 1
    admix_proportions <- admix_prop_1d_linear(n, k_subpops, sigma)
    F <- 1:k_subpops # scale doesn't matter right now...

    coancestry <- coanc_admix(admix_proportions, F) # in wrong scale but meh
    sWant <- mean(coancestry) / mean(diag(coancestry)) # this is the correct bias coeff, with uniform weights
    s <- bias_coeff_admix(admix_proportions, F) # calculation to compare to
    expect_equal(s, sWant)
    expect_true(s > 0) # other obvious properties...
    expect_true(s <= 1)

    # repeat with matrix F, missing (uniform) weights
    F_mat <- diag(F)
    s <- bias_coeff_admix(admix_proportions, F_mat) # calculation to compare to
    expect_equal(s, sWant)
    expect_true(s > 0) # other obvious properties...
    expect_true(s <= 1)
    
    # repeat with non-uniform weights...
    weights <- runif(n) # random weights for given number of individuals
    weights <- weights / sum(weights) # normalize to add up to 1! # NOTE: should check sum(weights)!= 0, meh...
    sWant <- drop(weights %*% coancestry %*% weights) / drop( diag(coancestry) %*% weights ) # this is the correct bias coeff, with uniform weights
    s <- bias_coeff_admix(admix_proportions, F, weights) # calculation to compare to
    expect_equal(s, sWant)
    expect_true(s > 0) # other obvious properties...
    expect_true(s <= 1)

    # repeat with matrix F and non-uniform weights
    s <- bias_coeff_admix(admix_proportions, F_mat, weights) # calculation to compare to
    expect_equal(s, sWant)
    expect_true(s > 0) # other obvious properties...
    expect_true(s <= 1)
})

test_that("admix_prop_indep_subpops returns valid admixture coefficients", {
    labs <- c(1, 2, 2, 3, 3, 3, 4, 4, 4, 4)
    n <- length(labs)
    k_subpops <- length(unique(labs))
    admix_proportions <- admix_prop_indep_subpops(labs)
    # general tests for admixture matrices
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k_subpops) # k_subpops columns
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n
    # specific tests for admix_prop_indep_subpops
    expect_true(all(admix_proportions %in% c(TRUE, FALSE)))
    expect_true(all(colnames(admix_proportions) == sort(unique(labs))))
    
    # test with provided subpops
    subpops <- 4:1
    admix_proportions <- admix_prop_indep_subpops(labs, subpops)
    # general tests for admixture matrices
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k_subpops) # k_subpops columns
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1,n)) # rows sum to 1, vector length n
    # specific tests for admix_prop_indep_subpops
    expect_true(all(admix_proportions %in% c(TRUE, FALSE)))
    expect_true(all(colnames(admix_proportions) == subpops))
    
    # test with provided subpops (additional labels)
    k_subpops <- 10
    subpops <- 1:k_subpops
    admix_proportions <- admix_prop_indep_subpops(labs, subpops)
    # general tests for admixture matrices
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k_subpops) # k_subpops columns
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n
    # specific tests for admix_prop_indep_subpops
    expect_true(all(admix_proportions %in% c(TRUE, FALSE)))
    expect_true(all(colnames(admix_proportions) == subpops))

    # test with provided subpops (missing labels, must die!)
    subpops <- 1:3 # missing 4!
    expect_error( admix_prop_indep_subpops(labs, subpops) )
})

test_that("admix_prop_1d_linear returns valid admixture coefficients", {
    n <- 10
    k_subpops <- 2
    admix_proportions <- admix_prop_1d_linear(n, k_subpops, sigma = 1)
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k_subpops) # k_subpops columns
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n

    # test with sigma == 0 (special case that makes usual formula break)
    admix_proportions <- admix_prop_1d_linear(n, k_subpops, sigma = 0)
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k_subpops) # k_subpops columns
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n
    # in this case it should equal independent subpopulations
    labs <- c( rep.int(1, 5), rep.int(2, 5) ) # two subpops
    admix_proportions2 <- admix_prop_indep_subpops(labs)
    dimnames(admix_proportions2) <- NULL # before comparing, must toss column names
    expect_equal(admix_proportions, admix_proportions2)

    # test bias_coeff version
    obj <- admix_prop_1d_linear(n, k_subpops, bias_coeff = 0.5, coanc_subpops = 1:k_subpops, fst = 0.1)
    admix_proportions <- obj$admix_proportions # returns many things in this case, get admix_proportions here
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k_subpops) # k_subpops columns
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n
})

test_that("admix_prop_1d_circular returns valid admixture coefficients", {
    n <- 10
    k_subpops <- 2
    admix_proportions <- admix_prop_1d_circular(n, k_subpops, sigma = 1)
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k_subpops) # k_subpops columns
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n

    # test with sigma == 0 (special case that makes usual formula break)
    admix_proportions <- admix_prop_1d_circular(n, k_subpops, sigma = 0)
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k_subpops) # k_subpops columns
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n
    # though the result is nearly island-like, there is an annoying shift I'd rather not try to figure out for this test...

    # test bias_coeff version
    obj <- admix_prop_1d_circular(n, k_subpops, bias_coeff = 0.5, coanc_subpops = 1:k_subpops, fst = 0.1)
    admix_proportions <- obj$admix_proportions # returns many things in this case, get admix_proportions here
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k_subpops) # k_subpops columns
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n
})

test_that("bias_coeff_admix_fit agrees with reverse func", {
    n <- 1000
    inbr_subpops <- c(0.1, 0.2, 0.3)
    k_subpops <- length(inbr_subpops)
    s_want <- 0.5
    coord_ind_first_linear <- 0.5
    coord_ind_last_linear <- k_subpops + 0.5
    coord_ind_first_circular <- 2 * pi / (2 * n)
    coord_ind_last_circular <- 2 * pi * (1 - 1 / (2 * n) )

    # test with admix_prop_1d_linear
    # get sigma
    sigma <- bias_coeff_admix_fit(
        bias_coeff = s_want,
        coanc_subpops = inbr_subpops,
        n_ind = n,
        k_subpops = k_subpops,
        func = admix_prop_1d_linear,
        coord_ind_first = coord_ind_first_linear,
        coord_ind_last = coord_ind_last_linear
    )
    # construct everything and verify s == s_want
    admix_proportions <- admix_prop_1d_linear(
        n_ind = n,
        k_subpops = k_subpops,
        sigma = sigma,
        coord_ind_first = coord_ind_first_linear,
        coord_ind_last = coord_ind_last_linear
    )
    coancestry <- coanc_admix(admix_proportions, inbr_subpops)
    s <- mean(coancestry) / mean(diag(coancestry)) # this is the correct bias coeff, with uniform weights
    expect_equal(s, s_want)
    # since we set 0 < s_want < 1, nothing else to test
    
    # test with admix_prop_1d_circular
    # get sigma
    sigma <- bias_coeff_admix_fit(
        bias_coeff = s_want,
        coanc_subpops = inbr_subpops,
        n_ind = n,
        k_subpops = k_subpops,
        func = admix_prop_1d_circular,
        coord_ind_first = coord_ind_first_circular,
        coord_ind_last = coord_ind_last_circular
    ) # get sigma
    # construct everything and verify s == s_want
    admix_proportions <- admix_prop_1d_circular(
        n_ind = n,
        k_subpops = k_subpops,
        sigma = sigma,
        coord_ind_first = coord_ind_first_circular,
        coord_ind_last = coord_ind_last_circular
    ) # now get admix_proportions from there
    coancestry <- coanc_admix(admix_proportions, inbr_subpops)
    s <- mean(coancestry) / mean(diag(coancestry)) # this is the correct bias coeff, with uniform weights
    expect_equal(s, s_want)
    # since we set 0 < s_want < 1, nothing else to test

    # test extreme case s_want = 1
    s_want <- 1

    # test with admix_prop_1d_linear
    # get sigma
    sigma <- bias_coeff_admix_fit(
        bias_coeff = s_want,
        coanc_subpops = inbr_subpops,
        n_ind = n,
        k_subpops = k_subpops,
        func = admix_prop_1d_linear,
        coord_ind_first = coord_ind_first_linear,
        coord_ind_last = coord_ind_last_linear
    )
    expect_true( is.infinite(sigma) ) # only `sigma = Inf` should achieve the max
    # construct everything and verify s == s_want
    admix_proportions <- admix_prop_1d_linear(
        n_ind = n,
        k_subpops = k_subpops,
        sigma = sigma,
        coord_ind_first = coord_ind_first_linear,
        coord_ind_last = coord_ind_last_linear
    ) # now get admix_proportions from there
    coancestry <- coanc_admix(admix_proportions, inbr_subpops)
    s <- mean(coancestry) / mean(diag(coancestry)) # this is the correct bias coeff, with uniform weights
    expect_equal(s, s_want)
    
    # test with admix_prop_1d_circular
    # get sigma
    sigma <- bias_coeff_admix_fit(
        bias_coeff = s_want,
        coanc_subpops = inbr_subpops,
        n_ind = n,
        k_subpops = k_subpops,
        func = admix_prop_1d_circular,
        coord_ind_first = coord_ind_first_circular,
        coord_ind_last = coord_ind_last_circular
    )
    expect_true( is.infinite(sigma) ) # only `sigma = Inf` should achieve the max
    # construct everything and verify s == s_want
    admix_proportions <- admix_prop_1d_circular(
        n_ind = n,
        k_subpops = k_subpops,
        sigma = sigma,
        coord_ind_first = coord_ind_first_circular,
        coord_ind_last = coord_ind_last_circular
    ) # now get admix_proportions from there
    coancestry <- coanc_admix(admix_proportions, inbr_subpops)
    s <- mean(coancestry) / mean(diag(coancestry)) # this is the correct bias coeff, with uniform weights
    expect_equal(s, s_want)

    # test extreme case s_want = minimum
    # test with admix_prop_1d_linear
    # construct directly
    admix_prop_bias_coeff_min <- admix_prop_1d_linear(
        n,
        k_subpops,
        sigma = 0,
        coord_ind_first = coord_ind_first_linear,
        coord_ind_last = coord_ind_last_linear
    )
    # this is the mminimum s_want
    s_want <- bias_coeff_admix(admix_prop_bias_coeff_min, inbr_subpops)
    # get sigma
    sigma <- bias_coeff_admix_fit(
        bias_coeff = s_want,
        coanc_subpops = inbr_subpops,
        n_ind = n,
        k_subpops = k_subpops,
        func = admix_prop_1d_linear,
        coord_ind_first = coord_ind_first_linear,
        coord_ind_last = coord_ind_last_linear
    )
    expect_equal( sigma, 0 ) # only `sigma = 0` should achieve the min
    # construct everything and verify s == s_want
    admix_proportions <- admix_prop_1d_linear(
        n_ind = n,
        k_subpops = k_subpops,
        sigma = sigma,
        coord_ind_first = coord_ind_first_linear,
        coord_ind_last = coord_ind_last_linear
    ) # now get admix_proportions from there
    coancestry <- coanc_admix(admix_proportions, inbr_subpops)
    s <- mean(coancestry) / mean(diag(coancestry)) # this is the correct bias coeff, with uniform weights
    expect_equal(s, s_want)
    
    # test with admix_prop_1d_circular
    # construct directly
    admix_prop_bias_coeff_min <- admix_prop_1d_circular(
        n,
        k_subpops,
        sigma = 0,
        coord_ind_first = coord_ind_first_circular,
        coord_ind_last = coord_ind_last_circular
    )
    # this is the mminimum s_want
    s_want <- bias_coeff_admix(admix_prop_bias_coeff_min, inbr_subpops)
    # get sigma
    sigma <- bias_coeff_admix_fit(
        bias_coeff = s_want,
        coanc_subpops = inbr_subpops,
        n_ind = n,
        k_subpops = k_subpops,
        func = admix_prop_1d_circular,
        coord_ind_first = coord_ind_first_circular,
        coord_ind_last = coord_ind_last_circular
    )
    expect_equal( sigma, 0 ) # only `sigma = 0` should achieve the min
    # construct everything and verify s == s_want
    admix_proportions <- admix_prop_1d_circular(
        n_ind = n,
        k_subpops = k_subpops,
        sigma = sigma,
        coord_ind_first = coord_ind_first_circular,
        coord_ind_last = coord_ind_last_circular
    ) # now get admix_proportions from there
    coancestry <- coanc_admix(admix_proportions, inbr_subpops)
    s <- mean(coancestry) / mean(diag(coancestry)) # this is the correct bias coeff, with uniform weights
    expect_equal(s, s_want)

    # test s_want with matrix coanc_subpops (all previous tests were vectors)
    # set up population structure
    n <- 1000
    k_subpops <- 3
    # non-trivial coancestry matrix for intermediate subpopulations!
    coanc_subpops <- matrix(
        c(
            0.1 , 0.05, 0   ,
            0.05, 0.2 , 0.05,
            0   , 0.05, 0.3
        ),
        nrow = k_subpops
    )
    s_want <- 0.5

    # test with admix_prop_1d_linear
    # get sigma
    sigma <- bias_coeff_admix_fit(
        bias_coeff = s_want,
        coanc_subpops = coanc_subpops,
        n_ind = n,
        k_subpops = k_subpops,
        func = admix_prop_1d_linear,
        coord_ind_first = coord_ind_first_linear,
        coord_ind_last = coord_ind_last_linear
    )
    # construct everything and verify s == s_want
    admix_proportions <- admix_prop_1d_linear(
        n_ind = n,
        k_subpops = k_subpops,
        sigma = sigma,
        coord_ind_first = coord_ind_first_linear,
        coord_ind_last = coord_ind_last_linear
    ) # now get admix_proportions from there
    coancestry <- coanc_admix(admix_proportions, coanc_subpops)
    s <- mean(coancestry) / mean(diag(coancestry)) # this is the correct bias coeff, with uniform weights
    expect_equal(s, s_want)
    # since we set 0 < s_want < 1, nothing else to test
    
    # test with admix_prop_1d_circular
    # get sigma
    sigma <- bias_coeff_admix_fit(
        bias_coeff = s_want,
        coanc_subpops = coanc_subpops,
        n_ind = n,
        k_subpops = k_subpops,
        func = admix_prop_1d_circular,
        coord_ind_first = coord_ind_first_circular,
        coord_ind_last = coord_ind_last_circular
    ) # get sigma
    # construct everything and verify s == s_want
    admix_proportions <- admix_prop_1d_circular(
        n_ind = n,
        k_subpops = k_subpops,
        sigma = sigma,
        coord_ind_first = coord_ind_first_circular,
        coord_ind_last = coord_ind_last_circular
    ) # now get admix_proportions from there
    coancestry <- coanc_admix(admix_proportions, coanc_subpops)
    s <- mean(coancestry) / mean(diag(coancestry)) # this is the correct bias coeff, with uniform weights
    expect_equal(s, s_want)
    # since we set 0 < s_want < 1, nothing else to test
})

test_that("rescale_coanc_subpops agrees with explicitly FST calculation", {
    n <- 5
    k_subpops <- 2
    sigma <- 1
    Fst <- 0.1
    admix_proportions <- admix_prop_1d_linear(n, k_subpops, sigma)
    F <- 1:k_subpops # scale doesn't matter right now...
    F2 <- rescale_coanc_subpops(admix_proportions, F, Fst) # calculation to compare to
    coancestry <- coanc_admix(admix_proportions, F2) # in wrong scale but meh
    Fst2 <- mean(diag(coancestry)) # this is the actual FST, with uniform weights
    expect_equal(Fst, Fst2)
    # since 0 < Fst=0.1 < 1, there's nothing else to test
})

test_that("draw_p_anc is in range", {
    m_loci <- 1000
    p_anc <- draw_p_anc(m_loci)
    expect_equal(length(p_anc), m_loci)
    expect_true(all(p_anc >= 0)) # all are non-negative
    expect_true(all(p_anc <= 1)) # all are smaller or equal than 1

    # repeat for Beta version
    beta <- 0.01
    p_anc <- draw_p_anc(m_loci, beta = beta)
    expect_equal(length(p_anc), m_loci)
    expect_true(all(p_anc >= 0)) # all are non-negative
    expect_true(all(p_anc <= 1)) # all are smaller or equal than 1
})

test_that("draw_p_subpops is in range", {
    
    # a typical example with vector p_anc and inbr_subpops
    m_loci <- 10
    inbr_subpops <- c(0.1, 0.2, 0.3)
    k_subpops <- length(inbr_subpops)
    p_anc <- draw_p_anc(m_loci)

    # make sure things die when parameters are missing
    expect_error( draw_p_subpops() )
    expect_error( draw_p_subpops(p_anc = p_anc) )
    expect_error( draw_p_subpops(inbr_subpops = inbr_subpops) )
    
    # test main use case
    p_subpops <- draw_p_subpops(p_anc, inbr_subpops)
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1

    # make sure this doesn't die or anything when extra parameters are passed but agree
    expect_silent( p_subpops <- draw_p_subpops(p_anc, inbr_subpops, m_loci = m_loci) )
    expect_silent( p_subpops <- draw_p_subpops(p_anc, inbr_subpops, k_subpops = k_subpops) )
    expect_silent( p_subpops <- draw_p_subpops(p_anc, inbr_subpops, m_loci = m_loci, k_subpops = k_subpops) )

    # special case of scalar p_anc
    p_subpops <- draw_p_subpops(p_anc = 0.5, inbr_subpops, m_loci = m_loci)
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    
    # special case of scalar inbr_subpops
    k_subpops <- 2
    p_subpops <- draw_p_subpops(p_anc, inbr_subpops = 0.2, k_subpops = k_subpops)
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    
    # both main parameters scalars but return value still matrix
    p_subpops <- draw_p_subpops(p_anc = 0.5, inbr_subpops = 0.2, m_loci = m_loci, k_subpops = k_subpops)
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    
    # passing scalar parameters without setting dimensions separately results in a 1x1 matrix
    p_subpops <- draw_p_subpops(p_anc = 0.5, inbr_subpops = 0.2)
    expect_equal(nrow(p_subpops), 1)
    expect_equal(ncol(p_subpops), 1)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
})

test_that("make_p_ind_admix is in range", {
    m_loci <- 10
    inbr_subpops <- c(0.1, 0.2, 0.3)
    k_subpops <- length(inbr_subpops)
    admix_proportions <- diag(rep.int(1, k_subpops)) # island model for test...
    # repeat so we have multiple people per island
    admix_proportions <- rbind(admix_proportions, admix_proportions, admix_proportions)
    n <- nrow(admix_proportions) # number of individuals (3 * k_subpops)
    p_anc <- draw_p_anc(m_loci)
    p_subpops <- draw_p_subpops(p_anc, inbr_subpops)
    p_ind <- make_p_ind_admix(p_subpops, admix_proportions)
    expect_equal(nrow(p_ind), m_loci)
    expect_equal(ncol(p_ind), n)
    expect_true(all(p_ind >= 0)) # all are non-negative
    expect_true(all(p_ind <= 1)) # all are smaller or equal than 1
})

test_that("draw_genotypes_admix is in range", {
    m_loci <- 10
    inbr_subpops <- c(0.1, 0.2, 0.3)
    k_subpops <- length(inbr_subpops)
    admix_proportions <- diag(rep.int(1, k_subpops)) # island model for test...
    # repeat so we have multiple people per island
    admix_proportions <- rbind(admix_proportions, admix_proportions, admix_proportions)
    n_ind <- nrow(admix_proportions) # number of individuals (3 * k_subpops)
    p_anc <- draw_p_anc(m_loci)
    p_subpops <- draw_p_subpops(p_anc, inbr_subpops)
    p_ind <- make_p_ind_admix(p_subpops, admix_proportions)
    
    # direct test
    X <- draw_genotypes_admix(p_ind)
    expect_equal(nrow(X), m_loci)
    expect_equal(ncol(X), n_ind)
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    
    # indirect draw test (default low_mem now!)
    X <- draw_genotypes_admix(p_subpops, admix_proportions)
    expect_equal(nrow(X), m_loci)
    expect_equal(ncol(X), n_ind)
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    
    # the default test has 10 loci and 9 individuals
    # create a test with more individuals to test (internal default) low_mem code when dimension ratio flips
    # easy way is to reduce m_loci only, redraw only parts as needed
    m_loci <- 5
    p_anc <- draw_p_anc(m_loci)
    p_subpops <- draw_p_subpops(p_anc, inbr_subpops)
    # indirect draw test (alternative low-mem algorithm triggered by dimension changes)
    X <- draw_genotypes_admix(p_subpops, admix_proportions)
    expect_equal(nrow(X), m_loci)
    expect_equal(ncol(X), n_ind)
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
})

# for recurring tests
draw_all_admix_names_ret_default <- c('X', 'p_anc')
draw_all_admix_names_ret_full <- c('X', 'p_anc', 'p_subpops', 'p_ind')
draw_all_admix_names_ret_low_mem <- c('X', 'p_anc', 'p_subpops') # excludes p_ind

test_that("draw_all_admix works well", {
    m_loci <- 10
    inbr_subpops <- c(0.1, 0.2, 0.3)
    k_subpops <- length(inbr_subpops)
    admix_proportions <- diag(rep.int(1, k_subpops)) # island model for test...
    # repeat so we have multiple people per island
    admix_proportions <- rbind(admix_proportions, admix_proportions, admix_proportions)
    n_ind <- nrow(admix_proportions) # number of individuals (3*k_subpops)

    # make sure things die when important parameters are missing
    # all missing
    expect_error( draw_all_admix() )
    # two missing
    expect_error( draw_all_admix(admix_proportions = admix_proportions) )
    expect_error( draw_all_admix(inbr_subpops = inbr_subpops) )
    expect_error( draw_all_admix(m_loci = m_loci) )
    # one missing
    expect_error( draw_all_admix(inbr_subpops = inbr_subpops, m_loci = m_loci) )
    expect_error( draw_all_admix(admix_proportions = admix_proportions, m_loci = m_loci) )
    expect_error( draw_all_admix(admix_proportions = admix_proportions, inbr_subpops = inbr_subpops) )
    
    # run draw_all_admix
    # first test default (p_ind and p_subpops not returned)
    out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci)
    expect_equal( names(out), draw_all_admix_names_ret_default )

    # now rerun with all outpts, so we can test them all
    out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci, want_p_ind = TRUE, want_p_subpops = TRUE)
    expect_equal( names(out), draw_all_admix_names_ret_full )
    X <- out$X # genotypes
    p_ind <- out$p_ind # IAFs
    p_subpops <- out$p_subpops # Intermediate AFs
    p_anc <- out$p_anc # Ancestral AFs
    
    # test X
    expect_equal(nrow(X), m_loci)
    expect_equal(ncol(X), n_ind)
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    expect_true(!any(fixed_loci(X))) # we don't expect any loci to be fixed
    
    # test p_ind
    expect_equal(nrow(p_ind), m_loci)
    expect_equal(ncol(p_ind), n_ind)
    expect_true(all(p_ind >= 0)) # all are non-negative
    expect_true(all(p_ind <= 1)) # all are smaller or equal than 1
    
    # test p_subpops
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    
    # test p_anc
    expect_equal(length(p_anc), m_loci)
    expect_true(all(p_anc >= 0)) # all are non-negative
    expect_true(all(p_anc <= 1)) # all are smaller or equal than 1
})

test_that("draw_all_admix beta works well", {
    m_loci <- 10
    beta <- 0.01
    inbr_subpops <- c(0.1, 0.2, 0.3)
    k_subpops <- length(inbr_subpops)
    admix_proportions <- diag(rep.int(1, k_subpops)) # island model for test...
    # repeat so we have multiple people per island
    admix_proportions <- rbind(admix_proportions, admix_proportions, admix_proportions)
    n_ind <- nrow(admix_proportions) # number of individuals (3*k_subpops)

    # run draw_all_admix
    # only test default (p_ind and p_subpops not returned)
    out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci, beta = beta)
    expect_equal( names(out), draw_all_admix_names_ret_default )

    X <- out$X # genotypes
    p_anc <- out$p_anc # Ancestral AFs
    
    # test X
    expect_equal(nrow(X), m_loci)
    expect_equal(ncol(X), n_ind)
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    expect_true(!any(fixed_loci(X))) # we don't expect any loci to be fixed
    
    # test p_anc
    expect_equal(length(p_anc), m_loci)
    expect_true(all(p_anc >= 0)) # all are non-negative
    expect_true(all(p_anc <= 1)) # all are smaller or equal than 1
})

test_that("draw_all_admix `require_polymorphic_loci = FALSE` works well", {
    # testing FALSE here since TRUE is default
    
    m_loci <- 1000
    inbr_subpops <- c(0.1, 0.2, 0.3)
    k_subpops <- length(inbr_subpops)
    admix_proportions <- diag(rep.int(1, k_subpops)) # island model for test...
    # repeat so we have multiple people per island
    admix_proportions <- rbind(admix_proportions, admix_proportions, admix_proportions)
    n_ind <- nrow(admix_proportions) # number of individuals (3*k_subpops)
    
    # run draw_all_admix
    out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci, want_p_ind = TRUE, want_p_subpops = TRUE, require_polymorphic_loci = FALSE)
    expect_equal( names(out), draw_all_admix_names_ret_full )
    X <- out$X # genotypes
    p_ind <- out$p_ind # IAFs
    p_subpops <- out$p_subpops # Intermediate AFs
    p_anc <- out$p_anc # Ancestral AFs
    
    # test X
    expect_equal(nrow(X), m_loci)
    expect_equal(ncol(X), n_ind)
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    # here loci may be fixed, don't require otherwise

    # test p_ind
    expect_equal(nrow(p_ind), m_loci)
    expect_equal(ncol(p_ind), n_ind)
    expect_true(all(p_ind >= 0)) # all are non-negative
    expect_true(all(p_ind <= 1)) # all are smaller or equal than 1

    # test p_subpops
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1

    # test p_anc
    expect_equal(length(p_anc), m_loci)
    expect_true(all(p_anc >= 0)) # all are non-negative
    expect_true(all(p_anc <= 1)) # all are smaller or equal than 1
})

test_that("draw_all_admix `want_p_ind = FALSE` works well", {
    # really tests low-memory scenario, which is the default `want_p_ind = FALSE` case (but originally we tested non-default `want_p_ind = TRUE` instead)
    m_loci <- 1000
    inbr_subpops <- c(0.1, 0.2, 0.3)
    k_subpops <- length(inbr_subpops)
    admix_proportions <- diag(rep.int(1, k_subpops)) # island model for test...
    # repeat so we have multiple people per island
    admix_proportions <- rbind(admix_proportions, admix_proportions, admix_proportions)
    n_ind <- nrow(admix_proportions) # number of individuals (3*k_subpops)
    
    # run draw_all_admix
    out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci, want_p_subpops = TRUE)
    expect_equal( names(out), draw_all_admix_names_ret_low_mem )
    X <- out$X # genotypes
    p_subpops <- out$p_subpops # Intermediate AFs
    p_anc <- out$p_anc # Ancestral AFs
    
    # test X
    expect_equal(nrow(X), m_loci)
    expect_equal(ncol(X), n_ind)
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    expect_true(!any(fixed_loci(X))) # we don't expect any loci to be fixed
    
    # test p_subpops
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    
    # test p_anc
    expect_equal(length(p_anc), m_loci)
    expect_true(all(p_anc >= 0)) # all are non-negative
    expect_true(all(p_anc <= 1)) # all are smaller or equal than 1
})

test_that("draw_all_admix with provided p_anc (scalar) works well", {
    m_loci <- 10
    inbr_subpops <- c(0.1, 0.2, 0.3)
    k_subpops <- length(inbr_subpops)
    admix_proportions <- diag(rep.int(1, k_subpops)) # island model for test...
    # repeat so we have multiple people per island
    admix_proportions <- rbind(admix_proportions, admix_proportions, admix_proportions)
    n_ind <- nrow(admix_proportions) # number of individuals (3*k_subpops)
    # this is key, to pass so all loci have p_anc 0.2, passed as scalar (not length m_loci)
    p_anc <- 0.2
    
    # run draw_all_admix
    # only test default (p_ind and p_subpops not returned)
    out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci, p_anc = p_anc)
    expect_equal( names(out), draw_all_admix_names_ret_default )

    X <- out$X # genotypes
    
    # test X
    expect_equal(nrow(X), m_loci)
    expect_equal(ncol(X), n_ind)
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    expect_true(!any(fixed_loci(X))) # we don't expect any loci to be fixed
    
    # test p_anc, should just match what we passed
    expect_equal( out$p_anc, p_anc )
})

test_that("draw_all_admix with provided p_anc (vector) works well", {
    m_loci <- 10
    inbr_subpops <- c(0.1, 0.2, 0.3)
    k_subpops <- length(inbr_subpops)
    admix_proportions <- diag(rep.int(1, k_subpops)) # island model for test...
    # repeat so we have multiple people per island
    admix_proportions <- rbind(admix_proportions, admix_proportions, admix_proportions)
    n_ind <- nrow(admix_proportions) # number of individuals (3*k_subpops)
    # construct p_anc separately here, but will demand that the code use it without changes 
    p_anc <- runif( m_loci )
    
    # run draw_all_admix
    # only test default (p_ind and p_subpops not returned)
    out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci, p_anc = p_anc)
    expect_equal( names(out), draw_all_admix_names_ret_default )

    X <- out$X # genotypes
    
    # test X
    expect_equal(nrow(X), m_loci)
    expect_equal(ncol(X), n_ind)
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    expect_true(!any(fixed_loci(X))) # we don't expect any loci to be fixed
    
    # test p_anc, should just match what we passed
    expect_equal( out$p_anc, p_anc )
})

