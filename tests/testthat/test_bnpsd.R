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
    k <- 2
    sigma <- 1
    admix_proportions <- admix_prop_1d_linear(n, k, sigma)
    F <- 1:k # scale doesn't matter right now...

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
    k <- length(unique(labs))
    admix_proportions <- admix_prop_indep_subpops(labs)
    # general tests for admixture matrices
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k) # k columns
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
    expect_equal(ncol(admix_proportions), k) # k columns
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1,n)) # rows sum to 1, vector length n
    # specific tests for admix_prop_indep_subpops
    expect_true(all(admix_proportions %in% c(TRUE, FALSE)))
    expect_true(all(colnames(admix_proportions) == subpops))
    
    # test with provided subpops (additional labels)
    k <- 10
    subpops <- 1:k
    admix_proportions <- admix_prop_indep_subpops(labs, subpops)
    # general tests for admixture matrices
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k) # k columns
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
    k <- 2
    admix_proportions <- admix_prop_1d_linear(n, k, sigma = 1)
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k) # k columns
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n

    # test with sigma == 0 (special case that makes usual formula break)
    admix_proportions <- admix_prop_1d_linear(n, k, sigma = 0)
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k) # k columns
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n
    # in this case it should equal independent subpopulations
    labs <- c( rep.int(1, 5), rep.int(2, 5) ) # two subpops
    admix_proportions2 <- admix_prop_indep_subpops(labs)
    dimnames(admix_proportions2) <- NULL # before comparing, must toss column names
    expect_equal(admix_proportions, admix_proportions2)

    # test bias_coeff version
    obj <- admix_prop_1d_linear(n, k, bias_coeff = 0.5, coanc_subpops = 1:k, fst = 0.1)
    admix_proportions <- obj$admix_proportions # returns many things in this case, get admix_proportions here
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k) # k columns
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n
})

test_that("admix_prop_1d_circular returns valid admixture coefficients", {
    n <- 10
    k <- 2
    admix_proportions <- admix_prop_1d_circular(n, k, sigma = 1)
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k) # k columns
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n

    # test with sigma == 0 (special case that makes usual formula break)
    admix_proportions <- admix_prop_1d_circular(n, k, sigma = 0)
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k) # k columns
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n
    # though the result is nearly island-like, there is an annoying shift I'd rather not try to figure out for this test...

    # test bias_coeff version
    obj <- admix_prop_1d_circular(n, k, bias_coeff = 0.5, coanc_subpops = 1:k, fst = 0.1)
    admix_proportions <- obj$admix_proportions # returns many things in this case, get admix_proportions here
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k) # k columns
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n
})

test_that("bias_coeff_admix_fit agrees with reverse func", {
    n <- 1000
    inbr_subpops <- c(0.1, 0.2, 0.3)
    k <- length(inbr_subpops)
    s_want <- 0.5
    coord_ind_first_linear <- 0.5
    coord_ind_last_linear <- k + 0.5
    coord_ind_first_circular <- 0
    coord_ind_last_circular <- 2 * pi

    # test with admix_prop_1d_linear
    # get sigma
    sigma <- bias_coeff_admix_fit(
        bias_coeff = s_want,
        coanc_subpops = inbr_subpops,
        n_ind = n,
        k_subpops = k,
        func = admix_prop_1d_linear,
        coord_ind_first = coord_ind_first_linear,
        coord_ind_last = coord_ind_last_linear
    )
    # construct everything and verify s == s_want
    admix_proportions <- admix_prop_1d_linear(
        n_ind = n,
        k_subpops = k,
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
        k_subpops = k,
        func = admix_prop_1d_circular,
        coord_ind_first = coord_ind_first_circular,
        coord_ind_last = coord_ind_last_circular
    ) # get sigma
    # construct everything and verify s == s_want
    admix_proportions <- admix_prop_1d_circular(
        n_ind = n,
        k_subpops = k,
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
        k_subpops = k,
        func = admix_prop_1d_linear,
        coord_ind_first = coord_ind_first_linear,
        coord_ind_last = coord_ind_last_linear
    )
    expect_true( is.infinite(sigma) ) # only `sigma = Inf` should achieve the max
    # construct everything and verify s == s_want
    admix_proportions <- admix_prop_1d_linear(
        n_ind = n,
        k_subpops = k,
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
        k_subpops = k,
        func = admix_prop_1d_circular,
        coord_ind_first = coord_ind_first_circular,
        coord_ind_last = coord_ind_last_circular
    )
    expect_true( is.infinite(sigma) ) # only `sigma = Inf` should achieve the max
    # construct everything and verify s == s_want
    admix_proportions <- admix_prop_1d_circular(
        n_ind = n,
        k_subpops = k,
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
    admix_prop_bias_coeff_min <- admix_prop_1d_linear(n, k, sigma = 0)
    # this is the mminimum s_want
    s_want <- bias_coeff_admix(admix_prop_bias_coeff_min, inbr_subpops)
    # get sigma
    sigma <- bias_coeff_admix_fit(
        bias_coeff = s_want,
        coanc_subpops = inbr_subpops,
        n_ind = n,
        k_subpops = k,
        func = admix_prop_1d_linear,
        coord_ind_first = coord_ind_first_linear,
        coord_ind_last = coord_ind_last_linear
    )
    expect_equal( sigma, 0 ) # only `sigma = 0` should achieve the min
    # construct everything and verify s == s_want
    admix_proportions <- admix_prop_1d_linear(
        n_ind = n,
        k_subpops = k,
        sigma = sigma,
        coord_ind_first = coord_ind_first_linear,
        coord_ind_last = coord_ind_last_linear
    ) # now get admix_proportions from there
    coancestry <- coanc_admix(admix_proportions, inbr_subpops)
    s <- mean(coancestry) / mean(diag(coancestry)) # this is the correct bias coeff, with uniform weights
    expect_equal(s, s_want)
    
    # test with admix_prop_1d_circular
    # construct directly
    admix_prop_bias_coeff_min <- admix_prop_1d_circular(n, k, sigma = 0)
    # this is the mminimum s_want
    s_want <- bias_coeff_admix(admix_prop_bias_coeff_min, inbr_subpops)
    # get sigma
    sigma <- bias_coeff_admix_fit(
        bias_coeff = s_want,
        coanc_subpops = inbr_subpops,
        n_ind = n,
        k_subpops = k,
        func = admix_prop_1d_circular,
        coord_ind_first = coord_ind_first_circular,
        coord_ind_last = coord_ind_last_circular
    )
    expect_equal( sigma, 0 ) # only `sigma = 0` should achieve the min
    # construct everything and verify s == s_want
    admix_proportions <- admix_prop_1d_circular(
        n_ind = n,
        k_subpops = k,
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
    k <- 3
    # non-trivial coancestry matrix for intermediate subpopulations!
    coanc_subpops <- matrix(
        c(
            0.1 , 0.05, 0   ,
            0.05, 0.2 , 0.05,
            0   , 0.05, 0.3
        ),
        nrow = k
    )
    s_want <- 0.5

    # test with admix_prop_1d_linear
    # get sigma
    sigma <- bias_coeff_admix_fit(
        bias_coeff = s_want,
        coanc_subpops = coanc_subpops,
        n_ind = n,
        k_subpops = k,
        func = admix_prop_1d_linear,
        coord_ind_first = coord_ind_first_linear,
        coord_ind_last = coord_ind_last_linear
    )
    # construct everything and verify s == s_want
    admix_proportions <- admix_prop_1d_linear(
        n_ind = n,
        k_subpops = k,
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
        k_subpops = k,
        func = admix_prop_1d_circular,
        coord_ind_first = coord_ind_first_circular,
        coord_ind_last = coord_ind_last_circular
    ) # get sigma
    # construct everything and verify s == s_want
    admix_proportions <- admix_prop_1d_circular(
        n_ind = n,
        k_subpops = k,
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
    k <- 2
    sigma <- 1
    Fst <- 0.1
    admix_proportions <- admix_prop_1d_linear(n, k, sigma)
    F <- 1:k # scale doesn't matter right now...
    F2 <- rescale_coanc_subpops(admix_proportions, F, Fst) # calculation to compare to
    coancestry <- coanc_admix(admix_proportions, F2) # in wrong scale but meh
    Fst2 <- mean(diag(coancestry)) # this is the actual FST, with uniform weights
    expect_equal(Fst, Fst2)
    # since 0 < Fst=0.1 < 1, there's nothing else to test
})

test_that("draw_p_anc is in range", {
    m <- 1000
    p_anc <- draw_p_anc(m)
    expect_equal(length(p_anc), m)
    expect_true(all(p_anc >= 0)) # all are non-negative
    expect_true(all(p_anc <= 1)) # all are smaller or equal than 1
})

test_that("draw_p_subpops is in range", {
    m <- 1000
    F <- c(0.1, 0.2, 0.3)
    k <- length(F)
    p_anc <- draw_p_anc(m)
    B <- draw_p_subpops(p_anc, F)
    expect_equal(nrow(B), m)
    expect_equal(ncol(B), k)
    expect_true(all(B >= 0)) # all are non-negative
    expect_true(all(B <= 1)) # all are smaller or equal than 1
})

test_that("rpiaf is in range", {
    m <- 1000
    F <- c(0.1, 0.2, 0.3)
    k <- length(F)
    admix_proportions <- diag(rep.int(1, k)) # island model for test...
    admix_proportions <- rbind(admix_proportions, admix_proportions, admix_proportions) # repeat so we have multiple people per island
    n <- nrow(admix_proportions) # number of individuals (3 * k)
    p_anc <- draw_p_anc(m)
    B <- draw_p_subpops(p_anc, F)
    P <- rpiaf(B, admix_proportions)
    expect_equal(nrow(P), m)
    expect_equal(ncol(P), n)
    expect_true(all(P >= 0)) # all are non-negative
    expect_true(all(P <= 1)) # all are smaller or equal than 1
})

test_that("rgeno is in range", {
    m <- 1000
    F <- c(0.1, 0.2, 0.3)
    k <- length(F)
    admix_proportions <- diag(rep.int(1, k)) # island model for test...
    admix_proportions <- rbind(admix_proportions, admix_proportions, admix_proportions) # repeat so we have multiple people per island
    n <- nrow(admix_proportions) # number of individuals (3 * k)
    p_anc <- draw_p_anc(m)
    B <- draw_p_subpops(p_anc, F)
    P <- rpiaf(B, admix_proportions)
    X <- rgeno(P) # direct test
    expect_equal(nrow(X), m)
    expect_equal(ncol(X), n)
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    X <- rgeno(B, admix_proportions) # indirect draw test
    expect_equal(nrow(X), m)
    expect_equal(ncol(X), n)
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    X <- rgeno(B, admix_proportions, lowMem = TRUE) # indirect draw test with low-mem algo
    expect_equal(nrow(X), m)
    expect_equal(ncol(X), n)
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
})

test_that("rbnpsd works well", {
    m <- 1000
    F <- c(0.1, 0.2, 0.3)
    k <- length(F)
    admix_proportions <- diag(rep.int(1, k)) # island model for test...
    admix_proportions <- rbind(admix_proportions, admix_proportions, admix_proportions) # repeat so we have multiple people per island
    n <- nrow(admix_proportions) # number of individuals (3*k)
    # run rbnpsd
    out <- rbnpsd(admix_proportions, F, m)
    X <- out$X # genotypes
    P <- out$P # IAFs
    B <- out$B # Intermediate AFs
    p_anc <- out$Pa # Ancestral AFs
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
    # test p_anc
    expect_equal(length(p_anc), m)
    expect_true(all(p_anc >= 0)) # all are non-negative
    expect_true(all(p_anc <= 1)) # all are smaller or equal than 1
})

test_that("rbnpsd `noFixed = TRUE` works well", {
    m <- 1000
    F <- c(0.1, 0.2, 0.3)
    k <- length(F)
    admix_proportions <- diag(rep.int(1, k)) # island model for test...
    admix_proportions <- rbind(admix_proportions, admix_proportions, admix_proportions) # repeat so we have multiple people per island
    n <- nrow(admix_proportions) # number of individuals (3*k)
    # run rbnpsd
    out <- rbnpsd(admix_proportions, F, m, noFixed = TRUE)
    X <- out$X # genotypes
    P <- out$P # IAFs
    B <- out$B # Intermediate AFs
    p_anc <- out$Pa # Ancestral AFs
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
    # test p_anc
    expect_equal(length(p_anc), m)
    expect_true(all(p_anc >= 0)) # all are non-negative
    expect_true(all(p_anc <= 1)) # all are smaller or equal than 1
})

test_that("rbnpsd `lowMem = TRUE` works well", {
    m <- 1000
    F <- c(0.1, 0.2, 0.3)
    k <- length(F)
    admix_proportions <- diag(rep.int(1, k)) # island model for test...
    admix_proportions <- rbind(admix_proportions, admix_proportions, admix_proportions) # repeat so we have multiple people per island
    n <- nrow(admix_proportions) # number of individuals (3*k)
    # run rbnpsd
    out <- rbnpsd(admix_proportions, F, m, lowMem = TRUE)
    X <- out$X # genotypes
    B <- out$B # Intermediate AFs
    p_anc <- out$Pa # Ancestral AFs
    # test X
    expect_equal(nrow(X), m)
    expect_equal(ncol(X), n)
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    # test B
    expect_equal(nrow(B), m)
    expect_equal(ncol(B), k)
    expect_true(all(B >= 0)) # all are non-negative
    expect_true(all(B <= 1)) # all are smaller or equal than 1
    # test p_anc
    expect_equal(length(p_anc), m)
    expect_true(all(p_anc >= 0)) # all are non-negative
    expect_true(all(p_anc <= 1)) # all are smaller or equal than 1
})
