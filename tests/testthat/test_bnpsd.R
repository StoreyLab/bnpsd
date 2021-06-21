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
            0, 0, 1, # a somewhat "rare" variant
            NA, NA, NA # completely missing locus (will be treated as fixed)
        ),
        ncol = 3,
        byrow = TRUE
    )
    indexes_expected <- c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE)
    expect_silent(
        indexes <- fixed_loci( X )   
    )
    # test that we get the desired values
    # includes names test implicitly
    expect_equal( indexes, indexes_expected )

    # a version with names
    rownames( X ) <- paste0( 'l', 1 : nrow( X ) )
    colnames( X ) <- paste0( 'i', 1 : ncol( X ) ) # should not be used, but meh
    # test fails unless expected indexes also have names
    names( indexes_expected ) <- rownames( X )
    # repeat test
    expect_silent(
        indexes <- fixed_loci( X )   
    )
    # test that we get the desired values
    # includes names test implicitly
    expect_equal( indexes, indexes_expected )

    # a version with non-trivial `maf_min`
    maf_min <- 1 / 6 # marks only one more locus as fixed, noting that equality marks things too
    indexes_expected[ 5 ] <- TRUE # this one
    # repeat test
    expect_silent(
        indexes <- fixed_loci( X, maf_min = maf_min )   
    )
    # test that we get the desired values
    # includes names test implicitly
    expect_equal( indexes, indexes_expected )

    # cause errors on purpose by passing bad values of `maf_min`
    expect_error( fixed_loci( X, maf_min = 'a' ) )
    expect_error( fixed_loci( X, maf_min = (1:10)/20 ) )
    expect_error( fixed_loci( X, maf_min = NA ) )
    expect_error( fixed_loci( X, maf_min = -1 ) )
    expect_error( fixed_loci( X, maf_min = 0.6 ) )
})

test_that( "names_coanc works", {
    # only has one required argument
    expect_error( names_coanc() )
    
    # scalar returns NULL, regardless of names
    coanc <- 1
    expect_true( is.null( names_coanc( coanc ) ) )
    names(coanc) <- 'a'
    expect_true( is.null( names_coanc( coanc ) ) )

    # vector without names returns NULL
    inbr_subpops <- c(0.1, 0.2, 0.3)
    expect_true( is.null( names_coanc( inbr_subpops ) ) )
    # add names!
    k_subpops <- length(inbr_subpops)
    names( inbr_subpops ) <- paste0( 'S', 1 : k_subpops )
    expect_equal( names_coanc( inbr_subpops ), names( inbr_subpops ) )

    # create a random symmetric matrix
    # first some random uniform data
    coanc_subpops <- matrix(
        runif( k_subpops^2 ),
        k_subpops,
        k_subpops
    )
    # make symmetric
    coanc_subpops <- ( coanc_subpops + t( coanc_subpops ) ) / 2
    # no names first
    expect_true( is.null( names_coanc( coanc_subpops ) ) )
    # add names
    rownames( coanc_subpops ) <- names( inbr_subpops )
    colnames( coanc_subpops ) <- names( inbr_subpops )
    expect_equal( names_coanc( coanc_subpops ), names( inbr_subpops ) )
    
    # data from a tree
    tree <- ape::rtree( k_subpops )
    # already has names, test that case first
    expect_equal( names_coanc( tree ), tree$tip.label )
    # copy without names
    tree_nameless <- tree
    tree_nameless$tip.label <- rep.int( '', k_subpops )
    expect_true( is.null( names_coanc( tree_nameless ) ) )
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
    expect_silent(
        coancestry <- coanc_admix(admix_proportions, F)
    )
    # equality implicitly tests names too
    expect_equal(coancestry, coancestry_expected)

    # add names to inputs and outputs
    rownames( admix_proportions ) <- c('i1', 'i2')
    colnames( admix_proportions ) <- c('S1', 'S2')
    names( F ) <- colnames( admix_proportions )
    colnames( coancestry_expected ) <- rownames( admix_proportions )
    rownames( coancestry_expected ) <- rownames( admix_proportions )
    expect_silent(
        coancestry <- coanc_admix(admix_proportions, F)
    )
    # equality implicitly tests names too
    expect_equal(coancestry, coancestry_expected)
    
    # same admix_proportions, scalar F
    F <- 0.2
    coancestry_expected <- diag(c(F, F)) # the coancestry we expect for this setup
    colnames( coancestry_expected ) <- rownames( admix_proportions )
    rownames( coancestry_expected ) <- rownames( admix_proportions )
    expect_silent(
        coancestry <- coanc_admix(admix_proportions, F)
    )
    # equality implicitly tests names too
    expect_equal(coancestry, coancestry_expected)

    # same admix_proportions, matrix F
    Fv <- c(0.1, 0.4) # vector version
    F <- diag(Fv) # matrix version
    colnames( F ) <- colnames( admix_proportions )
    rownames( F ) <- colnames( admix_proportions )
    coancestry_expected <- F # same as output here, but names should be different
    colnames( coancestry_expected ) <- rownames( admix_proportions )
    rownames( coancestry_expected ) <- rownames( admix_proportions )
    expect_silent(
        coancestry <- coanc_admix(admix_proportions, F)
    )
    # equality implicitly tests names too
    expect_equal(coancestry, coancestry_expected)

    # now create some errors
    # first labels agree in content but not order
    # easiest to scramble admix_proportions
    # there's only two subpopulations, so reversing gives us the only other order here
    colnames( admix_proportions ) <- rev( colnames( admix_proportions ) )
    expect_error( coanc_admix(admix_proportions, F) )
    # completely different names should also trigger error
    colnames( admix_proportions ) <- c('a', 'b')
    expect_error( coanc_admix(admix_proportions, F) )

    # most complex case, just a general math check
    admix_proportions <- matrix(c(0.7, 0.3, 0.2, 0.8), nrow = 2, byrow = TRUE)
    F <- matrix(c(0.3, 0.1, 0.1, 0.3), nrow = 2, byrow = TRUE)
    coancestry_expected <- admix_proportions %*% F %*% t(admix_proportions) # the coancestry we expect for this setup (slower but more explicit version)
    expect_silent(
        coancestry <- coanc_admix(admix_proportions, F)
    )
    # equality implicitly tests names too
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

# expected names for return object when sigma is fit
names_admix_prop_1d <- c('admix_proportions', 'coanc_subpops', 'sigma', 'coanc_factor')

test_that("admix_prop_1d_linear returns valid admixture coefficients", {
    n <- 10
    k_subpops <- 2
    admix_proportions <- admix_prop_1d_linear(n, k_subpops, sigma = 1)
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k_subpops) # k_subpops columns
    expect_true( !anyNA( admix_proportions ) ) # no missing data
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n
    # dimnames should be NULL in this case
    expect_true( is.null( dimnames( admix_proportions ) ) )

    # test with sigma == 0 (special case that makes usual formula break)
    admix_proportions <- admix_prop_1d_linear(n, k_subpops, sigma = 0)
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k_subpops) # k_subpops columns
    expect_true( !anyNA( admix_proportions ) ) # no missing data
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n
    # in this case it should equal independent subpopulations
    labs <- c( rep.int(1, 5), rep.int(2, 5) ) # two subpops
    admix_proportions2 <- admix_prop_indep_subpops(labs)
    dimnames(admix_proportions2) <- NULL # before comparing, must toss column names
    expect_equal(admix_proportions, admix_proportions2)
    # dimnames should be NULL in this case
    expect_true( is.null( dimnames( admix_proportions ) ) )

    # test bias_coeff version
    expect_silent(
        obj <- admix_prop_1d_linear(n, k_subpops, bias_coeff = 0.5, coanc_subpops = 1:k_subpops, fst = 0.1)
    )
    expect_true( is.list( obj ) )
    expect_equal( names( obj ), names_admix_prop_1d )
    admix_proportions <- obj$admix_proportions # returns many things in this case, get admix_proportions here
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k_subpops) # k_subpops columns
    expect_true( !anyNA( admix_proportions ) ) # no missing data
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n
    # dimnames should be NULL in this case
    expect_true( is.null( dimnames( admix_proportions ) ) )

    # test edge case that used to give NA coefficients (now dies if that happens)
    # these params were verified to give NAs in the previous version
    n <- 10
    k_subpops <- 10
    sigma <- 0.01
    expect_silent(
        admix_proportions <- admix_prop_1d_linear( n_ind = n, k_subpops = k_subpops, sigma = sigma)
    )
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k_subpops) # k_subpops columns
    expect_true( !anyNA( admix_proportions ) ) # no missing data
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n
    # dimnames should be NULL in this case
    expect_true( is.null( dimnames( admix_proportions ) ) )

    # test transfer of names for coanc_subpops vector version (must be bias_coeff version)
    # NOTE: case of names for matrix coanc_subpops is tested later where we test application for trees
    inbr_subpops <- 1 : k_subpops
    names( inbr_subpops ) <- letters[ inbr_subpops ]
    expect_silent(
        obj <- admix_prop_1d_linear(n, k_subpops, bias_coeff = 0.5, coanc_subpops = inbr_subpops, fst = 0.1)
    )
    expect_true( is.list( obj ) )
    expect_equal( names( obj ), names_admix_prop_1d )
    admix_proportions <- obj$admix_proportions # returns many things in this case, get admix_proportions here
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k_subpops) # k_subpops columns
    expect_true( !anyNA( admix_proportions ) ) # no missing data
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n
    # dimnames should be non-NULL in this case
    expect_true( !is.null( dimnames( admix_proportions ) ) )
    # no row names (individuals don't have natural names here)
    expect_true( is.null( rownames( admix_proportions ) ) )
    # but column names should be the inbr_subpops names
    expect_equal( colnames( admix_proportions ), names( inbr_subpops ) )
})

test_that("admix_prop_1d_circular returns valid admixture coefficients", {
    n <- 10
    k_subpops <- 2
    admix_proportions <- admix_prop_1d_circular(n, k_subpops, sigma = 1)
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k_subpops) # k_subpops columns
    expect_true( !anyNA( admix_proportions ) ) # no missing data
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n
    # dimnames should be NULL in this case
    expect_true( is.null( dimnames( admix_proportions ) ) )

    # test with sigma == 0 (special case that makes usual formula break)
    admix_proportions <- admix_prop_1d_circular(n, k_subpops, sigma = 0)
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k_subpops) # k_subpops columns
    expect_true( !anyNA( admix_proportions ) ) # no missing data
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n
    # dimnames should be NULL in this case
    expect_true( is.null( dimnames( admix_proportions ) ) )
    # though the result is nearly island-like, there is an annoying shift I'd rather not try to figure out for this test...

    # test bias_coeff version
    expect_silent(
        obj <- admix_prop_1d_circular(n, k_subpops, bias_coeff = 0.5, coanc_subpops = 1:k_subpops, fst = 0.1)
    )
    expect_true( is.list( obj ) )
    expect_equal( names( obj ), names_admix_prop_1d )
    admix_proportions <- obj$admix_proportions # returns many things in this case, get admix_proportions here
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k_subpops) # k_subpops columns
    expect_true( !anyNA( admix_proportions ) ) # no missing data
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n
    # dimnames should be NULL in this case
    expect_true( is.null( dimnames( admix_proportions ) ) )
    
    # test edge case that used to give NA coefficients (now dies if that happens)
    # these params were verified to give NAs in the previous version
    n <- 10
    k_subpops <- 10
    sigma <- 0.01
    expect_silent(
        admix_proportions <- admix_prop_1d_circular( n_ind = n, k_subpops = k_subpops, sigma = sigma)
    )
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k_subpops) # k_subpops columns
    expect_true( !anyNA( admix_proportions ) ) # no missing data
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n
    # dimnames should be NULL in this case
    expect_true( is.null( dimnames( admix_proportions ) ) )

    # test transfer of names for coanc_subpops vector version (must be bias_coeff version)
    # NOTE: case of names for matrix coanc_subpops is tested later where we test application for trees
    inbr_subpops <- 1 : k_subpops
    names( inbr_subpops ) <- letters[ inbr_subpops ]
    expect_silent(
        obj <- admix_prop_1d_circular(n, k_subpops, bias_coeff = 0.5, coanc_subpops = inbr_subpops, fst = 0.1)
    )
    expect_true( is.list( obj ) )
    expect_equal( names( obj ), names_admix_prop_1d )
    admix_proportions <- obj$admix_proportions # returns many things in this case, get admix_proportions here
    expect_equal(nrow(admix_proportions), n) # n rows
    expect_equal(ncol(admix_proportions), k_subpops) # k_subpops columns
    expect_true( !anyNA( admix_proportions ) ) # no missing data
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n)) # rows sum to 1, vector length n
    # dimnames should be non-NULL in this case
    expect_true( !is.null( dimnames( admix_proportions ) ) )
    # no row names (individuals don't have natural names here)
    expect_true( is.null( rownames( admix_proportions ) ) )
    # but column names should be the inbr_subpops names
    expect_equal( colnames( admix_proportions ), names( inbr_subpops ) )
    
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
    expect_silent( 
        obj <- rescale_coanc_subpops(admix_proportions, F, Fst) # calculation to compare to
    )
    F2 <- obj$coanc_subpops
    coancestry <- coanc_admix(admix_proportions, F2) # in wrong scale but meh
    Fst2 <- mean(diag(coancestry)) # this is the actual FST, with uniform weights
    expect_equal(Fst, Fst2)
    # since 0 < Fst=0.1 < 1, there's nothing else to test

    # test factor too
    factor <- obj$factor
    expect_equal( length( factor ), 1 )
    expect_equal( factor, mean( F2 ) / mean( F ) )
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
    # error if p_anc vector is out of range
    expect_error( draw_p_subpops( 1000 * p_anc, inbr_subpops ) )
    # or negative
    expect_error( draw_p_subpops( - p_anc, inbr_subpops ) )
    # error if inbr vector is out of range
    expect_error( draw_p_subpops( p_anc, 10 * inbr_subpops ) )
    # or negative
    expect_error( draw_p_subpops( p_anc, - inbr_subpops ) )
    
    # test main use case
    expect_silent(
        p_subpops <- draw_p_subpops(p_anc, inbr_subpops)
    )
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    # no names in this case
    expect_true( is.null( dimnames( p_subpops ) ) )

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
    # no names in this case
    expect_true( is.null( dimnames( p_subpops ) ) )
    
    # special case of scalar inbr_subpops
    p_subpops <- draw_p_subpops(p_anc, inbr_subpops = 0.2, k_subpops = k_subpops)
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    # no names in this case
    expect_true( is.null( dimnames( p_subpops ) ) )
    
    # both main parameters scalars but return value still matrix
    p_subpops <- draw_p_subpops(p_anc = 0.5, inbr_subpops = 0.2, m_loci = m_loci, k_subpops = k_subpops)
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    # no names in this case
    expect_true( is.null( dimnames( p_subpops ) ) )
    
    # passing scalar parameters without setting dimensions separately results in a 1x1 matrix
    p_subpops <- draw_p_subpops(p_anc = 0.5, inbr_subpops = 0.2)
    expect_equal(nrow(p_subpops), 1)
    expect_equal(ncol(p_subpops), 1)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    # no names in this case
    expect_true( is.null( dimnames( p_subpops ) ) )

    # now a case with names, in this case for both inputs
    names( p_anc ) <- paste0( 'l', 1 : m_loci )
    names( inbr_subpops ) <- paste0( 'S', 1 : k_subpops )
    expect_silent(
        p_subpops <- draw_p_subpops(p_anc, inbr_subpops)
    )
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    # test names
    expect_true( !is.null( dimnames( p_subpops ) ) )
    expect_equal( rownames( p_subpops ), names( p_anc ) )
    expect_equal( colnames( p_subpops ), names( inbr_subpops ) )
})

test_that("make_p_ind_admix is in range", {
    m_loci <- 10
    inbr_subpops <- c(0.1, 0.2, 0.3)
    k_subpops <- length(inbr_subpops)
    admix_proportions <- diag(rep.int(1, k_subpops)) # island model for test...
    # repeat so we have multiple people per island
    admix_proportions <- rbind(admix_proportions, admix_proportions, admix_proportions)
    n_ind <- nrow(admix_proportions) # number of individuals (3 * k_subpops)
    p_anc <- draw_p_anc(m_loci)
    p_subpops <- draw_p_subpops(p_anc, inbr_subpops)
    # start test
    expect_silent(
        p_ind <- make_p_ind_admix(p_subpops, admix_proportions)
    )
    expect_equal(nrow(p_ind), m_loci)
    expect_equal(ncol(p_ind), n_ind)
    expect_true(all(p_ind >= 0)) # all are non-negative
    expect_true(all(p_ind <= 1)) # all are smaller or equal than 1
    # no names in this case
    expect_true( is.null( dimnames( p_ind ) ) )

    # a case with full names for all inputs
    rownames( p_subpops ) <- paste0( 'l', 1 : m_loci )
    colnames( p_subpops ) <- paste0( 'S', 1 : k_subpops )
    rownames( admix_proportions ) <- paste0( 'i', 1 : n_ind )
    colnames( admix_proportions ) <- colnames( p_subpops ) # these two ought to match
    # repeat test
    expect_silent(
        p_ind <- make_p_ind_admix(p_subpops, admix_proportions)
    )
    expect_equal(nrow(p_ind), m_loci)
    expect_equal(ncol(p_ind), n_ind)
    expect_true(all(p_ind >= 0)) # all are non-negative
    expect_true(all(p_ind <= 1)) # all are smaller or equal than 1
    # test all names
    expect_equal( rownames( p_ind ), rownames( p_subpops ) )
    expect_equal( colnames( p_ind ), rownames( admix_proportions ) )

    # cause errors on purpose by having admix_proportions names that disagree with p_subpops
    # first just change order of subpopulations (in names only)
    colnames( admix_proportions ) <- rev( colnames( admix_proportions ) )
    expect_error( make_p_ind_admix(p_subpops, admix_proportions) )
    # change names entirely
    colnames( admix_proportions ) <- letters[ 1 : k_subpops ]
    expect_error( make_p_ind_admix(p_subpops, admix_proportions) )
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
    expect_true( !anyNA( X ) ) # no missing values
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    # dimnames should be NULL in this case
    expect_true( is.null( dimnames( X ) ) )
    
    # indirect draw test (default low_mem now!)
    X <- draw_genotypes_admix(p_subpops, admix_proportions)
    expect_equal(nrow(X), m_loci)
    expect_equal(ncol(X), n_ind)
    expect_true( !anyNA( X ) ) # no missing values
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    # dimnames should be NULL in this case
    expect_true( is.null( dimnames( X ) ) )
    
    # test names transfers
    # first add names to p_ind
    rownames( p_ind ) <- paste0( 'l', 1 : m_loci )
    colnames( p_ind ) <- paste0( 'i', 1 : n_ind )
    # now run code
    X <- draw_genotypes_admix(p_ind)
    expect_equal(nrow(X), m_loci)
    expect_equal(ncol(X), n_ind)
    expect_true( !anyNA( X ) ) # no missing values
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    # dimnames should be non-NULL in this case
    expect_true( !is.null( dimnames( X ) ) )
    # and test each dimension
    expect_equal( rownames( X ), rownames( p_ind ) )
    expect_equal( colnames( X ), colnames( p_ind ) )

    # and the same but for version with admix_proportions
    rownames( p_subpops ) <- rownames( p_ind )
    colnames( p_subpops ) <- paste0( 'S', 1 : k_subpops )
    rownames( admix_proportions ) <- colnames( p_ind )
    colnames( admix_proportions ) <- colnames( p_subpops )
    X <- draw_genotypes_admix(p_subpops, admix_proportions)
    expect_equal(nrow(X), m_loci)
    expect_equal(ncol(X), n_ind)
    expect_true( !anyNA( X ) ) # no missing values
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    # dimnames should be non-NULL in this case
    expect_true( !is.null( dimnames( X ) ) )
    # and test each dimension
    expect_equal( rownames( X ), rownames( p_subpops ) )
    expect_equal( colnames( X ), rownames( admix_proportions ) )

    # cause errors, suggesting misalignment or inconsistency of p_subpops and admix_proportions, via name disagreements
    # reverse names first
    colnames( admix_proportions ) <- rev( colnames( admix_proportions ) )
    expect_error( draw_genotypes_admix(p_subpops, admix_proportions) )
    # replace names completely
    colnames( admix_proportions ) <- letters[ 1 : k_subpops ]
    expect_error( draw_genotypes_admix(p_subpops, admix_proportions) )
    # fix names again for last test
    colnames( admix_proportions ) <- colnames( p_subpops )
        
    # do last because we change p_subpops and that messes up name transfer tests
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
    expect_true( !anyNA( X ) ) # no missing values
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    # this is a cool test because only one of the two dimnames is non-NULL
    # dimnames should be non-NULL in this case, because `admix_proportions` had names
    expect_true( !is.null( dimnames( X ) ) )
    # rownames should be null
    expect_true( is.null( rownames( X ) ) )
    # column names should match rownames of `admix_proportions`
    expect_equal( colnames( X ), rownames( admix_proportions ) )

})

# for recurring tests
draw_all_admix_names_ret_default <- c('X', 'p_anc')
draw_all_admix_names_ret_full <- c('X', 'p_anc', 'p_subpops', 'p_ind')
draw_all_admix_names_ret_low_mem <- c('X', 'p_anc', 'p_subpops') # excludes p_ind

test_that("draw_all_admix works", {
    # let's use names by default, these were tested in separate parts before but let's just add a test for the "all" version
    m_loci <- 10
    inbr_subpops <- c(0.1, 0.2, 0.3)
    k_subpops <- length(inbr_subpops)
    names( inbr_subpops ) <- paste0( 'S', 1 : k_subpops )
    admix_proportions <- diag(rep.int(1, k_subpops)) # island model for test...
    # repeat so we have multiple people per island
    admix_proportions <- rbind(admix_proportions, admix_proportions, admix_proportions)
    n_ind <- nrow(admix_proportions) # number of individuals (3*k_subpops)
    colnames( admix_proportions ) <- names( inbr_subpops )
    rownames( admix_proportions ) <- paste0( 'i', 1 : n_ind )
    
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
    expect_silent(
        out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci)
    )
    expect_equal( names(out), draw_all_admix_names_ret_default )

    # now rerun with all outputs, so we can test them all
    expect_silent(
        out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci, want_p_ind = TRUE, want_p_subpops = TRUE)
    )
    expect_equal( names(out), draw_all_admix_names_ret_full )
    X <- out$X # genotypes
    p_ind <- out$p_ind # IAFs
    p_subpops <- out$p_subpops # Intermediate AFs
    p_anc <- out$p_anc # Ancestral AFs
    
    # test X
    expect_equal(nrow(X), m_loci)
    expect_equal(ncol(X), n_ind)
    expect_true( !anyNA( X ) ) # no missing values
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    expect_true(!any(fixed_loci(X))) # we don't expect any loci to be fixed
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( rownames( X ) ) )
    # individuals do have names
    expect_equal( colnames( X ), rownames( admix_proportions ) )
    
    # test p_ind
    expect_equal(nrow(p_ind), m_loci)
    expect_equal(ncol(p_ind), n_ind)
    expect_true(all(p_ind >= 0)) # all are non-negative
    expect_true(all(p_ind <= 1)) # all are smaller or equal than 1
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( rownames( p_ind ) ) )
    # individuals do have names
    expect_equal( colnames( p_ind ), rownames( admix_proportions ) )
    
    # test p_subpops
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( rownames( p_subpops ) ) )
    # subpopulations do have names
    expect_equal( colnames( p_subpops ), colnames( admix_proportions ) )
    
    # test p_anc
    expect_equal(length(p_anc), m_loci)
    expect_true(all(p_anc >= 0)) # all are non-negative
    expect_true(all(p_anc <= 1)) # all are smaller or equal than 1
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( names( p_anc ) ) )

    # cause errors on purpose by making labels disagree
    # reverse names
    colnames( admix_proportions ) <- rev( colnames( admix_proportions ) )
    expect_error( draw_all_admix(admix_proportions, inbr_subpops, m_loci) )
    # completely replace names
    colnames( admix_proportions ) <- letters[ 1 : k_subpops ]
    expect_error( draw_all_admix(admix_proportions, inbr_subpops, m_loci) )

    # reset admix_proportions names
    colnames( admix_proportions ) <- names( inbr_subpops )
    # test case of scalar inbr_subpops
    inbr_subpops <- 0.1
    # ask for p_subpops as another means into testing that the correct k_subpops was used!
    expect_silent(
        out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci, want_p_subpops = TRUE)
    )
    expect_equal( names(out), draw_all_admix_names_ret_low_mem ) # exclude p_ind
    X <- out$X # genotypes
    p_ind <- out$p_ind # IAFs
    p_subpops <- out$p_subpops # Intermediate AFs
    p_anc <- out$p_anc # Ancestral AFs
    
    # test X
    expect_equal(nrow(X), m_loci)
    expect_equal(ncol(X), n_ind)
    expect_true( !anyNA( X ) ) # no missing values
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    expect_true(!any(fixed_loci(X))) # we don't expect any loci to be fixed
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( rownames( X ) ) )
    # individuals do have names
    expect_equal( colnames( X ), rownames( admix_proportions ) )
    
    # test p_subpops
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( rownames( p_subpops ) ) )
    # subpopulations do have names
    expect_equal( colnames( p_subpops ), colnames( admix_proportions ) )
    
    # test p_anc
    expect_equal(length(p_anc), m_loci)
    expect_true(all(p_anc >= 0)) # all are non-negative
    expect_true(all(p_anc <= 1)) # all are smaller or equal than 1
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( names( p_anc ) ) )
    
})

test_that("draw_all_admix with `maf_min > 0` works", {
    # let's use names by default, these were tested in separate parts before but let's just add a test for the "all" version
    m_loci <- 10
    inbr_subpops <- c(0.1, 0.2, 0.3)
    k_subpops <- length(inbr_subpops)
    names( inbr_subpops ) <- paste0( 'S', 1 : k_subpops )
    admix_proportions <- diag(rep.int(1, k_subpops)) # island model for test...
    # repeat so we have multiple people per island
    admix_proportions <- rbind(admix_proportions, admix_proportions, admix_proportions)
    n_ind <- nrow(admix_proportions) # number of individuals (3*k_subpops)
    colnames( admix_proportions ) <- names( inbr_subpops )
    rownames( admix_proportions ) <- paste0( 'i', 1 : n_ind )

    # here's the key parameter
    # there's only 9 individuals, so minimum non-zero freq is 1/18.
    # let's pick something larger to make test non-trivial
    maf_min <- 1/5
    
    # run draw_all_admix
    # only test default (p_ind and p_subpops not returned)
    # only X should be different anyway
    expect_silent(
        out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci, maf_min = maf_min)
    )
    expect_equal( names(out), draw_all_admix_names_ret_default )

    X <- out$X # genotypes
    p_anc <- out$p_anc # Ancestral AFs
    
    # test X
    expect_equal(nrow(X), m_loci)
    expect_equal(ncol(X), n_ind)
    expect_true( !anyNA( X ) ) # no missing values
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    expect_true(!any(fixed_loci( X, maf_min = maf_min ))) # we don't expect any loci to be fixed (use same MAF threshold here!)
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( rownames( X ) ) )
    # individuals do have names
    expect_equal( colnames( X ), rownames( admix_proportions ) )
    
    # test p_anc
    expect_equal(length(p_anc), m_loci)
    expect_true(all(p_anc >= 0)) # all are non-negative
    expect_true(all(p_anc <= 1)) # all are smaller or equal than 1
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( names( p_anc ) ) )

    # repeat test with more stringent threshold
    # here we force our code to redraw most loci lots of times, which in a specific real case (with very rare variants `p_anc`, though `maf_min = 0`, we just want to force large numbers of redraws) resulted in a "node stack overflow".
    # there's only 9 individuals, so minimum non-zero freq is 1/18.
    # pick something large but must be below 0.5 = 9/18
    # this combination reliably triggered the error in the original version, yet it has a good runtime in the fixed version (the crazy number of iterations is still super fast for toy example).
    maf_min <- 8/18
    p_anc <- 0.1
    # actual errors produced in this test prior to fix:
    # "Error: C stack usage  7971780 is too close to the limit"
    
    # run draw_all_admix
    # only test default (p_ind and p_subpops not returned)
    # only X should be different anyway
    expect_silent(
        out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci, maf_min = maf_min, p_anc = p_anc)
    )
    expect_equal( names(out), draw_all_admix_names_ret_default )

    X <- out$X # genotypes
    p_anc <- out$p_anc # Ancestral AFs
    
    # test X
    expect_equal(nrow(X), m_loci)
    expect_equal(ncol(X), n_ind)
    expect_true( !anyNA( X ) ) # no missing values
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    expect_true(!any(fixed_loci( X, maf_min = maf_min ))) # we don't expect any loci to be fixed (use same MAF threshold here!)
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( rownames( X ) ) )
    # individuals do have names
    expect_equal( colnames( X ), rownames( admix_proportions ) )
    
    # test p_anc
    expect_equal(length(p_anc), m_loci)
    expect_true(all(p_anc >= 0)) # all are non-negative
    expect_true(all(p_anc <= 1)) # all are smaller or equal than 1
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( names( p_anc ) ) )
})

test_that("draw_all_admix beta works", {
    # let's use names by default, these were tested in separate parts before but let's just add a test for the "all" version
    m_loci <- 10
    beta <- 0.01
    inbr_subpops <- c(0.1, 0.2, 0.3)
    k_subpops <- length(inbr_subpops)
    names( inbr_subpops ) <- paste0( 'S', 1 : k_subpops )
    admix_proportions <- diag(rep.int(1, k_subpops)) # island model for test...
    # repeat so we have multiple people per island
    admix_proportions <- rbind(admix_proportions, admix_proportions, admix_proportions)
    n_ind <- nrow(admix_proportions) # number of individuals (3*k_subpops)
    colnames( admix_proportions ) <- names( inbr_subpops )
    rownames( admix_proportions ) <- paste0( 'i', 1 : n_ind )

    # run draw_all_admix
    # only test default (p_ind and p_subpops not returned)
    out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci, beta = beta)
    expect_equal( names(out), draw_all_admix_names_ret_default )

    X <- out$X # genotypes
    p_anc <- out$p_anc # Ancestral AFs
    
    # test X
    expect_equal(nrow(X), m_loci)
    expect_equal(ncol(X), n_ind)
    expect_true( !anyNA( X ) ) # no missing values
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    expect_true(!any(fixed_loci(X))) # we don't expect any loci to be fixed
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( rownames( X ) ) )
    # individuals do have names
    expect_equal( colnames( X ), rownames( admix_proportions ) )
    
    # test p_anc
    expect_equal(length(p_anc), m_loci)
    expect_true(all(p_anc >= 0)) # all are non-negative
    expect_true(all(p_anc <= 1)) # all are smaller or equal than 1
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( names( p_anc ) ) )

    # NOTE: incompatible labels between admix_proportions and inbr_subpops not retested because that happens early and is independent of this and all subsequent feature tests (except tree case)
})

test_that("draw_all_admix `require_polymorphic_loci = FALSE` works", {
    # testing FALSE here since TRUE is default
    
    m_loci <- 1000
    inbr_subpops <- c(0.1, 0.2, 0.3)
    k_subpops <- length(inbr_subpops)
    names( inbr_subpops ) <- paste0( 'S', 1 : k_subpops )
    admix_proportions <- diag(rep.int(1, k_subpops)) # island model for test...
    # repeat so we have multiple people per island
    admix_proportions <- rbind(admix_proportions, admix_proportions, admix_proportions)
    n_ind <- nrow(admix_proportions) # number of individuals (3*k_subpops)
    colnames( admix_proportions ) <- names( inbr_subpops )
    rownames( admix_proportions ) <- paste0( 'i', 1 : n_ind )
    
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
    expect_true( !anyNA( X ) ) # no missing values
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    # here loci may be fixed, don't require otherwise
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( rownames( X ) ) )
    # individuals do have names
    expect_equal( colnames( X ), rownames( admix_proportions ) )

    # test p_ind
    expect_equal(nrow(p_ind), m_loci)
    expect_equal(ncol(p_ind), n_ind)
    expect_true(all(p_ind >= 0)) # all are non-negative
    expect_true(all(p_ind <= 1)) # all are smaller or equal than 1
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( rownames( p_ind ) ) )
    # individuals do have names
    expect_equal( colnames( p_ind ), rownames( admix_proportions ) )

    # test p_subpops
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( rownames( p_subpops ) ) )
    # subpopulations do have names
    expect_equal( colnames( p_subpops ), colnames( admix_proportions ) )

    # test p_anc
    expect_equal(length(p_anc), m_loci)
    expect_true(all(p_anc >= 0)) # all are non-negative
    expect_true(all(p_anc <= 1)) # all are smaller or equal than 1
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( names( p_anc ) ) )
})

test_that("draw_all_admix `want_p_ind = FALSE` works", {
    # really tests low-memory scenario, which is the default `want_p_ind = FALSE` case (but originally we tested non-default `want_p_ind = TRUE` instead)
    m_loci <- 1000
    inbr_subpops <- c(0.1, 0.2, 0.3)
    k_subpops <- length(inbr_subpops)
    names( inbr_subpops ) <- paste0( 'S', 1 : k_subpops )
    admix_proportions <- diag(rep.int(1, k_subpops)) # island model for test...
    # repeat so we have multiple people per island
    admix_proportions <- rbind(admix_proportions, admix_proportions, admix_proportions)
    n_ind <- nrow(admix_proportions) # number of individuals (3*k_subpops)
    colnames( admix_proportions ) <- names( inbr_subpops )
    rownames( admix_proportions ) <- paste0( 'i', 1 : n_ind )
    
    # run draw_all_admix
    out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci, want_p_subpops = TRUE)
    expect_equal( names(out), draw_all_admix_names_ret_low_mem )
    X <- out$X # genotypes
    p_subpops <- out$p_subpops # Intermediate AFs
    p_anc <- out$p_anc # Ancestral AFs
    
    # test X
    expect_equal(nrow(X), m_loci)
    expect_equal(ncol(X), n_ind)
    expect_true( !anyNA( X ) ) # no missing values
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    expect_true(!any(fixed_loci(X))) # we don't expect any loci to be fixed
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( rownames( X ) ) )
    # individuals do have names
    expect_equal( colnames( X ), rownames( admix_proportions ) )
    
    # test p_subpops
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( rownames( p_subpops ) ) )
    # subpopulations do have names
    expect_equal( colnames( p_subpops ), colnames( admix_proportions ) )
    
    # test p_anc
    expect_equal(length(p_anc), m_loci)
    expect_true(all(p_anc >= 0)) # all are non-negative
    expect_true(all(p_anc <= 1)) # all are smaller or equal than 1
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( names( p_anc ) ) )
})

test_that("draw_all_admix with provided p_anc (scalar) works", {
    m_loci <- 10
    inbr_subpops <- c(0.1, 0.2, 0.3)
    k_subpops <- length(inbr_subpops)
    names( inbr_subpops ) <- paste0( 'S', 1 : k_subpops )
    admix_proportions <- diag(rep.int(1, k_subpops)) # island model for test...
    # repeat so we have multiple people per island
    admix_proportions <- rbind(admix_proportions, admix_proportions, admix_proportions)
    n_ind <- nrow(admix_proportions) # number of individuals (3*k_subpops)
    colnames( admix_proportions ) <- names( inbr_subpops )
    rownames( admix_proportions ) <- paste0( 'i', 1 : n_ind )
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
    expect_true( !anyNA( X ) ) # no missing values
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    expect_true(!any(fixed_loci(X))) # we don't expect any loci to be fixed
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( rownames( X ) ) )
    # individuals do have names
    expect_equal( colnames( X ), rownames( admix_proportions ) )
    
    # test p_anc, should just match what we passed
    # but always returns vector, compare as vector
    expect_equal( out$p_anc, rep.int( p_anc, m_loci ) )
    # this case shouldn't have names
    expect_true( is.null( names( p_anc ) ) )
})

test_that("draw_all_admix with provided p_anc (vector) works", {
    m_loci <- 10
    inbr_subpops <- c(0.1, 0.2, 0.3)
    k_subpops <- length(inbr_subpops)
    names( inbr_subpops ) <- paste0( 'S', 1 : k_subpops )
    admix_proportions <- diag(rep.int(1, k_subpops)) # island model for test...
    # repeat so we have multiple people per island
    admix_proportions <- rbind(admix_proportions, admix_proportions, admix_proportions)
    n_ind <- nrow(admix_proportions) # number of individuals (3*k_subpops)
    colnames( admix_proportions ) <- names( inbr_subpops )
    rownames( admix_proportions ) <- paste0( 'i', 1 : n_ind )
    # construct p_anc separately here, but will demand that the code use it without changes 
    p_anc <- runif( m_loci )
    # loci should have names in this case!!!
    names( p_anc ) <- paste0( 'l', 1 : m_loci )
    
    # run draw_all_admix
    # only test default (p_ind and p_subpops not returned)
    out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci, p_anc = p_anc)
    expect_equal( names(out), draw_all_admix_names_ret_default )

    X <- out$X # genotypes
    
    # test X
    expect_equal(nrow(X), m_loci)
    expect_equal(ncol(X), n_ind)
    expect_true( !anyNA( X ) ) # no missing values
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    expect_true(!any(fixed_loci(X))) # we don't expect any loci to be fixed
    # loci and individuals have names!
    expect_equal( rownames( X ), names( p_anc ) )
    expect_equal( colnames( X ), rownames( admix_proportions ) )
    
    # test p_anc, should just match what we passed
    # (implicitly tests names too)
    expect_equal( out$p_anc, p_anc )
})

test_that( "validate_coanc_tree works", {
    # draw a random tree that is valid
    k_subpops <- 5
    # this one automatically draws edge lengths between 0 and 1, as desired
    tree <- ape::rtree( k_subpops )
    # let's test a tree without names, still valid
    tree_no_names <- ape::rtree( k_subpops, tip.label = rep.int( '', k_subpops ) )
    # unrooted tree (bad example)
    tree_unrooted <- ape::rtree( k_subpops, rooted = FALSE )

    # errors due to missing arguments
    # only `tree` is required
    expect_error( validate_coanc_tree() )
    # pass things that are not `phylo` objects
    expect_error( validate_coanc_tree( 'a' ) )
    expect_error( validate_coanc_tree( 1 ) )
    # a correct tree except without class
    tree_bad <- unclass( tree )
    expect_error( validate_coanc_tree( tree_bad ) )
    # a correct tree except class name is different
    tree_bad <- tree
    class( tree_bad ) <- 'phylo_coanc'
    expect_error( validate_coanc_tree( tree_bad ) )
    # a list with the right class but nothing else is right
    tree_bad <- list( a = 1 )
    class( tree_bad ) <- 'phylo'
    expect_error( validate_coanc_tree( tree_bad ) )
    # unrooted trees are not acceptable
    expect_error( validate_coanc_tree( tree_unrooted ) )
    
    # successful runs
    expect_silent( validate_coanc_tree( tree ) )
    expect_silent( validate_coanc_tree( tree, name = 'treebeard' ) )
    expect_silent( validate_coanc_tree( tree_no_names ) )
    expect_silent( validate_coanc_tree( tree_no_names, name = 'treebeard' ) )

    # valid trees except small details are wrong
    # here delete each of the 4 mandatory elements in turn
    tree_bad <- tree
    tree_bad$edge <- NULL # delete this element
    expect_error( validate_coanc_tree( tree_bad ) )
    tree_bad <- tree
    tree_bad$edge.length <- NULL # delete this element
    expect_error( validate_coanc_tree( tree_bad ) )
    tree_bad <- tree
    tree_bad$tip.label <- NULL # delete this element
    expect_error( validate_coanc_tree( tree_bad ) )
    tree_bad <- tree
    tree_bad$Nnode <- NULL # delete this element
    expect_error( validate_coanc_tree( tree_bad ) )
    # make an edge have value greater than 1
    tree_bad <- tree
    tree_bad$edge.length[1] <- 10
    expect_error( validate_coanc_tree( tree_bad ) )
    # negative edges
    tree_bad <- tree
    tree_bad$edge.length[1] <- -0.1
    expect_error( validate_coanc_tree( tree_bad ) )
    # mess with dimensions
    tree_bad <- tree
    tree_bad$edge.length <- tree_bad$edge.length[ -1 ]
    expect_error( validate_coanc_tree( tree_bad ) )
    tree_bad <- tree
    tree_bad$edge <- tree_bad$edge[ , -1 ]
    expect_error( validate_coanc_tree( tree_bad ) )
    tree_bad <- tree
    tree_bad$edge <- tree_bad$edge[ -1, ]
    expect_error( validate_coanc_tree( tree_bad ) )
    tree_bad <- tree
    tree_bad$Nnode <- tree_bad$Nnode + 1
    expect_error( validate_coanc_tree( tree_bad ) )
    tree_bad <- tree
    tree_bad$tip.label <- tree_bad$tip.label[ -1 ]
    expect_error( validate_coanc_tree( tree_bad ) )
    # mess with `tree$edge` table
    # max index is wrong
    tree_bad <- tree
    i_max <- max( tree_bad$edge )
    tree_bad$edge[ tree_bad$edge == i_max ] <- i_max + 1
    expect_error( validate_coanc_tree( tree_bad ) )
    # some other index is missing
    tree_bad <- tree
    tree_bad$edge[ tree_bad$edge == i_max - 1 ] <- i_max
    expect_error( validate_coanc_tree( tree_bad ) )
    
    # test that additional elements don't cause problems
    tree_ext <- tree
    tree_ext$extra <- 'a'
    expect_silent( validate_coanc_tree( tree_ext ) )
    # test that extended class doesn't cause problems
    tree_ext <- tree
    class( tree_ext ) <- c( class( tree_ext ), 'phylo_coanc' )
    expect_silent( validate_coanc_tree( tree_ext ) )

    # test root edge cases
    # first a good case
    tree_ext <- tree
    tree_ext$root.edge <- 0.9 # an extreme but acceptable case
    expect_silent( validate_coanc_tree( tree_ext ) )
    # now bad (out of range) cases
    tree_ext$root.edge <- -0.1
    expect_error( validate_coanc_tree( tree_ext ) )
    tree_ext$root.edge <- 1.1
    expect_error( validate_coanc_tree( tree_ext ) )
})

test_that( "draw_p_subpops_tree works", {
    k_subpops <- 5
    m_loci <- 10

    # simulate a random tree
    # this one automatically draws edge lengths between 0 and 1, as desired
    tree_subpops <- ape::rtree( k_subpops )
    # let's test a tree without names
    tree_subpops_no_names <- ape::rtree( k_subpops, tip.label = rep.int( '', k_subpops ) )
    # draw allele frequencies too
    p_anc <- runif( m_loci )
    # allele frequency vector with names
    p_anc_names <- p_anc
    names( p_anc_names ) <- paste0( 'l', 1 : m_loci )
    
    # create a version where internal nodes have names too
    tree_subpops_full_names <- tree_subpops
    tree_subpops_full_names$node.label <- paste0( 'n', 1 : tree_subpops_full_names$Nnode )
    
    # cause errors on purpose
    # first two arguments are required
    expect_error( draw_p_subpops_tree( ) )
    expect_error( draw_p_subpops_tree( p_anc = p_anc ) )
    expect_error( draw_p_subpops_tree( tree_subpops = tree_subpops ) )

    # now a successful run
    expect_silent(
        p_subpops <- draw_p_subpops_tree(
            p_anc = p_anc,
            tree_subpops = tree_subpops
        )
    )
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    # names are inherited from tree
    expect_equal( colnames( p_subpops ), tree_subpops$tip.label )
    # no names from p_anc
    expect_true( is.null( rownames( p_subpops ) ) )

    # version with names from p_anc
    expect_silent(
        p_subpops <- draw_p_subpops_tree(
            p_anc = p_anc_names,
            tree_subpops = tree_subpops
        )
    )
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    # names are inherited from tree
    expect_equal( colnames( p_subpops ), tree_subpops$tip.label )
    # match names from p_anc_names!
    expect_equal( rownames( p_subpops ), names( p_anc_names ) )
    
    # expect a warning if there's a root edge
    tree_subpops_warn <- tree_subpops
    tree_subpops_warn$root.edge <- 0.5
    expect_warning(
        p_subpops <- draw_p_subpops_tree(
            p_anc = p_anc,
            tree_subpops = tree_subpops_warn
        )
    )
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    # names are inherited from tree
    expect_equal( colnames( p_subpops ), tree_subpops_warn$tip.label )
    # no names from p_anc
    expect_true( is.null( rownames( p_subpops ) ) )

    # request AFs for all internal nodes too
    expect_silent(
        p_subpops <- draw_p_subpops_tree(
            p_anc = p_anc,
            tree_subpops = tree_subpops,
            nodes = TRUE
        )
    )
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), 2 * k_subpops - 1 ) # there should be these many columns instead (assumes bifurcating tree, which rtree does return)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    # names are inherited from tree (this only has tip names)
    expect_equal( colnames( p_subpops )[ 1 : k_subpops ], tree_subpops$tip.label )
    # test that non-tips names are blank
    expect_true( all( colnames( p_subpops )[ ( k_subpops + 1 ) : ( 2 * k_subpops - 1 ) ] == '' ) )
    # no names from p_anc
    expect_true( is.null( rownames( p_subpops ) ) )
    
    # pass scalar p_anc, specify m_loci separately
    expect_silent(
        p_subpops <- draw_p_subpops_tree(
            p_anc = 0.5,
            tree_subpops = tree_subpops,
            m_loci = m_loci
        )
    )
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    # names are inherited from tree
    expect_equal( colnames( p_subpops ), tree_subpops$tip.label )
    # no names from p_anc
    expect_true( is.null( rownames( p_subpops ) ) )
    
    # pass scalar p_anc, but don't specify m_loci
    expect_silent(
        p_subpops <- draw_p_subpops_tree(
            p_anc = 0.5,
            tree_subpops = tree_subpops
        )
    )
    expect_equal(nrow(p_subpops), 1 ) # expect a single locus
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    # names are inherited from tree
    expect_equal( colnames( p_subpops ), tree_subpops$tip.label )
    # no names from p_anc
    expect_true( is.null( rownames( p_subpops ) ) )

    # run a case with a tree without names
    expect_silent(
        p_subpops <- draw_p_subpops_tree(
            p_anc = p_anc,
            tree_subpops = tree_subpops_no_names
        )
    )
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    # names should be null
    expect_true( is.null( colnames( p_subpops ) ) )
    # no names from p_anc
    expect_true( is.null( rownames( p_subpops ) ) )

    # run a case with full names, but only tips are returned
    expect_silent(
        p_subpops <- draw_p_subpops_tree(
            p_anc = p_anc,
            tree_subpops = tree_subpops_full_names
        )
    )
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    # names are inherited from tree
    expect_equal( colnames( p_subpops ), tree_subpops_full_names$tip.label )
    # no names from p_anc
    expect_true( is.null( rownames( p_subpops ) ) )

    # run a case with full names, but all nodes are returned
    expect_silent(
        p_subpops <- draw_p_subpops_tree(
            p_anc = p_anc,
            tree_subpops = tree_subpops_full_names,
            nodes = TRUE
        )
    )
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), 2 * k_subpops - 1 ) # there should be these many columns instead (assumes bifurcating tree, which rtree does return)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    # names are inherited from tree (this only has tip names)
    expect_equal( colnames( p_subpops )[ 1 : k_subpops ], tree_subpops_full_names$tip.label )
    # test that non-tips names are blank
    expect_equal( colnames( p_subpops )[ ( k_subpops + 1 ) : ( 2 * k_subpops - 1 ) ], tree_subpops_full_names$node.label )
    # no names from p_anc
    expect_true( is.null( rownames( p_subpops ) ) )

    # a successful run with randomized edges (used to cause problems due to assumption that this wasn't allowed)
    tree_subpops_rand <- tree_subpops
    tree_subpops_rand$edge <- tree_subpops_rand$edge[ sample( ape::Nedge( tree_subpops_rand ) ), ]
    expect_silent(
        p_subpops <- draw_p_subpops_tree(
            p_anc = p_anc,
            tree_subpops = tree_subpops_rand
        )
    )
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    # names are inherited from tree
    expect_equal( colnames( p_subpops ), tree_subpops_rand$tip.label )
    # no names from p_anc
    expect_true( is.null( rownames( p_subpops ) ) )
    
})

test_that( "tree_additive works", {
    # draw a random tree that is valid
    k_subpops <- 5
    # this one automatically draws edge lengths between 0 and 1, as desired
    tree <- ape::rtree( k_subpops )

    # only errors are is if tree is missing or is not a tree
    expect_error( tree_additive() )
    expect_error( tree_additive( 1 : k_subpops ) )

    # now a successful run
    expect_silent(
        tree2 <- tree_additive( tree )
    )
    # validate tree
    expect_silent(
        validate_coanc_tree( tree2 )
    )
    # check for data that extends beyond regular tree
    expect_true( !is.null( tree2$edge.length.add ) )
    expect_equal( length( tree2$edge.length.add ), length( tree2$edge.length ) )
    expect_true( all( tree2$edge.length.add >= 0 ) )
    expect_true( all( tree2$edge.length.add <= 1 ) )

    # a successful run with an extreme but valid root edge
    tree$root.edge <- 0.9
    expect_silent(
        tree2 <- tree_additive( tree )
    )
    # validate tree
    expect_silent(
        validate_coanc_tree( tree2 )
    )
    # check for data that extends beyond regular tree
    expect_true( !is.null( tree2$edge.length.add ) )
    expect_equal( length( tree2$edge.length.add ), length( tree2$edge.length ) )
    expect_true( all( tree2$edge.length.add >= 0 ) )
    expect_true( all( tree2$edge.length.add <= 1 ) )
    # root edge doesn't change
    expect_equal( tree$root.edge, tree2$root.edge )

    # a successful run with randomized edges (used to cause problems due to assumption that this wasn't allowed)
    tree_rand <- tree
    # reorder edges fully so output tree matches last tree2 (upon the same reordering)
    indexes <- sample( ape::Nedge( tree_rand ) )
    tree_rand$edge <- tree_rand$edge[ indexes, ]
    tree_rand$edge.length <- tree_rand$edge.length[ indexes ]
    expect_silent(
        tree_rand2 <- tree_additive( tree_rand )
    )
    # if we reorder previous output, we expect to have the same exact tree
    tree_rand2_exp <- tree2
    tree_rand2_exp$edge <- tree_rand2_exp$edge[ indexes, ]
    tree_rand2_exp$edge.length <- tree_rand2_exp$edge.length[ indexes ]
    tree_rand2_exp$edge.length.add <- tree_rand2_exp$edge.length.add[ indexes ] # additive edges too!
    expect_equal( tree_rand2, tree_rand2_exp )
})


test_that( "tree_additive rev works", {
    # must construct an additive tree
    # best way is to use inverse function `tree_additive`, makes test easy

    # draw a random tree that is valid
    k_subpops <- 5
    # this one automatically draws edge lengths between 0 and 1, as desired
    tree_prob <- ape::rtree( k_subpops )
    # this calculates additive edges, but has them separate (still overall a probabilistic-edge tree)
    tree_prob <- tree_additive( tree_prob )
    # desired input to function to test, have to overwrite edges and remove extra entry (it is checked to not exist)
    tree_add <- tree_prob
    tree_add$edge.length <- tree_add$edge.length.add
    # causes error if passed a tree that already has additive edges calculated/stored
    expect_error( tree_additive( tree_add, rev = TRUE ) )
    # but goes on if `force = TRUE` (only test of this option)
    expect_silent( tree_additive( tree_add, rev = TRUE, force = TRUE ) )
    # now actually remove this element, for rest of tests to not need `force`
    tree_add$edge.length.add <- NULL
    
    # causes error if passed a tree that already has additive edges calculated/stored
    expect_error( tree_additive( tree_prob, rev = TRUE ) )
    
    # now a successful run
    expect_silent(
        tree_prob2 <- tree_additive( tree_add, rev = TRUE )
    )
    # confirms that `rev = TRUE` is inverse function of `rev = FALSE`!
    expect_equal( tree_prob, tree_prob2 )

    # expect error when a non-additive tree is passed
    tree_bad <- tree_add
    # force max edge to equal 1, so its sum to anything is sure to exceed 1
    tree_bad$edge.length <- 1.1 * tree_bad$edge.length / max( tree_bad$edge.length )
    expect_error( tree_additive( tree_bad, rev = TRUE ) )

    # construct test with extreme but valid root edge
    # redraw tree to have a clean slate
    tree_prob <- ape::rtree( k_subpops )
    # add root edge
    tree_prob$root.edge <- 0.9
    # repeat other steps
    tree_prob <- tree_additive( tree_prob )
    tree_add <- tree_prob
    tree_add$edge.length <- tree_add$edge.length.add
    tree_add$edge.length.add <- NULL
    # actual test
    expect_silent(
        tree_prob2 <- tree_additive( tree_add, rev = TRUE )
    )
    expect_equal( tree_prob, tree_prob2 )

    # a successful run with randomized edges (used to cause problems due to assumption that this wasn't allowed)
    # redraw tree to have a clean slate
    tree_prob <- ape::rtree( k_subpops )
    # randomize right away
    tree_prob$edge <- tree_prob$edge[ sample( ape::Nedge( tree_prob ) ), ]
    # repeat other steps
    tree_prob <- tree_additive( tree_prob )
    tree_add <- tree_prob
    tree_add$edge.length <- tree_add$edge.length.add
    tree_add$edge.length.add <- NULL
    # actual test
    expect_silent(
        tree_prob2 <- tree_additive( tree_add, rev = TRUE )
    )
    expect_equal( tree_prob, tree_prob2 )
})

test_that( "coanc_tree works", {
    # draw a random tree that is valid
    k_subpops <- 5
    # this one automatically draws edge lengths between 0 and 1, as desired
    tree <- ape::rtree( k_subpops )

    # only errors are is if tree is missing or is not a tree
    expect_error( coanc_tree() )
    expect_error( coanc_tree( 1 : k_subpops ) )

    # now a successful run
    expect_silent(
        coanc_mat <- coanc_tree( tree )
    )
    expect_true( is.numeric( coanc_mat ) )
    expect_true( !anyNA( coanc_mat ) )
    expect_true( is.matrix( coanc_mat ) )
    expect_equal( nrow( coanc_mat ), k_subpops )
    expect_equal( ncol( coanc_mat ), k_subpops )
    expect_true( isSymmetric( coanc_mat ) )
    expect_true( all( coanc_mat >= 0 ) )
    expect_true( all( coanc_mat <= 1 ) ) # because in additive scale everything should be below 1
    expect_equal( rownames( coanc_mat ), tree$tip.label ) # earlier isSymmetric ensures colnames is the same, not testing again
    # answer should also be positive-semidefinite at least, but meh we don't test that here

    # repeat with a tree with an extreme but valid root edge
    tree$root.edge <- 0.9
    expect_silent(
        coanc_mat <- coanc_tree( tree )
    )
    expect_true( is.numeric( coanc_mat ) )
    expect_true( !anyNA( coanc_mat ) )
    expect_true( is.matrix( coanc_mat ) )
    expect_equal( nrow( coanc_mat ), k_subpops )
    expect_equal( ncol( coanc_mat ), k_subpops )
    expect_true( isSymmetric( coanc_mat ) )
    expect_true( all( coanc_mat >= 0 ) )
    expect_true( all( coanc_mat <= 1 ) )
    expect_equal( rownames( coanc_mat ), tree$tip.label ) # earlier isSymmetric ensures colnames is the same, not testing again

    # repeat with randomized edges (don't expect errors, but let's just check)
    # keep root edge, meh
    tree_rand <- tree
    indexes <- sample( ape::Nedge( tree_rand ) )
    # to expect equality to previous run, `$edge` and `$edge.length` have to be reordered together
    tree_rand$edge <- tree_rand$edge[ indexes, ]
    tree_rand$edge.length <- tree_rand$edge.length[ indexes ]
    expect_silent(
        coanc_mat_rand <- coanc_tree( tree_rand )
    )
    # should match last output (makes things easy!)
    expect_equal( coanc_mat_rand, coanc_mat )
})

test_that("draw_all_admix works with a tree", {
    # draw a random tree that is valid
    k_subpops <- 3
    # this one automatically draws edge lengths between 0 and 1, as desired
    # these always have names!
    tree_subpops <- ape::rtree( k_subpops )
    # other info
    m_loci <- 10
    admix_proportions <- diag(rep.int(1, k_subpops)) # island model for test...
    # repeat so we have multiple people per island
    admix_proportions <- rbind(admix_proportions, admix_proportions, admix_proportions)
    n_ind <- nrow(admix_proportions) # number of individuals (3*k_subpops)
    colnames( admix_proportions ) <- tree_subpops$tip.label
    rownames( admix_proportions ) <- paste0( 'i', 1 : n_ind )

    # make sure things die when important parameters are missing
    # all missing
    expect_error( draw_all_admix() )
    # two missing
    expect_error( draw_all_admix(admix_proportions = admix_proportions) )
    expect_error( draw_all_admix(tree_subpops = tree_subpops) )
    expect_error( draw_all_admix(m_loci = m_loci) )
    # one missing
    expect_error( draw_all_admix(tree_subpops = tree_subpops, m_loci = m_loci) )
    expect_error( draw_all_admix(admix_proportions = admix_proportions, m_loci = m_loci) )
    expect_error( draw_all_admix(admix_proportions = admix_proportions, tree_subpops = tree_subpops) )

    # expect error if both inbr_subpops and tree_subpops are passed
    expect_error( draw_all_admix(admix_proportions, inbr_subpops = c(0.1, 0.2, 0.3), tree_subpops = tree_subpops, m_loci = m_loci) )
    
    # run draw_all_admix
    # first test default (p_ind and p_subpops not returned)
    expect_silent( 
        out <- draw_all_admix(admix_proportions, tree_subpops = tree_subpops, m_loci = m_loci)
    )
    expect_equal( names(out), draw_all_admix_names_ret_default )

    # now rerun with all outputs, so we can test them all
    expect_silent( 
        out <- draw_all_admix(admix_proportions, tree_subpops = tree_subpops, m_loci = m_loci, want_p_ind = TRUE, want_p_subpops = TRUE)
    )
    expect_equal( names(out), draw_all_admix_names_ret_full )
    X <- out$X # genotypes
    p_ind <- out$p_ind # IAFs
    p_subpops <- out$p_subpops # Intermediate AFs
    p_anc <- out$p_anc # Ancestral AFs
    
    # test X
    expect_equal(nrow(X), m_loci)
    expect_equal(ncol(X), n_ind)
    expect_true( !anyNA( X ) ) # no missing values
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    expect_true(!any(fixed_loci(X))) # we don't expect any loci to be fixed
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( rownames( X ) ) )
    # individuals do have names
    expect_equal( colnames( X ), rownames( admix_proportions ) )
    
    # test p_ind
    expect_equal(nrow(p_ind), m_loci)
    expect_equal(ncol(p_ind), n_ind)
    expect_true(all(p_ind >= 0)) # all are non-negative
    expect_true(all(p_ind <= 1)) # all are smaller or equal than 1
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( rownames( p_ind ) ) )
    # individuals do have names
    expect_equal( colnames( p_ind ), rownames( admix_proportions ) )
    
    # test p_subpops
    expect_equal(nrow(p_subpops), m_loci)
    expect_equal(ncol(p_subpops), k_subpops)
    expect_true(all(p_subpops >= 0)) # all are non-negative
    expect_true(all(p_subpops <= 1)) # all are smaller or equal than 1
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( rownames( p_subpops ) ) )
    # subpopulations do have names
    expect_equal( colnames( p_subpops ), colnames( admix_proportions ) )
    
    # test p_anc
    expect_equal(length(p_anc), m_loci)
    expect_true(all(p_anc >= 0)) # all are non-negative
    expect_true(all(p_anc <= 1)) # all are smaller or equal than 1
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( names( p_anc ) ) )

    # repeat tests with randomized edges
    tree_rand <- tree_subpops
    tree_rand$edge <- tree_rand$edge[ sample( ape::Nedge( tree_rand )), ]
    # test default outputs only
    expect_silent( 
        out <- draw_all_admix(admix_proportions, tree_subpops = tree_rand, m_loci = m_loci)
    )
    expect_equal( names(out), draw_all_admix_names_ret_default )
    X <- out$X # genotypes
    p_anc <- out$p_anc # Ancestral AFs
    # test X
    expect_equal(nrow(X), m_loci)
    expect_equal(ncol(X), n_ind)
    expect_true( !anyNA( X ) ) # no missing values
    expect_true(all(X %in% c(0, 1, 2))) # only three values allowed!
    expect_true(!any(fixed_loci(X))) # we don't expect any loci to be fixed
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( rownames( X ) ) )
    # individuals do have names
    expect_equal( colnames( X ), rownames( admix_proportions ) )
    # test p_anc
    expect_equal(length(p_anc), m_loci)
    expect_true(all(p_anc >= 0)) # all are non-negative
    expect_true(all(p_anc <= 1)) # all are smaller or equal than 1
    # p_anc was simulated inside function (not passed) so loci don't have names
    expect_true( is.null( names( p_anc ) ) )
    
    # expect a warning if there's a root edge
    tree_subpops_warn <- tree_subpops
    tree_subpops_warn$root.edge <- 0.5
    expect_warning( 
        out <- draw_all_admix(admix_proportions, tree_subpops = tree_subpops_warn, m_loci = m_loci)
    )
    # won't test outputs otherwise...
    expect_equal( names(out), draw_all_admix_names_ret_default )

    # cause errors on purpose by making labels disagree
    # reverse names
    colnames( admix_proportions ) <- rev( colnames( admix_proportions ) )
    expect_error( draw_all_admix(admix_proportions, tree_subpops = tree_subpops, m_loci = m_loci) )
    # completely replace names
    colnames( admix_proportions ) <- letters[ 1 : k_subpops ]
    expect_error( draw_all_admix(admix_proportions, tree_subpops = tree_subpops, m_loci = m_loci) )
})

test_that("admix_prop_1d_linear/circular bias_coeff work with tree", {
    # we don't pass a tree, but rather pass it's coancestry matrix, make sure that works
    # NOTE: no need to test trees with randomized edge orders, since that was already tested for `coanc_tree` and it doesn't come up otherwise

    # random trees are dangerous here, because minimum bias coeff can vary wildly
    # so instead pass a fixed, reasonable tree, which has a minimum bias_coeff that has been pretested to be less than the value we ask for below
    tree_str <- '(S1:0.1,(S2:0.1,(S3:0.1,(S4:0.1,S5:0.1)N3:0.1)N2:0.1)N1:0.1)T;'
    tree_subpops <- ape::read.tree( text = tree_str )
    # its true coancestry
    coanc_subpops <- coanc_tree( tree_subpops )
    k_subpops <- nrow( coanc_subpops )
    # other parameters required
    n_ind <- 10
    bias_coeff <- 0.6
    fst <- 0.1
    
    # try a successful run
    expect_silent(
        obj <- admix_prop_1d_linear(
            n_ind,
            k_subpops,
            bias_coeff = bias_coeff,
            coanc_subpops = coanc_subpops,
            fst = fst
        )
    )
    expect_true( is.list( obj ) )
    expect_equal( names( obj ), names_admix_prop_1d )
    admix_proportions <- obj$admix_proportions # returns many things in this case, get admix_proportions here
    expect_equal(nrow(admix_proportions), n_ind) # n rows
    expect_equal(ncol(admix_proportions), k_subpops) # k_subpops columns
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n_ind)) # rows sum to 1, vector length n
    # dimnames should be non-NULL in this case
    expect_true( !is.null( dimnames( admix_proportions ) ) )
    # no row names (individuals don't have natural names here)
    expect_true( is.null( rownames( admix_proportions ) ) )
    # but column names should be tree labels!
    expect_equal( colnames( admix_proportions ), tree_subpops$tip.label )

    # repeat with circular rather than linear version
    expect_silent(
        obj <- admix_prop_1d_circular(
            n_ind,
            k_subpops,
            bias_coeff = bias_coeff,
            coanc_subpops = coanc_subpops,
            fst = fst
        )
    )
    expect_true( is.list( obj ) )
    expect_equal( names( obj ), names_admix_prop_1d )
    admix_proportions <- obj$admix_proportions # returns many things in this case, get admix_proportions here
    expect_equal(nrow(admix_proportions), n_ind) # n rows
    expect_equal(ncol(admix_proportions), k_subpops) # k_subpops columns
    expect_true(all(admix_proportions >= 0)) # all are non-negative
    expect_true(all(admix_proportions <= 1)) # all are smaller or equal than 1
    expect_equal(rowSums(admix_proportions), rep.int(1, n_ind)) # rows sum to 1, vector length n
    # dimnames should be non-NULL in this case
    expect_true( !is.null( dimnames( admix_proportions ) ) )
    # no row names (individuals don't have natural names here)
    expect_true( is.null( rownames( admix_proportions ) ) )
    # but column names should be tree labels!
    expect_equal( colnames( admix_proportions ), tree_subpops$tip.label )

})

test_that( "is_ordered_tips_edges, tree_reindex_tips work", {
    # create a random tree
    # these always have ordered tips by construction
    k_subpops <- 5
    tree <- ape::rtree( k_subpops )

    # only one mandatory argument (both functions)
    expect_error( is_ordered_tips_edges() )
    expect_error( tree_reindex_tips() )

    # a succesful run
    expect_silent(
        test <- is_ordered_tips_edges( tree )
    )
    expect_true( test )

    # tree reindexing should not change the tree in this case
    expect_silent(
        tree2 <- tree_reindex_tips( tree )
    )
    expect_equal( tree2, tree )

    # purposefully reverse order
    tree_rev <- tree
    tree_rev$edge <- tree_rev$edge[ ape::Nedge( tree_rev ) : 1, ]
    # this should be false
    expect_silent(
        test <- is_ordered_tips_edges( tree_rev )
    )
    expect_false( test )

    # reindex reversed tree
    expect_silent(
        tree_rev <- tree_reindex_tips( tree_rev )
    )
    # now this should pass "ordered" test
    expect_silent(
        test <- is_ordered_tips_edges( tree_rev )
    )
    expect_true( test )
    # and because the reordering was exact reversal, then we should see that in the tip labels
    expect_equal( tree_rev$tip.label, rev( tree$tip.label ) )

    # test reordering under additive edges
    tree <- tree_additive( tree )
    # tree reindexing should not change the tree in this case
    expect_silent(
        tree2 <- tree_reindex_tips( tree )
    )
    expect_equal( tree2, tree )
    # purposefully reverse order
    tree_rev <- tree
    tree_rev$edge <- tree_rev$edge[ ape::Nedge( tree_rev ) : 1, ]
    # not sure what values to expect exactly for additive edges in this case, we just want non-failure
    expect_silent(
        tree_rev <- tree_reindex_tips( tree_rev )
    )
    expect_silent(
        test <- is_ordered_tips_edges( tree_rev )
    )
    expect_true( test )
    expect_equal( tree_rev$tip.label, rev( tree$tip.label ) )

    # make sure code doesn't die under complete random edge reorderings
    tree_rand <- tree
    indexes <- sample( ape::Nedge( tree_rand ) )
    tree_rand$edge <- tree_rand$edge[ indexes, ]
    # this should work, but by random chance we might get tips in the right order, so this could be true or false
    expect_silent(
        test <- is_ordered_tips_edges( tree_rand )
    )
    expect_true( is.logical( test ) )
    # reindexing should also work
    expect_silent(
        tree_rand <- tree_reindex_tips( tree_rand )
    )
    # this time tips must be ordered
    expect_silent(
        test <- is_ordered_tips_edges( tree_rand )
    )
    expect_true( test )
})

test_that( "tree_reorder works", {
    # create a random tree
    # these always have ordered tips by construction
    k_subpops <- 10
    tree <- ape::rtree( k_subpops )
    # get true labels, set as desired ones
    labels <- tree$tip.label
    # scramble tree
    tree_scrambled <- tree
    indexes <- sample( ape::Nedge( tree ) )
    tree_scrambled$edge <- tree_scrambled$edge[ indexes, ]
    tree_scrambled$edge.length <- tree_scrambled$edge.length[ indexes ]

    # check for missing arguments
    expect_error( tree_reorder( ) )
    expect_error( tree_reorder( tree ) )
    expect_error( tree_reorder( labels = labels ) )
    # pass invalid tree
    expect_error( tree_reorder( 1:10, labels ) )
    # invalid labels (wrong length)
    expect_error( tree_reorder( tree, labels[-1] ) )
    # invalid labels (wrong content)
    labels_bad <- labels
    labels_bad[1] <- 'WRONG!'
    expect_error( tree_reorder( tree, labels_bad ) )
    
    # a successful run
    # this one ought to not change tree
    expect_silent(
        tree2 <- tree_reorder( tree, labels )
    )
    expect_equal( tree2, tree )

    # now see if we can unscramble tree!
    expect_silent(
        tree2 <- tree_reorder( tree_scrambled, labels )
    )
    expect_equal( tree2, tree )
    
})

test_that( "edges_to_tips works", {
    # error if tree is missing
    expect_error( edges_to_tips() )
    
    # start with a tree where we know what the answer should be
    tree_str <- '(S1:0.1,(S2:0.1,(S3:0.1,(S4:0.1,S5:0.1)N3:0.1)N2:0.1)N1:0.1)T;'
    tree <- ape::read.tree( text = tree_str )
    # stared at `tree$edge` structure and determined that this is what I expect
    edge_to_tips_exp <- list(
        1,
        2:5,
        2,
        3:5,
        3,
        4:5,
        4,
        5
    )

    # a successful run
    expect_silent(
        edge_to_tips_obs <- edges_to_tips( tree )
    )
    # compare outputs
    expect_equal(
        edge_to_tips_obs,
        edge_to_tips_exp
    )

    # randomly reorder edges, which messes with calculations in older versions
    indexes <- sample( ape::Nedge( tree ) )
    tree_rand <- tree
    tree_rand$edge <- tree_rand$edge[ indexes, ]
    # NOTE: no edge lengths here to worry about (they're all the same values anyway), plus the output ignores these values
    expect_silent(
        edge_to_tips_rand_obs <- edges_to_tips( tree_rand )
    )
    # compare outputs
    expect_equal(
        edge_to_tips_rand_obs,
        edge_to_tips_exp[ indexes ] # this reorders list!
    )
    
    # now try the same on a random tree
    # we won't know exactly what to expect but there are some patterns we do expect
    k_subpops <- 5
    tree <- ape::rtree( k_subpops )
    expect_silent(
        edge_to_tips <- edges_to_tips( tree )
    )
    expect_true( is.list( edge_to_tips ) )
    expect_equal( length( edge_to_tips ), nrow( tree$edge ) )
    n_tips <- length( tree$tip.label )
    # only tip nodes are present, and each is present at least once
    expect_true( all( sort( unique( unlist( edge_to_tips ) ) ) == 1 : n_tips ) )
})

test_that( "fit_tree_single works", {
    # create a random tree
    k_subpops <- 5 # 100
    tree <- ape::rtree( k_subpops )
    # add a non-trivial root edge
    tree$root.edge <- runif( 1 )
    # and form its true, linear scale coancestry matrix (also from ape, not `coanc_tree`)
    # however, root edge is ignored here so add it (as it is on the root, it applies to all elements)
    coancestry <- ape::vcv( tree ) + tree$root.edge

    # expect errors when arguments are missing
    expect_error( fit_tree_single( coancestry ) )
    expect_error( fit_tree_single( tree = tree ) )
    # invalid coancestry
    expect_error( fit_tree_single( 1:10, tree ) )
    # invalid tree
    expect_error( fit_tree_single( coancestry, 1:10 ) )
    
    # successful run
    expect_silent(
        tree_fit <- fit_tree_single( coancestry, tree )
    )
    # in this case the fit tree should be basically the same as the input because the coancestry is noiseless
    # check RSS first, should be near zero
    expect_equal( tree_fit$rss, 0 )
    # remove that element, the rest should be the same tree
    tree_fit$rss <- NULL
    expect_equal( tree_fit, tree )

    # make sure this works with scrambled tree edges
    indexes <- sample( ape::Nedge( tree ) )
    tree_rand <- tree
    tree_rand$edge <- tree_rand$edge[ indexes, ]
    # reorder edge lengths too for validation
    tree_rand$edge.length <- tree_rand$edge.length[ indexes ]
    expect_silent(
        tree_rand_fit <- fit_tree_single( coancestry, tree_rand )
    )
    # expect again a perfect fit!
    expect_equal( tree_rand_fit$rss, 0 )
    tree_rand_fit$rss <- NULL
    expect_equal( tree_rand_fit, tree_rand )

    # cause a warning when tree and coancestry labels don't agree, even though matrices are aligned by construction
    tree2 <- tree
    tree2$tip.label <- 1 : k_subpops
    expect_warning( fit_tree_single( coancestry, tree2 ) )
    # repeat but here coancestry has no names (but tree2 has indexes as names, so this combination works without warnings)
    coancestry2 <- coancestry
    dimnames( coancestry2 ) <- NULL
    expect_silent( fit_tree_single( coancestry2, tree2 ) )
    # and now combine nameless coancestry with trees with non-index names, which does cause a warning
    expect_warning( fit_tree_single( coancestry2, tree ) )

    # now a successful run under permutated inputs, which should be realigned successfully!
    # permutation
    indexes <- sample( k_subpops )
    # permute coancestry only (tree objects are harder to manipulate)
    coancestry2 <- coancestry[ indexes, indexes ]
    # shouldn't get warning here
    expect_silent(
        tree_fit <- fit_tree_single( coancestry2, tree )
    )
    # in this case we should find the exact answer again!
    expect_equal( tree_fit$rss, 0 )
    tree_fit$rss <- NULL
    expect_equal( tree_fit, tree )

    # now try a more aggressive example, where we try to fit a topology that is completely wrong
    # here we don't expect a good fit but there shouldn't be any errors at least
    # draw a second random tree for this
    tree2 <- ape::rtree( k_subpops )
    # for this test, we have to ensure the values to fit don't exceed one
    if ( max( coancestry ) > 1 )
        coancestry <- coancestry / max( coancestry )
    # successful run
    expect_silent(
        tree2_fit <- fit_tree_single( coancestry, tree2 )
    )
    # tree is good if it passes our tests
    # remember that a linear-scale tree always passes the coancestry tree tests
    expect_silent(
        validate_coanc_tree( tree2_fit )
    )
    expect_true( tree2_fit$rss >= 0 )
})

test_that( "fit_tree works", {
    # create a random tree
    k_subpops <- 5
    tree <- ape::rtree( k_subpops )
    # add a non-trivial root edge
    tree$root.edge <- runif( 1 )
    # and form its true coancestry matrix
    coancestry <- coanc_tree( tree )
    
    # expect errors when arguments are missing
    expect_error( fit_tree( ) )
    # invalid coancestry
    expect_error( fit_tree( 1:10 ) )
    
    # successful run
    expect_silent(
        tree_fit <- fit_tree( coancestry )
    )
    # in this case the fit tree should be basically the same as the input because the coancestry is noiseless
    # check RSS first, should be near zero
    expect_equal( tree_fit$rss, 0 )
    # remove that element, the rest should be the same tree
    tree_fit$rss <- NULL
    expect_equal( tree_fit, tree )

    # repeat with a permuted coancestry input!
    # permutation
    indexes <- sample( k_subpops )
    # permuted coancestry
    coancestry2 <- coancestry[ indexes, indexes ]
    # expect to recover true tree again!
    expect_silent(
        tree_fit <- fit_tree( coancestry2 )
    )
    expect_equal( tree_fit$rss, 0 )
    tree_fit$rss <- NULL
    expect_equal( tree_fit, tree )
    
})

test_that( "scale_tree works", {
    # create a random tree
    k_subpops <- 5
    tree <- ape::rtree( k_subpops )
    # factors smaller than 1 always work!
    factor <- 0.5
    # this is the answer we expect
    tree_exp <- tree
    tree_exp$edge.length <- tree_exp$edge.length * factor

    # cause errors on purpose
    # missing params
    expect_error( scale_tree( ) )
    expect_error( scale_tree( tree ) )
    expect_error( scale_tree( factor = factor ) )
    # pass a bad tree
    expect_error( scale_tree( 1:10, factor ) )
    # bad factors
    expect_error( scale_tree( tree, -factor ) )
    expect_error( scale_tree( tree, 'a' ) )
    expect_error( scale_tree( tree, c( factor, factor ) ) )
    # select a factor that will exceed prob scale
    expect_error( scale_tree( tree, 1.1 / max( tree$edge.length ) ) )
    
    # now the successful example
    expect_silent(
        tree_obs <- scale_tree( tree, factor )
    )
    # see if we got back the tree we expected
    expect_equal( tree_obs, tree_exp )

    # add root edge, an extreme but valid case
    tree$root.edge <- 0.9
    tree_exp$root.edge <- tree$root.edge * factor
    # repeat test
    expect_silent(
        tree_obs <- scale_tree( tree, factor )
    )
    expect_equal( tree_obs, tree_exp )
    
    # add additive edges separately for both input and expected output
    tree <- tree_additive( tree )
    tree_exp <- tree_additive( tree_exp )
    # repeat test
    expect_silent(
        tree_obs <- scale_tree( tree, factor )
    )
    expect_equal( tree_obs, tree_exp )
})

test_that( "diff_af works", {
    # an internal function for these tests only
    m <- 100
    p <- runif( m )
    F <- 0.3 # choose something large for a large effect
    # differentiate distribution!
    expect_silent(
        p2 <- diff_af( p, F )
    )
    # boring requirements
    expect_equal( length( p2 ), m )
    expect_true( !anyNA( p2 ) )
    expect_true( min( p2 ) >= 0 )
    expect_true( max( p2 ) <= 1 )
    # the real test is that variance has indeed increased
    V1 <- mean( ( p - 0.5 )^2 )
    V2 <- mean( ( p2 - 0.5 )^2 )
    expect_true( V2 >= V1 )
})

test_that( "undiff_af works", {
    # start by differentiating some data, to know we're starting from something reasonable
    m <- 100
    p <- runif( m ) # original
    V1 <- mean( ( p - 0.5 )^2 ) # a statistic we use for tests
    F <- 0.2 # choose something large for a large effect
    p2 <- diff_af( p, F ) # differentiated
    V2 <- mean( ( p2 - 0.5 )^2 ) # a statistic we use a lot for tests
    F2 <- 0.05 # something smaller for avoiding edge cases where uniform mixing variance is too high

    # cause errors on purpose
    expect_error( undiff_af( p2 ) )
    expect_error( undiff_af( F = F ) )
    expect_error( undiff_af( p2, F, distr = 'madeup-nomatch' ) )

    # uniform works well when F2 isn't too large
    expect_silent(
        obj <- undiff_af( p2, F2, distr = 'uniform' )
    )
    p3 <- obj$p
    # boring requirements
    expect_equal( length( p3 ), m )
    expect_true( !anyNA( p3 ) )
    expect_true( min( p3 ) >= 0 )
    expect_true( max( p3 ) <= 1 )
    # the real test is that variance has indeed decreased
    V3 <- mean( ( p3 - 0.5 )^2 )
    expect_true( V2 >= V3 )
    # test weights too
    w <- obj$w
    expect_equal( length( w ), 1 )
    expect_true( !is.na( w ) )
    expect_true( w >= 0 )
    expect_true( w <= 1 )

    # point always succeeds, so use original (large) F here
    expect_silent(
        obj <- undiff_af( p2, F, distr = 'point' )
    )
    p3 <- obj$p
    # boring requirements
    expect_equal( length( p3 ), m )
    expect_true( !anyNA( p3 ) )
    expect_true( min( p3 ) >= 0 )
    expect_true( max( p3 ) <= 1 )
    # the real test is that variance has indeed decreased
    V3 <- mean( ( p3 - 0.5 )^2 )
    expect_true( V2 >= V3 )
    # test weights too
    w <- obj$w
    expect_equal( length( w ), 1 )
    expect_true( !is.na( w ) )
    expect_true( w >= 0 )
    expect_true( w <= 1 )

    # beta also more likely than uniform to succeed (when more concentrated than unif), so use original (large) F here
    expect_silent(
        obj <- undiff_af( p2, F, distr = 'beta', alpha = 2 )
    )
    p3 <- obj$p
    # boring requirements
    expect_equal( length( p3 ), m )
    expect_true( !anyNA( p3 ) )
    expect_true( min( p3 ) >= 0 )
    expect_true( max( p3 ) <= 1 )
    # the real test is that variance has indeed decreased
    V3 <- mean( ( p3 - 0.5 )^2 )
    expect_true( V2 >= V3 )
    # test weights too
    w <- obj$w
    expect_equal( length( w ), 1 )
    expect_true( !is.na( w ) )
    expect_true( w >= 0 )
    expect_true( w <= 1 )

    # finally, "auto" hacks a beta, which always succeeds, so use original (large) F here
    expect_silent(
        obj <- undiff_af( p2, F, distr = 'auto' )
    )
    p3 <- obj$p
    # boring requirements
    expect_equal( length( p3 ), m )
    expect_true( !anyNA( p3 ) )
    expect_true( min( p3 ) >= 0 )
    expect_true( max( p3 ) <= 1 )
    # the real test is that variance has indeed decreased
    V3 <- mean( ( p3 - 0.5 )^2 )
    expect_true( V2 >= V3 )
    # test weights too
    w <- obj$w
    expect_equal( length( w ), 1 )
    expect_true( !is.na( w ) )
    expect_true( w >= 0 )
    expect_true( w <= 1 )

    # finally, "auto" hacks a beta, which always succeeds, so use original (large) F here
    # here a more extrmeme case: undifferentiate uniform inputs!
    expect_silent(
        obj <- undiff_af( p, F, distr = 'auto' )
    )
    p3 <- obj$p
    # boring requirements
    expect_equal( length( p3 ), m )
    expect_true( !anyNA( p3 ) )
    expect_true( min( p3 ) >= 0 )
    expect_true( max( p3 ) <= 1 )
    # the real test is that variance has indeed decreased
    V3 <- mean( ( p3 - 0.5 )^2 )
    expect_true( V1 >= V3 )
    # test weights too
    w <- obj$w
    expect_equal( length( w ), 1 )
    expect_true( !is.na( w ) )
    expect_true( w >= 0 )
    expect_true( w <= 1 )
})
