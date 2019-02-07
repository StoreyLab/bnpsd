#' Simulate random allele frequencies and genotypes from the BN-PSD admixture model
#'
#' This function returns simulated ancestral, intermediate, and individual-specific allele frequencies and genotypes given the admixture structure, as determined by the admixture proportions and the vector of intermediate subpopulation \eqn{F_{ST}}{FST} values.
#' The function is a wrapper around \code{\link{rpanc}}, \code{\link{rpint}}, \code{\link{rpiaf}}, and \code{\link{rgeno}}.
#' Below \eqn{m} is the number of loci, \eqn{n} is the number of individuals, and \eqn{k} is the number of intermediate subpopulations.
#'
#' @param Q The \eqn{n \times k}{n-by-k} matrix of admixture proportions
#' @param F The length-\eqn{k} vector of intermediate subpopulation \eqn{F_{ST}}{FST} values
#' @param m The number of loci to draw
#' @param wantX If TRUE (default), calculates and includes the random genotype matrix in the return list
#' @param wantP If TRUE (default), includes the random IAF matrix in the return list
#' @param wantB If TRUE (default), includes the random intermediate pop allele freq matrix in the return list
#' @param wantPa If TRUE (default), includes the random ancestral allele freq matrix in the return list
#' @param lowMem If TRUE, uses a low-memory algorithm to raw genotypes without storing or returning the corresponding IAF matrix.
#' @param verbose If TRUE, prints messages for every stage in the algorithm
#' @param noFixed If TRUE, returned matrix will not include any fixed loci (loci that happened to be fixed are drawn again, starting from the ancestral allele frequency, and checked iteratively until no fixed loci remain, so that the final number of loci is exactly \eqn{m} as specified).
#'
#' @return A named list that includes the following random matrices: X=genotypes, P=IAFs, B=intermediate pop allele freqs, Pa=vector of ancestral allele frequencies.  Items may be omitted depending on the values of wantX, wantP, wantB, or wantPa above.
#'
#' @examples
#' # dimensions
#' m <- 10 # number of loci
#' n <- 5 # number of individuals
#' k <- 2 # number of intermediate subpops
#'
#' # define population structure
#' F <- c(0.1, 0.3) # FST values for k=2 subpopulations
#' sigma <- 1 # dispersion parameter of intermediate subpops
#' Q <- q1d(n, k, sigma) # admixture proportions from 1D geography
#'
#' # draw all random allele freqs and genotypes
#' out <- rbnpsd(Q, F, m)
#' X <- out$X # genotypes
#' P <- out$P # IAFs
#' B <- out$B # Intermediate AFs
#' pAnc <- out$Pa # Ancestral AFs
#' 
#' @export
rbnpsd <- function(Q, F, m, wantX=TRUE, wantP=TRUE, wantB=TRUE, wantPa=TRUE, lowMem=FALSE, verbose=FALSE, noFixed=FALSE) {
    ## simulation of Pritchard-Stephens-Donnelly "admixture" model for allele frequencies/genotypes.
    
    ## always ask for an admixture matrix (no default)
    if (missing(Q)) stop('Fatal: must provide Q (admixture proportion matrix)')
    if (missing(F)) stop('Fatal: must provide F (intermediate subpopulation Fst vector)')
    if (missing(m)) stop('Fatal: must provide m (number of loci to draw)')
    
    ## get dimensions, test coherence
    n <- nrow(Q)
    k <- ncol(Q)
    k2 <- length(F) # don't allow scalar F here!
    if (k != k2) stop('Fatal: Q and F are not compatible: ncol(Q) == ', k, ' != ', k2, ' == length(F)')
    
    ## generate the random ancestral allele frequencies, in usual range and with minimum threshold for simplicity
    ## don't do this if a Pa was already provided (a way to provide arbitrary distributions)
    if (verbose) message('rbnpsd: drawing Pa')
    Pa <- rpanc(m)
    
    ## draw intermediate allele frequencies from Balding-Nichols
    if (verbose) message('rbnpsd: drawing B')
    B <- rpint(Pa, F)
    
    if (lowMem) {
        ## draw genotypes!
        if (wantX) {
            if (verbose) message('rbnpsd: drawing X')
            X <- rgeno(B, Q, lowMem=lowMem)
        }
    } else {
        if (wantP || wantX) {
            if (verbose) message('rbnpsd: drawing P')
            P <- rpiaf(B, Q)
        }
        
        ## draw genotypes from P
        if (wantX) {
            if (verbose) message('rbnpsd: drawing X')
            X <- rgeno(P)
        }
    }
    
    if (noFixed && wantX) {
        # check for fixed loci, draw them again if needed
        # Note this only applies to wantX==TRUE, since P and B are continuous and therefore practically never truly fixed
        fixedLoci <- fixed_loci(X) # boolean vector identifies fixed loci
        m2 <- sum(fixedLoci) # number of cases
        if (m2 > 0) {
            # call self with desired number of loci, all the same parameters otherwise
            # note that since this is also called asking for no fixed loci, it will work recursively within itself and return when all loci desired are not fixed!
            obj <- rbnpsd(Q, F, m2, wantX=wantX, wantP=wantP, wantB=wantB, wantPa=wantPa, lowMem=lowMem, verbose=verbose, noFixed=noFixed)
            # overwrite fixed loci with redrawn polymorphic data
            X[fixedLoci,] <- obj$X # guaranteed to be there
            if (wantP && !lowMem) P[fixedLoci,] <- obj$P
            if (wantB) B[fixedLoci,] <- obj$B
            if (wantPa) Pa[fixedLoci] <- obj$Pa
        }
    }
    
    ## now prepare output list
    out <- list()
    if (wantX) out$X <- X
    if (wantP && !lowMem) out$P <- P # don't have when low-mem!
    if (wantB) out$B <- B
    if (wantPa) out$Pa <- Pa
    return(out)
}

