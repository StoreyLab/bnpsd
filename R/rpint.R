#' Draw intermediate subpopulation allele frequencies
#'
#' Intermediate subpopulation allele frequencies \eqn{p_i^{S_u}} for subpopulation \eqn{S_u} at locus \eqn{i} are drawn from the Balding-Nichols distribution with ancestral allele frequency \eqn{p_i^T} and \eqn{F_{ST}}{FST} parameter \eqn{f^T_{S_u}} as
#' \deqn{p_i^{S_u} \sim \mbox{Beta}(\nu_u p_i^T, \nu_u (1-p_i^T)),}
#' where \eqn{\nu_u = 1/f^T_{S_u} - 1}.
#' Below \eqn{m} is the number of loci and \eqn{k} is the number of subpopulations.
#'
#' @param pAnc The length-\eqn{m} vector of ancestral allele frequencies per locus
#' @param F The length-\eqn{k} vector of subpopulation \eqn{F_{ST}}{FST} values
#'
#' @return The \eqn{m \times k}{m-by-k} matrix of intermediate subpopulation allele frequencies
#'
#' @examples
#' m <- 10 # number of loci
#' pAnc <- rpanc(m) # random vector of ancestral allele frequencies
#' F <- c(0.1, 0.3) # FST values for two subpopulations
#' B <- rpint(pAnc, F) # matrix of intermediate subpop allele freqs
#'
#' @export
rpint <- function(pAnc, F) {
    ## generate at each locus random intermediate allele frequencies from BN model
    ## inputs are
    ## - pAnc: ancestral allele frequencies (must be defined here!). Vector of length m (number of loci)
    ## - F: Fst vector (different for each intermediate population)

    ## basic param checking
    if (missing(pAnc)) stop('Fatal in rpint: ancestral allele frequencies missing!')
    if (missing(F)) stop('Fatal in rpint: F=Fst scalar or vector missing!')
    ## number of loci and subpopulations to simulate
    m <- length(pAnc)
    k <- length(F)
    
    ## let's translate parameters for Balding-Nichols case
    nu <- 1/F-1 # nu is a vector or a scalar, same as F (whatever that is)
    pAncM <- 1-pAnc # precompute vector of "minus" pAnc's shared by all subpopulations below
    ## vectorization makes a lot of sense for each population... (doing all SNPs together)
    B <- matrix(nrow=m, ncol=k) # matrix of intermediate allele frequencies we want...
    for (j in 1:k) {
        nuj <- nu[j]
        if (is.infinite(nuj)) {
            ## there is no drift from pAnc in this special case (F==0)
            ## (coded separately because rbeta incorrectly returns 0.5 instead)
            B[,j] <- pAnc
        } else {
            B[,j] <- stats::rbeta(m, nuj*pAnc, nuj*pAncM) # draw all SNPs for this population, store immediately
        }
    }
    
    B # the only thing we want out of this
}
