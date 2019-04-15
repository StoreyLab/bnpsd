#' A package for modeling and simulating an admixed population
#'
#' The underlying model is called the BN-PSD admixture model, which combines the Balding-Nichols (BN) allele frequency model for the intermediate subpopulations with the Pritchard-Stephens-Donnelly (PSD) model of individual-specific admixture proportions.
#' The BN-PSD model enables the simulation of complex population structures, ideal for illustrating challenges in kinship coefficient and \eqn{F_{ST}}{FST} estimation.
#' Note that simulated loci are drawn independently (in linkage equilibrium).
#'
#' @examples
#' # dimensions of data/model
#' m_loci <- 10 # number of loci
#' n_ind <- 5 # number of individuals
#' k_subpops <- 2 # number of intermediate subpops
#'
#' # define population structure
#' inbr_subpops <- c(0.1, 0.3) # FST values for k=2 subpopulations
#' sigma <- 1 # dispersion parameter of intermediate subpops
#' # admixture proportions from 1D geography
#' admix_proportions <- admix_prop_1d_linear(n_ind, k_subpops, sigma)
#'
#' # get pop structure parameters of the admixed individuals
#' # the coancestry matrix
#' coancestry <- coanc_admix(admix_proportions, inbr_subpops)
#' # FST of admixed individuals
#' Fst <- fst_admix(admix_proportions, inbr_subpops)
#'
#' # draw all random allele freqs and genotypes
#' out <- rbnpsd(admix_proportions, inbr_subpops, m_loci)
#' X <- out$X # genotypes
#' p_ind <- out$P # individual-specific AFs
#' p_subops <- out$B # independent subpops (intermediate) AFs
#' p_anc <- out$Pa # ancestral AFs
#'
#' # OR... draw each vector or matrix separately
#' # provided for additional flexibility
#' # ancestral AFs
#' p_anc <- draw_p_anc(m_loci)
#' # independent subpops (intermediate) AFs
#' p_subpops <- draw_p_subpops(p_anc, inbr_subpops)
#' # individual-specific AFs
#' p_ind <- make_p_ind_admix(p_subpops, admix_proportions)
#' # genotypes
#' X <- rgeno(p_ind)
#' 
#' @docType package
#' @name bnpsd
"_PACKAGE"
