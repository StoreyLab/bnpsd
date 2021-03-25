#' A package for modeling and simulating an admixed population
#'
#' The underlying model is called the BN-PSD admixture model, which combines the Balding-Nichols (BN) allele frequency model for the intermediate subpopulations with the Pritchard-Stephens-Donnelly (PSD) model of individual-specific admixture proportions.
#' The BN-PSD model enables the simulation of complex population structures, ideal for illustrating challenges in kinship coefficient and FST estimation.
#' Simulated loci are drawn independently (in linkage equilibrium).
#'
#' @examples
#' # dimensions of data/model
#' # number of loci
#' m_loci <- 10
#' # number of individuals
#' n_ind <- 5
#' # number of intermediate subpops
#' k_subpops <- 2
#'
#' # define population structure
#' # FST values for k = 2 subpopulations
#' inbr_subpops <- c( 0.1, 0.3 )
#' # admixture proportions from 1D geography
#' admix_proportions <- admix_prop_1d_linear( n_ind, k_subpops, sigma = 1 )
#' # also available:
#' # - admix_prop_1d_circular
#' # - admix_prop_indep_subpops
#'
#' # get pop structure parameters of the admixed individuals
#' # the coancestry matrix
#' coancestry <- coanc_admix( admix_proportions, inbr_subpops )
#' # FST of admixed individuals
#' Fst <- fst_admix( admix_proportions, inbr_subpops )
#'
#' # draw all random allele freqs and genotypes
#' out <- draw_all_admix( admix_proportions, inbr_subpops, m_loci )
#' # genotypes
#' X <- out$X
#' # ancestral allele frequencies (AFs)
#' p_anc <- out$p_anc
#'
#' # OR... draw each vector or matrix separately
#' # provided for additional flexibility
#' # ancestral AFs
#' p_anc <- draw_p_anc( m_loci )
#' # independent subpops (intermediate) AFs
#' p_subpops <- draw_p_subpops( p_anc, inbr_subpops )
#' # individual-specific AFs
#' p_ind <- make_p_ind_admix( p_subpops, admix_proportions )
#' # genotypes
#' X <- draw_genotypes_admix( p_ind )
#'
#' ### Examples with a tree for intermediate subpopulations
#'
#' # tree allows for correlated subpopulations
#' # (prev examples had independent subpopulations)
#' 
#' # best to start by specifying tree in Newick string format
#' tree_str <- '(S1:0.1,(S2:0.1,S3:0.1)N1:0.1)T;'
#' # and turn it into `phylo` object using the `ape` package
#' library(ape)
#' tree_subpops <- read.tree( text = tree_str )
#' # true coancestry matrix corresponding to this tree
#' coanc_subpops <- coanc_tree( tree_subpops )
#' 
#' # admixture proportions from 1D geography
#' # (constructed again but for k=3 tree)
#' k_subpops <- nrow( coanc_subpops )
#' admix_proportions <- admix_prop_1d_linear( n_ind, k_subpops, sigma = 0.5 )
#' 
#' # get pop structure parameters of the admixed individuals
#' # the coancestry matrix
#' coancestry <- coanc_admix( admix_proportions, coanc_subpops )
#' # FST of admixed individuals
#' Fst <- fst_admix( admix_proportions, coanc_subpops )
#'
#' # draw all random allele freqs and genotypes, tree version
#' out <- draw_all_admix( admix_proportions, tree_subpops = tree_subpops, m_loci = m_loci )
#' # genotypes
#' X <- out$X
#' # ancestral allele frequencies (AFs)
#' p_anc <- out$p_anc
#'
#' # OR... draw tree subpops (intermediate) AFs separately
#' p_subpops_tree <- draw_p_subpops_tree( p_anc, tree_subpops )
#' 
#' @docType package
#' @name bnpsd
"_PACKAGE"
