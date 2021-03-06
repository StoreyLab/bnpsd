#' Construct admixture proportion matrix for independent subpopulations
#'
#' This function constructs an admixture proportion matrix where every individuals is actually unadmixed (draws its full ancestry from a single intermediate subpopulation).
#' The inputs are the vector of subpopulation labels `labs` for every individual (length `n`), and
#' the length-`k` vector of unique subpopulations `subpops` in the desired order.
#' If `subpops` is missing, the sorted unique subpopulations observed in `labs` is used.
#' This function returns the admixture proportion matrix, for each individual 1 for the column corresponding to its subpopulation, 0 otherwise.
#'
#' @param labs Length-`n` vector of subpopulation labels
#' @param subpops Optional length-`k` vector of unique subpopulations in desired order.
#' Stops if `subpops` does not contain all unique labels in `labs` (no error if `subpops` contains additional labels).
#'
#' @return The `n`-by-`k` admixture proportion matrix.
#' The unique subpopulation labels are given in the column names.
#'
#' @examples
#' # vector of subpopulation memberships
#' labs <- c(1, 1, 1, 2, 2, 3, 1)
#' # admixture matrix with subpopulations (along columns) sorted
#' admix_proportions <- admix_prop_indep_subpops(labs)
#'
#' # declare subpopulations in custom order
#' subpops <- c(3, 1, 2)
#' # columns will be reordered to match subpops as provided
#' admix_proportions <- admix_prop_indep_subpops(labs, subpops)
#'
#' # declare subpopulations with unobserved labels
#' subpops <- 1:5
#' # note columns 4 and 5 will be false for all individuals
#' admix_proportions <- admix_prop_indep_subpops(labs, subpops)
#'
#' @export
admix_prop_indep_subpops <- function(
                                     labs,
                                     subpops = sort( unique( labs ) )
                                     ) {
    # this argument is mandatory
    if ( missing( labs ) )
        stop( '`labs` is required!' )
    
    # check if subpops is provided but was inconsistent with labs
    if ( !all( labs %in% subpops ) )
        stop("provided `subpops` does not contain all labels in `labs`!")

    # get dimensions from inputs
    n_ind <- length(labs) # number of individuals
    k_subpops <- length(subpops) # number of subpopulations
    
    admix_proportions <- matrix(0, nrow = n_ind, ncol = k_subpops) # initialize the desired admixture matrix
    colnames(admix_proportions) <- subpops # return subpopulation names in admixture matrix column names
    
    # fill in admix_proportions, navigating each subpopulation
    for (i in 1 : k_subpops) {
        admix_proportions[, i] <- labs == subpops[i] # indicator for whether the labels were of subpopulation i or not
    }
    
    return( admix_proportions )
}
