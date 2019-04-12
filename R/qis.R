#' Construct admixture proportion matrix for independent subpopulations
#'
#' This function constructs an admixture proportion matrix where every individuals is actually unadmixed (draws its full ancestry from a single intermediate subpopulation).
#' The inputs are the vector of subpopulation labels \code{labs} for every individual (length \eqn{n}), and
#' the length-\eqn{k} vector of unique subpopulations \code{subpops} in the desired order.
#' If \code{subpops} is missing, the sorted unique subpopulations observed in \code{labs} is used.
#' This function returns the admixture proportion matrix \eqn{Q} with individuals along the rows and subpopulations along the columns, marking for each individual TRUE for the column corresponding to its subpopulation, FALSE otherwise.
#' Treating the entries of \eqn{Q} as numerical, every individual has an admixture proportion of 1 for its subpopulation and 0 for all other subpopulations.
#'
#' @param labs Length-\eqn{n} vector of subpopulation labels
#' @param subpops Optional length-\eqn{k} vector of unique subpopulations in desired order.
#' Stops if \code{subpops} does not contain all unique labels in \code{labs} (no error if \code{subpops} contains additional labels).
#'
#' @return The \eqn{n \times k}{n-by-k} admixture proportion matrix \eqn{Q}.
#' The unique subpopulation labels are given in \code{colnames(Q)}.
#'
#' @examples
#' # vector of subpopulation memberships
#' labs <- c(1, 1, 1, 2, 2, 3, 1)
#' # admixture matrix with subpopulations (along columns) sorted
#' Q <- qis(labs)
#'
#' # declare subpopulations in custom order
#' subpops <- c(3, 1, 2)
#' # columns will be reordered to match subpops as provided
#' Q <- qis(labs, subpops)
#'
#' # declare subpopulations with unobserved labels
#' subpops <- 1:5
#' # note columns 4 and 5 will be false for all individuals
#' Q <- qis(labs, subpops)
#'
#' @export
qis <- function(labs, subpops) {
    if (missing(subpops)) {
        # get unique list of subpopulations, sorted for maximum stability
        subpops <- sort(unique(labs))
    } else {
        if (!all(labs %in% subpops))
            stop("provided 'subpops' does not contain all labels in 'labs'!")
    }
    n <- length(labs) # number of individuals
    k <- length(subpops) # number of subpopulations
    Q <- matrix(0, nrow = n, ncol = k) # initialize the desired admixture matrix
    colnames(Q) <- subpops # return subpopulation names in admixture matrix column names
    # fill in Q, navigating each subpopulation
    for (i in 1:k) {
        Q[,i] <- labs == subpops[i] # booleas for whether the labels were of subpopulation i or not
    }
    Q # return!
}
