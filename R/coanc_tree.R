#' Calculate coancestry matrix corresponding to a tree
#'
#' This function calculates the coancestry matrix of the subpopulations which are tip nodes in the input tree.
#' The edges of this tree are coancestry values relative to the parent node of each child node.
#' 
#' The calculation takes into account that total coancestries are non-linear functions of the per-edge coancestries.
#' Interestingly, the calculation can be simplified by a simple transformation performed by [tree_additive()], see that for more information.
#' The self-coancestry (diagonal values) are the total coancestries of the tip nodes.
#' The coancestry between different subpopulations is the total coancestry of their last common ancestor node.
#'
#' @param tree The coancestry tree relating the `k_subpops` subpopulations.
#' Must be a `phylo` object from the `ape` package (see [ape::read.tree()]).
#' This tree may have a valid root edge (non-NULL `tree$root.edge` between 0 and 1), which is incorporated in the output calculations.
#'
#' @return The `k_subpops`-by-`k_subpops` coancestry matrix.
#'
#' @examples
#' # for simulating a tree with `rtree`
#' library(ape)
#' 
#' # a typical, non-trivial example
#' # number of tip subpopulations
#' k_subpops <- 3
#' # simulate a random tree
#' tree_subpops <- rtree( k_subpops )
#' # coancestry matrix of subpopulations
#' coancestry <- coanc_tree( tree_subpops )
#'
#' @seealso
#' [tree_additive()] for calculating the additive edges.
#' This function is called internally by `coanc_tree` but the additive edges are not returned here, so call [tree_additive()] if you desired them.
#' 
#' @export
coanc_tree <- function( tree ) {
    # check mandatory data
    if ( missing( tree ) )
        stop( '`tree` is required!' )
    
    # calculate branch lengths on additive scale if needed
    # if the data we want is already there, use it
    if ( is.null( tree$edge.length.add ) ) {
        # `tree_additive` runs `validate_coanc_tree` (so don't repeat check)
        tree <- tree_additive( tree )
    } else {
        # run through validator for additional checks
        # (needed since we didn't run `tree_additive`)
        validate_coanc_tree( tree )
    }

    # replace branches here so the `ape` code just works
    tree$edge.length <- tree$edge.length.add

    # calculate coancestry matrix from the tree with the modified edges
    coanc <- ape::vcv( tree )

    # add root edge to all elements, if present
    # NOTE: root edge is the same in both additive and probabilistic scales, so there isn't a different version to consider here
    if ( !is.null( tree$root.edge ) )
        coanc <- coanc + tree$root.edge
    
    # return this coancestry matrix
    return( coanc )
}
