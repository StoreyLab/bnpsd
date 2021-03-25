#' Calculate coancestry matrix corresponding to a tree
#'
#' This function calculates the coancestry matrix of the subpopulations which are tip nodes in the input tree.
#' The edges of this tree are coancestry values relative to the parent node of each child node.
#' 
#' The calculation takes into account that total coancestries are non-linear combinations of the per-edge coancestries.
#' For example, if the root node is `A`, and subpopulation `C` is connected to `A` only through an internal node `B`, then its total self-coancestry `coanc_A_C` relative to `A` is given by `coanc_A_B` (the coancestry between `A` and `B`) and `coanc_B_C` (the coancestry between `B` and `C`) by
#' `coanc_A_C = coanc_A_B + coanc_B_C * ( 1 - coanc_A_B )`.
#' This transformation ensures that the total coancestry is a probability that does not exceed one even when the per-edge coancestries sum to a value greater than one.
#' This transformation is iterated for all nodes in the tree.
#' The self-coancestry (diagonal values) are the total coancestries of the tip nodes.
#' The coancestry between different subpopulations is the total coancestry of their last common ancestor node.
#'
#' @param tree The coancestry tree relating the `k_subpops` subpopulations.
#' Must be a `phylo` object from the `ape` package (see [ape::read.tree()]).
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
#' @export
coanc_tree <- function( tree ) {
    # check mandatory data
    if ( missing( tree ) )
        stop( '`tree` is required!' )
    
    # calculate branch lengths on additive scale if needed
    # if the data we want is already there, use it
    if ( is.null( tree$edge.length.add ) ) {
        # `tree_coanc_to_linear` runs `validate_coanc_tree` (so don't repeat check)
        tree <- tree_coanc_to_linear( tree )
    } else {
        # run through validator for additional checks
        # (needed since we didn't run `tree_coanc_to_linear`)
        validate_coanc_tree( tree )
    }

    # replace branches here so the `ape` code just works
    tree$edge.length <- tree$edge.length.add

    # calculate coancestry matrix from the tree with the modified edges
    coanc <- ape::vcv( tree )

    # return this coancestry matrix
    return( coanc )
}
