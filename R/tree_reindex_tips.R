#' Reindex tree tips in order of appearance in edges
#'
#' The `phylo` object from the `ape` package (see [ape::read.tree()]) permits mismatches between the order of tips as present in the tip labels vector (`tree$tip.label`) and in the edge matrix (`tree$edge`), which can cause mismatches in plots and other settings.
#' This function reindexes edges and tips so that they agree in tip order.
#'
#' Order mismatches between tip labels vector and edge matrix can lead to disagreeing order downstream, in variables that order tips as in the labels vector (such as the coancestry matrices output by our [coanc_tree()]) and tree plots (whose order is guided by the edge matrix, particularly when the edge matrix is ordered by clade as in [ape::reorder.phylo()]).
#'
#' This function first reorders the edges of the input tree using [ape::reorder.phylo()] with default option `order = 'cladewise'`, which will list edges and tip nodes in plotting order.
#' Then tips are indexed so that the first tip to appear has index 1 in the edge matrix (and consequently appears first in the tip labels vector), the second has index 2, and so on.
#' Thus, the output tree has both edges and tips reordered, to have a consistent tip order and best experience visualizing trees and their coancestry matrices.
#'
#' @param tree A `phylo` object from the `ape` package (see [ape::read.tree()]).
#' Works with standard `phylo` objects, and also with our extended trees (in that additive edges `tree$edge.length.add` are recalculated after reordering if they were present).
#'
#' @return The modified `tree` (`phylo` object) with reordered edges and tips.
#'
#' @examples
#' # create a random tree
#' library(ape)
#' k <- 5
#' tree <- rtree( k )
#'
#' # trees constructed this way are already ordered as desired, so this function doesn't change it:
#' tree2 <- tree_reindex_tips( tree )
#' # tree2 equals tree!
#'
#' # let's scramble the edges on purpose
#' # to create an example where reindexing is needed
#' 
#' tree_rand <- tree
#' # new order of edges
#' indexes <- sample( Nedge( tree_rand ) )
#' # reorder all edge values
#' tree_rand$edge <- tree_rand$edge[ indexes, ]
#' tree_rand$edge.length <- tree_rand$edge.length[ indexes ]
#' # now let's reorder edges slightly so tree is more reasonable-looking
#' # (otherwise plot looks tangled)
#' tree_rand <- reorder( tree_rand, order = 'postorder' )
#' # you can now see that, unless permutation got lucky,
#' # the order of the tip labels in the vector and on the plot disagree:
#' tree_rand$tip.label
#' plot( tree_rand )
#'
#' # now reorder tips to match plotting order (as defined in the `edge` matrix)
#' tree_rand <- tree_reindex_tips( tree_rand )
#' # now tip labels vector and plot should agree in order:
#' # (plot was not changed)
#' tree_rand$tip.label
#' plot( tree_rand )
#'
#' @seealso
#' [tree_reorder()] for reordering tree structure so that tips appear in a particular desired order.
#' 
#' @export
tree_reindex_tips <- function( tree ) {
    if ( missing( tree ) )
        stop( '`tree` is required!' )

    # NOTE: we skip `validate_coanc_tree` because `tree_additive` will do that anyway, and in it we handle the case that additive edges will be wrong and have to be recalculated (otherwise the validation step trows a fatal error)
    
    # if we don't reorder by edges clade, this doesn't work
    tree <- ape::reorder.phylo( tree )
    
    # from here on, edges are not reordered again, only tips are reordered/reindexed

    # recalculate additive edges if they were present (these weren't reordered correctly in the last step)
    if ( !is.null( tree$edge.length.add ) )
        tree <- tree_additive( tree, force = TRUE )
    
    # get number of tips
    n_tips <- ape::Ntip( tree )
    # extract second column of edge matrix
    # all tips appear here, and we are interested in this order
    indexes <- tree$edge[, 2]
    # exclude non-tip indexes
    indexes <- indexes[ indexes <= n_tips ]
    
    # reorder only if needed
    if ( !all( indexes == 1 : n_tips ) ) {
        
        # reordering tips is easy!
        tree$tip.label <- tree$tip.label[ indexes ]

        # need to walk through edge matrix, can't see an easier way
        # number of edges
        n_edges <- ape::Nedge( tree )
        # navigate edges
        for ( e in 1 : n_edges ) {
            # tips appear only in second column
            j_child <- tree$edge[ e , 2 ]
            if ( j_child <= n_tips ) {
                # we've got a tip!
                # overwrite with desired index
                tree$edge[ e , 2 ] <- which( indexes == j_child )
            }
        }
        
    }
    return( tree )
}

# for internal tests of `tree_reindex_tips` above
# for correspondence between tree and coancestry, need to make sure tips are listed in the same order in $edge and in $tip.label (implicitly by index)
is_ordered_tips_edges <- function( tree ) {
    if ( missing( tree ) )
        stop( '`tree` is required!' )

    # if we don't reorder edges by clade, this doesn't work
    tree <- ape::reorder.phylo( tree )
    
    # get number of tips
    n_tips <- ape::Ntip( tree )
    # extract second column of edge matrix
    # all tips appear here, and we are interested in this order
    indexes <- tree$edge[, 2]
    # exclude non-tip indexes
    indexes <- indexes[ indexes <= n_tips ]
    # this is the test
    test <- all( indexes == 1 : n_tips )
    return( test )
}

