#' Reorder tree tips to best match a desired order
#'
#' This functions reorganizes the tree structure so that its tips appear in a desired order if possible, or in a reasonably close order when an exact solution is impossible.
#' This tip order in the output tree is the same in both the tip labels vector (`tree$tip.label`) and edge matrix (`tree$edge`), ensured by using [tree_reindex_tips()] internally.
#'
#' This function has the same goal as [ape::rotateConstr()], which implements a different heuristic algorithm that did not perform well in our experience.
#'
#' @param tree A `phylo` object from the `ape` package (see [ape::read.tree()]).
#' Works with standard `phylo` objects, and also with our extended trees (in that additive edges `tree$edge.length.add` are recalculated after reordering if they were present).
#' @param labels A character vector with all tip labels in the desired order.
#' Must contain each tip label in `tree` exactly once.
#'
#' @return The modified `tree` (`phylo` object) with reordered edges and tips.
#'
#' @examples
#' # create a random tree
#' library(ape)
#' k <- 5
#' tree <- rtree( k )
#' # let's set the current labels as the desired order
#' labels <- tree$tip.label
#'
#' # now let's scramble the edges on purpose
#' # to create an example where reordering is needed
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
#' # the order of the tip labels in the vector and on the plot disagree with each other:
#' tree_rand$tip.label
#' plot( tree_rand )
#'
#' # now reorder tree object so tips are in the desired order:
#' tree_rand <- tree_reorder( tree_rand, labels )
#' # now tip labels vector and plot should agree in order:
#' # (the original tree was recovered!)
#' tree_rand$tip.label
#' plot( tree_rand )
#'
#' # order the tree in a different way than the original order
#' labels <- paste0( 't', 1 : k )
#' # in this case, it's often impossible to get a perfect output order
#' # (because the tree structure constrains the possible plot orders),
#' # but this function tries its best to get close to the desired order
#' tree2 <- tree_reorder( tree, labels )
#' plot(tree2)
#' 
#' @seealso
#' [tree_reindex_tips()] to reorder tips in the labels vector to match the edge matrix order, which ensures agreement in plots (assuming plot show desired order already).
#'
#' @export
tree_reorder <- function( tree, labels ) {
    if ( missing( tree ) )
        stop( '`tree` is required!' )
    if ( missing( labels ) )
        stop( '`labels` is required!' )

    # validate tree
    validate_coanc_tree( tree )
    
    # number of edges
    n_edges <- ape::Nedge( tree )
    # number of tips
    n_tips <- ape::Ntip( tree )
    # number of internal nodes
    n_nodes <- ape::Nnode( tree )

    # validate labels as being same set as tree$tip.label
    if ( length( labels ) != n_tips )
        stop( '`labels` must have length equal to number of tips in `tree`!' )
    if ( !all( labels %in% tree$tip.label ) )
        stop( '`labels` must equal in content the tip labels in `tree`!' )
    
    # code assumes names (tip labels/indexes) are in order of appearance on tree (edge matrix)
    # ensure that this is correct
    tree <- tree_reindex_tips( tree )
    
    # correspondence between desired and existing labels
    indexes <- match( tree$tip.label, labels )
    
    # reordering edges by postorder helps this setup a lot!
    # (this doesn't change the actual plot, just the internal organization)
    # (should not change the tip order!!!)
    tree <- ape::reorder.phylo( tree, order = 'postorder' )
    
    # expand indexes to hold values for internal nodes too (only average matters for these)
    indexes <- c( indexes, rep.int( NA, n_nodes ) )
    # actually will keep numerators and denominators separate since branches can be unbalanced, this way the numbers make sense
    indexes_n <- c( rep.int( 1, n_tips ), rep.int( NA, n_nodes ) )
    
    # navigate edges
    j_node <- 1 # current/last node, to know when to skip stuff that has already been processed
    for ( e in 1 : n_edges ) {
        # get parent node for this edge
        # decide whether to process or not
        j_parent <- tree$edge[ e, 1 ]
        if ( j_parent == j_node ) {
            next
        } else
            j_node <- j_parent

        # edges are continuous by parent node here (because we reordered by "postorder")
        # so we need to find extreme of range
        e2 <- e + 1 # next edge
        while ( e2 <= n_edges && tree$edge[ e2, 1 ] == j_node )
            e2 <- e2 + 1
        # so when we stop, e2 is the next edge that belongs to another case
        # decrement so it's an inclusive range
        e2 <- e2 - 1
        # there's always at least two edges per node in proper data, but let's check
        if ( e2 == e )
            stop( 'Internal node ', j_node, ' only had one edge!' )

        # this should work for 2 or more edges!
        # collect desired indexes for each child node now
        j_children <- tree$edge[ e : e2 , 2 ]
        # average indexes per children
        indexes_children <- indexes[ j_children ] / indexes_n[ j_children ]
        # in "postorder", all children should have values, but let's check
        # this is an unexpected error, not sure it's worthy of an error message
        stopifnot( !anyNA( indexes_children ) )
        # this is the desired order (since we want indexes to be increasing)
        order_children <- order( indexes_children )
        # reordering is needed if this isn't 1:n_children
        if ( !all( order_children == 1 : length( order_children ) ) ) {
            # change order in edge matrix
            # NOTE: first column was chosen so all elements were `j_node` (in these rows), so that one doesn't need reordering
            tree$edge[ e : e2 , 2 ] <- j_children[ order_children ]
            # also need to change order in edge lengths vector
            tree$edge.length[ e : e2 ] <- ( tree$edge.length[ e : e2 ] )[ order_children ]
        }
        # done reordering!
        # for next round, so this parent node has an average position, save it as the sum of the children values (which turns into an average upon dividing the two numbers)
        indexes[ j_parent ] <- sum( indexes[ j_children ] )
        indexes_n[ j_parent ] <- sum( indexes_n[ j_children ] )
    }

    # now we should reindex tips in order of appearance!
    # if there were additive edges, they are recalculated now
    tree <- tree_reindex_tips( tree )
    
    # return edited tree!
    return( tree )
}
