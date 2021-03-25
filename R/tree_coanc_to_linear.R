# let's edit the tree so the edges are scaled to be in the linear (additive) scale
# if the edge has coancestry t1 and the parent has edge t2 (already in the total/linear scale, i.e. condensing earlier edges iteratively as here), then the total branch length is
# 1 - ( 1 - t1 ) * ( 1 - t2 )
# = t1 + t2 - t1 * t2
# = t2 + t1 * ( 1 - t2 )
# so the edge beyond the t2 contribution, t1 * ( 1 - t2 ), is not just t1 but it is scaled down by ( 1 - t2 ).
# this code performs that scaling

# transforms the edges to they are additive in the linear scale
# most useful for predicting the covariance matrix
# lets "extend" `phylo` object by adding more elements that are useful for our operations
tree_coanc_to_linear <- function( tree, name = 'tree' ) {
    if ( missing( tree ) )
        stop( '`tree` is required!' )
    # run through validator for additional checks
    validate_coanc_tree( tree, name = name )
    
    # number of edges
    n_edges <- nrow( tree$edge )
    # determine root node
    j_root <- tree$edge[ 1, 1 ]
    # keep track of additive edges for each child from the root as it becomes a parent
    # a list in parallel to the normal list of edges for each node to its parent only
    edge_length_add <- vector( 'numeric', n_edges )
    # a list with the same info but indexed by nodes instead, easier to retrieve in this case
    # initialize this one to NA to catch issues
    node_edge_from_root <- rep.int( NA, max( tree$edge ) )
    # only root node is pre-set to zero
    node_edge_from_root[ j_root ] <- 0
    # navigate all edges
    for ( e in 1 : n_edges ) {
        # get parent and child nodes for this edge
        j_parent <- tree$edge[ e, 1 ]
        j_child <- tree$edge[ e, 2 ]
        # edge length is FST we need
        fst_child_from_parent <- tree$edge.length[ e ]
        # get edge length of parent from root
        fst_parent_from_root <- node_edge_from_root[ j_parent ]
        # check
        if ( is.na( fst_parent_from_root ) )
            stop( 'Node index ', j_parent, ' was not assigned coancestry from root! (unexpected)' )
        # shrink edge now, store in new edge vector
        edge_length_add[ e ] <- fst_child_from_parent * ( 1 - fst_parent_from_root )
        # and store total edge for child (as it may become a parent in later iterations)
        node_edge_from_root[ j_child ] <- fst_parent_from_root + edge_length_add[ e ]
    }
    # we could store `node_edge_from_root`, but we won't for now
    # we include additive edges vector, without erasing original for now
    tree$edge.length.add <- edge_length_add
    # return this extended tree
    return( tree )
}
