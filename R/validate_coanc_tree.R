validate_coanc_tree <- function( tree, name = 'tree' ) {
    if ( missing( tree ) )
        stop( '`', name, '` is required!' )
    
    # for now support only class "ape::phylo"
    # (can have additional classes)
    if ( !( 'phylo' %in% class( tree ) ) )
        stop( '`', name, '` must be of class "phylo"!' )

    # use ape to test if the tree is rooted (it must be)
    if ( ! ape::is.rooted( tree ) )
        stop( '`', name, '` must be rooted!' )
    
    # check that the four named elements are there
    # ignore additional names, but these 4 can't be missing
    names_exp <- c('edge', 'tip.label', 'Nnode', 'edge.length')
    names_obs <- names( tree )
    names_mis <- setdiff( names_exp, names_obs )
    if ( length( names_mis ) > 0 )
        stop( '`', name, '` is missing phylo class-required named elements: ', toString( names_mis ) )
    
    # in addition to the normal requirements of this class, we require each edge to be a coancestry value (an input to BN), so it must be a probability in the correct range
    if ( any( tree$edge.length > 1 ) )
        stop( '`', name, '` edge lengths must all be smaller than 1!  Bad values observed: ', toString( tree$edge.length[ tree$edge.length > 1 ] ) )
    # might as well check for negatives
    if ( any( tree$edge.length < 0 ) )
        stop( '`', name, '` edge lengths must all be non-negative!  Bad values observed: ', toString( tree$edge.length[ tree$edge.length < 0 ] ) )

    # if root edge is present, it must also satisfy bounds
    if ( !is.null( tree$root.edge ) ) {
        if ( tree$root.edge > 1 )
            stop( '`', name, '` root edge length must be smaller than 1!  Observed: ', tree$root.edge )
        if ( tree$root.edge < 0 )
            stop( '`', name, '` root edge length must be non-negative!  Observed: ', tree$root.edge )
    }
    
    # check overall dimensions consistency
    # edges should always have two columns only, check that immediately to catch most egregious errors
    if ( ncol( tree$edge ) != 2 )
        stop( '`tree$edge` does not have 2 columns (actual: ', ncol( tree$edge ), ')!' )
    n_tips <- length( tree$tip.label )
    n_nodes <- tree$Nnode
    # edges according to the first two counts
    # subtract one because root node is counted but has no edge
    n_edges <- n_tips + n_nodes - 1
    # now the other matrices/vectors must match this count
    if ( length( tree$edge.length ) != n_edges )
        stop( 'Number of edges in `', name, '$edge.length` (', length( tree$edge.length ), ') disagrees with expected number of edges (', n_edges, ') based on number of tips (', n_tips, ') and internal nodes excluding root (', n_nodes - 1, ')!' )
    if ( nrow( tree$edge ) != n_edges )
        stop( 'Number of edges in `', name, '$edge` (', nrow( tree$edge ), ' rows) and `', name, '$edge.length` (', n_edges, ') disagrees!' )
    # test contents of tree$edge
    # these are node indexes
    i_max <- max( tree$edge )
    if ( i_max != n_tips + n_nodes )
        stop( 'Number of node indexes in `', name, '$edge` (', i_max, ') does not match expected number of nodes from counting tips and nodes (', n_tips + n_nodes, ')!' )
    # all nodes must be present
    indexes_all <- 1 : i_max
    indexes_missing <- !( indexes_all %in% tree$edge )
    if ( any( indexes_missing ) )
        stop( 'These expected node indexes are missing in `', name, '$edge`: ', toString( indexes_missing ) )
    
    # if we didn't stop for any reason, we're good! (nothing to return)
}

## > names(tree)
## [1] "edge"        "tip.label"   "Nnode"       "edge.length"
## > tree$edge
##      [,1] [,2]
## [1,]    6    1
## [2,]    6    7
## [3,]    7    8
## [4,]    8    2
## [5,]    8    9
## [6,]    9    3
## [7,]    9    4
## [8,]    7    5
## > tree$tip.label
## [1] "t1" "t2" "t3" "t5" "t4"
## > tree$Nnode
## [1] 4
## > tree$edge.length
## [1] 0.2058876 0.7738696 0.2760511 0.7556557 0.3127623 0.6355037 0.9236313
## [8] 0.6934122
