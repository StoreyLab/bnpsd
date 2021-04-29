#' Calculate additive edges for a coancestry tree, or viceversa
#'
#' A coancestry tree has IBD probabilities as edge values, which describe how each child and parent subpopulation in the tree is related.
#' This means that each parameter is relative to its parent subpopulation (varies per edge), and they are not in general IBD probabilities from the root.
#' This function computes "additive" edges that corresponds more closely with the coancestry matrix of this tree, which gives parameters relative to the root (ancestral) population (see details below).
#' The additive edges are computed on a new element of the tree `phylo` object, so they do not overwrite the probabilistic edges.
#' The reverse option assumes that the main edges of the `phylo` object are additive edges, then calculates the probabilistic edges and stores as the main edges and moves the original additive edges on the new element.
#' 
#' The calculation takes into account that total coancestries are non-linear combinations of the per-edge coancestries.
#' For example, if the root node is `A`, and subpopulation `C` is connected to `A` only through an internal node `B`, then its total self-coancestry `coanc_A_C` relative to `A` is given by `coanc_A_B` (the coancestry between `A` and `B`) and `coanc_B_C` (the coancestry between `B` and `C`) by
#' `coanc_A_C = coanc_A_B + coanc_B_C * ( 1 - coanc_A_B )`.
#' This transformation ensures that the total coancestry is a probability that does not exceed one even when the per-edge coancestries sum to a value greater than one.
#' The "additive" edge for `B` and `C` is `coanc_B_C * ( 1 - coanc_A_B )`, so it is the probabilistic edge `coanc_B_C` shrunk by `1 - coanc_A_B`, which can then just be added to the parent edge `coanc_A_B` to give the total coancestry `coanc_A_C`.
#' This transformation is iterated for all nodes in the tree, noting that roots `B` connected to the root node `A` have equal probabilistic and additive edges `coanc_A_B` (unless the tree has a root edge, in which case that one is used to transform as above), and the edge of a node `C` not directly connected to a root uses the calculated edge `coanc_A_C` as above to shrink its children's edges into the additive scale.
#'
#' @param tree The coancestry tree with either probabilistic edges (if `rev = FALSE`) or additive edges (if `rev = TRUE`) as the main edges (stored in `tree$edge.length`).
#' Must be a `phylo` object from the `ape` package (see [ape::read.tree()]).
#' This tree may have a valid root edge (non-NULL `tree$root.edge` between 0 and 1), which is incorporated in the output calculations.
#' Function stops if the input data is not valid.
#' Probabilistic edges are valid if and only if they are all between zero and one.
#' Additive edges are valid if and only if they are all non-negative and the sum of edges between the root and each tip (leaf node) does not exceed 1.
#' @param rev If `FALSE` (default), assumes the main edges are probabilistic values, and calculates additive values.
#' If `TRUE`, assumes main edges are additive values, and calculates probabilistic values.
#' @param force If `FALSE` (default), function stops if input tree already has additive edges (if `tree$edge.length.add` is not `NULL`).
#' If `TRUE`, these values are ignored and overwritten.
#'
#' @return The input `phylo` object extended so that the main edges (`tree$edge.length`) are probabilistic edges, and the additive edges are stored in `tree$edge.length.add`.
#' This is so for both values of `rev`
#'
#' @examples
#' # for simulating a tree with `rtree`
#' library(ape)
#' 
#' # SETUP: number of tip subpopulations
#' k <- 5
#' # simulate a random tree
#' # edges are drawn from Uniform(0, 1), so these are valid probabilistic edges
#' tree <- rtree( k )
#' # inspect edges
#' tree$edge.length
#' 
#' # RUN calculate additive edges (safe to overwrite object)
#' tree <- tree_additive( tree )
#' # inspect edges again
#' # probabilistic edges are still main edges:
#' tree$edge.length
#' # additive edges are here
#' tree$edge.length.add
#'
#' 
#' # let's go backwards now, starting from the additive edges
#' # SETUP
#' # these are harder to simulate, so let's copy the previous value to the main edges
#' tree$edge.length <- tree$edge.length.add
#' # and delete the extra entry (if it's present, function stops)
#' tree$edge.length.add <- NULL
#' # inspect before
#' tree$edge.length
#' 
#' # RUN reverse version (safe to overwrite object again)
#' tree <- tree_additive( tree, rev = TRUE )
#' # inspect after
#' # probabilistic edges are main edges:
#' tree$edge.length
#' # additive edges (previously main edges) were moved here
#' tree$edge.length.add
#'
#' @seealso
#' [coanc_tree()], the key application facilitated by additive edges.
#' 
#' @export
tree_additive <- function( tree, rev = FALSE, force = FALSE ) {
    if ( missing( tree ) )
        stop( '`tree` is required!' )
    # run through validator for additional checks
    # NOTE: additive-scale trees always pass this because their edges are non-negative and always below 1 (so input trees in both `rev` cases should pass this)
    validate_coanc_tree( tree )
    
    # input tree shouldn't have additive data calculated already (either way it means there's nothing to do)
    if ( !is.null( tree$edge.length.add ) ) {
        if ( force ) {
            # just pretend we didn't have them
            # set to NULL for safety, though code isn't supposed to use this anyway (it ought to just get ovewritten later)
            tree$edge.length.add <- NULL
        } else 
            stop( '`tree$edge.length.add` is not `NULL`, suggests this tree already has probabilistic and additive edges!  Tip: if you want to ignore these values, use option `force = TRUE`!' )
    }
    
    # algorithm is sensitive to edge ordering
    # assumption is that we move from the root up, which is the reverse of postorder
    order_edges <- rev( ape::postorder( tree ) )
    
    # number of edges
    n_edges <- ape::Nedge( tree )
    # determine root node
    # it is very first parent node (in reverse postorder)
    j_root <- tree$edge[ order_edges[ 1 ], 1 ]
    # keep track of additive or probability edges for each child from the root as it becomes a parent
    # a list in parallel to the normal list of edges for each node to its parent only
    if ( rev ) {
        edge_length_prob <- vector( 'numeric', n_edges )
    } else {
        edge_length_add <- vector( 'numeric', n_edges )
    }
    # a list with the same info but indexed by nodes instead, easier to retrieve in this case
    # initialize this one to NA to catch issues
    node_edge_from_root <- rep.int( NA, max( tree$edge ) )
    # only root node is pre-set to zero, or to value of root edge if present
    node_edge_from_root[ j_root ] <- if ( is.null( tree$root.edge ) ) 0 else tree$root.edge

    # navigate all edges in reverse postorder!
    for ( e in order_edges ) {
        # get parent and child nodes for this edge
        j_parent <- tree$edge[ e, 1 ]
        j_child <- tree$edge[ e, 2 ]
        if (rev) {
            # edge length is in additive scale
            edge_child_add <- tree$edge.length[ e ]
        } else {
            # edge length is FST we need
            edge_child_prob <- tree$edge.length[ e ]
        }
        # get edge length of parent from root
        fst_parent_from_root <- node_edge_from_root[ j_parent ]
        # check
        if ( is.na( fst_parent_from_root ) )
            stop( 'Node index ', j_parent, ' was not assigned coancestry from root! (unexpected)' )
        if (rev) {
            # stretch edge now, store in new edge vector
            edge_length_prob[ e ] <- edge_child_add / ( 1 - fst_parent_from_root )
        } else {
            # shrink edge now, store in new edge vector
            edge_child_add <- edge_child_prob * ( 1 - fst_parent_from_root )
            edge_length_add[ e ] <- edge_child_add
        }
        # and store total edge for child to root (as it may become a parent in later iterations)
        # the additive edge is what we want to add here (both `rev` cases)
        node_edge_from_root[ j_child ] <- fst_parent_from_root + edge_child_add
    }

    if (rev) {
        # we could store `node_edge_from_root`, but we won't for now
        # we include additive edges vector, without erasing original for now
        # let's do a switcheroo so standard edge is probabilistic one and additive one is the extra info
        tree$edge.length.add <- tree$edge.length
        tree$edge.length <- edge_length_prob
        
        # a final check, if the input tree was not a valid additive tree, several parameters will be out of bounds ("probabilities" greater than 1 or negative)
        # let's look for those, stop if encountered
        # the most obvious tell-tale sign is that total edge lengths from roots, which are sums of the original edges, exceed 1 (as `validate_coanc_tree` guaranteed non-negative edges, their sums are also non-negative)
        # because of limited machine precision, let's add an epsilon here
        indexes_bad <- node_edge_from_root > 1 + sqrt( .Machine$double.eps )
        if ( any( indexes_bad ) )
            stop( 'Input tree was not coancestry tree in "additive" form, since total edge (sum) for some tip nodes exceeded 1!  These node indexes were bad: ', toString( which( indexes_bad ) ), '.  These were their sum of edge values from root: ', toString( node_edge_from_root[ indexes_bad ] ), '; excess: ', toString( node_edge_from_root[ indexes_bad ] - 1 ) )
        # if no total edges exceeded 1, then we definitely won't get negative prob edges (division by `1 - fst_parent_from_root` above is always positive).
        # also, total edge for node j_child above not exceeding 1 means
        # `fst_parent_from_root + edge_child_add <= 1`,
        # which guarantees that
        # `edge_child_add / ( 1 - fst_parent_from_root ) <= 1`,
        # so its probabilistic edge is also guaranteed to be in range
        # so there's nothing else to test for
    } else {
        # we could store `node_edge_from_root`, but we won't for now
        # we include additive edges vector, without erasing original for now
        tree$edge.length.add <- edge_length_add
    }
    
    # return this extended tree
    return( tree )
}
