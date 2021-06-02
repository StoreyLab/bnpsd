#' Fit a tree structure to a coancestry matrix
#'
#' Implements a heuristic algorithm to find the optimal tree topology based on joining pairs of subpopulations with the highest between-coancestry values, and averaging parent coancestries for the merged nodes (a variant of WPGMA hierarchical clustering [stats::hclust()]).
#' Branch lengths are optimized using the non-negative least squares approach [nnls::nnls()], which minimize the residual square error between the given coancestry matrix and the model coancestry implied by the tree.
#' This algorithm recovers the true tree when the input coancestry truly belongs to this tree and there is no estimation error (i.e. it is the inverse function of [coanc_tree()]), and performs well for coancestry estimates (i.e. the result of simulating genotypes from the true tree, i.e. from [draw_all_admix()], and estimating the coancestry using [popkin::popkin()] followed by [popkin::inbr_diag()]).
#' 
#' The tree is bifurcating by construction, but edges may equal zero if needed, effectively resulting in multifurcation, although this code makes no attempt to merge nodes with zero edges.
#' For best fit to arbitrary data, the root edge is always fit to the data (may also be zero).
#' Data fit may be poor if the coancestry does not correspond to a tree, particularly if there is admixture between subpopulations.
#'
#' @param coancestry The coancestry matrix to fit a tree to.
#' @param method The hierarchical clustering method (passed to [stats::hclust()]).
#' While all [stats::hclust()] methods are supported, only two really make sense for this application: "mcquitty" (i.e. WPGMA, default) and "average" (UPGMA).
#'
#' @return A `phylo` object from the `ape` package (see [ape::read.tree()]), with two additional list elements at the end:
#' - `edge`: (standard `phylo`.)  A matrix describing the tree topology, listing parent and child node indexes.
#' - `Nnode`: (standard `phylo`.)  Number of internal (non-leaf) nodes.
#' - `tip.label`: (standard `phylo`.)  Labels for tips (leaf nodes), in order of index as in `edge` matrix above.
#'   These match the row names of the input `coancestry` matrix, or if names are missing, the row indexes of this matrix are used (in both cases, labels may be reordered compared to `rownames( coancestry )`).
#'   Tips are ordered as they appear in the above `edge` matrix (ensuring visual agreement in plots between the tree and its resulting coancestry matrix via [coanc_tree()]), and in an order that matches the input `coancestry`'s subpopulation order as much as possible (tree constraints do not permit arbitrary tip orders, but a heuristic implemented in [tree_reorder()] is used to determine a reasonable order when an exact match is not possible).
#' - `edge.length`: (standard `phylo`.)  Values of edge lengths in same order as rows of `edge` matrix above.
#' - `root.edge`: (standard `phylo`.)  Value of root edge length.
#' - `rss`:  The residual sum of squares of the model coancestry versus the input coancestry.
#' - `edge.length.add`: Additive edge values (regarding their effect on coancestry, as opposed to `edge.length` which are probabilities, see [coanc_tree()])
#'
#' @examples
#' # create a random tree
#' library(ape)
#' k <- 5
#' tree <- rtree( k )
#' # true coancestry matrix for this tree
#' coanc <- coanc_tree( tree )
#'
#' # now recover tree from this coancestry matrix:
#' tree2 <- fit_tree( coanc )
#' # tree2 equals tree!
#'
#' @seealso
#' [coanc_tree()]) for the inverse function (when the coancestry corresponds exactly to a tree).
#' 
#' [tree_reorder()] for reordering tree structure so that tips appear in a particular desired order.
#' 
#' @export
fit_tree <- function( coancestry, method = 'mcquitty' ) {
    if ( missing( coancestry ) )
        stop( '`coancestry` is required!' )
    # coancestry should be symmetric, etc
    if ( !isSymmetric( coancestry ) )
        stop( '`coancestry` must be symmetric!' )

    # data dimensions
    k_subpops <- nrow( coancestry )

    # identify original entries by name
    # if there are no names, use indexes as names
    if ( is.null( dimnames( coancestry ) ) ) {
        names <- 1 : k_subpops
        rownames( coancestry ) <- names
        colnames( coancestry ) <- names
    }
    # NOTE: isSymmetric tested earlier guarantees that names are the same in both dimensions
    
    # hclust requires a distance matrix
    # hack distance that agrees with our algorithm
    coanc_dist <- max( coancestry ) - coancestry
    # normal distances are zero for self, hack that into this too (won't matter otherwise)
    diag( coanc_dist ) <- 0
    # turn into proper distance object (doesn't change values)
    coanc_dist <- stats::as.dist( coanc_dist )

    # this corresponds exactly to my original algorithm!
    tree <- stats::hclust( coanc_dist, method = method )
    # turn hclust obtect into a phylo object
    tree <- ape::as.phylo( tree )
    
    # fit edges that minimize squared error
    tree <- fit_tree_single( coancestry, tree )
    # NOTE: do before `tree_reorder` because that one, in validation, requires `tree$edge.length` to be present!
    
    # now perform a more elaborate reordering of edges and tips
    # trees in the previous step are not guaranteed to be in the same order as the input coancestry
    # in general, it's not be possible for trees to be reordered to match arbitrary tip orders, but this code tries its best to give a best-match fit, in a sense, by greedily rotating nodes in "postorder"
    # the output tree has tip labels in matching order with plot (clade) order
    labels <- rownames( coancestry )
    tree <- tree_reorder( tree, labels )
    
    # expand from linear to probabilistic scale!
    tree <- tree_additive( tree, rev = TRUE )
    
    # done, return!
    return( tree )
}

# fits the tree (with a fixed topology) that best first a given coancestry matrix
# works on additive scale (input and outpus are so, the probabilistic scale is not used at any point)
# optimizes branch lengths of a given tree (input)
# so input is starting point, but also topology is fixed in this version
fit_tree_single <- function( coancestry, tree ) {
    if ( missing( coancestry ) )
        stop( '`coancestry` is required!' )
    if ( missing( tree ) )
        stop( '`tree` is required!' )

    # coancestry should be symmetric, etc
    if ( !isSymmetric( coancestry ) )
        stop( '`coancestry` must be symmetric!' )
    
    # run through validator for additional checks
    # NOTE: additive-scale trees always pass this because their edges are non-negative and always below 1 (so input trees in both `rev` cases should pass this)
    # NOTE: actually doesn't make sense to test edges since we just want topology, edges can be potentially bad (if they came from hclust) and this should still succeed
    # validate_coanc_tree( tree )

    # tree and coancestry data may be misaligned, it's important to get that right!
    # rely on names, if they already match (happens automatically for fit_tree and fit_tree_hclust, both when coancestry has labels as well as when they are absent, in which case tree labels are indexes!)
    # tree tips always have names, but maybe the coancestry matrix doesn't, so start checking there
    names_coanc <- rownames( coancestry )
    # NOTE: tree labels are always defined, though they can be empty strings (vector is always of the right length)
    names_tree <- tree$tip.label
    # again, the way things work internally, names_tree is always non-NULL, and names_coanc have an obvious default, which we apply now if needed
    if ( is.null( names_coanc ) ) {
        k_subpops <- nrow( coancestry )
        names_coanc <- 1 : k_subpops
        rownames( coancestry ) <- names_coanc
        colnames( coancestry ) <- names_coanc
    }
    # now both vectors are defined always, let's warn if they don't agree, reorder successfully otherwise
    # suffices to see if the sets agree (what we require ultimately)
    if ( all( names_coanc %in% names_tree ) ) {
        # now apply necessary reordering for coancestry
        indexes <- match( names_tree, names_coanc )
        coancestry <- coancestry[ indexes, indexes ]
        # check that the names actually agree now
        stopifnot( rownames( coancestry ) == names_tree )
    } else
        warning( 'Could not align `coancestry` to `tree` because names disagree!  Names in `coancestry`: ', toString( names_coanc ), '.  Names in `tree`: ', toString( names_tree ) )
    
    # optimization requires each branch to have an indicator block
    # construct then starting now

    # first step is really to know which are the tips and which ones are children to each internal edge
    # separated into another function
    edge_to_tips <- edges_to_tips( tree )
    # create indicator block matrices now
    n_tips <- ape::Ntip( tree )
    ## edge_to_blocks <- lapply( edge_to_tips, tips_to_block, n_tips ) # list-of-matrices version
    X <- sapply( edge_to_tips, tips_to_block, n_tips ) # direct to matrix-of-vectors version

    # there are several alternatives to do NNLS
    # https://www.r-bloggers.com/2019/11/non-negative-least-squares/

    # frame as linear regression with non-negative coefficients
    # have to linearize the coancestry matrix, which is the data to fit
    y <- sym_mat_to_vec( coancestry )
    # same for the edges, but this results in a (design) matrix
    ## X <- sapply( edge_to_blocks, sym_mat_to_vec ) # matrix-of-vectors from list-of-matrices version
    # hackily add root node to beginning
    X <- cbind( 1, X )
    # all the hard work is done by NNLS
    obj <- nnls::nnls( X, y )
    # extract fit coefficients and residual sum of squares
    edge_length_fit <- obj$x
    rss <- sum( obj$residuals^2 )

    # validate model fits
    if ( any( edge_length_fit < 0 ) )
        warning( 'Optimal tree edges were negative!' )
    if ( any( edge_length_fit > 1 ) )
        warning( 'Optimal tree edges were > 1!' )
    
    # when we're done, overwrite edges on tree and return that!
    # always fit root edge, but have to separate it here
    tree$edge.length <- edge_length_fit[ -1 ]
    tree$root.edge <- edge_length_fit[ 1 ]
    # use RSS from linear model fit
    tree$rss <- rss
    
    return( tree )
}

sym_mat_to_vec <- function( x ) {
    # flatten a symmetric matrix, in the most R way possible
    # NOTE: no need in this application for a reverse function
    x <- x[ lower.tri( x, diag = TRUE ) ]
    return( x )
}

tips_to_block <- function( tips, n_tips ) {
    # create a blank square matrix of dimensions n_tips
    B <- matrix( 0, n_tips, n_tips )
    # fill in the desired tips with ones (all pairs between tips)
    B[ tips, tips ] <- 1
    # flatten matrix into vector
    B <- sym_mat_to_vec( B )
    # all done!
    return( B )
}

# get tips covered by each edge, very useful for tree fitting
# surprisingly I didn't find this functionality in ape/phytools
# assumes tree has been validated already
edges_to_tips <- function( tree ) {
    if ( missing( tree ) )
        stop( '`tree` is required!' )

    # optimization requires each branch to have an indicator block
    # construct then starting now

    # first step is really to know which are the tips and which ones are children to each internal edge
    
    # number of edges
    n_edges <- ape::Nedge( tree )
    # number of tips
    n_tips <- ape::Ntip( tree )
    # number of internal nodes
    n_nodes <- ape::Nnode( tree )
    # edge to tips
    edge_to_tips <- vector( 'list', n_edges )
    # and a separate list of internal nodes to tips (similar info but it is organized differently)
    # this is an intermediate structure, ultimately not of interest
    node_to_tips <- vector( 'list', n_nodes + n_tips )
    
    # algorithm is sensitive to edge ordering
    # assumption is that we move from the tips down, which is postorder
    order_edges <- ape::postorder( tree )
    
    # navigate all edges in postorder!
    for ( e in order_edges ) {
        # get parent and child nodes for this edge
        j_parent <- tree$edge[ e, 1 ]
        j_child <- tree$edge[ e, 2 ]

        # extract total descendants of child
        # actually we're only interested in tips, so internal nodes are excluded
        if ( j_child <= n_tips ) {
            # child is a tip node, so this is complete list
            js_tips <- j_child
        } else {
            # child is not a tip, so it has its own children
            # retrieve those values, this gives complete list
            # don't include j_child because is an internal node!
            js_tips <- node_to_tips[[ j_child ]]
            #js_tips <- c( j_child, node_to_tips[[ j_child ]] ) # this includes internal nodes too
        }

        # add list of tip nodes to descendants of current edge
        edge_to_tips[[ e ]] = js_tips
        
        # add list of tip nodes to descendants of current parent
        # needed for next round
        # NOTE: parents can appear multiple times (as parents of multiple edges), so always concatenate
        # sort for extra niceness, though this is probably not necessary
        node_to_tips[[ j_parent ]] = sort( c( node_to_tips[[ j_parent ]], js_tips ) )
    }

    # done, return this list!
    return( edge_to_tips )
}
