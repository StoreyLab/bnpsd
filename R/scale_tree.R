#' Scale a coancestry tree
#'
#' Scale a tree in the additive scale by `factor`.
#' This results in the coancestry matrix of this tree being scaled by the same factor.
#' The (non-additive) IBD edges are transformed in a non-linear fashion.
#' 
#' Internally, additive edges (`$edge.length.add`) are calculated using [tree_additive()] if not already present, then they are all scaled, including the root edge (`$root.edge`) if present, and lastly IBD edges (`$edge.length`) are recalculated from the additive edges using [tree_additive()] with option `rev = TRUE`.
#' Stops if any of the edges exceed 1 before or after scaling (since these edges are IBD probabilities).
#'
#' @param tree The coancestry tree to edit.
#' Must be a `phylo` object from the `ape` package.
#' @param factor The scalar factor to multiply all edges in additive scale.
#' Must be non-negative, and not be so large that any edge exceeds 1 after scaling.
#'
#' @return The edited tree with all edges scaled as desired.
#' This tree will always contain (correctly scaled) additive edges in addition to the default non-additive edges.
#'
#' @examples
#' # create a random tree
#' library(ape)
#' k <- 5
#' tree <- rtree( k )
#'
#' # scale this tree
#' tree_scaled <- scale_tree( tree, 0.5 )
#'
#' @seealso
#' [ape::read.tree()] for reading `phylo` objects and their definition.
#'
#' [tree_additive()] for difference between IBD and additive edges.
#'
#' @export
scale_tree <- function( tree, factor ) {
    # require both inputs
    if ( missing( tree ) )
        stop( '`tree` is required!' )
    if ( missing( factor ) )
        stop( '`factor` is required!' )

    # tree should be valid already
    validate_coanc_tree( tree )
    # validate factor
    if ( length( factor ) != 1 )
        stop( '`factor` must be scalar!  Observed length: ', length( factor ) )
    if ( !is.numeric( factor ) )
        stop( '`factor` must be numeric!  Observed value: ', factor )
    if ( factor < 0 )
        stop( '`factor` must be non-negative!  Observed value: ', factor )

    # start editing values, but will check each carefully to make sure we don't go out of bounds!
    # edge lengths are the majority of the work usually:
    # easiest way to do it is to modify additive edges first, then convert back to non-additive
    # calculate additive edges if not present
    if ( is.null( tree$edge.length.add ) )
        tree <- tree_additive( tree )
    # apply factor!
    tree$edge.length.add <- tree$edge.length.add * factor # correct!
    # die if any edge exceeds one now
    # (negatives don't occur given what we've checked already)
    if ( any( tree$edge.length.add > 1 ) )
        stop( 'At least one `tree` additive edge length exceeds 1 after scaling, max: ', max( tree$edge.length.add ) )
    # edit root edge the same way if present
    if ( !is.null( tree$root.edge ) ) {
        # this also gets checked
        tree$root.edge <- tree$root.edge * factor
        # die if any edge exceeds one now
        # (negatives don't occur given what we've checked already)
        if ( any( tree$root.edge > 1 ) )
            stop( 'Root edge length exceeds 1 after scaling: ', tree$root.edge )
    }
    
    # now recalculate and overwrite non-additive edges
    # (there's no easy formula for scaling these non-additive values, this is easiest)
    # have to do it in this awkward way, first overwrite main edges with additive ones, then delete additive edges, then calculate non-additive ones with the second command
    tree$edge.length <- tree$edge.length.add
    tree <- tree_additive( tree, rev = TRUE, force = TRUE )
    # check non-additive edges too
    if ( any( tree$edge.length > 1 ) )
        stop( 'At least one `tree` edge length exceeds 1 after scaling, max: ', max( tree$edge.length ) )
    
    # return now if all was scaled and validated
    return( tree )
}
