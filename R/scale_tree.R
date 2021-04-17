#' Scale a coancestry tree
#'
#' Scale by a scalar `factor` all the edges (`$edge.length`) of a `phylo` object from the `ape` package, including the root edge (`$root.edge`) if present, and additive edges (`$edge.length.add`, present in trees returned by [fit_tree()]).
#' Stops if any of the edges exceed 1 before or after scaling (since these edges are IBD probabilities).
#'
#' @param tree The coancestry tree to edit.
#' @param factor The scalar factor to multiply all edges.
#' Must be non-negative, and not be so large that any edge exceeds 1 after scaling.
#'
#' @return The edited tree with all edges scaled as desired.
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
    # apply factor
    tree$edge.length <- tree$edge.length * factor
    # die if any edge exceeds one now
    # (negatives don't occur given what we've checked already)
    if ( any( tree$edge.length > 1 ) )
        stop( 'At least one `tree` edge length exceeds 1 after scaling, max: ', max( tree$edge.length ) )
    # edit additive edges if present
    # this case doesn't need any checking (if the previous step worked, all is good)
    if ( !is.null( tree$edge.length.add ) )
        tree$edge.length.add <- tree$edge.length.add * factor
    
    # edit root edge if present
    if ( !is.null( tree$root.edge ) ) {
        # this also gets checked
        tree$root.edge <- tree$root.edge * factor
        # die if any edge exceeds one now
        # (negatives don't occur given what we've checked already)
        if ( any( tree$root.edge > 1 ) )
            stop( 'Root edge length exceeds 1 after scaling: ', tree$root.edge )
    }

    # return now if all was scaled and validated
    return( tree )
}
