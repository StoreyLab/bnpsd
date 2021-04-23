# a flexible function that handles inputs that are vectors, matrices, and trees
names_coanc <- function( coanc ) {
    if ( missing( coanc ) )
        stop( '`coanc` is required!' )
    my_names <- NULL
    if ( is.matrix( coanc ) ) {
        # get rownames directly
        my_names <- rownames( coanc )
        if ( !is.null( my_names ) ) {
            # should ensure they're the same on both dimensions of matrix
            if ( !all( colnames( coanc ) == my_names ) )
                stop( 'Coancestry matrix does not have the same names along rows and columns!' )
        }
        # else we return the names we got (which may be NULL)
    } else if ( 'phylo' %in% class( coanc ) ) {
        # have a tree!
        # get tip labels
        my_names <- coanc$tip.label
        # these names are always non-NULL/NA, but if they're all blank we should return NULL
        if ( all( my_names == '' ) )
            my_names <- NULL
    } else {
        # should be left with vector/scalar case
        if ( length( coanc ) > 1 ) {
            # leave scalar case as NULL (regardless of name on that scalar)
            # for proper vectors, this returns desired names (which may be NULL)
            my_names <- names( coanc )
        }
    }
    return( my_names )
}
