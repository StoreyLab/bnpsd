# stops if both names exist but disagree (with descriptive messages), nothing happens otherwise
compare_names <- function( names1, names2, var1, var2 ) {
    # die if things are missing
    if ( missing( names1 ) )
        stop( '`names1` is required!' )
    if ( missing( names2 ) )
        stop( '`names2` is required!' )
    if ( missing( var1 ) )
        stop( '`var1` is required!' )
    if ( missing( var2 ) )
        stop( '`var2` is required!' )
    
    if ( !is.null( names1 ) && !is.null( names2 ) ) {
        if ( !all( names1 == names2 ) ) {
            # if names agree, then we're good
            # here we investigate why they don't agree even though both exist
            # could they be permuted?
            if ( all( names1 %in% names2 ) ) {
                # in this case yes, the labels are the same except they are out of order
                # best to stop so user can clean data up
                # (could reorder automatically, but that could lead to messiness elsewhere)
                stop( 'The names of `', var1, '` and `', var2, '` are both defined, contain the same values, but are not in the same order!  Please correct ordering and check your data!' )
            } else {
                # in this case they disagree in content as well, stop too but with a different message
                stop( 'The names of `', var1, '` and `', var2, '` are both defined but disagree in their contents!  Please check your data!' )
            }
        }
    }
}
