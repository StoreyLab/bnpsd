#' Draw allele frequencies for subpopulations related by a tree
#'
#' The allele frequency matrix `P` for `m_loci` loci (rows) and `k_subpops` subpopulations (columns) are drawn from the Balding-Nichols distribution.
#' If the allele frequency in the parent node is `p` and the FST parameter of the child node from the parent node is `F`, then the allele frequency in the child node is drawn from
#' `rbeta( 1, nu * p, nu * ( 1 - p ) )`,
#' where `nu <- 1 / F - 1`.
#' This function iterates drawing allele frequencies through the tree structure, therefore allowing covariance between subpopulations that share branches.
#'
#' @param p_anc The scalar or length-`m_loci` vector of ancestral allele frequencies per locus.
#' @param tree_subpops The coancestry tree relating the `k_subpops` subpopulations.
#' Must be a `phylo` object from the `ape` package (see [ape::read.tree()]).
#' The edge lengths of this tree must be the FST parameters relating parent and child subpopulations.
#' @param m_loci If `p_anc` is scalar, optionally provide the desired number of loci (lest only one locus be simulated).
#' Stops if both `length(p_anc) > 1` and `m_loci` is not `NA` and they disagree.
#' @param nodes  If `FALSE` (default), returns allele frequencies of "tips" subpopulations only; otherwise returns all allele frequencies, including internal nodes
#'
#' @return The `m_loci`-by-`k_subpops` matrix of independent subpopulation allele frequencies.
#' If `nodes = FALSE`, columsn include only tip subpopulations.
#' If `nodes = TRUE`, internal node subpopulations are also included (in this case the input `p_anc` is returned in the column corresponding to the root node).
#' In all cases subpopulations are ordered as indexes in the tree, which normally implies the tip nodes are listed first, followed by the internal nodes (as in `tree_subpops$edge` matrix, see [ape::read.tree()] for details).
#'
#' @examples
#' # for simulating a tree with `rtree`
#' library(ape)
#' 
#' # a typical, non-trivial example
#' # number of tip subpopulations
#' k_subpops <- 3
#' # number of loci
#' m_loci <- 10
#' # random vector of ancestral allele frequencies
#' p_anc <- draw_p_anc(m_loci)
#' # simulate tree
#' tree_subpops <- rtree( k_subpops )
#' # matrix of intermediate subpop allele freqs
#' p_subpops <- draw_p_subpops_tree(p_anc, tree_subpops)
#'
#' @seealso
#' [draw_p_subpops()] for version for independent subpopulations.
#' 
#' For "phylo" tree class, see [ape::read.tree()]
#' 
#' @export
draw_p_subpops_tree <- function(
                                p_anc,
                                tree_subpops,
                                m_loci = NA,
                                nodes = FALSE
                                ) {
    # basic param checking
    if ( missing( p_anc ) )
        stop('ancestral allele frequencies `p_anc` are required!')
    if ( missing( tree_subpops ) )
        stop('`tree_subpops` (FST) scalar or vector are required!')
    
    # run through special tree validator
    validate_coanc_tree( tree_subpops, name = 'tree_subpops' )

    if ( !is.null( tree_subpops$root.edge ) )
        warning( 'Input `tree_subpops` has a root edge (`tree_subpops$root.edge = ', tree_subpops$root.edge, '`) that will be ignored by function `draw_p_subpops_tree`!' )
    
    # look at data ranges
    # all of these are probabilities so problems happen when they're out of range
    if ( any( p_anc < 0 ) )
        stop( '`p_anc` cannot be negative!' )
    if ( any( p_anc > 1 ) )
        stop( '`p_anc` cannot exceed 1!' )
    
    # number of loci to simulate
    # allow p_anc to be a scalar, m_loci can be passed separately
    # actual length
    m_loci_p <- length(p_anc)
    # make sure both things were not set and contradict each other
    if (
        m_loci_p > 1 &&    # length(p_anc) > 1
        !is.na(m_loci) &&  # and m_loci was also set
        m_loci != m_loci_p # and they disagree
    )
        stop('length of `p_anc` (', m_loci_p, ') disagrees with passed parameter `m_loci` (', m_loci, ')')
    # if this is missing, always set to actual value (even if it is 1)
    if ( is.na( m_loci ) )
        m_loci <- m_loci_p
    
    # number of subpopulations to simulate
    k_subpops <- length( tree_subpops$tip.label )

    # and total number of nodes (initially we keep them all in the output matrix
    k_subpops_all <- max( tree_subpops$edge )
    
    # matrix of intermediate allele frequencies we want
    p_subpops <- matrix( nrow = m_loci, ncol = k_subpops_all )

    # it'd be nice to keep names if they're available
    # they're separate, unfortunately
    # need to join them before assigning
    # 1) tip labels
    # unfortunately this is always defined, even when there are no labels
    names_subpops <- tree_subpops$tip.label
    if ( any( names_subpops != '' ) ) {
        # we have non-trivial labels!
        # NOTE: we won't incorporate node labels if there aren't any tip labels (that'd be weird)
        # 2) node labels, these are NULL if missing
        names_node <- tree_subpops$node.label
        # fill with empty strings, since the final names must have length `k_subpops_all`
        if ( is.null( names_node ) ) 
            names_node <- rep.int( '', tree_subpops$Nnode )
        # concatenate what we have
        names_subpops <- c( names_subpops, names_node )
        # check final length
        if ( length( names_subpops ) != k_subpops_all ) 
            stop( 'Names of tips and nodes have length ', length( names_subpops ), ', expected ', k_subpops_all, '!  Is the `tree_subpops` object malformed?' )
        # add names to matrix
        colnames( p_subpops ) <- names_subpops
    }
    
    # initialize the root AFs
    # this is its column
    j_root <- tree_subpops$edge[ 1, 1 ]
    # make sure we got it right, should be 1 + the number of leaf nodes
    # the code doesn't really depend on this assumption exactly, although it'd be very bad if j_root <= k_subpops
    if ( j_root != k_subpops + 1 )
        stop( 'The root node index in `tree_subpops$edge` (', j_root, ') does not match `k_subpops + 1` (', k_subpops + 1, ') where `k_subpops` is the number of tips!  Is the `tree_subpops` object malformed?' )
    # now copy input AFs into that "root" column
    # this works if p_anc is a scalar as well as a length-m vector
    p_subpops[ , j_root ] <- p_anc

    # for extra safety, keep track of nodes that have been initialized with AFs
    nodes_done <- j_root # just this one for now
    
    # now we navigate the edges, which in the current structure move away from the root, so it's the perfect order for our operations
    n_edges <- nrow( tree_subpops$edge )
    # make sure "edge lengths" vector has right length
    if ( length( tree_subpops$edge.length ) != n_edges )
        stop( 'Length of `tree_subpops$edge.length` (', length( tree_subpops$edge.length ), ') does not match number of rows of `tree_subpops$edge` (', n_edges, ')!  Is `tree_subpops` object malformed?' )
    # start navigating
    for ( e in 1 : n_edges ) {
        # get parent and child nodes for this edge
        j_parent <- tree_subpops$edge[ e, 1 ]
        j_child <- tree_subpops$edge[ e, 2 ]
        # edge length is FST we need
        fst <- tree_subpops$edge.length[ e ]
        
        # make sure parent has been processed already
        if ( !( j_parent %in% nodes_done ) )
            stop( 'Parent node index ', j_parent, ' has not yet been processed (nodes processed: ', toString( nodes_done ), ')!  Is `tree_subpops` object malformed?' )
        # make sure child hasn't been processed already
        if ( j_child %in% nodes_done )
            stop( 'Child node index ', j_child, ' has already been processed (nodes processed: ', toString( nodes_done ), ')!  Is `tree_subpops` object malformed?' )
        
        # get allele frequencies from parent column
        p_parent <- p_subpops[ , j_parent ]
        # transform FST into this other factor
        nu <- 1 / fst - 1

        # handle infinite "nu" special case
        if ( is.infinite( nu ) ) {
            # there is no drift from p_anc in this special case (inbr_subpops == 0)
            # (coded separately because rbeta incorrectly returns 0.5 instead)
            p_subpops[ , j_child ] <- p_parent
        } else {
            # draw all SNPs for this population, store immediately
            p_subpops[ , j_child ] <- stats::rbeta( m_loci, nu * p_parent, nu * ( 1 - p_parent ) )
        }
        
        # add child node to processed nodes now
        nodes_done <- c( nodes_done, j_child )
    }

    # check that all nodes were processed
    nodes_exp <- 1 : k_subpops_all
    nodes_mis <- setdiff( nodes_done, nodes_exp )
    if ( length( nodes_mis ) > 0 )
        stop( 'Some nodes were not processed: ', toString( nodes_mis ), '!  Is `tree_subpops` object malformed?' )
    
    # subset to return "tips" (leaf nodes) only if required
    if ( !nodes )
        p_subpops <- p_subpops[ , 1 : k_subpops, drop = FALSE ]

    # done!
    return( p_subpops )
}
