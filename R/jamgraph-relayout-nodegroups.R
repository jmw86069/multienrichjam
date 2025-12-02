

## Design idea:
# take Cnet layout
# - iterate each nodegroup, relayout just that group
# - goal is for nodegroups to be preferentially co-located


#' Relayout each nodegroup in a bipartite (Cnet) graph, experimental
#' 
#' Relayout each nodegroup in a bipartite (Cnet) graph, experimental
#' 
#' This function iteratively re-applies a layout function
#' to each nodegroup in a graph, constraining the position of
#' all other nodes for each iteration.
#' Currently the layout uses `relayout_with_qfr()`, in future it
#' may use any layout function.
#' 
#' The purpose is to "encourage" nodes in a nodegroup to become
#' bundled together.
#' 
#' Strategy:
#' 
#' * Each nodegroup is isolated, and `relayout_with_qfr()` is called
#' on nodes in each nodegroup, while constraining all other nodes so
#' they cannot move.
#' * In theory, using the same `repulse`, the nodes would not move at all.
#' Changing the repulse force could encourage nodes to stay together.
#' * By default `add_edges=TRUE` which adds phantom edges to connect all
#' nodes in a node group.
#' 
#'    * The edge weight is scaled down by the number
#'    of nodes.
#'    * These phantom edges are intended to help 'encourage' the nodegroup
#'    nodes to group together.
#'    * Otherwise, most layout algorithms are only focused on specific
#'    edge forces, and not secondary forces which are common in Cnet plots.
#'    * For example, nodes in a nodegroup all share
#'    the same network connections, however they are not otherwise
#'    attracted to each other in a network layout. In absence of any
#'    repulsive force, they would all be co-located. But with some
#'    repulsive force, they are repelled from each other, and sometimes
#'    end up radially positioned around the plot, and not grouped together.
#'    * The phantom edges add some minimal force for nodes to be
#'    grouped closer together, and are removed once the layout is complete.
#' 
#' * As a final polishing step, `do_final_relayout=TRUE` enables a final
#' round of global node layout, with `final_repulse`.
#' 
#'    * We observed that sometimes the nodegroups are too clumped,
#'    in a big circular "ball" due to the phantom edge process above.
#'    The final relayout is helpful to allow nodes to space out somewhat.
#'    * It may be helpful to pass `niter` to control the number of
#'    layout iterations in this final step. The default is 500.
#' 
#' ## Todo
#' 
#' Bonus points:
#' 
#' * Determine if each nodegroup is "split" around another nodegroup.
#' If so, iterate that nodegroup using varying 'repulse' or other
#' strategies to attempt to minimize the issue.
#' 
#' @family jam cnet utilities
#' 
#' @returns `igraph` object with updated layout.
#' 
#' @param cnet `igraph` object with node layout already defined
#' @param repulse `numeric` default 3.5, passed to `relayout_with_qfr()`
#' @param fix_set_nodes `logical` default TRUE, whether to fix all
#'    nodes with nodeType=='Set' to prevent them from moving.
#' @param spread_labels `logical` default TRUE, whether to apply
#'    `spread_igraph_labels()` after the relayout iterations are complete.
#' @param nodegroups `list` default NULL, intended to pass a custom
#'    set of nodegroups. When NULL it uses `get_cnet_nodeset()`.
#' @param ... additional arguments are passed to
#'    `spread_igraph_labels()`
#' 
#' @examples
#' cnet <- make_cnet_test();
#' ns <- get_cnet_nodeset(cnet)
#' # mark.groups: highlight just one nodegroup
#' # nodegroups: enables the edge_bundling to work properly
#' jam_igraph(cnet, mark.groups=ns["SetA,SetB"], nodegroups=ns,
#'    main="One nodegroup is split")
#' 
#' cnet2 <- relayout_nodegroups(cnet, nodegroups=ns["SetA,SetB"])
#' jam_igraph(cnet2, mark.groups=ns["SetA,SetB"], nodegroups=ns,
#'    main="This nodegroup is 'encouraged' to re-group")
#' 
#' cnet2 <- relayout_nodegroups(cnet, nodegroups=ns)
#' jam_igraph(cnet2, mark.groups=ns["SetA,SetB"], nodegroups=ns,
#'    main="Re-layout across all nodegroups")
#' 
#' @export
relayout_nodegroups <- function
(cnet,
 repulse=3.5,
 fix_set_nodes=TRUE,
 spread_labels=TRUE,
 nodegroups=NULL,
 add_edges=TRUE,
 edge_factor=2,
 do_final_relayout=NULL,
 final_repulse=3.5,
 ...)
{
	#
	cnet2 <- cnet;
	if (length(nodegroups) == 0) {
		nodegroups <- get_cnet_nodeset(cnet)
		# sort by size
		ngnum <- lengths(strsplit(names(nodegroups), ","))
		ngdf <- data.frame(ng=names(nodegroups),
			ngnum=ngnum,
			ngsize=lengths(nodegroups))
		ngdf <- jamba::mixedSortDF(ngdf, byCols=c("-ngnum", "-ngsize", "ng"))
		nodegroups <- nodegroups[ngdf$ng];
		print(ngdf);# debug
	}
	
	if (length(do_final_relayout) == 0) {
		if (length(nodegroups) > 1) {
			do_final_relayout <- TRUE;
		} else {
			do_final_relayout <- FALSE;
		}
	}
	# Fix the set nodes
	fixed_nodes <- NULL;
	unfixed_nodes <- NULL;
	if (length(fix_set_nodes) == 0) {
		# do nothing
	} else if (isTRUE(fix_set_nodes)) {
		fixed_nodes <- igraph::V(cnet2)[
			igraph::V(cnet2)$nodeType %in% "Set"]$name;
	} else {
		unfixed_nodes <- igraph::V(cnet2)[
			igraph::V(cnet2)$nodeType %in% "Set"]$name;
	}
	# iterate each nodegroup
	for (i in seq_along(nodegroups)) {
		move_nodes <- nodegroups[[i]];
		use_constrain <- unique(c(
			setdiff(igraph::V(cnet2)$name,
				c(move_nodes, unfixed_nodes)),
			fixed_nodes))
		
		# optionally add edges between nodes
		if (isTRUE(add_edges)) {
			inodes <- sort(match(move_nodes, igraph::V(cnet2)$name));
			if (length(inodes) > 1) {
				iedges <- (combn(inodes, 2));
				print(length(inodes));# debug
				print(ncol(iedges));# debug
				cnet2e <- cnet2;
				if (!"weight" %in% igraph::edge_attr_names(cnet2)) {
					igraph::E(cnet2e)$weight <- 1;
				}
				# add complete edges
				cnet2e <- igraph::add_edges(cnet2e, iedges,
					attr=list(weight=edge_factor/sqrt(length(inodes))))
					# attr=list(weight=edge_factor))
				# print(tail(igraph::as_data_frame(cnet2e, "edges")), ncol(iedges) + 1);# debug
				# relayout
				cnet2e <- relayout_with_qfr(cnet2e,
					repulse=repulse,
					constrain=use_constrain,
					spread_labels=FALSE,
					...)
				cnet2 <- set_igraph_layout(cnet2,
					layout=get_igraph_layout(cnet2e));
				## for debug, keep the intra-nodegroup edges
				# cnet2 <- cnet2e;
			}
		} else {
			# relayout
			cnet2 <- relayout_with_qfr(cnet2,
				repulse=repulse,
				constrain=use_constrain,
				spread_labels=FALSE,
				...)
		}
	}
	# final relayout
	if (isTRUE(do_final_relayout)) {
		cnet2 <- relayout_with_qfr(cnet2,
			repulse=final_repulse,
			spread_labels=spread_labels,
			...);
	} else {
		# spread
		if (isTRUE(spread_labels)) {
			cnet2 <- spread_igraph_labels(cnet2,
				...)
		}
	}
	return(cnet2)
}



#' Get Cnet nodesets as a named vector
#' 
#' Get Cnet nodesets as a named vector, with nodegroups named by node
#' 
#' This function is a simple conversion of `list` output from
#' `get_cnet_nodeset()` that returns a `character` vector alternative.
#' 
#' It is easier to answer: "What nodegroup is this node in?"
#' 
#' @returns `character` vector named by node, whose value is the
#'    name of the nodegroup. All nodes in `g` are returned in order,
#'    with `NA` values used when there is no nodegroup assigned to
#'    a particular node.
#' 
#' @family jam cnet utilities
#' 
#' @param g `igraph` Cnet object
#' @param ... additional arguments passed to `get_cnet_nodeset()`
#' 
#' @examples
#' cnet <- make_cnet_test()
#' nsv <- get_cnet_nodeset_vector(cnet)
#' nsv
#' 
#' @export
get_cnet_nodeset_vector <- function
(g,
 ...)
{
	#
	nsl <- get_cnet_nodeset(g, ...)
	nsv <- jamba::nameVector(rep(names(nsl), lengths(nsl)), unlist(nsl))
	if ("name" %in% igraph::vertex_attr_names(g)) {
		vnames <- igraph::V(g)$name;
		nsv <- nsv[vnames];
		names(nsv) <- vnames;
	}
	return(nsv)
}
