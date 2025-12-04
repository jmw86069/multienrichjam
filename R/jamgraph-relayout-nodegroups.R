

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
#' @param nodegroups `list` of node names, or `communities` object,
#'    passed to `communities2nodegroups()`.
#'    When NULL, it checks for supporting data in this order:
#'    1. If graph attribute 'mark.groups' is defined, it is used.
#'    2. If vertex attribute 'nodeType' exists, it calls
#'    `get_cnet_nodeset()`.
#'    3. Finally, it calls `igraph::cluster_optimal()`, hoping this
#'    method will be appropriate for the graph size.
#' @param repulse `numeric` default 3.5, passed to `relayout_with_qfr()`
#' @param fix_set_nodes `logical` default TRUE, whether to fix all
#'    nodes with nodeType=='Set' to prevent them from moving.
#' @param spread_labels `logical` default TRUE, whether to apply
#'    `spread_igraph_labels()` after the relayout iterations are complete.
#' @param add_edges `logical` default TRUE, whether to add edges within
#'    nodes of each nodegroup. The same is accomplished by setting
#'    `edge_factor=0`.
#' @param edge_factor `numeric` default 2, used as the numerator in
#'    new edge weights with equation `edge_factor/sqrt(n)` where 'n'
#'    is the number of nodes in the nodegroup.
#'    * Set `edge_factor=0` or `add_edges=FALSE` to skip this step.
#' @param do_final_relayout `logical` default NULL, whether to apply
#'    one more `relayout_with_qfr()` after each nodegroup is adjusted.
#'    It uses `final_repulse`.
#' @param final_repulse `numeric` used when `do_final_relayout` is TRUE,
#'    used as the 'repulse' argument in `relayout_with_qfr()`.
#' @param apply_by_size `logical` default TRUE, whether to apply the
#'    relayout to nodegroups ordered by size, using `byCols` to sort.
#'    * The default applies layout such that nodegroups with the most terms
#'    in `names(nodegroups)` are applied first, then largest to smallest
#'    nodegroups.
#'    * For Cnet plots, Gene nodes in nodegroups connected
#'    to the most Set nodes are adjusted first, then largest nodegroups,
#'    then sorted by nodegroup name.
#'    * It is unclear if the order is useful, future iterations of
#'    this approach may "move" other nodegroups aside first, then
#'    re-introduce each nodegroup into the layout one by one.
#'    Otherwise nodes could become "tangled" in the center, with
#'    no ideal method to optimize separation by nodegroup.
#' @param byCols `character` vector used when `apply_by_size` is TRUE.
#'    * '-num_terms': reverse-order by the number of terms in each
#'    nodegroup name, assuming comma-delimited terms.
#'    * '-num_nodes': reverse-order by the number of nodes in each nodegroup.
#'    * 'nodegroup': alphnumeric sort of the `names(nodegroups)`.
#' @param verbose `logical` whether to print verbose output.
#' @param ... additional arguments are passed to
#'    `communities2nodegroups()` for argument 'sep', or to
#'    `spread_igraph_labels()` for arguments regarding node ordering, etc.
#' 
#' @examples
#' cnet <- make_cnet_test();
#' ns <- get_cnet_nodeset(cnet)
#' # mark.groups: highlight just one nodegroup
#' # nodegroups: enables the edge_bundling to work properly
#' jam_igraph(cnet, mark.groups=ns["SetA,SetB"], nodegroups=ns,
#'    main="One nodegroup is split")
#' 
#' cnet2 <- relayout_nodegroups(cnet, nodegroups=ns["SetA,SetB"], do_final_layout=TRUE, verbose=TRUE)
#' jam_igraph(cnet2, mark.groups=ns["SetA,SetB"], nodegroups=ns,
#'    main="This nodegroup is 'encouraged' to re-group")
#' 
#' # by default, it is applied to all nodes
#' cnet3 <- relayout_nodegroups(cnet)
#' jam_igraph(cnet3, mark.groups=ns, nodegroups=ns,
#'    main="Re-layout across all nodegroups")
#' 
#' @export
relayout_nodegroups <- function
(cnet,
 nodegroups=NULL,
 repulse=3.5,
 fix_set_nodes=TRUE,
 spread_labels=TRUE,
 add_edges=TRUE,
 edge_factor=2,
 do_final_relayout=NULL,
 final_repulse=3.5,
 apply_by_size=TRUE,
 byCols=c("-num_terms",
 	"-num_nodes",
 	"nodegroup"),
 verbose=FALSE,
 ...)
{
	#
	cnet2 <- cnet;
	if (length(nodegroups) == 0) {
		# obtain from graph attributes if present
		if ("mark.groups" %in% igraph::graph_attr_names(cnet)) {
			if (verbose) {
				jamba::printDebug("relayout_nodegroups(): ",
					"Using graph attribute 'mark.groups'.");
			}
			nodegroups <- igraph::graph_attr(cnet, "mark.groups")
			nodegroups <- communities2nodegroups(nodegroups,
				...)
		} else if ("nodeType" %in% igraph::vertex_attr_names(cnet)) {
			# determine Cnet nodegroups
			if (verbose) {
				jamba::printDebug("relayout_nodegroups(): ",
					"Calling get_cnet_nodeset().");
			}
			nodegroups <- get_cnet_nodeset(cnet)
		} else {
			# Todo: Community-detection?
			if (verbose) {
				jamba::printDebug("relayout_nodegroups(): ",
					"Calling igraph::cluster_optimal().");
			}
			nodegroups <- igraph::cluster_optimal(cnet)
		}
	} else {
		# convert to list, accepting list or communities input
		if (verbose) {
			jamba::printDebug("relayout_nodegroups(): ",
				"Using nodegroups as supplied.");
		}
		nodegroups <- communities2nodegroups(nodegroups,
			...)
	}
	
	# sort by terms, nodegroup size, nodegroup name	
	if (isTRUE(apply_by_size) && length(nodegroups) > 1) {
		# sort by size
		ngnum <- lengths(strsplit(names(nodegroups), ","))
		ngdf <- data.frame(
			nodegroup=names(nodegroups),
			num_terms=ngnum,
			num_nodes=lengths(nodegroups))
		ngdf <- jamba::mixedSortDF(ngdf, byCols=byCols)
		nodegroups <- nodegroups[ngdf$nodegroup];
		if (verbose) {
			jamba::printDebug("relayout_nodegroups(): ",
				"Sorted nodegroups using by: ",
				paste0("'", byCols, "'"));
		}
	}

	# only apply final relayout if there is more than one nodegroup?
	# unless forced with do_final_relayout=TRUE
	if (length(do_final_relayout) == 0) {
		if (length(nodegroups) > 1) {
			do_final_relayout <- TRUE;
		} else {
			do_final_relayout <- FALSE;
		}
	}
	
	# Fix the set node coordinates
	# (Un-necessary for Cnet plots, they will always be fixed.)
	fixed_nodes <- NULL;
	unfixed_nodes <- NULL;
	if (!"nodeType" %in% igraph::vertex_attr_names(cnet2)) {
		fix_set_nodes <- NULL;
	}
	
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
		if (verbose) {
			jamba::printDebug("relayout_nodegroups(): ",
				"Applying to nodegroup '",
				names(nodegroups)[i], "' with members: ",
				nodegroups[[i]]);
		}
		# optionally add edges between nodes
		if (isTRUE(add_edges) && length(edge_factor) == 1 && edge_factor > 0) {
			inodes <- sort(match(move_nodes, igraph::V(cnet2)$name));
			if (length(inodes) > 1) {
				iedges <- (combn(inodes, 2));
				cnet2e <- cnet2;
				if (!"weight" %in% igraph::edge_attr_names(cnet2)) {
					igraph::E(cnet2e)$weight <- 1;
				}
				# add complete edges
				cnet2e <- igraph::add_edges(cnet2e, iedges,
					attr=list(weight=edge_factor/sqrt(length(inodes))))
					# attr=list(weight=edge_factor))
				
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
	nsl <- communities2nodegroups(nodegroups);
	nsv <- jamba::nameVector(rep(names(nsl), lengths(nsl)), unlist(nsl))
	if ("name" %in% igraph::vertex_attr_names(g)) {
		vnames <- igraph::V(g)$name;
		nsv <- nsv[vnames];
		names(nsv) <- vnames;
	}
	return(nsv)
}
