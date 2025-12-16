
# get_cnet_nodeset functions
#
# Todo:
# - DONE. make it work with unknown bipartite graphs
# - detect bipartite graphs and "do something" when it does
#   not make sense to return node groupings, for example
#   if there are no multi-node nodegroups.

#' Get Cnet node set by connected Sets
#'
#' Get Cnet node set by connected Sets, using nodeType to
#' delineate the two types of nodes in the graph.
#'
#' This function operates on a Cnet `igraph` object,
#' distinguished by node attribute `'nodeType'` with
#' values: `'Gene'` for Gene nodes, and `'Set'` for Set nodes.
#' This function returns Gene nodes that are connected
#' to the Sets given by argument `set_nodes`. It is
#' useful to identify which genes are connected to
#' Set `"A"` and `"B"` for example, in other words the
#' Gene nodes in a specific subcluster of the Cnet
#' `igraph`.
#'
#' When `set_nodes` is `NULL`, this function returns a `list`
#' containing `character` vectors, where each vector represents
#' Gene nodes that are connected to each combination of
#' Set nodes. This option is useful for obtaining all
#' possible Gene subclusters.
#'
#' This function is also called by `adjust_cnet_nodeset()` in
#' order to manipulate all nodes in a subcluster as a group,
#' for example calling `nudge_igraph_node()` on the whole set
#' of nodes.
#' 
#' Note that when input data do not contain vertex attribute 'nodeType',
#' data are still returned but cannot be filtered by 'Set' and 'Gene',
#' and so will contain all possible node-connection groupings.
#' 
#' @family jam igraph utilities
#' @family jam cnet utilities
#'
#' @param g `igraph` object specifically containing Cnet plot data,
#'    with node attribute `"nodeType"` that has values `c("Gene", "Set")`.
#' @param set_nodes `character` vector, default NULL, used to subset
#'    the returned data to nodes which are connected to these 'Set' nodes.
#'    of one or more set names, which are
#'    defined here as the names of the nodes with `nodeType="Set"`
#'    to which each node with `nodeType="Gene"` is connected.
#'    For example `set_nodes=c("A", "B")` will return all Gene nodes
#'    that are connected only to Set nodes `"A"` and `"B"`. If
#'    `set_nodes=NULL` then this function returns a `list` with
#'    all set_nodes observed.
#' @param sep `character` string used as a delimited when combining
#'    `set_nodes` into a name.
#' @param filter_set_only `logical` default TRUE, whether to exclude
#'    Set nodes, with vertex attribute nodeType='Set', from being
#'    members of nodegroups. In principle, 'Gene' nodes are typically
#'    members of Cnet nodesets, and 'Set' nodes do not usually
#'    share the same 'Gene' nodes with other 'Set' nodes, and therefore
#'    are not included.
#' @param ... additional arguments are ignored.
#'
#' @return `list` named by nodeset, with `character` vectors of
#'    node names in each nodeset. When `filter_set_only=TRUE` by default,
#'    the nodes will not contain 'Set' nodes.
#'    The `list` is named by nodes connected by nodes in each node group,
#'    using `sep` as delimiter.
#' 
#' @examples
#' cnet1 <- make_cnet_test()
#' get_cnet_nodeset(cnet1)
#' 
#' get_cnet_nodeset(cnet1, c("SetB,SetD"))
#' 
#' @export
get_cnet_nodeset <- function
(g,
 set_nodes=NULL,
 sep=",",
 filter_set_only=TRUE,
 min_size=1,
 ...)
{
	# alternative approach for speed
	# 1. assemble edgelist data.frame
	gel1 <- igraph::as_edgelist(g);
	if (length(gel1) == 0) {
		# no edges!
		return(NULL);
	}
	gel2 <- gel1[,2:1, drop=FALSE];
	gel <- data.frame(check.names=FALSE,
		unique(rbind(gel1, gel2)));
	colnames(gel) <- c("A", "B")
	
	# optionally limit output to nodes connected to nodeType="Set" nodes
	if (!"nodeType" %in% igraph::vertex_attr_names(g)) {
		filter_set_only <- FALSE;
	}
	if (TRUE %in% filter_set_only) {
		use_set_nodes <- igraph::V(g)[igraph::V(g)$nodeType %in% "Set"]$name;
		gel <- subset(gel, gel[,2] %in% use_set_nodes)
	}
	
	gel_sorted <- jamba::mixedSortDF(gel)
	
	# 2. cPaste() because unique and sorted are already accomplished above
	neighbor_list <- split(gel_sorted$B, gel_sorted$A);
	neighborG1 <- jamba::cPaste(neighbor_list,
		sep=sep)
	cnet_nodesets <- split(names(neighborG1), neighborG1)
	
	if (length(min_size) == 1 && is.numeric(min_size) && min_size > 1) {
		cnet_nodesets <- cnet_nodesets[lengths(cnet_nodesets) >= min_size];
	}
	if (length(set_nodes) > 0) {
		set_nodes_v <- jamba::cPasteS(set_nodes,
			sep=sep)
		cnet_nodesets <- cnet_nodesets[names(cnet_nodesets) %in% set_nodes_v];
	}
	return(cnet_nodesets);
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
#' @param nodegroups default is NULL, which calls `get_cnet_nodeset()`.
#'    Also accepts `list` or `communities`.
#' @param filter_set_only `logical` default TRUE, passed to
#'    `get_cnet_nodeset()` to determine whether 'Set' nodes are
#'    assigned a nodeset.
#'    * When TRUE, 'Set' nodes are filtered out,
#'    therefore their nodeset becomes NA.
#'    * When FALSE, 'Set' nodes are also assigned to a nodeset.
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
 nodegroups=NULL,
 filter_set_only=TRUE,
 ...)
{
	# fill nodegroups
	if (length(nodegroups) == 0) {
		nodegroups <- get_cnet_nodeset(g, ...)
	}
	# confirm proper format, tolerate: communities,clusters,list
	nsl <- communities2nodegroups(nodegroups);
	
	nsv <- jamba::nameVector(rep(names(nsl), lengths(nsl)), unlist(nsl))
	if ("name" %in% igraph::vertex_attr_names(g)) {
		vnames <- igraph::V(g)$name;
		nsv <- nsv[vnames];
		names(nsv) <- vnames;
	}
	return(nsv)
}
