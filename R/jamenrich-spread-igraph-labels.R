#' Spread igraph node labels by angle from node center
#'
#' Spread igraph node labels by angle from node center
#'
#' This function uses the igraph vertex attribute
#' `"label.degree"`, which describes the angular offset for
#' each vertex label. (The `"label.degree"` values are
#' represented as radians, not degrees, starting at 0 for
#' right, and proceeding clockwise starting from the right,
#' down, left, top, right.)
#'
#' This function requires a network layout, which must be fixed
#' in order for the vertex labels to be properly oriented.
#' Labels are oriented opposite the most dominant angular mean
#' of edges from each network node. Typically the side of a node
#' with the fewest edges has the most space to place a label.
#' No further checks are performed for overlapping labels.
#'
#' Note that this function only modifies the other important
#' attribute `"label.dist"` when `label_min_dist`` is defined,
#' in order to enforce a minimum label distance from the center
#' of each node. There is no other logic to position small or
#' large labels to avoid overlapping labels.
#'
#' @family jam igraph layouts
#'
#' @param g igraph object
#' @param layout numeric matrix representing the x and y
#'    coordinates of each node in `g`, in the same order as `V(g)`.
#'    When `layout` is not supplied, nodes are checked for
#'    attributes `c("x", "y")` which define a fixed internal
#'    layout. When `force_layout=TRUE` these coordinates are ignored.
#'    If that is not supplied, then `layout_with_qfr()`
#'    is called along with the `repulse` argument. Subsequent
#'    coordinates are stored in `V(g)$x` and `V(g)$y` when
#'    argument `update_g_coords=TRUE`.
#' @param y_bias numeric value indicating the tendency to spread
#'    labels on the y-axis rather than symmetrically around each node.
#'    This argument elongates the circle surrounding a node into
#'    an ellipse with this ratio.
#' @param update_g_coords logical indicating whether the layout
#'    coordinates will be stored in `graph_attr(g, "layout")`.
#' @param do_reorder logical indicating whether to call
#'    `reorderIgraphNodes()` which re-distributes equivalent nodes
#'    based upon the node color(s). A node is "equivalent" to another
#'    node if both nodes have identical edges.
#' @param sortAttributes,nodeSortBy arguments passed to
#'    `reorderIgraphNodes()` when `do_reorder=TRUE`.
#' @param repulse argument passed to `layout_with_qfr()` only
#'    when `layout` is not supplied, and the layout is not stored
#'    in `c(V(g)$x, V(g)$y)`.
#' @param force_relayout logical indicating whether the `igraph` layout
#'    should be recalculated, in order to override coordinates that
#'    may be previously stored in the `igraph` object itself.
#'    Note that when `layout` is supplied, it is always used.
#' @param label_min_dist `numeric` minimum distance, default 0.5,
#'    used to ensure all labels are at least some distance from the center.
#'    Note that this argument is not used when `label_dist_l` is non-NULL.
#'    The units are roughly in units of one line height of text.
#' @param label_dist_l `list` with default to apply distance by 'nodeType'
#'    node vertex attribute: Gene distance 5, Set distance 0.
#'    When supplied as NULL, or when there is no node matched by
#'    this argument, the label distance uses `label_min_dist`.
#' @param ... additional arguments are passed to `layout_with_qfr()`
#'    when needed.
#'
#' @export
spread_igraph_labels <- function
(g,
 layout=NULL,
 y_bias=1,
 update_g_coords=TRUE,
 do_reorder=TRUE,
 sortAttributes=NULL,
 nodeSortBy=c("x", "-y"),
 repulse=3.5,
 force_relayout=FALSE,
 label_min_dist=0.5,
 label_dist_l=list(nodeType=c(Gene=2.5, Set=0)),
 verbose=FALSE,
 ...)
{
	##
	if (verbose) {
		jamba::printDebug("spread_igraph_labels(): ",
			"vcount:", igraph::vcount(g));
	}
	if (length(layout) == 0) {
		if (!force_relayout) {
			if ("layout" %in% igraph::graph_attr_names(g)) {
				if (verbose) {
					jamba::printDebug("spread_igraph_labels(): ",
						"Using ","layout"," from graph attributes.");
				}
				layout <- igraph::graph_attr(g, "layout");
			} else if (all(c("x", "y") %in% igraph::vertex_attr_names(g))) {
				if (verbose) {
					jamba::printDebug("spread_igraph_labels(): ",
						"Using ","x,y"," from vertex attributes.");
				}
				layout <- cbind(x=igraph::V(g)$x, y=V(g)$y);
			} else {
				layout <- layout_with_qfr(g,
					repulse=repulse,
					...);
			}
		} else {
			if (verbose) {
				jamba::printDebug("spread_igraph_labels(): ",
					"Calling ","layout_with_qfr()"," for node coordinates.");
			}
			layout <- layout_with_qfr(g,
				repulse=repulse,
				verbose=verbose,
				...);
		}
	} else if (is.function(layout)) {
		if (verbose) {
			jamba::printDebug("spread_igraph_labels(): ",
				"Calling ","layout()"," for node coordinates.");
		}
		layout <- layout(g);
	}
	
	if (length(rownames(layout)) == 0) {
		rownames(layout) <- igraph::V(g)$name;
	} else if (!all(rownames(layout) == igraph::V(g)$name)) {
		# ensure the rownames(layout) matches properly to V(g)$name
		if (!all(rownames(layout) %in% igraph::V(g)$name)) {
			stop("Not all V(g)$name node names are in rownames(layout).")
		}
		rowmatch <- match(igraph::V(g)$name, rownames(layout));
		layout <- layout[rowmatch, , drop=FALSE];
	}
	
	if (do_reorder) {
		if (verbose) {
			jamba::printDebug("spread_igraph_labels(): ",
				"Calling multienrichjam::reorderIgraphNodes()");
			jamba::printDebug("spread_igraph_labels(): ",
				"head(layout) before:");
			print(head(layout));
		}
		# if sortAttributes is empty, use defaults from reorderIgraphNodes()
		if (length(sortAttributes) == 0) {
			sortAttributes <- eval(formals(reorderIgraphNodes)$sortAttributes);
		}
		# apply node re-ordering step
		g <- reorderIgraphNodes(g,
			layout=layout,
			nodeSortBy=nodeSortBy,
			sortAttributes=sortAttributes,
			verbose=verbose,
			...);
		layout <- igraph::graph_attr(g, "layout");
		if (verbose) {
			jamba::printDebug("spread_igraph_labels(): ",
				"head(layout) after:");
			print(head(layout));
		}
	}
	g_angle <- jamba::nameVector(sapply(seq_len(igraph::vcount(g)), function(i){
		xy1 <- layout[i,1:2,drop=FALSE];
		xy2 <- layout[as.numeric(igraph::ego(g, nodes=i, mindist=1)[[1]]),1:2,drop=FALSE];
		if (length(xy2) == 0) {
			xy2 <- matrix(ncol=2, c(0,0));
		}
		xymean <- colMeans(xy1[rep(1, nrow(xy2)),,drop=FALSE] - xy2);
		-(xyAngle(xymean[1], xymean[2]*y_bias, directed=TRUE) + 0) %% 360
	}), igraph::V(g)$name);
	if (isTRUE(update_g_coords)) {
		g <- set_igraph_layout(g, layout=layout)
	}
	
	# store label angle (in radians)
	igraph::V(g)$label.degree <- jamba::deg2rad(g_angle);
	
	# establish minimum label distance
	igraph::V(g)$label.dist <- label_min_dist;
	if (length(label_dist_l) > 0 &&
			is.list(label_dist_l) &&
			length(unlist(label_dist_l)) > 0) {
		#
		vattrnames <- igraph::vertex_attr_names(g);
		for (iname in names(label_dist_l)) {
			if (iname %in% vattrnames) {
				vattrvalues <- igraph::vertex_attr(g, name=iname);
				lnames <- names(label_dist_l[[iname]]);
				k <- which(vattrvalues %in% lnames);
				# assign attributes to nodes with matching value
				if (length(k) > 0) {
					igraph::vertex_attr(g, name="label.dist", index=k) <- (
						label_dist_l[[iname]][vattrvalues[k]]);
				}
			}
		}
	}
	g;
}
