
#' Get partite/connected graph nodesets
#'
#' Get partite/connected graph nodesets defined by shared connections
#'
#' This method is under development, the intent is to bundle
#' edges where a large subset of nodes are all connected to
#' the same node neighbors. A typical graph may not have any
#' two nodes with the same neighbors, but this situation tends
#' to happen much more often with bipartite graphs,
#' where nodes of one type are only permitted to have node
#' neighbors of the other type. It is not required for this
#' method to work, however.
#'
#' The driving scenario is with Cnet (concept network) plots,
#' which is a bipartite network with `"Gene"` and `"Set"` nodes.
#' It is fairly common to have multiple genes present in the
#' same one or few pathways. As a result, these nodes are
#' most often positioned near each other as a natural
#' by-product of having the same connected neighbor nodes.
#'
#' Identifying a nodeset with identical node neighbors enables
#' some other useful operations:
#'
#' * re-positioning, rotating, expanding, compressing the
#' whole nodeset layout to improve network graph aesthetics,
#' node label readability, reducing overlaps
#' * edge bundling to improve visual distinctiveness between
#' multiple nodesets
#'
#' @family jam igraph functions
#'
#' @param g `igraph` object that contains one attribute column with
#'    node type.
#' @param type `character` string of the node/vertex attribute that
#'    represents the node type.
#' @param set_nodes `character` or `NULL`, which contains the set
#'    of node neighbors for the requested nodeset. For example,
#'    one might want all nodes that connect with `c("A", "B", "C")`.
#'    When `set_nodes=NULL` then all nodesets are returned.
#' @param sep `character` string used as a delimiter between
#'    node names when defining a nodeset name
#' @param return_type `character` string indicating the type of
#'    data to return:
#'    * `"list"` returns a list of nodesets, each element in the `list`
#'    is a `character` vector with node names.
#'    * `"df"` returns a `data.frame` with more detailed annotation
#'    for each node, including nodesets, neighbor nodes, etc.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
get_bipartite_nodeset <- function
(g,
   type="nodeType",
   set_nodes=NULL,
   sep=",",
   return_type=c("list", "df"),
   verbose=FALSE,
   ...)
{
   return_type <- match.arg(return_type);

   ## Enforce a "name" for each vertex node
   if (!"name" %in% igraph::list.vertex.attributes(g)) {
      igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)));
   }

   ## comma-delimited neighboring nodes for each node
   neighbor_list <- lapply(seq_len(igraph::vcount(g)), function(v){
      names(igraph::neighbors(g, v, mode="all"));
   });
   names(neighbor_list) <- igraph::V(g)$name;
   neighbor_v <- jamba::cPasteS(neighbor_list,
      sep=sep);

   neighbor_df <- data.frame(stringsAsFactors=FALSE,
      node=igraph::V(g)$name,
      num_neighbors=lengths(neighbor_list),
      neighbor_group_size=as.vector(table(neighbor_v)[neighbor_v]),
      neighbors=as.character(neighbor_v),
      neighbor=as.character(neighbor_v));
   if (type %in% igraph::list.vertex.attributes(g)) {
      neighbor_df$type <- igraph::vertex_attr(g, type);
   }

   neighbor_tall <- deconcat_df2(neighbor_df,
      column="neighbor",
      split=",")
   neighbor_tall$edge <- jamba::cPasteS(
      strsplit(
         jamba::pasteByRow(neighbor_tall[,c("node", "neighbor"),drop=FALSE],
            sep=":!:"),
         split=":!:"));
   neighbor_tall <- jamba::mixedSortDF(neighbor_tall,
      byCols=c(
         "edge",
         "-neighbor_group_size",
         "num_neighbors",
         "neighbors"));

   # get unique edges
   umatch <- match(unique(neighbor_tall$edge),
      neighbor_tall$edge);
   neighbor_tall_unique <- jamba::mixedSortDF(
      neighbor_tall[umatch,,drop=FALSE],
      byCols=c("neighbors"));
   if ("df" %in% return_type) {
      return(neighbor_tall_unique);
   }

   neighbor_tall_unique2 <- unique(neighbor_tall_unique[,c("node", "neighbors")]);
   nodesets <- split(neighbor_tall_unique2$node,
      neighbor_tall_unique2$neighbors);
   # add singlet nodes to their own solo group
   missing_nodes <- setdiff(neighbor_df$node,
      unlist(nodesets));
   if (length(missing_nodes) > 0) {
      names(missing_nodes) <- paste0("singlet_", missing_nodes);
      nodesets <- c(nodesets,
         lapply(missing_nodes, function(i){i}));
   }
   if (length(set_nodes) == 0) {
      if (verbose) {
         jamba::printDebug("get_bipartite_nodeset(): ",
            "Returning all nodesets.");
      }
      return(nodesets);
   }
   set_nodes_v <- jamba::cPasteSU(set_nodes,
      sep=sep);
   if (all(set_nodes_v %in% names(nodesets))) {
      if (verbose) {
         jamba::printDebug("get_bipartite_nodeset(): ",
            "Returning nodesets that match set_nodes: ",
            set_nodes_v);
      }
      return(nodesets[unique(match(set_nodes_v, names(nodesets)))]);
   }
   if (verbose) {
      jamba::printDebug("get_bipartite_nodeset(): ",
         "Returning nodesets that contain: ",
         set_nodes_v);
   }
   use_set_nodes <- ifelse(set_nodes_v %in% names(nodesets),
      set_nodes_v,
      neighbor_tall_unique2[match(set_nodes_v, neighbor_tall_unique2$node),"neighbors"]);
   return(nodesets[unique(match(use_set_nodes, names(nodesets)))]);
}


#' Bundle edges in a bipartite graph
#'
#'
#' @export
edge_bundle_bipartite <- function
(g,
   type="nodeType",
   ...)
{
   # get bipartite groups
   neighbor_tall_unique <- get_bipartite_nodeset(g,
      type=type)

   # define nodegroups
   nodegroups <- get_bipartite_nodeset(g,
      type=type,
      return_type="list");
   names(nodegroups) <- paste0("to_",
      names(nodegroups));
   # add singlet nodes
   singlet_nodes <- setdiff(igraph::V(g)$name,
      unlist(nodegroups));
   if (length(singlet_nodes) > 0) {
      nodegroups <- c(nodegroups,
         split(singlet_nodes, singlet_nodes));
   }
   names(nodegroups) <- jamba::makeNames(names(nodegroups));
   return(nodegroups);
}


#' Bundle edges using node groups
#'
#' Bundle edges using node groups
#'
#' This edge bundling technique relies upon some form of
#' node grouping, usually derived from network community
#' detection, or from bipartite nodesets (see
#' `get_bipartite_nodeset()` for details.)
#'
#' Given a set of node groups, edges are bundled entering
#' and exiting each node group, along the linear path between
#' the two node group center positions, using a spline
#' function and intermediate control points.
#'
#' The default spline uses the initial node positions, and the
#' midpoint along the line between the two respective node groups.
#' The midpoints can be adjusted with the argument `midpoint`
#' and a vector of one or more fractional positions between `0` and `1`.
#' A useful alternative is `midpoint=c(0.3, 0.7)` which adds
#' two control points along the linear path between node group
#' centers, and tends to make the edges bundle closer together
#' for a longer distance.
#'
#' When used with bipartite nodesets, edges are bundled between
#' each nodeset and individual nodes. The edge bundling rules are
#' the same, with the default `midpoint=0.5` being exactly half
#' the distance between the nodeset center, and the single node.
#' In this case, the `midpoint` is directional, always pointing
#' from the nodeset to the single node, therefore can be adjusted
#' closer to the nodeset center with `midpoint=0.2` or closer to
#' the single node with `midpoint=0.8`.
#'
#' @family jam igraph functions
#'
#' @examples
#' # using community detection
#' karate <- igraph::make_graph("Zachary")
#' igraph::V(karate)$name <- as.character(seq_len(igraph::vcount(karate)))
#'
#' # run any igraph::cluster_*()
#' wc <- igraph::cluster_louvain(karate)
#' # define list
#' nodegroups_wc <- split(igraph::V(karate)$name, wc$membership)
#'
#' # bonus points for colorizing nodes and edges by community
#' igraph::V(karate)$color <- colorjam::group2colors(igraph::membership(wc));
#' igraph::V(karate)$label.color <- jamba::setTextContrastColor(igraph::V(karate)$color);
#' igraph::V(karate)$frame.color <- jamba::makeColorDarker(igraph::V(karate)$color);
#' karate <- color_edges_by_nodes(karate);
#'
#' # update graph layout
#' layout_xy <- igraph::layout_with_graphopt(karate);
#' igraph::graph_attr(karate, "layout") <- layout_xy;
#'
#' jam_igraph(karate,
#'    edge_bundling="nodegroups",
#'    nodegroups=nodegroups_wc,
#'    use_shadowText=TRUE);
#'
#' @param g `igraph` that contains layout coordinates in
#'    graph attributes, stored as `igraph::graph_attr(g, "layout")`.
#' @param nodegroups `list` of node names, or object with
#'    class `"communities"` as produced by `igraph::cluster_*`
#'    methods such as `igraph::cluster_walktrap()`. Note that
#'    every node must be represented.
#' @param midpoint `numeric` vector of one or two values ranging
#'    from `0` to `1` that define control point positions between
#'    each node group center, used to create each edge spline.
#' @param detail `integer` number of intermediate points along
#'    the spline to render for each edge.
#' @param draw_lines `logical` indicating whether to render the edge
#'    splines after calculating them.
#' @param ... additional arguments are ignored.
#'
#' @return `data.frame` with each edge spline point represented
#'    on its own row, with breaks in edges defined by `NA` coordinates.
#'
#' @export
edge_bundle_nodegroups <- function
(g,
   nodegroups,
   midpoint=0.5,
   detail=10,
   draw_lines=TRUE,
   ...)
{
   # node layout
   layout_xy <- igraph::graph_attr(g, "layout");
   colnames(layout_xy)[1:2] <- c("x", "y");
   if (!"name" %in% igraph::list.vertex.attributes(g)) {
      igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)));
   }
   rownames(layout_xy) <- igraph::V(g)$name;

   # accept class "communities"
   if ("communities" %in% class(nodegroups)) {
      nodegroups <- split(igraph::V(g)$name,
         igraph::membership(nodegroups))
   }

   # nodegroup_df
   # todo: require that every node is contained in a nodegroup
   nodegroup_df <- data.frame(check.names=FALSE,
      stringsAsFactors=FALSE,
      node=unlist(nodegroups),
      nodegroup=rep(names(nodegroups),
         lengths(nodegroups)))

   # center point for each nodegroup
   # todo: replace with "node hull" logic, with center inside the hull
   nodegroup_centers <- lapply(nodegroups, function(i){
      colMeans(layout_xy[i,c("x", "y"),drop=FALSE])
   });
   nodegroup_centers_df <- data.frame(check.names=FALSE,
      stringsAsFactors=FALSE,
      jamba::rbindList(nodegroup_centers));
   nodegroup_centers_df$nodegroup <- names(nodegroups);

   # get edge data.frame
   edge_df <- data.frame(check.names=FALSE,
      stringsAsFactors=FALSE,
      igraph::as_edgelist(g));
   edge_df$nodegroup1 <- nodegroup_df$nodegroup[match(edge_df[,1], nodegroup_df$node)];
   edge_df$nodegroup2 <- nodegroup_df$nodegroup[match(edge_df[,2], nodegroup_df$node)];
   edge_df$nodegroup1_2 <- jamba::pasteByRow(
      edge_df[,c("nodegroup1", "nodegroup2"),
         drop=FALSE]);

   # midpoint between node and nodeset_centers
   midpoint <- jamba::noiseFloor(midpoint,
      minimum=0,
      ceiling=1);
   midpoint_df <- unique(edge_df[,c("nodegroup1", "nodegroup2"), drop=FALSE]);
   midpoint_df$x1 <- nodegroup_centers_df$x[match(midpoint_df$nodegroup1, nodegroup_centers_df$nodegroup)]
   midpoint_df$y1 <- nodegroup_centers_df$y[match(midpoint_df$nodegroup1, nodegroup_centers_df$nodegroup)]
   midpoint_df$x3 <- nodegroup_centers_df$x[match(midpoint_df$nodegroup2, nodegroup_centers_df$nodegroup)]
   midpoint_df$y3 <- nodegroup_centers_df$y[match(midpoint_df$nodegroup2, nodegroup_centers_df$nodegroup)]
   xmids <- do.call(cbind, lapply(midpoint, function(i){
      midpoint_df$x1 * (1 - i) + midpoint_df$x3 * (i)
   }))
   colnames(xmids) <- jamba::makeNames(rep("x2", length(midpoint)),
      suffix="_");
   ymids <- do.call(cbind, lapply(midpoint, function(i){
      midpoint_df$y1 * (1 - i) + midpoint_df$y3 * (i)
   }))
   colnames(ymids) <- jamba::makeNames(rep("y2", length(midpoint)),
      suffix="_");
   midpoints_df <- data.frame(check.names=FALSE,
      stringsAsFactors=FALSE,
      midpoint_df[,c("nodegroup1", "nodegroup2", "x1"), drop=FALSE],
      xmids,
      midpoint_df[,c("x3", "y1"), drop=FALSE],
      ymids,
      midpoint_df[,c("y3"), drop=FALSE]);
   midpoints_df$nodegroup1_2 <- jamba::pasteByRow(
      midpoints_df[,c("nodegroup1", "nodegroup2"),
         drop=FALSE]);
   xcols <- colnames(xmids);
   ycols <- colnames(ymids);

   # control points for each edge spline
   edge_df$x1 <- layout_xy[match(edge_df[,1], rownames(layout_xy)),"x"];
   ematch <- match(edge_df$nodegroup1_2, midpoints_df$nodegroup1_2);
   edge_df[,xcols] <- midpoints_df[ematch, xcols, drop=FALSE];
   edge_df$x3 <- layout_xy[match(edge_df[,2], rownames(layout_xy)),"x"];
   edge_df$y1 <- layout_xy[match(edge_df[,1], rownames(layout_xy)), "y"];
   edge_df[,ycols] <- midpoints_df[ematch, ycols, drop=FALSE];
   edge_df$y3 <- layout_xy[match(edge_df[,2], rownames(layout_xy)), "y"];
   use_xcols <- jamba::vigrep("^x[123]", colnames(edge_df));
   use_ycols <- jamba::vigrep("^y[123]", colnames(edge_df));

   # calculate each spline
   edge_splines <- lapply(seq_len(nrow(edge_df)), function(n){
      x1 <- unlist(edge_df[n, use_xcols]);
      y1 <- unlist(edge_df[n, use_ycols]);
      path_xy <- rbind(
         ggforce:::bezierPath(
            x=x1,
            y=y1,
            detail=detail),
         c(NA, NA));
      path_df <- data.frame(stringsAsFactors=FALSE,
         path_xy,
         edge_row=n);
   })
   edge_spline_df <- jamba::rbindList(edge_splines);

   edge_attr_names <- igraph::list.edge.attributes(g);
   if ("color" %in% edge_attr_names) {
      edge_spline_df$color <- igraph::E(g)$color[edge_spline_df$edge_row];
   } else {
      edge_spline_df$color <- "grey60";
   }
   if ("width" %in% edge_attr_names) {
      edge_spline_df$width <- igraph::E(g)$width[edge_spline_df$edge_row];
   } else {
      edge_spline_df$width <- 1;
   }

   if (draw_lines) {
      # note that lines() is not vectorized for multiple col, lwd, lty values
      # and should be re-run for each unique set
      #print(head(edge_spline_df, 30));
      row_split <- split(seq_len(nrow(edge_spline_df)),
         jamba::pasteByRow(edge_spline_df[,c("color", "width"),drop=FALSE]));
      for (i in row_split) {
         lines(edge_spline_df[i,1:2,drop=FALSE],
            col=edge_spline_df$color[i],
            lwd=edge_spline_df$width[i],
            ...);
      }
   }
   return(invisible(edge_spline_df));
}
