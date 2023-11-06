
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
#' Bundle edges in a bipartite graph
#'
#' This function performs edge bundling in bipartite network graphs,
#' which are expected to contain two classes of nodes. In general
#' this situation lends itself well to bundling edges by shared
#' connections, where a subset of nodes of one class all bind to
#' the same set of nodes in the other class. These nodes are
#' typically co-located in the network layout, which works well
#' with this style of edge bundling.
#'
#' @family jam igraph functions
#'
#' @param g `igraph` object
#' @param type `character` name of vertex attribute that defines the
#'    node type.
#' @param ... additional arguments are ignored.
#'
#' @export
edge_bundle_bipartite <- function
(g,
 type="nodeType",
 ...)
{
   # get bipartite groups
   # neighbor_tall_unique <- get_bipartite_nodeset(g,
   #    type=type)

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
#' the same, with the default `midpoint=c(0.4, 0.6)` being centered at half
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
#' @param shape `character` (optional) used to override the `vertex.shape`
#'    passed in `params`. It is recycled to the number of nodes,
#'    for example by `igraph::vcount(g)`.
#' @param params `function` representing `igraph` plotting parameters
#'    used at rendering time. The output is also produced by
#'    `parse_igraph_plot_params()` for use in `jam_igraph()`
#'    plotting, and is passed to other node and edge rendering
#'    functions.
#' @param midpoint `numeric` vector of one or more values ranging
#'    from `0` to `1` that define control point positions along the
#'    line between two nodegroup center coordinates. When one nodegroup
#'    contains only one node, this line segment is shortened to end
#'    at that node border after clipping the corresponding node shape.
#'    The position along the line is defined relative to the first node
#'    in the edge, toward the second node in the edge.
#'    Using `midpoint=0.5` guarantees the control point is the exact middle,
#'    while `midpoint=c(0.2, 0.8)` will use two control points at 20% and
#'    80% distance along the line segment resulting in an edge that more
#'    closely follows the line segment.
#' @param detail `integer` number of intermediate points along
#'    the spline to render for each edge.
#' @param draw_lines `logical` indicating whether to render the edge
#'    splines after calculating them.
#' @param nodegroup_midpoints `list` experimental support for defining
#'    specific control points used by bundled edges. Not fully implemented
#'    as yet. In future, it will require two nodegroups to be defined
#'    for each set of control point coordinates, with no requirement
#'    for the location of control points.
#' @param linear_cor_threshold `numeric` value between 0 and 1.
#'    Coordinates for each edge, and intermediate control point
#'    coordinates are used in `xspline()` to create a curved spline
#'    from node to node. However, when the nodes and intermediate
#'    control points are already linear, the edge will be treated
#'    as a linear edge. To test for linearity, `cor()` correlation
#'    is calculated, and values at or above `linear_cor_threshold`
#'    are considered linear.
#'
#'    The driving problem is when the control point is colinear with
#'    two nodes, and the control point is positioned outside the two
#'    nodes. Without this modification, the line would appear to pass
#'    from one node beyond the other node, with an arrow (if directed)
#'    pointing back to the other node from the opposite direction.
#' @param bundle_style `character` string describing the type of curvature
#'    to use for edge bundles:
#'    * `"bezier"`: (default) calls `bezier::bezier()` to define a bezier
#'    curve using the edge control points.
#'    * `"xspline"`: calls `graphics::xspline()` to define an XSpline curve
#'    using the edge control points, however the method is  customized
#'    to include each edge endpoint twice, which makes the intermediate
#'    curve much rounder than normal.
#'    * `"angular"`: calls `graphics::xspline()` to define an XSpline curve
#'    using the edge control points. This shape tends to appear angular,
#'    thus the name.
#'    * `"bezierPath"`: calls `ggforce:::bezierPath()` when `ggforce` is
#'    available, producing a bezier curve using the edge control points.
#'    Note this method appears identical to `"bezier"` above, and will
#'    likely be removed in a future release.
#'    * `"subway"`: experimental method that uses the `"angular"` appearance,
#'    with more repeated intermediate control points intended to group
#'    all bundled edges to the same linear segment. The intent is to "dodge"
#'    edges along the line segment, similar to the appearance of subway maps,
#'    however it is not fully implemented.
#' @param bundle_self `logical` to indicate whether edges that begin and
#'    end in the same nodegroup should be bundled through the nodegroup
#'    center coordinate.
#'    * `bundle_self=FALSE` forces all edges within a nodegroup to be
#'    rendered as straight lines, therefore not using the nodegroup center
#'    as the control point.
#'    * `bundle_self=TRUE` overrides the validation check that
#'    requires the distance between center points of two nodegroups to
#'    have distance at least 0.5% the layout coordinate span. It can
#'    be a visual aid to have connections bundle through the center
#'    of the nodegroup, especially when the nodegroup is almost fully
#'    connected.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param debug `logical` indicating whether to plot debug output that
#'    may be helpful in understanding the edge bundle control points.
#'    To specify debug only for edge bundling, use the substring "bundl",
#'    for example `options("debug"="bundling")`.
#' @param ... additional arguments are ignored.
#'
#' @return `data.frame` with each edge spline point represented
#'    on its own row, with breaks in edges defined by `NA` coordinates.
#'
#' @export
edge_bundle_nodegroups <- function
(g,
 nodegroups,
 shape=NULL,
 params=NULL,
 midpoint=0.5,
 detail=10,
 draw_lines=TRUE,
 nodegroup_midpoints=NULL,
 linear_cor_threshold=1,
 bundle_style=getOption("jam.bundle_style", "bezier"),
 bundle_self=FALSE,
 verbose=FALSE,
 debug=getOption("debug", FALSE),
 ...)
{
   # Todo:
   # - 1. DONE: skip edge bundling when two node centers are "identical"
   #         (no control line)
   # - 2. DONE: skip edge bundling when nodes and control points are co-linear
   #         (no curve)
   # - 3. DONE: clip bundled edges using control points instead of linear edge
   #         (edge begins at node boundary point facing nearest control point)
   # - 4. DONE: mark bundling valid=FALSE when both nodegroups are singlets
   #         (nothing to bundle)
   # - 5. DONE: when one nodegroup has one node, clip the nodegroup line
   #         using the node border, such that midpoint is calculated from there
   #         (midpoint 0.99 should not be inside a node boundary)
   # - 6. DONE: when one nodegroup has one node, mark valid=FALSE when
   #         one nodegroup center is inside the other single-node-nodegroup border.
   #         This test could be if the from_to distance is shorted than from-clipped_to
   #         which would indicate the clipped_to is farther away.
   # - 7. For scenario 2 above, when control points are not perfectly
   #         co-linear, consider only invalidating the bundling if the
   #         control point goes past either node. If the control point lies
   #         inside the two nodes, leave it as-is -- problem happens when
   #         one edge in a bundle is not included in the group. Unsure how
   #         best to handle this exception, it is similar to scenario 6.
   #         When control point is not inside the final node edge, then
   #         valid=TRUE when the co-linearity test is not perfect.

   # bundle_style
   if (length(bundle_style) == 0) {
      bundle_style <- "xspline";
   }
   bundle_style <- head(bundle_style, 1);
   bundle_styles <- c("xspline",
      "bezierPath",
      "bezier",
      "angular",
      "subway");
   if (!any(bundle_style %in% bundle_styles)) {
      bundle_style <- "xspline";
   }
   # node layout
   layout_xy <- igraph::graph_attr(g, "layout");
   if (length(layout_xy) == 0 || nrow(layout_xy) == 0) {
      stop('Edge bundling requires layout in graph_attr(g, "layout"), none was found.');
   }
   colnames(layout_xy)[1:2] <- c("x", "y");
   vct <- igraph::vcount(g);
   if (!"name" %in% igraph::list.vertex.attributes(g)) {
      igraph::V(g)$name <- as.character(seq_len(vct));
   }
   rownames(layout_xy) <- igraph::V(g)$name;

   # define params() if not provided
   if (length(params) == 0) {
      params <- parse_igraph_plot_params(g, list(...));
   }
   # arrow mode if relevant
   arrow.mode <- params("edge", "arrow.mode")
   arrow.mode <- get_igraph_arrow_mode(g, arrow.mode)
   arrow.size <- params("edge", "arrow.size")
   arrow.width <- params("edge", "arrow.width")
   elen <- igraph::ecount(g);
   arrow.mode <- rep(arrow.mode, length.out=elen);
   arrow.size <- rep(arrow.size, length.out=elen);
   arrow.width <- rep(arrow.width, length.out=elen);

   # edge label (if present)
   edge.label <- params("edge", "label")
   edge.label.color <- params("edge", "label.color")
   edge.label.family <- params("edge", "label.family")
   edge.label.font <- params("edge", "label.font")
   edge.label.cex <- params("edge", "label.cex")

   # validate shape
   if (length(shape) == 0) {
      shape <- params("vertex", "shape");
   }
   shape <- rep(shape,
      length.out=vct);
   if (verbose) {
      jamba::printDebug("edge_bundle_nodegroups(): ",
         "shape: ", head(shape, 10));
   }

   # accept class "communities"
   if ("communities" %in% class(nodegroups)) {
      nodegroups <- split(igraph::V(g)$name,
         igraph::membership(nodegroups))
   }

   # nodegroup_df
   # complete: require that every node is contained in a nodegroup
   if (length(names(nodegroups)) == 0) {
      names(nodegroups) <- paste0("nodegroup_",
         jamba::colNum2excelName(seq_along(nodegroups)))
   }
   if (!all(igraph::V(g)$name %in% unlist(nodegroups))) {
      add_nodes <- setdiff(igraph::V(g)$name,
         unlist(nodegroups));
      add_nodegroups <- as.list(jamba::nameVector(
         add_nodes));
      nodegroups <- c(nodegroups, add_nodegroups)
   }
   nodegroup_df <- data.frame(check.names=FALSE,
      stringsAsFactors=FALSE,
      node=unlist(nodegroups),
      nodegroup=rep(names(nodegroups),
         lengths(nodegroups)))

   # center point for each nodegroup
   # todo: replace with "node hull" logic, with center inside the hull
   if (length(nodegroup_midpoints) > 0) {
      # TODO: implement
   }
   nodegroup_centers <- lapply(nodegroups, function(i){
      colMeans(layout_xy[i,c("x", "y"),drop=FALSE])
   });
   nodegroup_centers_df <- data.frame(check.names=FALSE,
      stringsAsFactors=FALSE,
      jamba::rbindList(nodegroup_centers));
   nodegroup_centers_df$nodegroup <- names(nodegroups);

   # get edge data.frame
   el <- igraph::as_edgelist(g, names=FALSE);
   edge_df <- data.frame(check.names=FALSE,
      stringsAsFactors=FALSE,
      igraph::as_edgelist(g));
   edge_df$nodegroup1 <- nodegroup_df$nodegroup[match(edge_df[,1], nodegroup_df$node)];
   edge_df$nodegroup2 <- nodegroup_df$nodegroup[match(edge_df[,2], nodegroup_df$node)];
   edge_df$nodegroup1_2 <- jamba::pasteByRow(
      edge_df[,c("nodegroup1", "nodegroup2"),
         drop=FALSE]);
   if (length(edge.label) == 0) {
      edge.label <- ""
   }
   if (length(edge.label) < nrow(edge_df)) {
      edge.label <- rep(edge.label,
         length.out=nrow(edge_df))
   }
   edge_df$label <- edge.label;

   # define edge.coords
   edge.coords <- matrix(0,
      nrow=nrow(edge_df),
      ncol=4)
   edge.coords[, 1] <- layout_xy[match(edge_df[, 1], rownames(layout_xy)), 1]
   edge.coords[, 2] <- layout_xy[match(edge_df[, 1], rownames(layout_xy)), 2]
   edge.coords[, 3] <- layout_xy[match(edge_df[, 2], rownames(layout_xy)), 1]
   edge.coords[, 4] <- layout_xy[match(edge_df[, 2], rownames(layout_xy)), 2]
   layout_scale <- max(c(diff(range(layout_xy[,1])),
      diff(range(layout_xy[,2]))));

   # midpoint defined by nodegroup centers
   # which defines a straight line from one nodegroup center to another
   midpoint <- jamba::noiseFloor(midpoint,
      minimum=0,
      ceiling=1);
   midpoint_df <- unique(edge_df[,c("nodegroup1", "nodegroup2"), drop=FALSE]);
   match1 <- match(midpoint_df$nodegroup1, nodegroup_centers_df$nodegroup);
   match2 <- match(midpoint_df$nodegroup2, nodegroup_centers_df$nodegroup);
   midpoint_df$x1 <- nodegroup_centers_df$x[match1]
   midpoint_df$y1 <- nodegroup_centers_df$y[match1]
   midpoint_df$x3 <- nodegroup_centers_df$x[match2]
   midpoint_df$y3 <- nodegroup_centers_df$y[match2]

   # first assert that distance between nodegroup centers should be non-zero
   # (and add a little buffer based upon overall plot scale)
   # must be at least greater than 0.5% the layout width
   dist13 <- sqrt((midpoint_df$x1 - midpoint_df$x3)^2 + (midpoint_df$y1 - midpoint_df$y3)^2)
   if (verbose) {
      jamba::printDebug("edge_bundle_nodegroups(): ",
         "nodeset from-to distance:", dist13,
         ", layout_scale:", layout_scale,
         ", layout_scale*0.005:", layout_scale * 0.005,
         ", result=", ifelse(dist13 >= (layout_scale * 0.005),
            "valid", "not valid"))
      print(midpoint_df);
   }
   midpoint_df$valid <- (dist13 >= (layout_scale * 0.005) |
      midpoint_df$nodegroup1 == midpoint_df$nodegroup2)

   # optionally do not bundle edges in the same nodegroup
   if (!TRUE %in% bundle_self) {
      same_nodegroup <- (midpoint_df$nodegroup1 == midpoint_df$nodegroup2)
      if (any(same_nodegroup)) {
         midpoint_df$valid[same_nodegroup] <- FALSE;
      }
   }

   # when one nodegroup is a single node
   # clip this line to the edge of that node boundary
   # define singlet nodegroups
   singlet_nodegroups <- names(nodegroups)[lengths(nodegroups) == 1];
   # mark entries invalid if both nodegroups are singlets
   both_singlets <- (midpoint_df$nodegroup1 %in% singlet_nodegroups &
      midpoint_df$nodegroup2 %in% singlet_nodegroups)
   if (any(both_singlets)) {
      midpoint_df$valid[both_singlets] <- FALSE;
   }
   if (length(singlet_nodegroups) > 0) {
      doclip1 <- (midpoint_df$valid %in% TRUE &
         midpoint_df$nodegroup1 %in% singlet_nodegroups);
      doclip2 <- (midpoint_df$valid %in% TRUE &
            midpoint_df$nodegroup2 %in% singlet_nodegroups);
      if (length(debug) > 0 &&
            (TRUE %in% debug || any(grepl("bundl", debug)) ) ) {
         jamba::printDebug("edge_bundle_nodegroups(): ",
            "doclip1:", doclip1);
         jamba::printDebug("edge_bundle_nodegroups(): ",
            "doclip2:", doclip2);
         jamba::printDebug("edge_bundle_nodegroups(): ",
            "singlet_nodegroups:");
         print(singlet_nodegroups);
         jamba::printDebug("edge_bundle_nodegroups(): ",
            "shape:");
         print(shape);
         jamba::printDebug("edge_bundle_nodegroups(): ",
            "midpoint_df:");
         print(midpoint_df);
      }
      if (any(doclip2)) {
         # must supply coords for only the nodes required
         # and make edgelist consistent with those nodes
         midclip_ec <- as.matrix(
            midpoint_df[doclip2, c("x1", "y1", "x3", "y3"), drop=FALSE]);
         # calculate distance from x1,y1 to x3,y3
         midclip_ec_dist <- sqrt(
            (midpoint_df$x1[doclip2] - midpoint_df$x3[doclip2])^2 +
            (midpoint_df$y1[doclip2] - midpoint_df$y3[doclip2])^2);
         #
         singlets <- unname(unlist(nodegroups[midpoint_df$nodegroup2[doclip2]]));
         singlets_num <- match(singlets, igraph::V(g)$name);
         midclip_el <- cbind(from=singlets_num,
            to=singlets_num);
         if (length(debug) > 0 &&
               (TRUE %in% debug || any(grepl("bundl", debug)) ) ) {
            jamba::printDebug("edge_bundle_nodegroups(): ",
               "midclip_ec2:");
            print(midclip_ec);
            jamba::printDebug("edge_bundle_nodegroups(): ",
               "midclip_el2:");
            print(midclip_el);
         }
         # iterate each unique shape
         clip_xy <- jamba::rbindList(lapply(seq_len(sum(doclip2)), function(ix){
            igraph::shapes(shape[midclip_el[ix, 2]])$clip(
               midclip_ec[ix, , drop=FALSE],
               midclip_el[ix, , drop=FALSE],
               params=params,
               end="to")
         }))
         if (length(debug) > 0 &&
               (TRUE %in% debug || any(grepl("bundl", debug)) ) ) {
            jamba::printDebug("edge_bundle_nodegroups(): ",
               "clip_xy2:");
            print(clip_xy);
         }
         newclip_ec_dist <- sqrt(
            (clip_xy[,1] - midclip_ec[,"x3"])^2 +
            (clip_xy[,2] - midclip_ec[,"y3"])^2);
         doclip2_invalid <- (newclip_ec_dist > midclip_ec_dist);
         if (any(doclip2_invalid)) {
            midpoint_df$valid[doclip2 & doclip2_invalid] <- FALSE;
         }
      }
      if (any(doclip1)) {
         # must supply coords for only the nodes required
         # and make edgelist consistent with those nodes
         midclip_ec <- as.matrix(
            midpoint_df[doclip1, c("x1", "y1", "x3", "y3"), drop=FALSE])
         # calculate distance from x1,y1 to x3,y3
         midclip_ec_dist <- sqrt(
            (midpoint_df$x1[doclip1] - midpoint_df$x3[doclip1])^2 +
            (midpoint_df$y1[doclip1] - midpoint_df$y3[doclip1])^2);
         #
         singlets <- unname(unlist(nodegroups[midpoint_df$nodegroup1[doclip1]]));
         singlets_num <- match(singlets, igraph::V(g)$name);
         midclip_el <- cbind(from=singlets_num,
            to=singlets_num);
         if (length(debug) > 0 &&
               (TRUE %in% debug || any(grepl("bundl", debug)) ) ) {
            jamba::printDebug("edge_bundle_nodegroups(): ",
               "midclip_ec1:");
            print(midclip_ec);
            jamba::printDebug("edge_bundle_nodegroups(): ",
               "midclip_el1:");
            print(midclip_el);
         }
         # iterate each unique shape
         clip_xy <- jamba::rbindList(lapply(seq_len(sum(doclip1)), function(ix){
            igraph::shapes(shape[midclip_el[ix, 1]])$clip(
               midclip_ec[ix, , drop=FALSE],
               midclip_el[ix, , drop=FALSE],
               params=params,
               end="from")
         }))
         if (length(debug) > 0 &&
               (TRUE %in% debug || any(grepl("bundl", debug)) ) ) {
            jamba::printDebug("edge_bundle_nodegroups(): ",
               "clip_xy1:");
            print(clip_xy);
         }
         midpoint_df[doclip1, c("x1", "y1")] <- clip_xy;
         newclip_ec_dist <- sqrt(
            (clip_xy[,1] - midclip_ec[,"x1"])^2 +
            (clip_xy[,2] - midclip_ec[,"y1"])^2);
         doclip1_invalid <- (newclip_ec_dist >= midclip_ec_dist);
         if (any(doclip1_invalid)) {
            midpoint_df$valid[doclip1 & doclip1_invalid] <- FALSE;
         }
      }
   }
   if (length(debug) > 0 &&
         (TRUE %in% debug || any(grepl("bundl", debug)) ) ) {
      jamba::printDebug("edge_bundle_nodegroups(): ",
         "midpoint_df:");
      print(midpoint_df);
   }
   # calculate each midpoint position along the line as a fraction
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
      midpoint_df[,c("y3", "valid"), drop=FALSE]);
   midpoints_df$nodegroup1_2 <- jamba::pasteByRow(
      midpoints_df[,c("nodegroup1", "nodegroup2"),
         drop=FALSE]);
   xcols <- colnames(xmids);
   ycols <- colnames(ymids);

   # Validate edgepoint-controlpoints-edgepoint
   # to require absolute correlation less than linear_cor_threshold
   # otherwise it is handled like a straight line
   if (any(midpoints_df$valid %in% TRUE)) {
      # only test entries that are still valid
      x2_colnames <- jamba::vigrep("^x2", colnames(midpoints_df))
      y2_colnames <- gsub("^x2", "y2", x2_colnames);
      edgetest_cor_values <- sapply(seq_len(nrow(el)), function(ix1) {
         ematch <- match(edge_df$nodegroup1_2[ix1], midpoints_df$nodegroup1_2);
         etest_xy <- rbind(
            edge.coords[ix1, 1:2, drop=FALSE],
            cbind(x=unlist(midpoints_df[ematch, x2_colnames]),
               y=unlist(midpoints_df[ematch, y2_colnames])),
            edge.coords[ix1, 3:4, drop=FALSE])
         # jamba::printDebug("etest_xy:");print(etest_xy);
         # test by correlation
         cor_xy <- cor(etest_xy[,1], etest_xy[,2]);
         ifelse(is.na(cor_xy), 1, cor_xy)
      })
      if (length(debug) > 0 &&
            (TRUE %in% debug || any(grepl("bundl", debug)) ) ) {
         jamba::printDebug("edge_bundle_nodegroups(): ",
            "edgetest_cor_values: ", edgetest_cor_values)
      }
      edgetest_cor_valid <- abs(edgetest_cor_values) < linear_cor_threshold;
      # jamba::printDebug("edge_bundle_nodegroups(): ",
      #    "edgetest_cor_valid: ", edgetest_cor_valid)
   } else {
      edgetest_cor_valid <- rep(FALSE, nrow(el))
   }

   # Determine edge clipping using control points
   # 1. node1 and first control point
   # 2. node2 and last control point
   if (length(shape) > 0) {
      xcp1 <- head(jamba::vigrep("^x2", colnames(midpoints_df)), 1);
      ycp1 <- gsub("^x2", "y2", xcp1);
      xcp3 <- tail(jamba::vigrep("^x2", colnames(midpoints_df)), 1);
      ycp3 <- gsub("^x2", "y2", xcp3);
      xcp1invalid <- "x3";
      ycp1invalid <- "y3";
      xcp3invalid <- "x1";
      ycp3invalid <- "y1";

      # apply clipping only when all shapes exist in igraph framework
      # Todo: rewrite for vectorized calculations by "from" shape, "to" shape
      valid_shapes <- igraph::shapes();
      if (all(shape %in% valid_shapes)) {
         ## clip edge coordinates using node shape functions
         shape <- rep(shape, length=igraph::vcount(g))
         ec <- edge.coords;
         ec[, 1:2] <- t(sapply(seq_len(nrow(el)), function(ix1) {
            ematch <- match(edge_df$nodegroup1_2[ix1], midpoints_df$nodegroup1_2);
            use_ec <- edge.coords[ix1, 1:4, drop=FALSE];
            use_ec[,3] <- ifelse(midpoints_df$valid[ematch] %in% TRUE &
                  edgetest_cor_valid[ix1] %in% TRUE,
               midpoints_df[ematch, xcp1],
               use_ec[,3])
            use_ec[,4] <- ifelse(midpoints_df$valid[ematch] %in% TRUE &
                  edgetest_cor_valid[ix1] %in% TRUE,
               midpoints_df[ematch, ycp1],
               use_ec[,4])
            igraph::shapes(shape[el[ix1, 1]])$clip(
               use_ec,
               el[ix1, , drop=FALSE],
               params=params,
               end="from")
         }))
         ec[, 3:4] <- t(sapply(seq_len(nrow(el)), function(ix1) {
            ematch <- match(edge_df$nodegroup1_2[ix1], midpoints_df$nodegroup1_2);
            use_ec <- edge.coords[ix1, 1:4, drop=FALSE];
            use_ec[,1] <- ifelse(midpoints_df[ematch, "valid"] %in% TRUE &
                  edgetest_cor_valid[ix1] %in% TRUE,
               midpoints_df[ematch, xcp3],
               use_ec[,1])
            use_ec[,2] <- ifelse(midpoints_df[ematch, "valid"] %in% TRUE &
                  edgetest_cor_valid[ix1] %in% TRUE,
               midpoints_df[ematch, ycp3],
               use_ec[,2])
            igraph::shapes(shape[el[ix1, 2]])$clip(
               use_ec,
               el[ix1, , drop=FALSE],
               params=params,
               end="to")
         }))
         # jamba::printDebug("edge.coords:");print(edge.coords);
         # jamba::printDebug("clipped ec:");print(ec);
         # assign into edge.coords
         edge.coords <- ec;
      } else {
         # Edge clipping is skipped because one or more igraph shapes
         # are not recognized by igraph::shapes().
      }
   }

   # control points for each edge spline
   edge_df$x1 <- edge.coords[, 1];

   ematch <- match(edge_df$nodegroup1_2, midpoints_df$nodegroup1_2);
   edge_df[,xcols] <- midpoints_df[ematch, xcols, drop=FALSE];

   edge_df$x3 <- edge.coords[, 3];

   edge_df$y1 <- edge.coords[, 2];

   edge_df[,ycols] <- midpoints_df[ematch, ycols, drop=FALSE];

   edge_df$y3 <- edge.coords[, 4];

   use_xcols <- jamba::vigrep("^x[123]", colnames(edge_df));
   use_ycols <- jamba::vigrep("^y[123]", colnames(edge_df));
   # propagate valid from midpoints_df
   edge_df$valid <- midpoints_df$valid[ematch];

   # print("head(edge_df):");print(head(edge_df));
   # calculate each spline
   edge_splines <- lapply(seq_len(nrow(edge_df)), function(n){
      x1 <- unlist(edge_df[n, use_xcols]);
      y1 <- unlist(edge_df[n, use_ycols]);

      # test correlation for perfect linearity
      if (FALSE %in% edgetest_cor_valid[n] ||
            FALSE %in% edge_df$valid[n]) {
         # linear segment
         x1c <- c(head(x1, 1),
            tail(x1, 1));
         y1c <- c(head(y1, 1),
            tail(y1, 1));
         path_xy <- rbind(
            cbind(x=c(head(x1c, 1),
               mean(x1c),
               tail(x1c, 1)),
               y=c(head(y1c, 1),
                  mean(y1c),
                  tail(y1c, 1))),
            c(NA, NA));
      } else {
         # bezierPath
         if ("bezierPath" %in% bundle_style) {
            path_xy <- rbind(
               ggforce:::bezierPath(
                  x=x1,
                  y=y1,
                  detail=10),
               c(NA, NA));
            colnames(path_xy) <- c("x", "y")
         } else if ("bezier" %in% bundle_style) {
            pts <- cbind(x=x1, y=y1);
            path_xy <- rbind(
               bezier::bezier(
                  t=seq(from=0, to=1, length.out=25),
                  p=tail(head(pts, -1), -1),
                  start=head(pts, 1),
                  end=tail(pts, 1)),
               c(NA, NA));
            colnames(path_xy) <- c("x", "y")
         } else if ("angular" %in% bundle_style) {
            # xspline without added control points, tends to look
            # curved but with angular turns
            path_xy <- rbind(do.call(cbind, xspline(
               x=c(x1),
               y=c(y1),
               shape=c(0, rep(1, length.out=length(x1) - 2), 0),
               # option below adds weight to endpoints
               # x=c(head(x1, 1), x1, tail(x1, 1)),
               # y=c(head(y1, 1), y1, tail(y1, 1)),
               # shape=c(0, rep(1, length.out=length(x1)), 0),
               open=TRUE,
               draw=FALSE)),
               c(NA, NA));
         } else if ("subway" %in% bundle_style) {
            # experimental method to emulate subway parallelism
            # currently works best with 4 midpoints, the first two
            # define soft curvature, the inner two define the line
            #
            # design idea: For each shared edge, encode "offset"
            # with integer values indicating the line width offset
            # above (+) or below (-) the line.
            #
            # challenge: If there were 75 edges shared, the line would
            # become width=75. In that case, apply max line width,
            # and scale the offset values accordingly.
            subway_x <- c(
               rep(head(x1, 2), c(1, 1)),
               tail(head(x1, -2), -2),
               rep(tail(x1, 2), c(1, 1))
            );
            subway_y <- c(
               rep(head(y1, 2), c(1, 1)),
               tail(head(y1, -2), -2),
               rep(tail(y1, 2), c(1, 1))
            );
            subway_shape <- c(
               c(0, 1),
               rep(0, length.out=length(tail(head(y1, -2), -2))),
               c(1, 0)
            )
            # print(data.frame(subway_x, subway_y, subway_shape));
            # Todo: apply perpendicular offset for each edge
            # that shares the same internal control points
            path_xy <- rbind(do.call(cbind, xspline(
               x=subway_x,
               y=subway_y,
               shape=subway_shape,
               open=TRUE,
               draw=FALSE)),
               c(NA, NA));
         } else {#if ("xspline" %in% bundle_style) {
            # xspline with added control points to make them rounder,
            # resembling bezierPath() output above
            path_xy <- rbind(do.call(cbind, xspline(
               # x=c(x1),
               # y=c(y1),
               # shape=c(0, rep(1, length.out=length(x1) - 2), 0),
               # option below adds weight to endpoints
               x=c(head(x1, 1), x1, tail(x1, 1)),
               y=c(head(y1, 1), y1, tail(y1, 1)),
               shape=c(0, rep(1, length.out=length(x1)), 0),
               open=TRUE,
               draw=FALSE)),
               c(NA, NA));
         }
      }

      # edge label by middle coordinate, not necessarily centered
      # path_xy_is_edge_label <- (seq_len(nrow(path_xy)) == floor(nrow(path_xy)/2))
      # edge label by most centered coordinate
      dist123 <- as.matrix(dist(as.matrix(head(path_xy[,c("x", "y"), drop=FALSE], -1))))
      dist12 <- unname(pmin(na.rm=TRUE,
         dist123[1,],
         dist123[nrow(dist123),]))
      path_xy_is_edge_label <- rep(FALSE, nrow(path_xy));
      path_xy_is_edge_label[which.max(dist12)] <- TRUE;
      path_df <- data.frame(stringsAsFactors=FALSE,
         path_xy,
         edge_row=n,
         is_edge_label=path_xy_is_edge_label,
         label=ifelse(path_xy_is_edge_label,
            edge.label[n],
            ""));
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
   if ("lty" %in% edge_attr_names) {
      edge_spline_df$lty <- igraph::E(g)$lty[edge_spline_df$edge_row];
   } else {
      edge_spline_df$lty <- 1;
   }
   # correct any NA entries
   if (any(is.na(edge_spline_df$width))) {
      edge_spline_df$width <- jamba::rmNA(edge_spline_df$width,
         naValue=1);
   }
   if (any(is.na(edge_spline_df$color))) {
      edge_spline_df$color <- jamba::rmNA(edge_spline_df$color,
         naValue="darkgrey");
   }
   if (any(is.na(edge_spline_df$lty))) {
      if (is.character(edge_spline_df$lty)) {
         edge_spline_df$lty <- jamba::rmNA(edge_spline_df$lty,
            naValue="solid");
      } else {
         edge_spline_df$lty <- jamba::rmNA(edge_spline_df$lty,
            naValue=1);
      }
   }
   edge_spline_df$arrow.mode <- arrow.mode[edge_spline_df$edge_row];
   edge_spline_df$arrow.size <- arrow.size[edge_spline_df$edge_row];
   edge_spline_df$arrow.width <- arrow.width[edge_spline_df$edge_row];
   if (FALSE) {
      jamba::printDebug("head(midpoint_df, 5):");print(head(midpoint_df, 5));
      jamba::printDebug("head(midpoints_df, 5):");print(head(midpoints_df, 5));
      jamba::printDebug("edge_spline_df:");print(subset(edge_spline_df, edge_row %in% 1));
   }

   if (draw_lines) {
      # note that lines() is not vectorized for multiple col, lwd, lty values
      # and should be re-run for each unique set
      row_split <- split(seq_len(nrow(edge_spline_df)),
         jamba::pasteByRow(edge_spline_df[,c("color", "width", "lty"),drop=FALSE]));
      for (i in row_split) {
         # render edge labels if relevant
         if (!all(edge_spline_df$label %in% c("", NA))) {
            edge_label_df <- subset(edge_spline_df, !label %in% c("", NA) &
                  is_edge_label %in% TRUE);
            if (length(edge.label.family) == 1) {
               edge.label.color <- rep(edge.label.color,
                  length.out=max(edge_label_df$edge_row))
            }
            if (length(edge.label.family) == 1) {
               edge.label.family <- rep(edge.label.family,
                  length.out=max(edge_label_df$edge_row))
            }
            if (length(edge.label.family) == 1) {
               edge.label.font <- rep(edge.label.font,
                  length.out=max(edge_label_df$edge_row))
            }
            if (length(edge.label.family) == 1) {
               edge.label.cex <- rep(edge.label.cex,
                  length.out=max(edge_label_df$edge_row))
            }
            # edge_label_df$edge_row
            edge_label_df$label_color <- edge.label.color[edge_label_df$edge_row];
            edge_label_df$label_family <- edge.label.family[edge_label_df$edge_row];
            edge_label_df$label_font <- edge.label.font[edge_label_df$edge_row];
            edge_label_df$label_cex <- edge.label.cex[edge_label_df$edge_row];
            text_subsets <- jamba::pasteByRow(edge_label_df[,c("label_family"), drop=FALSE])
            for (k in split(seq_along(text_subsets), text_subsets)) {
               text(
                  x=edge_label_df$x[k],
                  y=edge_label_df$y[k],
                  labels=edge_label_df$label[k],
                  col=edge_label_df$label_color[k],
                  family=head(edge_label_df$label_family[k], 1),
                  font=edge_label_df$label_font[k],
                  cex=edge_label_df$label_cex[k])
            }
            # points(
            #    x=edge_label_df$x,
            #    y=edge_label_df$y,
            #    col="red", cex=2);
            # points(
            #    x=edge_spline_df$x,
            #    y=edge_spline_df$y,
            #    col="red", cex=0.5);
            # text
         }

         # render the edge lines
         lines(edge_spline_df[i, 1:2, drop=FALSE],
            col=jamba::rmNA(edge_spline_df$color[i],
               naValue="darkgrey"),
            lty=jamba::rmNA(edge_spline_df$lty[i],
               naValue=1),
            lwd=jamba::rmNA(edge_spline_df$width[i],
               naValue=1))

         # next render arrow heads when necessary
         if (any(edge_spline_df$arrow.mode[i] %in% c(1, 2, 3))) {
            # first arrows on the end node
            code <- c(2, 3)
            # define valid subset of rows
            valid <- (edge_spline_df$arrow.mode[i] %in% code)
            if (any(valid)) {
               # determine first and second points
               row1 <- which(!duplicated(edge_spline_df$edge_row[i][valid]));
               row2 <- row1 + 2;
               # determine last and penultimate points
               row4 <- rev(length(edge_spline_df$edge_row[i][valid]) -
                     which(!duplicated(rev(edge_spline_df$edge_row[i][valid]))));
               row3 <- row4 - 2;
               xx <- edge_spline_df$x[i][valid]
               yy <- edge_spline_df$y[i][valid]
               if (FALSE) {
                  jamba::printDebug("");
                  print(data.frame(
                     xx_row1=xx[row1], xx_row2=xx[row2], xx_row3=xx[row3], xx_row4=xx[row4],
                     yy_row1=yy[row1], yy_row2=yy[row2], yy_row3=yy[row3], yy_row4=yy[row4],
                     row1, row2, row3, row4))
               }
               # sometimes when edges are edited the new edges have NA
               ec <- jamba::rmNA(edge_spline_df$color[i][valid],
                  naValue="darkgrey");
               ew <- jamba::rmNA(edge_spline_df$width[i][valid],
                  naValue=1)
               elty <- jamba::rmNA(edge_spline_df$lty[i][valid],
                  naValue=1)
               # acode <- edge_spline_df$arrow.code[i][valid]
               asize <- edge_spline_df$arrow.size[i][valid]
               awidth <- edge_spline_df$arrow.width[i][valid]
               # render end point arrows where necessary
               lc <- jam_igraph_arrows(
                  x1=xx[row3],
                  y1=yy[row3],
                  x2=xx[row4],
                  y2=yy[row4],
                  code=2,
                  sh.col=ec[row4],
                  # sh.col="navy",# debug
                  sh.lwd=ew[row4],
                  sh.lty=elty[row4],
                  sh.adj=1,
                  h.col=ec[row4],
                  # h.col="#FF000055",# debug
                  h.lwd=ew[row4],
                  h.lty=elty[row4],
                  open=FALSE,
                  size=asize[row4],
                  width=awidth[row4],
                  arrows_only=TRUE,
                  curved=FALSE)
            }
            #
            # next, arrows on the start node
            code <- c(1, 3)
            # define valid subset of rows
            valid <- (edge_spline_df$arrow.mode[i] %in% code)
            if (any(valid)) {
               # determine first and second points
               row1 <- which(!duplicated(edge_spline_df$edge_row[i][valid]));
               row2 <- row1 + 1;
               # determine last and penultimate points
               row4 <- rev(length(edge_spline_df$edge_row[i][valid]) -
                     which(!duplicated(rev(edge_spline_df$edge_row[i][valid]))));
               row3 <- row4 - 1;
               xx <- edge_spline_df$x[i][valid]
               yy <- edge_spline_df$y[i][valid]
               # sometimes when edges are edited the new edges have NA
               ec <- jamba::rmNA(edge_spline_df$color[i][valid],
                  naValue="darkgrey");
               ew <- jamba::rmNA(edge_spline_df$width[i][valid],
                  naValue=1)
               elty <- jamba::rmNA(edge_spline_df$lty[i][valid],
                  naValue=1)
               # acode <- edge_spline_df$arrow.code[i][valid]
               asize <- edge_spline_df$arrow.size[i][valid]
               awidth <- edge_spline_df$arrow.width[i][valid]
               # render start point arrows where necessary
               lc <- jam_igraph_arrows(
                  xx[row2],
                  yy[row2],
                  xx[row1],
                  yy[row1],
                  code=2,
                  sh.col=ec[row1],
                  sh.lwd=ew[row1],
                  sh.lty=elty[row1],
                  h.col=ec[row1],
                  h.lwd=ew[row1],
                  h.lty=elty[row1],
                  open=FALSE,
                  size=asize[row1],
                  width=awidth[row1],
                  arrows_only=TRUE,
                  curved=FALSE)
            }
         }
      }
      if (verbose) {
         jamba::printDebug("edge_bundle_nodegroups(): ",
            "midpoint_df:");print(midpoint_df);
         jamba::printDebug("edge_bundle_nodegroups(): ",
            "edge_df:");print(edge_df);
      }
      # optional debug showing control points
      if (length(debug) > 0 &&
            (TRUE %in% debug || any(grepl("bundl", debug)) ) ) {
         # indicate linear region between nodegroups
         jamba::printDebug("edge_bundle_nodegroups(): ",
            "midpoints_df:");
         print(midpoints_df);
         for (ix in seq_len(nrow(midpoints_df))) {
            xcols <- jamba::vigrep("^x[123]", colnames(midpoints_df));
            ycols <- gsub("^x", "y", xcols);
            xvals <- unlist(midpoints_df[ix, xcols]);
            yvals <- unlist(midpoints_df[ix, ycols]);
            segments(col="darkorange3",
               x0=head(xvals, 1),
               x1=tail(xvals, 1),
               y0=head(yvals, 1),
               y1=tail(yvals, 1));
            points(x=xvals, y=yvals,
               pch=rep(c(5, 1, 5), c(1, length(xcols)-2, 1)),
               cex=2.5,
               col="darkorange3")
            points(midpoints_df[,c("x1", "y1")],
               col="darkorange4",
               pch="1",
               cex=0.7);
         }

         # indicate control points
         points(edge_df[,c("x1", "y1")],
            col="darkorange4",
            pch="1",
            cex=0.7)
         points(edge_df[,c("x1", "y1")],
            col="darkorange4",
            pch=1,
            cex=2.5)
         for (x2name in jamba::vigrep("^x2", colnames(edge_df))) {
            y2name <- gsub("^x2", "y2", x2name);
            text(x=edge_df[,c(x2name, y2name)],
               col="darkorange4",
               labels=gsub("_", ".",
                  gsub("^x2", "2", x2name)),
               cex=0.7)
         }
         points(unique(edge_df[,c("x3", "y3")]),
            col="darkorange4",
            pch="3",
            cex=0.7);
         points(unique(edge_df[,c("x3", "y3")]),
            pch=1,
            cex=2.5,
            col="#00000099")
         # bottom legend
         legend("bottom",
            title="Debug legend",
            cex=0.7,
            pt.cex=2,
            y.intersp=1.5, x.intersp=1.5,
            box.col=NA,
            bg=NA,
            legend=c("<1> = start nodegroup center",
               "1 = start clipped edge",
               "2.1 .. 2.n = midpoint nodegrouop control points",
               "3 = end clipped edge",
               "<3> = end nodegroup center"),
            pch=c(5, 1, 1, 1, 5))
         legend("bottom",
            title="Debug legend",
            cex=0.7,
            pt.cex=0.6,
            y.intersp=1.5, x.intersp=1.5,
            box.col=NA,
            bg=NA,
            legend=c("<1> = start nodegroup center",
               "1 = start clipped edge",
               "2.1 .. 2.n = midpoint nodegrouop control points",
               "3 = end clipped edge",
               "<3> = end nodegroup center"),
            pch=c("1", "1", "2", "3", "3"))
      }
   }
   return(invisible(edge_spline_df));
}
