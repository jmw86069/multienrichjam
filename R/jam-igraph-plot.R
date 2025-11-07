#
# jam_igraph() specific plot functions for igraph objects

#' Jam custom function to plot an igraph network
#'
#' Jam custom function to plot an igraph network
#'
#' This function is a drop-in replacement of `igraph::plot.igraph()`,
#' intended to provide substantially faster vectorized plotting,
#' to render bundled edges when requested, and
#' to handle `rescale=FALSE` without requiring further adjustments.
#' Note that this function focuses on recognizing graph options and
#' settings, then passes the work off to `jam_plot_igraph()`
#' which performs the heavy work of rendering the graph.
#'
#' It also provides some convenient methods to adjust node
#' size, label font size, and label distance from node center,
#' based upon node attributes.
#'
#' ## vectorized plots
#'
#' This function calls `jam_plot_igraph()` as a replacement for
#' `igraph::plot.igraph()`, and that function implements vectorized
#' plot features when `vectorized_node_shapes=TRUE` by default:
#'
#' 1. When there are multiple different vertex `"shape"` attributes,
#' this function renders nodes vectorized, one shape at a time.
#' In this scenario, the original `igraph::plot.igraph()` draws
#' each individual vertex,
#' which is substantially slower (minutes compared to 1-2 seconds)
#' for large `igraph` objects.
#' 2. When there are multiple font families, labels are rendered in
#' groups by font family, in order to comply with limitations
#' in `graphics::text()`. This situation is fairly rare, however
#' the speed improvement is substantial, again roughly minutes down
#' to 1-2 seconds. The `igraph::plot.igraph()` renders each node label
#' individually when there are multiple font families.
#'
#' ## rescale=FALSE
#'
#' The default `igraph::plot.igraph()` uses `rescale=TRUE`, which
#' can distort the layout coordinates to fit within a fixed
#' x- and y-axis range `c(-1, 1)`. When using `rescale=FALSE` the
#' `xlim` and `ylim` values are not adjusted to the actual coordinate
#' range. The desired effect from this function `jam_igraph()` is
#' to apply `aspect=1` (`asp=1`) which fixes the aspect ratio so
#' the coordinates represent visual Euclidean distance between nodes,
#' and to define `xlim` and `ylim` to accomodate the full layout.
#' This function also adjusts node `vertex.size` and `vertex.label.dist`
#' proportionally.
#'
#' ## edge bundling
#'
#' When `edge_bundling` is something other than `edge_bundling="none"`,
#' edge connections between nodes are rendered using a specific
#' function by drawing curved splines for each bundle of edges.
#' The approach in `igraph::plot.igraph()` only draws
#' straight edges between nodes. The recommended method is
#' `edge_bundling="default"` which will try to detect an appropriate
#' method to bundle edges. When `mark.groups` and `nodegroups` are not
#' defined, the default method is `"connections"` which bundles edges
#' only among nodes that share the same connections. The assumption is
#' that nodes that share the same connections usually have very similar
#' layout coordinates, so edge bundling is usually intuitive. In fact,
#' for a very large set of nodes, they are often in a round cluster
#  in the layout, so the edges bundle together and result in a much
#  cleaner overall visualization.
#'
#' ## Adjust node size, label size, label distance
#'
#' The following arguments apply scaling to all nodes or edges:
#'
#' * `node_factor` - `numeric` multiplied by node size
#' * `edge_factor` - `numeric` multiplied by edge size
#' * `label_factor` - `numeric` multiplied by label font size
#' * `label_dist_factor` - `numeric` multiplied by label distance from
#' node center
#'
#' New attributes `vertex.label.fontsize` and `edge.label.fontsize`
#' which define fixed fontsize in points for nodes and edges, respectively.
#' These values are not modified by `vertex.label.cex` nor `edge.label.cex`
#' and are intended to allow control over specific fonts used in the final
#' figure. Note their calculations are based upon `par("ps")` which should
#' represent device-dependent point size. If this value is inappropriate,
#' it should be adjusted to control the font sizing.
#'
#' The following arguments apply scale factor based upon node attribute:
#'
#' * `node_factor_l` - `list` of named vectors applied to node size
#' * `label_factor_l` - `list` of named vectors applied as `label_factor`
#'.   as a multiplier to label font size.
#' * `label_fontsize_l` - `list` of named vectors applied to define a
#'    specific, fixed label fontsize in points, which is not modified
#'.   by `vertex.label.cex` nor `label_factor`.
#' * `label_dist_factor_l` - `list` of named vectors applied to label
#' distance from node center
#'
#' The `factor_l` technique is as follows:
#'
#' ```R
#' node_factor_l = list(node_attr_name = c(
#'    attr_value1 = factor1,
#'    attr_value2 = factor2))
#' ```
#'
#' A specific example:
#' ```R
#' node_factor_l = list(nodeType = c(
#'    Gene=1.5,
#'    Set=2))
#' ```
#'
#' In this case, nodes with attribute
#' `igraph::V(g)$nodeType == "Gene"` will use factor `1.5`
#' Nodes with attribute
#' `igraph::V(g)$nodeType == "Set"` will use factor `2`
#' All other nodes will not be adjusted.
#'
#' ## Other features
#'
#' The plot layout by default is **not** rescaled to `c(-1, 1)`, therefore
#' allowing direct control over plot dimensions and node sizes.
#' The plot aspect ratio is fixed at 1, which renders many network layouts
#' in their intended form, as opposed to scaling each axis to `c(-1, 1)`,
#' which can impose distortion of intended layout node distances.
#'
#' When `use_shadowText=TRUE` node labels call `jamba::shadowText()`
#' which draws a small partly transparent outline around labels, making
#' them more legible when they overlap colored nodes. This step
#' effectively draws each label `n` times, which can slightly slow
#' the rendering of the overall figure.
#'
#' When `pie_to_jampie=TRUE`, any nodes with `shape="pie"` are
#' changed to `shape="jampie"` for the purpose of rendering pie
#' shapes in vectorized fashion, instead of being drawn for each
#' node separately. This change is a substantial improvement in
#' rendering time. In addition, optional node attributes are available:
#' * `pie.border` to control individual pie wedge borders, which are
#' drawn as inner borders so each pie wedge border is visible without
#  overlap.
#' * `pie.lwd` to control line width of pie wedge borders.
#' * `pie.lty` to control the line type of pie wedge borders
#' * `frame.color` to control the frame border color drawn around the
#' full circular pie node. This border is drawn as an outer border, so
#' it will not overlap any internal pie wedge border colors.
#'
#' Default colors for marked node groups `mark.col` and `mark.border`
#' when not defined upfront, will call `colorjam::rainbowJam()`
#' and not `grDevices::rainbow()`. The `colorjam::rainbowJam()`
#' produces more visually distinct categorical colors.
#' This behavior can be controlled by supplying a `character`
#' vector with specific colors for `mark.col` and `mark.border`. Note
#' that the border should match the colors, or it can be set to `"grey45"`
#' for a generally visible border.
#'
#' When `names(mark.groups)` is defined, the values are used as
#' labels, positioned at the outer edge of each polygon. The label
#' text size is adjusted with `label.cex`, and the position can
#' be adjusted with `mark.x.nudge`, `mark.y.nudge`, in units of
#' fraction of the maximum x- or y-axis range (effectively fraction
#' of the layout size).
#'
#' Optional argument `nodegroups` can be supplied, which is a `list`
#' of vectors, where each vector represents a group of nodes. The
#' `nodegroups` can be used with `edge_bundling="nodegroups"` to
#' define custom edge bundling. This option is useful for defining a
#' group of nodes for edge bundling, when those nodes should not be used
#' to render group borders as with `mark.groups`.
#'
#' Finally, individual plot components can be individually disabled:
#'
#' * `render_nodes=FALSE`
#' * `render_edges=FALSE`
#' * `render_groups=FALSE`
#' * `render_nodelabels=FALSE`
#'
#' @family multienrichjam core functions
#' @family jam igraph utilities
#'
#' @inheritParams jam_plot_igraph
#' @param x `igraph` object to be plotted
#' @param ... additional arguments are passed to `igraph::plot.igraph()`
#' @param expand `numeric` value used to expand the x and y axis ranges,
#'    where `0.03` expands each size `3%`.
#' @param rescale `logical` indicating whether to rescale the layout
#'    coordinates to `c(-1, 1)`. When `rescale=FALSE` the original
#'    layout coordinates are used as-is without change.
#' @param mark.groups `list` or `igraph::components` object from one of
#'    the many igraph cluster functions.
#'    * When not explicitly supplied, it will use graph attributes
#'    "mark.groups" as a way of re-using stored sub-clusters.
#'    Graph attributes are defined by `mem2emap()` for example,
#'    `igraph::graph_attr(x, "mark.groups")`. Note 'mark.colors' is also
#'    used when defined together with 'mark.groups'.
#'    * When 'mark.groups' is `FALSE`, it will not apply any 'mark.groups'
#'    even when graph attributes contain 'mark.groups'.
#' @param node_factor `numeric` value multiplied by `igraph::V(x)$size` to adjust
#'    the relative size of all nodes by a common numeric scalar value.
#' @param edge_factor `numeric` value multiplied by `igraph::E(x)$width` to adjust
#'    the relative width of all edges by a common numeric scalar value.
#' @param label_factor `numeric` value multiplied by `igraph::V(x)$label.cex`
#'    and `igraph::E(x)$label.cex` to adjust the relative size of all labels on
#'    nodes and edges by a common numeric scalar value.
#' @param label_dist_factor `numeric` value multiplied by `igraph::V(x)$label.dist`
#'    to adjust the relative distance of all nodes labels from the node center
#'    by a common numeric scalar value.
#' @param node_factor_l,label_factor_l,label_dist_factor_l `list`
#'    of vectors, where the names of the `list` are attribute
#'    names, and the names of each vector are attributes values.
#'    These values are applied in addition to `node_factor`, `label_factor`,
#'    `label_dist_factor`, respectively.
#'    The vector values are used as scalar multipliers, analogous to
#'    `node_factor`. The purpose is to apply scalar values to different
#'    subsets of nodes. For example, consider:
#'    `node_factor_l=list(nodeType=c(Gene=1, Set=2)`. The list name
#'    `"nodeType"` says to look at `igraph::vertex_attr(x, "nodeType")`.
#'    Nodes with `nodeType="Gene"` will use `1`, and `nodeType="Set"`
#'    will use `2` as the scalar value.
#' @param plot_function `function` that renders the graph, not intended to
#'    be changed except for very customized uses. By default
#'    `plot_function=jam_plot_igraph()` which calls a modified variant of
#'    `igraph:::plot.igraph()`.
#'
#' @examples
#' # Example with karate
#' karate <- igraph::make_graph("Zachary");
#' cl <- igraph::cluster_louvain(karate);
#' jam_igraph(karate,
#'    layout=layout_with_qfrf(repulse=3.5),
#'    mark.groups=cl,
#'   mark.lwd=c(1:4),
#'   mark.lty=1:4, mark.shape=1,
#'   edge_bundling="default")
#'
#' # create example cnet data
#' cnet <- make_cnet_test(num_sets=3)
#'
#' ## example showing how to use the list form
#' ## This form resizes nodes where igraph::V(g)$nodeType %in% "Gene" by 2x,
#' ## and resizes nodes where igraph::V(g)$nodeType %in% "Set" by 3x.
#' node_factor_l <- list(nodeType=c(Gene=1.2, Set=2));
#'
#' ## This form multiplies label.dist for nodeType="Gene" nodes by 1,
#' ## and multiplies label.dist for nodeType="Set" nodes by 0.5
#' label_dist_factor_l <- list(nodeType=c(Gene=1, Set=0.5))
#'
#' par("mar"=c(0, 0, 0, 0) + 0.5);
#' jam_igraph(cnet,
#'    use_shadowText=TRUE,
#'    node_factor_l=node_factor_l,
#'    label_factor=0.6,
#'    label_factor_l=list(nodeType=c(Gene=1, Set=2)))
#' par("mar"=c(2, 2, 2, 2));
#'
#' # Example using edge bundling by community detection
#' g <- igraph::make_graph("Zachary");
#' gcom <- igraph::cluster_leading_eigen(g);
#'
#' jam_igraph(g,
#'    layout=layout_with_qfr,
#'    edge_bundling="nodegroups",
#'    mark.groups=gcom,
#'    nodegroups=gcom,
#'    vertex.color=colorjam::group2colors(gcom$membership))
#'
#' cfuncs <- list(cluster_leading_eigen=igraph::cluster_leading_eigen,
#'    cluster_edge_betweenness=igraph::cluster_edge_betweenness,
#'    cluster_fast_greedy=igraph::cluster_fast_greedy,
#'    cluster_spinglass=igraph::cluster_spinglass)
#' for (i in seq_along(cfuncs)) {
#'    cfunc <- cfuncs[[i]];
#'    gcom <- cfunc(g);
#'    igraph::V(g)$color <- colorjam::group2colors(gcom$membership);
#'    g <- color_edges_by_nodes(g);
#'    set.seed(123);
#'    jam_igraph(g,
#'       layout=layout_with_qfr,
#'       edge_bundling="nodegroups",
#'       nodegroups=gcom,
#'       mark.groups=gcom)
#'    title(main=names(cfuncs)[i]);
#' }
#'
#' # fancy example showing mark.groups and colorizing
#' # edges using node colors
#' gcom <- igraph::cluster_spinglass(g);
#' igraph::V(g)$color <- colorjam::group2colors(gcom$membership);
#' g <- color_edges_by_nodes(g);
#' jam_igraph(g,
#'    layout=layout_with_qfrf(repulse=3.2),
#'    edge_bundling="nodegroups",
#'    nodegroups=gcom,
#'    mark.groups=gcom)
#' title(main=paste0("cluster_spinglass()\n",
#'    "edge_bundling='nodegroups'"))
#'
#' # same but different edge_style
#' jam_igraph(g,
#'    layout=layout_with_qfrf(repulse=3.2),
#'    edge_bundling="nodegroups",
#'    nodegroups=gcom,
#'    mark.groups=gcom,
#'    bundle_style="xspline",
#'    detail=14)
#' title(main="bundle_style='xspline'")
#'
#' # same but using node connections
#' jam_igraph(g,
#'    layout=layout_with_qfrf(repulse=3.2),
#'    edge_bundling="connections",
#'    nodegroups=gcom,
#'    mark.groups=gcom)
#' title(main="edge_bundling='connections'")
#'
#' @export
jam_igraph <- function
(x,
 ...,
 xlim=NULL,
 ylim=NULL,
 expand=0.03,
 rescale=FALSE,
 node_factor=1,
 node_factor_l=NULL,
 edge_factor=1,
 edge_factor_l=NULL,
 label_factor=1,
 label_factor_l=NULL,
 label_fontsize_l=NULL,
 label_dist_factor=1,
 label_dist_factor_l=1,
 use_shadowText=FALSE,
 edge_bundling=c(
    "default",
    "connections",
    "none",
    "mark.groups",
    "nodegroups"),
 bundle_self=FALSE,
 nodegroups=NULL,
 render_nodes=TRUE,
 render_edges=TRUE,
 render_nodelabels=TRUE,
 render_groups=TRUE,
 vectorized_node_shapes=TRUE,
 plot_grid=FALSE,
 plot_function=jam_plot_igraph,
 mark.groups=list(),
 mark.shape=1/2,
 mark.col=NULL,
 mark.alpha=0.2,
 mark.border=NULL,
 mark.expand=NULL,
 mark.lwd=2,
 mark.lty=1,
 mark.smooth=TRUE,
 mark.cex=1,
 mark.x.nudge=0,
 mark.y.nudge=0,
 verbose=FALSE,
 debug=NULL)
{
   ##
   # validate edge_bundling input
   if (length(edge_bundling) == 0) {
      edge_bundling <- "none";
   } else if (is.atomic(edge_bundling)) {
      edge_bundling <- head(intersect(edge_bundling,
         eval(formals(jam_igraph)$edge_bundling)), 1);
      if (length(edge_bundling) == 0) {
         stop("edge_bundling method not recognized")
      }
   } else if (is.function(edge_bundling)) {
      edge_function <- edge_bundling;
      edge_bundling <- "function";
   }

   # params <- igraph:::i.parse.plot.params(x, list(...));
   params <- parse_igraph_plot_params(x, list(...));

   # create layout coordinates
   layout <- params("plot", "layout");
   # store back in params object environment for persistence
   environment(params)$p$plot$layout <- layout;

   vertex.size <- params("vertex", "size");
   # jamba::printDebug("jam_igraph(): ", "params('vertex', 'size'): ", vertex.size);# debug
   vertex.size2 <- params("vertex", "size2");

   vertex.label.cex <- params("vertex", "label.cex");
   edge.label.cex <- params("edge", "label.cex");
   if (is.function(label_factor)) {
      vertex.label.cex <- label_factor(vertex.label.cex);
      edge.label.cex <- label_factor(edge.label.cex);
      if (verbose) {
         jamba::printDebug("jam_igraph(): ",
            "Applying: ", "label_factor(label.cex)");
         jamba::printDebug("jam_igraph(): ",
            "vertex.label.cex: ", head(vertex.label.cex, 10));
         jamba::printDebug("jam_igraph(): ",
            "edge.label.cex: ", head(edge.label.cex, 10));
      }
   } else {
      if (verbose) {
         jamba::printDebug("jam_igraph(): ",
            "Applying: ", "label.cex * label_factor");
      }
      vertex.label.cex <- vertex.label.cex * label_factor;
      edge.label.cex <- edge.label.cex * label_factor;
   }

   if (is.function(node_factor)) {
      vertex.size <- node_factor(vertex.size);
      vertex.size2 <- node_factor(vertex.size2);
   } else {
      vertex.size <- vertex.size * node_factor;
      vertex.size2 <- vertex.size2 * node_factor;
   }

   vertex.label.dist <- params("vertex", "label.dist");
   if (is.function(label_dist_factor)) {
      vertex.label.dist <- label_dist_factor(vertex.label.dist);
   } else {
      vertex.label.dist <- vertex.label.dist * label_dist_factor;
   }

   ## label_fontsize_l=list(nodeType=c(Gene=1, Set=2))
   vertex.label.fontsize <- NULL;
   if (length(label_fontsize_l) > 0) {
      vertex.label.fontsize <- handle_igraph_param_list(x,
         attr="label.fontsize",
         factor_l=label_fontsize_l,
         i_values=1,
         attr_type="node")
   }

   ## label_factor_l=list(nodeType=c(Gene=1, Set=2))
   if (length(label_factor_l) > 0) {
      vertex.label.cex <- handle_igraph_param_list(x,
         attr="label.cex",
         factor_l=label_factor_l,
         i_values=vertex.label.cex,
         attr_type="node")
   }
   if (length(label_dist_factor_l) > 0) {
      vertex.label.dist <- handle_igraph_param_list(x,
         attr="label.dist",
         factor_l=label_dist_factor_l,
         i_values=vertex.label.dist,
         attr_type="node");
   }
   if (length(node_factor_l) > 0) {
      if (verbose) {
         jamba::printDebug("jam_igraph(): ",
            "Applying: ", "node_factor_l to vertex.size");
         jamba::printDebug("jam_igraph(): ",
            "vertex.size: ");print(head(vertex.size, 10));
      }
      vertex.size <- handle_igraph_param_list(x,
         attr="size",
         factor_l=node_factor_l,
         i_values=vertex.size,
         attr_type="node");
      vertex.size2 <- handle_igraph_param_list(x,
         attr="size2",
         factor_l=node_factor_l,
         i_values=vertex.size2,
         attr_type="node");
   }

   edge.width <- params("edge", "width");
   if (is.function(edge_factor)) {
      edge.width <- edge_factor(edge.width);
   } else {
      edge.width <- edge.width * edge_factor;
   }

   if (length(edge_factor_l) > 0) {
      edge.width <- handle_igraph_param_list(x,
         attr="width",
         edge_factor_l,
         i_values=edge.width,
         attr_type="edge");
   }

   ## Optional axis range scaling
   if (!TRUE %in% rescale) {
      dist_factor <- 4;
      x_range <- range(layout[,1], na.rm=TRUE);
      y_range <- range(layout[,2], na.rm=TRUE);
      max_xy_range <- max(c(
         diff(x_range),
         diff(y_range)));
      xlim_asis <- TRUE;
      if (length(xlim) == 0) {
         xlim <- x_range;
         xlim_asis <- FALSE;
      }
      ylim_asis <- TRUE;
      if (length(ylim) == 0) {
         ylim <- y_range;
         ylim_asis <- FALSE;
      }
      if (length(debug) > 0 && any(c("vertex.label.dist","label.dist","labels") %in% debug)) {
         jamba::printDebug("jam_igraph(): ",
            "xlim before:",
            xlim);
         jamba::printDebug("jam_igraph(): ",
            "head(vertex.size, 20) before:",
            head(vertex.size, 20));
      }
      vertex.size <- vertex.size * (max_xy_range) / 1.8;
      vertex.size2 <- vertex.size2 * (max_xy_range) / 1.8;
      vertex.label.dist <- vertex.label.dist * (max_xy_range) / dist_factor / 1.8;
      if (!xlim_asis) {
         xlim <- xlim + diff(xlim) * c(-1, 1) * expand;
      }
      if (!ylim_asis) {
         ylim <- ylim + diff(ylim) * c(-1, 1) * expand;
      }
   } else {
      if (length(xlim) == 0) {
         xlim <- c(-1, 1);
      }
      if (length(ylim) == 0) {
         ylim <- c(-1, 1);
      }
   }

   # store back in params object environment for persistence
   # vertex attributes
   # size, size2, label.cex, label.dist
   # 0.0.101.900: divide by 200 and stop dividing by 200 in jampie shape plot
   vertex.size <- vertex.size / 200;
   vertex.size2 <- vertex.size2 / 200;
   #
   # assign to params environment for use by downstream functions
   environment(params)$p$vertex$size <- vertex.size;
   # jamba::printDebug("assigned vertex.size: ", vertex.size);# debug
   
   environment(params)$p$vertex$size2 <- vertex.size2;
   environment(params)$p$vertex$label.cex <- vertex.label.cex;
   environment(params)$p$vertex$label.dist <- vertex.label.dist;
   if (length(vertex.label.fontsize) > 0) {
      environment(params)$p$vertex$label.fontsize <- vertex.label.fontsize;
   }
   # edge attributes
   # label.cex, width
   environment(params)$p$edge$width <- edge.width;
   environment(params)$p$edge$label.cex <- edge.label.cex;
   # if (length(edge.label.fontsize) > 0) {
   #    environment(params)$p$edge$label.fontsize <- edge.label.fontsize;
   # }

   if (length(debug) > 0 &&
         any(c("vertex.label.dist","label.dist","labels") %in% debug)) {
      jamba::printDebug("jam_igraph(): ",
         "xlim after:",
         xlim);
      jamba::printDebug("jam_igraph(): ",
         "head(vertex.label.dist, 20):",
         head(vertex.label.dist, 20));
      jamba::printDebug("jam_igraph(): ",
         "head(vertex.size, 20) after:",
         head(vertex.size, 20));
   }
   
   # determine whether mark.groups is encoded in graph_attr
   if (length(mark.groups) == 0 &&
         "mark.groups" %in% igraph::graph_attr_names(x)) {
      mark.groups <- igraph::graph_attr(x, "mark.groups");
      # convert to list
      mark.groups <- setNames(
         split(mark.groups$names, mark.groups$membership),
         jamba::cPaste(mark.groups$cluster_names, ",\n"))
      # confirm all elements of mark.groups are present as node names
      if (!all(unlist(mark.groups) %in% igraph::V(x)$name)) {
         mark.groups <- NULL;
      } else {
         if (length(mark.col) == 0 &&
               "mark.colors" %in% igraph::graph_attr_names(x)) {
            mark.col <- igraph::graph_attr(x, "mark.colors");
         }
      }
   }

   plot_function(x=x,
      ...,
      rescale=rescale,
      ## the following parameters are now included inside params
      # vertex.size=vertex.size,
      # vertex.size2=vertex.size2,
      # vertex.label.dist=vertex.label.dist,
      # vertex.label.cex=vertex.label.cex,
      # edge.label.cex=edge.label.cex,
      # edge.width=edge.width,
      mark.groups=mark.groups,
      mark.shape=mark.shape,
      mark.col=mark.col,
      mark.alpha=mark.alpha,
      mark.border=mark.border,
      mark.expand=mark.expand,
      mark.lwd=mark.lwd,
      mark.lty=mark.lty,
      mark.smooth=mark.smooth,
      mark.cex=mark.cex,
      mark.x.nudge=mark.x.nudge,
      mark.y.nudge=mark.y.nudge,
      use_shadowText=use_shadowText,
      xlim=xlim,
      ylim=ylim,
      render_nodes=render_nodes,
      render_edges=render_edges,
      render_nodelabels=render_nodelabels,
      render_groups=render_groups,
      edge_bundling=edge_bundling,
      bundle_self=bundle_self,
      nodegroups=nodegroups,
      vectorized_node_shapes=vectorized_node_shapes,
      plot_grid=plot_grid,
      params=params,
      debug=debug);
   return(invisible(params));
}



#' Handle igraph attribute parameter list
#'
#' Handle igraph attribute parameter list, internal function for `jam_igraph()`
#'
#' This mechanism is intended to help update `igraph` attributes
#' in bulk operations by the attribute values associated with
#' nodes or edges. Most commonly, the argument `factor_l` is
#' multiplied by numeric attributes to scale attribute values,
#' for example label font size, or node size.
#'
#' For example:
#'
#' ```
#' handle_igraph_param_list(x,
#'    attr="size",
#'    factor_l=list(nodeType=c(Gene=1, Set=2)),
#'    i_values=rep(1, igraph::vcount(x)),
#'    attr_type="node")
#' ```
#'
#' This function call will match node attribute `nodeType`,
#' the size of nodes with attribute value
#' `nodeType="Set"` are multiplied `size * 2`,
#' `nodeType="Gene"` are multiplied `size * 1`.
#'
#' @return `vector` of attribute values representing `attr`.
#'
#' @family jam utility functions
#'
#' @param x `igraph` object
#' @param attr `character` name of the attribute to update in `x`.
#' @param factor_l `list` or `numeric` vector or `function`:
#'    * `list` of `numeric` vectors where `names(factor_l)`
#'    correspond to attribute names, and the names of numeric vectors
#'    are attribute values The attribute names and attribute values are
#'    used to match relevant entities of type `attr_type`.
#'    For matching entities, attribute values are used as defined by
#'    attribute name `attr`, and are multiplied by the matching
#'    numeric value in `factor_l`.
#'    * `numeric` vector which is directly multiplied by `i_values`
#'    to produce an adjusted output vector `i_values`.
#'    * `function` which is used to modify `i_values` by calling
#'    `factor_l(i_values)` to produce adjusted output `i_values`.
#' @param i_values `vector` of attribute values that represent the current
#'    attribute values in `x` for the attribute `attr`.
#' @param attr_type `character` string indicating the type of entity
#'    being adjusted in `x`:
#'    * `"node"` or `"vertex"` refers to `igraph::vertex_attr()`
#'    * `"edge"` refers to `igraph::edge_attr()`
#' @param ... additional arguments are ignored.
#'
#'
#' @export
handle_igraph_param_list <- function
(x,
 attr,
 factor_l,
 i_values=NULL,
 attr_type=c("node", "vertex", "edge"),
 verbose=FALSE,
 ...)
{
   attr_type <- match.arg(attr_type);
   if (attr_type %in% c("node", "vertex")) {
      vct <- igraph::vcount(x);
      x_attr_names <- igraph::vertex_attr_names(x);
      if (length(i_values) == 0) {
         i_values <- igraph::vertex_attr(x, attr);
      }
      if (length(i_values) > 0 && length(i_values) < vct) {
         i_values <- rep(i_values,
            length.out=vct);
      }
   } else if (attr_type %in% c("edge")) {
      ect <- igraph::ecount(x);
      x_attr_names <- igraph::edge_attr_names(x);
      if (length(i_values) == 0) {
         i_values <- igraph::edge_attr(x, attr);
      }
      if (length(i_values) > 0 && length(i_values) < ect) {
         i_values <- rep(i_values,
            length.out=ect);
      }
   }
   if (is.numeric(factor_l)) {
      if (verbose) {
         jamba::printDebug("jam_igraph(): ",
            "Applying ['", attr, "'] as i_values * (factor_l)");
      }
      i_values <- factor_l * i_values;
   } else if (is.function(factor_l)) {
      if (verbose) {
         jamba::printDebug("jam_igraph(): ",
            "Applying ['", attr, "'] using factor_l(i_values)");
      }
      i_values <- factor_l(i_values);
   } else if (length(factor_l) > 0 && is.list(factor_l)) {
      # iterate each name in factor_l which corresponds to node/vertex attribute names
      for (i in names(factor_l)) {
         j <- factor_l[[i]];
         # iterate the name of each vector element in factor_l[[i]]
         # to identify nodes whose attribute value matches the name
         for (k in names(j)) {
            adjust_value <- factor_l[[i]][[k]];
            if (i %in% x_attr_names) {
               if (attr_type %in% c("node", "vertex")) {
                  i_update <- (igraph::vertex_attr(x, i) %in% k);
               } else {
                  i_update <- (igraph::edge_attr(x, i) %in% k);
               }
               if (any(i_update)) {
                  if (!is.function(adjust_value)) {
                     if (is.numeric(adjust_value)) {
                        adj_label <- paste0(" ", attr_type,": ['", attr, "'] * (", adjust_value, ")");
                        i_values[i_update] <- i_values[i_update] * adjust_value;
                     } else {
                        adj_label <- paste0(" ", attr_type,": ['", attr, "'] <- '", adjust_value, "'");
                        i_values[i_update] <- adjust_value;
                     }
                  } else if (is.function(adjust_value)) {
                     adj_label <- paste0(" ", attr_type,": adjust_function(['", attr, "'])");
                     i_values[i_update] <- adjust_value(i_values[i_update]);
                  }
                  # optional verbose output
                  if (verbose) {
                     jamba::printDebug(sep="",
                        c("jam_igraph(): ",
                           "Applying ",
                           " factor_l[['",i,"']][['",k,"']] to ",
                           jamba::formatInt(sum(i_update)),
                           adj_label));
                  }
               }
            }
         }
      }
   }
   return(i_values);
}


#' Plot layout scale by percentage of coordinate range
#'
#' Plot layout scale by percentage of coordinate range
#'
#' The `layout` argument supplied coordinates, and largest numeric
#' range of any column is used to define 100 percent scale.
#' A grey grid is drawn on a base R plot to indicate the big
#' and small steps across the range, using `big_tick` and `small_tick`,
#' respectively.
#'
#' @family jam plot functions
#'
#' @param layout `matrix` or `data.frame` with at least two columns,
#'    only the first two columns are used for the grid.
#' @param grid_colors `character` colors used for the small and big grid
#'    lines, respectively.
#' @param grid_lty `integer` or `character` line type used for the
#'    small and big grid lines, respectively.
#' @param big_tick,small_tick `numeric` values in percent, the step size
#'    between big grid lines, and small grid lines, respectively.
#' @param ... additional arguments are ignored.
#'
#' @export
plot_layout_scale <- function
(layout,
 grid_colors=c("grey80", "grey70"),
 grid_lty=c("dotted", "solid"),
 big_tick=10,
 small_tick=2.5,
 ...)
{
   #
   if (!any(c("numeric", "matrix", "data.frame", "tbl_df") %in% class(layout))) {
      stop("layout must contain numeric coordinates as matrix or data.frame")
   }
   if (length(grid_colors) == 0) {
      grid_colors <- c("grey80", "grey70")
   }
   grid_colors <- rep(grid_colors, length.out=2)
   if (length(grid_lty) == 0) {
      grid_lty <- c("dotted", "solid")
   }
   grid_lty <- rep(grid_lty, length.out=2)
   xy_ranges <- apply(layout, 2, function(i){
      jamba::nameVector(range(i, na.rm=TRUE),
         c("min", "max"))})
   xy_mids <- apply(xy_ranges, 2, mean)
   xy_spans <- apply(xy_ranges, 2, diff)
   xy_max <- max(xy_spans)
   # tick mark steps
   big_step <- big_tick;
   smol_step <- small_tick;
   seq50_big <- seq(from=0, to=100, by=big_step)
   seq50_smol <- setdiff(seq(from=0, to=100, by=smol_step), seq50_big)
   # tick mark positions
   seq50 <- sort(c(seq50_big, seq50_smol))
   seq50_issmol <- (seq50 %in% seq50_smol)
   ticks_x <- (seq50 - 50) * xy_max / 100 + xy_mids[1]
   ticks_y <- (seq50 - 50) * xy_max / 100 + xy_mids[2]
   # apply ablines
   abline(v=ticks_x,
      lty=ifelse(seq50_issmol, grid_lty[[1]], grid_lty[[2]]),
      col=ifelse(seq50_issmol, grid_colors[[1]], grid_colors[[2]]))
   abline(h=ticks_y,
      lty=ifelse(seq50_issmol, grid_lty[[1]], grid_lty[[2]]),
      col=ifelse(seq50_issmol, grid_colors[[1]], grid_colors[[2]]))
   rect(xleft=min(ticks_x), xright=max(ticks_x),
      ybottom=min(ticks_y), ytop=max(ticks_y),
      lwd=2,
      lty=grid_lty[[2]],
      col=NA, border=grid_colors[[2]])

}

