
# Jam Cnet plot functions

#' Nudge igraph layout by node
#'
#' Nudge igraph layout by node
#'
#' @family jam cnet igraph functions
#' @family jam igraph functions
#'
#' This function takes an `igraph` object that contains a layout
#' stored in graph attributes, `graph_attr(g, "layout")`, and adjusts
#' the position of one or more nodes based upon relative x- and y-axis
#' size of the layout.
#'
#' @param g `igraph` object that contains layout coordinates stored as
#'    graph attribute `"layout"`, for example `graph_attr(g, "layout")`.
#' @param nodes `character` vector indicating which nodes in `g` should
#'    be nudged.
#' @param x,y `numeric` values indicating the amount to move each node
#'    defined by `nodes`. These values are relative to the x- and y-axis
#'    ranges of the layout coordinates, and based upon `aspect` below.
#' @param use_grep `logical` indicating whether to match values in `nodes`
#'    to `V(g)$name` and `V(g)$label` using `jamba::provigrep()`, which
#'    follows case-insensitive `grep()`. When `use_grep=FALSE` the
#'    values in either `V(g)$name` or `V(g)$label` must be identical to
#'    `nodes`.
#' @param aspect `numeric` indicating the aspect ratio of the output
#'    plot. Any value other than `aspect=1` will use the observed x-axis
#'    and y-axis range of the layout coordinates. When `aspect=1` then
#'    the highest of x-axis and y-axis ranges is used for the relative
#'    `x` and `y` adjustment. Note that `igraph::plot()` does not maintain
#'    aspect ratio 1 by default, but `jam_igraph()` does maintain
#'    aspect ratio 1 and is preferred.
#' @param debug `logical` indicating whether to plot the layout before
#'    and after adjustment, drawing arrows to indicate the movement
#'    of particular nodes. This plot is very basic, using base R `plot()`,
#'    and is only intended as a quick visual review.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#'
#' @export
nudge_igraph_node <- function
(g,
 nodes,
 x=0,
 y=0,
 use_grep=TRUE,
 aspect=1,
 debug=FALSE,
 verbose=FALSE,
 ...)
{
   x <- rep(x, length.out=length(nodes));
   y <- rep(y, length.out=length(nodes));
   if (!"layout" %in% igraph::list.graph.attributes(g)) {
      stop("g must have graph attribute 'layout'.");
   }
   layout <- igraph::graph_attr(g, "layout");
   if (aspect == 1) {
      xrange <- max(c(
         diff(range(layout[,1])),
         diff(range(layout[,2]))));
      yrange <- xrange;
   } else {
      xrange <- diff(range(layout[,1]));
      yrange <- diff(range(layout[,2]));
   }
   if (is.numeric(nodes)) {
      idx <- nodes;
   } else {
      if (!use_grep || all(nodes %in% igraph::V(g)$name | nodes %in% igraph::V(g)$label)) {
         matchlist <- list(match(nodes, igraph::V(g)$name),
            match(nodes, igraph::V(g)$label));
         pmin_narm <- function(...){pmin(..., na.rm=TRUE)};
         idx <- Reduce("pmin_narm", matchlist);
      } else if (use_grep) {
         idx1 <- jamba::provigrep(nodes, igraph::V(g)$name, returnType="list", value=FALSE)
         idx2 <- jamba::provigrep(nodes, igraph::V(g)$label, returnType="list", value=FALSE)
         idx <- lapply(seq_along(idx1), function(i){
            sort(unique(c(idx1[[i]], idx2[[i]])));
         });
         x <- rep(x, lengths(idx));
         y <- rep(y, lengths(idx));
         if (verbose) {
            jamba::printDebug("nudge_igraph_node(): ",
               "idx1:");
            print(idx1);
            jamba::printDebug("nudge_igraph_node(): ",
               "idx2:");
            print(idx2);
            jamba::printDebug("nudge_igraph_node(): ",
               "idx:");
            print(idx);
            jamba::printDebug("nudge_igraph_node(): ",
               "x:");
            print(x);
            jamba::printDebug("nudge_igraph_node(): ",
               "y:");
            print(y);
         }
         idx <- unlist(idx);
      }
   }
   idx_na <- is.na(idx);
   if (any(idx_na)) {
      x <- x[!idx_na];
      y <- y[!idx_na];
      idx <- idx[!idx_na];
   }
   xadj <- x * xrange;
   yadj <- y * yrange;
   if (verbose) {
      nudge_labels <- format(justify="right",
         c(igraph::V(g)$name[idx],
            round(digits=2, xadj),
            round(digits=2, yadj)));
      n <- length(idx);
      jamba::printDebug("nudge_igraph_node(): ",
         "node names: ",
         head(nudge_labels, n));
      jamba::printDebug("nudge_igraph_node(): ",
         "nudge x:    ",
         tail(head(nudge_labels, n*2), n));
      jamba::printDebug("nudge_igraph_node(): ",
         "nudge y:    ",
         tail(nudge_labels, n));
   }
   if (debug) {
      plot(layout,
         asp=aspect);
      layout1 <- layout[idx,,drop=FALSE];
      points(layout1, col="red3", pch=20)
   }
   layout[idx,1] <- layout[idx,1] + xadj;
   layout[idx,2] <- layout[idx,2] + yadj;
   if (debug) {
      layout2 <- layout[idx,,drop=FALSE];
      points(layout2, col="green3", pch=20, cex=3)
      arrows(x0=layout1[,1],
         x1=layout2[,1],
         y0=layout1[,2],
         y1=layout2[,2],
         length=min(par("pin"))/65);
   }
   g <- igraph::set_graph_attr(g,
      name="layout",
      value=layout);
   return(g);
}


#' Get Cnet node set by connected Sets
#'
#' Get Cnet node set by connected Sets
#'
#' This function operates on a Cnet `igraph` object,
#' distinguished by node attribute `"nodeType"` with
#' value `"Gene"` for Gene nodes, and `"Set"` for Set nodes.
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
#' @family jam cnet igraph functions
#'
#' @param g `igraph` object specifically containing Cnet plot data,
#'    with node attribute `"nodeType"` that has values `c("Gene", "Set")`.
#' @param set_nodes `character` vector of one or more set names, which are
#'    defined here as the names of the nodes with `nodeType="Set"`
#'    to which each node with `nodeType="Gene"` is connected.
#'    For example `set_nodes=c("A", "B")` will return all Gene nodes
#'    that are connected only to Set nodes `"A"` and `"B"`. If
#'    `set_nodes=NULL` then this function returns a `list` with
#'    all set_nodes observed.
#' @param sep `character` string used as a delimited when combining
#'    `set_nodes` into a name.
#' @param ... additional arguments are ignored.
#'
#' @return `character` vector when argument `set_nodes` is supplied
#'    as a `character` vector; when `set_nodes` is `NULL`,
#'    this function returns a `list` of `character` vectors whose
#'    names represent the Set node names connected to each
#'    Gene node, where multiple Set node names are
#'    delimited by `sep`.
#'
#' @export
get_cnet_nodeset <- function
(g,
 set_nodes=NULL,
 sep=",",
 ...)
{
   ## comma-delimited neighboring nodes for each node
   neighborG <- jamba::cPasteS(sep=sep,
      lapply(seq_len(igraph::vcount(g)), function(v){
         names(igraph::neighbors(g, v, mode="all"));
      }));
   names(neighborG) <- igraph::V(g)$name;
   if (length(set_nodes) == 0) {
      useG <- neighborG[igraph::V(g)$nodeType %in% "Gene"];
      useG <- split(names(useG), useG);
   } else {
      set_nodes_v <- jamba::cPasteS(set_nodes,
         sep=sep);
      useG <- names(neighborG)[which(neighborG %in% set_nodes_v)];
   }
   return(useG);
}


#' Adjust Cnet node set
#'
#' Adjust Cnet node set
#'
#' This function is intended to help adjust a node set,
#' defined as a set of `"Gene"` nodes that all connect
#' to the same `"Set"` nodes. These nodes are usually clustered
#' together in a Cnet plot, and can be manipulated as a
#' group using this function. Use `get_cnet_nodeset(g)`
#' to get all node sets from a Cnet `igraph` object.
#'
#' Specifically, this function can be used to expand,
#' rotate, and move the cluster of nodes in a Cnet
#' network layout. The purpose is usually to help the
#' visual clarity of the plot, and to reduce node label
#' overlaps.
#'
#' @family jam cnet igraph functions
#'
#' @examples
#' # the examples below are for internal data
#' # and will be rewritten to work with example data
#' if (1 == 2) {
#' g <- mem_dmjdm_gp_20_plots_jikl$cnet_collapsed_set;
#' set_nodes <- c("J", "L")
#' g2 <- adjust_cnet_nodeset(g, set_nodes=c("J", "L"), expand=0.5)
#'
#' jam_igraph(g)
#' g2 <- adjust_cnet_nodeset(g, set_nodes=c("I", "K", "L"), x=0.01, y=-0.01, expand=0.2)
#' g2 <- adjust_cnet_nodeset(g2, set_nodes=c("I", "J", "K", "L"), x=-0.08, y=-0.05, expand=0.2)
#' g2 <- adjust_cnet_nodeset(g2, set_nodes=c("I", "L"), x=-0.05, y=0.05, expand=0.3)
#' g2 <- adjust_cnet_nodeset(g2, set_nodes=c("I", "J", "L"), x=-0.02, y=0.03, expand=0.3)
#' g2 <- adjust_cnet_nodeset(g2, set_nodes=c("I", "K"), x=0.03, y=0.09, expand=0.)
#' g2 <- adjust_cnet_nodeset(g2, set_nodes=c("I", "J", "K"), x=0.02, y=-0.01, expand=0.)
#' g2 <- adjust_cnet_nodeset(g2, set_nodes=c("I", "J"), x=0.01, y=-0.00, expand=0.1)
#' g2 <- adjust_cnet_nodeset(g2, set_nodes=c("J"), x=-0.03, y=-0.02, expand=0.)
#' jam_igraph(g2)
#' }
#'
#' @param g `igraph` Cnet plot, with vertex attribute `nodeType` that
#'    contains values `"Gene"` and `"Set"`, as produced by `memIM2cnet()`.
#' @param set_nodes `character` vector of one or more set names, which are
#'    defined here as the names of the nodes with `nodeType="Set"`
#'    to which each node with `nodeType="Gene"` is connected.
#'    For example `set_nodes=c("A", "B")` will return all Gene nodes
#'    that are connected only to Set nodes `"A"` and `"B"`. Alternatively,
#'    one can refer to a cluster by the node name of one member of the
#'    cluster, as a convenience option.
#' @param x,y `numeric` indicating the relative amount to move nodes
#'    in the node set, based upon the overall range of x and y values
#'    in the layout coordinates.
#' @param expand `numeric` the relative amount to expand node coordinates,
#'    where `expand=0` does no expansion, `expand=1` will expand coordinates
#'    100% of the distance to the center of the cluster, twice the original
#'    distance; and `expand=-1` will expand coordinates to half the
#'    distance to the center.
#' @param rotate_degrees `numeric` value indicating the rotation in
#'    degrees to adjust the node coordinates, relative to the center of
#'    the coordinates.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
adjust_cnet_nodeset <- function
(g,
 set_nodes=NULL,
 x=0,
 y=0,
 expand=0,
 rotate_degrees=0,
 verbose=FALSE,
 ...)
{
   ## List input for bulk operations through one function call
   if (is.list(set_nodes)) {
      if (verbose) {
         jamba::printDebug("adjust_cnet_nodeset(): ",
            "Bulk operation for ",
            length(set_nodes),
            "nodesets.");
      }
      n <- length(set_nodes);
      x <- rep(x, length.out=n)
      y <- rep(y, length.out=n)
      expand <- rep(expand, length.out=n)
      rotate_degrees <- rep(rotate_degrees, length.out=n)
      for (i in seq_len(n)) {
         if (verbose) {
            jamba::printDebug("adjust_cnet_nodeset(): ",
               "iteration i:", i)
         }
         g <- adjust_cnet_nodeset(g=g,
            set_nodes=set_nodes[[i]],
            x=x[[i]],
            y=y[[i]],
            expand=expand[[i]],
            rotate_degrees=rotate_degrees[[i]],
            verbose=verbose,
            ...)
      }
      return(g)
   }

   ## comma-delimited neighboring nodes for each node
   set_nodes_v <- jamba::cPasteS(set_nodes);

   # get nodes in the defined nodeset
   cnet_nodesets <- get_cnet_nodeset(g);
   if (set_nodes_v %in% names(cnet_nodesets)) {
      useG <- cnet_nodesets[[set_nodes_v]];
      if (verbose) {
         jamba::printDebug("adjust_cnet_nodeset(): ",
            "Matched nodeset ",
            set_nodes_v);
      }
   } else if (any(set_nodes_v %in% unlist(cnet_nodesets))) {
      set_nodes_v2 <- rep(names(cnet_nodesets),
         lengths(cnet_nodesets))[match(set_nodes_v, unlist(cnet_nodesets))];
      if (verbose) {
         jamba::printDebug("adjust_cnet_nodeset(): ",
            "Matched nodeset ",
            set_nodes_v2,
            " using nodeset member ",
            set_nodes_v);
      }
      useG <- cnet_nodesets[[set_nodes_v2]];
   } else {
      stop(paste0("The set_nodes provided do not match nodeset names, nor node members."))
   }
   if (length(useG) == 0) {
      warning(paste("No nodes found with set_nodes:", set_nodes_v));
      return(g);
   }

   #
   if (verbose) {
      jamba::printDebug("adjust_cnet_nodeset(): ",
         "Adjusting nodes connected to ", set_nodes_v,
         ":",
         useG);
   }

   # nudge nodes if given x,y non-zero
   if (any(x != 0) || any(y != 0)) {
      if (verbose) {
         jamba::printDebug("adjust_cnet_nodeset(): ",
            "Applying nudge_igraph_node() on x,y coordinates.");
      }
      g <- nudge_igraph_node(g,
         nodes=useG,
         x=rep(x, length.out=length(useG)),
         y=rep(y, length.out=length(useG)),
         use_grep=FALSE,
         verbose=verbose);
   }

   # optionally expand
   if (any(expand != 0) && length(useG) > 1) {
      if (verbose) {
         jamba::printDebug("adjust_cnet_nodeset(): ",
            "Applying expand on x,y coordinates.");
      }
      layout <- igraph::graph_attr(g, "layout");
      rownames(layout) <- igraph::V(g)$name;
      use_xy <- layout[useG,,drop=FALSE];
      center_xy <- matrix(ncol=2,
         byrow=TRUE,
         colMeans(use_xy),
         nrow=nrow(use_xy));
      diff_xy <- use_xy - center_xy;

      # expand by a numeric factor
      expand <- ifelse(expand < 0,
         1 / (abs(expand[expand < 0]) + 1),
         expand + 1)
      expand <- rep(rep(expand, length.out=ncol(use_xy)),
         each=nrow(diff_xy));
      new_xy <- diff_xy * (expand) + center_xy;

      layout[useG,] <- new_xy;
      g <- igraph::set_graph_attr(g,
         name="layout",
         value=layout);
      #xrange <- max(c(
      #   diff(range(layout[,1])),
      #   diff(range(layout[,2]))));
      #yrange <- xrange;
   }

   # optionally rotate
   rotate_degrees <- head(rotate_degrees, 1);
   if (rotate_degrees != 0) {
      if (verbose) {
         jamba::printDebug("adjust_cnet_nodeset(): ",
            "Applying rotate_degrees:",
            rotate_degrees);
      }
      layout <- igraph::graph_attr(g, "layout");
      rownames(layout) <- igraph::V(g)$name;
      use_xy <- layout[useG,,drop=FALSE];
      center_xy <- matrix(ncol=2,
         byrow=TRUE,
         colMeans(use_xy),
         nrow=nrow(use_xy));
      diff_xy <- use_xy - center_xy;
      co <- cos(-jamba::deg2rad(rotate_degrees));
      si <- sin(-jamba::deg2rad(rotate_degrees));
      new_xy <- cbind(
         co * diff_xy[,1] - si * diff_xy[,2],
         si * diff_xy[,1] + co * diff_xy[,2]) + center_xy;
      layout[useG,] <- new_xy;
      g <- igraph::set_graph_attr(g,
         name="layout",
         value=layout);
   }

   return(g)
}


#' Adjust Set nodes then relayout Gene nodes
#'
#' Adjust Cnet Set nodes then relayout Gene nodes
#'
#' This function operates on a Cnet `igraph` object,
#' distinguished by node attribute `"nodeType"` with
#' value `"Gene"` for Gene nodes, and `"Set"` for Set nodes.
#'
#' This function is intended to help move a `Set` node
#' to improve visual spacing between nodes, then it
#' re-positions only the `Gene` nodes using `layout_with_qfr()`,
#' keeping the `Set` nodes in fixed positions.
#'
#' It calls `nudge_igraph_node()` using arguments `nodes,x,y`,
#' then defines the layout coordinates of Set nodes
#' as `constraints` that therefore are not allowed to change
#' when calling `layout_with_qfr()`. Note that Gene node
#' coordinates are allowed to change, even if Gene nodes
#' were included in `nodes`.
#'
#' @family jam cnet igraph functions
#'
#' @param g `igraph` Cnet object
#' @param nodes `character` vector of one or more Set nodes.
#' @param x,y `numeric` values passed to `nudge_igraph_node()`.
#' @param use_grep `logical` passed to `nudge_igraph_node()`.
#' @param do_reorder `logical` indicating whether nodes should
#'    be re-positioned within each subcluster by calling
#'    `reorderIgraphNodes()`.
#' @param spread_labels `logical` indicating whether labels
#'    should be re-oriented around each node by calling
#'    `spread_igraph_labels()`.
#' @param repulse `numeric` value passed to `layout_with_qfr()`
#'    to indicate the repulsion force. Typical values range between
#'    `3` for loosely-packed nodes, and `5` or higher for more
#'    closely-packed nodes.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to functions
#'    `nudge_igraph_node()`, `layout_with_qfr()`, `spread_igraph_labels()`,
#'    `reorderIgraphNodes()`.
#'
#' @export
adjust_cnet_set_relayout_gene <- function
(g,
 nodes=NULL,
 x=0,
 y=0,
 use_grep=TRUE,
 do_reorder=TRUE,
 spread_labels=TRUE,
 repulse=4,
 verbose=FALSE,
 ...)
{
   # call nudge_igraph_node()
   if (verbose) {
      jamba::printDebug("adjust_cnet_set_relayout_gene(): ",
         "Calling nudge_igraph_node()");
   }
   xb <- nudge_igraph_node(g,
      nodes=nodes,
      x=x,
      y=y,
      use_grep=use_grep,
      ...)
   g_types <- igraph::V(g)$nodeType;

   # get layout
   layout4b <- igraph::graph_attr(xb, "layout");
   if (length(layout4b) == 0) {
      stop("This function assumes a node layout already exists.");
   }
   rownames(layout4b) <- igraph::V(xb)$name;

   # give Gene nodes coordinates NA so they will be repositioned
   layout4b[g_types %in% "Gene",] <- NA;

   # call layout_with_qfr
   if (verbose) {
      jamba::printDebug("adjust_cnet_set_relayout_gene(): ",
         "Calling layout_with_qfr()");
   }
   layout4c <- layout_with_qfr(xb,
      constraints=layout4b,
      repulse=repulse,
      ...)
   rownames(layout4c) <- rownames(layout4b)
   # apply to igraph object
   igraph::graph_attr(xb, "layout") <- layout4c;

   # optionally call spread_igraph_labels()
   if (spread_labels) {
      if (verbose) {
         jamba::printDebug("adjust_cnet_set_relayout_gene(): ",
            "Calling spread_igraph_labels()");
      }
      xb <- spread_igraph_labels(xb,
         do_reorder=do_reorder,
         force_relayout=FALSE,
         verbose=verbose,
         ...)
   } else if (do_reorder) {
      if (verbose) {
         jamba::printDebug("adjust_cnet_set_relayout_gene(): ",
            "Calling reorderIgraphNodes()");
      }
      xb <- reorderIgraphNodes(xb,
         ...)
   }

   return(xb)
}

#' Apply minimum node spacing for each Cnet node set
#'
#' Apply minimum node spacing for each Cnet node set
#'
#' This function automates the process of spacing nodesets
#' in a Cnet plot by calling `adjust_cnet_nodeset()` on
#' all nodesets whose median nearest node distance is
#' above or below a threshold.
#'
#' By default nodes must be at least `percent_spacing`
#' distance apart, relative to the maximum x-axis and
#' y-axis layout coordinate range. Roughly speaking,
#' nodes are clearly visible with distance between nodes
#' of at least 3 percent of the coordinate range.
#'
#' Nodes are not compressed down to this distance unless
#' using the argument `apply_negative=TRUE`.
#'
#' Note that this function does not detect when two
#' clusters overlap. For other specific layout adjustments
#' use `adjust_cnet_nodeset()`. In future, this function
#' may attempt to re-position nodesets relative to
#' each other, but that process tends to be iterative
#' and is not guaranteed to have a solution.
#'
#' @family jam cnet igraph functions
#'
#' @return `igraph` Cnet object, whose nodes are adjusted
#'    according to the arguments given.
#'
#' @param cnet `igraph` object that contains attribute
#'    `"nodeType"` and values `nodeType=c("Gene", "Set")`,
#'    consistent with use of `get_cnet_nodesets()` and
#'    `adjust_cnet_nodeset()`.
#'
#' @export
apply_nodeset_spacing <- function
(cnet,
 percent_spacing=3,
 apply_negative=FALSE,
 ...)
{
   # get node layout
   layout_xy <- igraph::graph_attr(cnet, "layout")
   if (length(layout_xy) == 0) {
      stop("apply_nodeset_spacing() requires layout stored in graph_attr(cnet, 'layout')");
   }
   rownames(layout_xy) <- igraph::V(cnet)$name;

   # get nodesets
   cnet_nodesets <- get_cnet_nodeset(cnet)

   # get x- and y-range for overall scale of the layout
   max_xy <- max(c(
      diff(range(layout_xy[,1])),
      diff(range(layout_xy[,2]))));

   # determine nodeset spacing relative to layout size
   cnet_spacing <- sapply(cnet_nodesets, function(i){
      use_xy <- layout_xy[i,,drop=FALSE]
      i1 <- dist(use_xy / max_xy * 100)
      i2 <- as.matrix(i1);
      diag(i2) <- NA;
      median(matrixStats::colMins(i2, na.rm=TRUE))
   })

   # determine nodesets to expand
   adjust_nodesets <- names(which(cnet_spacing < (percent_spacing/1.05)));
   if (apply_negative) {
      adjust_nodesets_neg <- names(which(!is.infinite(cnet_spacing) &
            cnet_spacing > (percent_spacing * 1.05)))
      adjust_nodesets <- c(adjust_nodesets,
         adjust_nodesets_neg);
   }
   if (length(adjust_nodesets) == 0) {
      return(cnet)
   }
   adjust_nodesets_expand <- percent_spacing / cnet_spacing[adjust_nodesets];
   adjust_nodesets_expand <- ifelse(adjust_nodesets_expand >= 1,
      adjust_nodesets_expand - 1,
      -1 * (1/adjust_nodesets_expand - 1));

   # apply adjustments
   cnet_new <- adjust_cnet_nodeset(cnet,
      set_nodes=strsplit(adjust_nodesets, ","),
      x=0,
      y=0,
      expand=adjust_nodesets_expand)
   return(cnet_new)
}

