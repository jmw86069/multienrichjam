
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
#' @param nodes_xy `list` alternative to using arguments `nodes,x,y`.
#'    This argument assumes a `list` of `numeric` adjustments to `x` and `y`,
#'    where the `names(nodes_xy)` are node names. For example:
#'    * `nodes_xy=list(APOE=c(x=0.05, y=-0.01), GAPDH=c(x=0.1, y=0.0)`
#'    however the `numeric` vector does not need to contain names.
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
 nodes=NULL,
 x=0,
 y=0,
 nodes_xy=NULL,
 use_grep=FALSE,
 aspect=1,
 debug=FALSE,
 verbose=FALSE,
 ...)
{
   # expand x,y to length(nodes)
   x <- rep(x, length.out=length(nodes));
   y <- rep(y, length.out=length(nodes));

   # add optional nodes_xy
   if (length(nodes_xy) > 0) {
      nodes <- c(nodes, names(nodes_xy));
      x <- c(x, sapply(nodes_xy, function(ixy){ixy[1]}));
      y <- c(y, sapply(nodes_xy, function(ixy){ixy[2]}));
   }

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
 filter_set_only=TRUE,
 ...)
{
   # alternative approach for speed
   # 1. assemble edgelist data.frame
   gel1 <- igraph::as_edgelist(g);
   gel2 <- gel1[,2:1, drop=FALSE];
   gel <- data.frame(check.names=FALSE,
      unique(rbind(gel1, gel2)));
   colnames(gel) <- c("A", "B")

   # optionally limit output to nodes connected to nodeType="Set" nodes
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

   if (length(set_nodes) > 0) {
      set_nodes_v <- jamba::cPasteS(set_nodes,
         sep=sep)
      cnet_nodesets <- cnet_nodesets[names(cnet_nodesets) %in% set_nodes_v];
   }
   return(cnet_nodesets);

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
#' @param cnet_nodesets `list` optional to provide the known nodesets
#'    upfront, in order to avoid re-determining the nodesets.
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
 percent_spacing=NULL,
 rotate_degrees=0,
 cnet_nodesets=NULL,
 verbose=FALSE,
 ...)
{
   # define cnet_nodesets if not supplied
   if (length(cnet_nodesets) == 0) {
      cnet_nodesets <- get_cnet_nodeset(g);
   }

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
      if (length(percent_spacing) > 0) {
         percent_spacing <- rep(percent_spacing, length.out=n)
      }

      rotate_degrees <- rep(rotate_degrees, length.out=n)
      for (i in seq_len(n)) {
         if (verbose) {
            jamba::printDebug("adjust_cnet_nodeset(): ",
               "nodeset number i:", i)
         }
         g <- adjust_cnet_nodeset(g=g,
            set_nodes=set_nodes[[i]],
            x=x[[i]],
            y=y[[i]],
            expand=expand[[i]],
            percent_spacing=percent_spacing[[i]],
            rotate_degrees=rotate_degrees[[i]],
            cnet_nodesets=cnet_nodesets,
            verbose=verbose,
            ...)
      }
      return(g)
   }

   ## comma-delimited neighboring nodes for each node
   set_nodes_v <- jamba::cPasteS(set_nodes);

   # get nodes in the defined nodeset
   if (set_nodes_v %in% names(cnet_nodesets)) {
      useG <- cnet_nodesets[[set_nodes_v]];
      if (verbose) {
         jamba::printDebug("adjust_cnet_nodeset(): ",
            "Matched nodeset: ",
            set_nodes_v);
      }
   } else if (any(set_nodes_v %in% unlist(cnet_nodesets))) {
      cnet_nodesets_df <- data.frame(
         nodeset=rep(names(cnet_nodesets), lengths(cnet_nodesets)),
         node=unname(unlist(cnet_nodesets)))
      use_nodeset <- unique(subset(cnet_nodesets_df, node %in% set_nodes_v)$nodeset);
      if (length(use_nodeset) > 1) {
         stop("Handle multiple nodesets with list input.")
      }
      useG <- cnet_nodesets[[use_nodeset]];
      if (verbose) {
         jamba::printDebug("adjust_cnet_nodeset(): ",
            "Matched nodeset: ",
            use_nodeset,
            " using nodeset member: ",
            set_nodes_v);
      }
   } else {
      # jamba::printDebug("cnet_nodesets:");
      # print(cnet_nodesets);
      warning(paste0("The set_nodes '",
         set_nodes_v,
         "' provided do not match nodeset names, nor node members."));
      return(g);
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

   # optional percent_spacing, convert to expand
   if (length(percent_spacing) > 0) {
      if (verbose) {
         jamba::printDebug("adjust_cnet_nodeset(): ",
            "Determining percent_spacing via apply_nodeset_spacing()");
      }
      current_spacing <- apply_nodeset_spacing(g,
         nodesets=set_nodes_v,
         percent_spacing=percent_spacing,
         debug=TRUE,
         verbose=FALSE,
         ...)
      # define expand
      if (verbose) {
         jamba::printDebug("adjust_cnet_nodeset(): ",
            "current_spacing:");
         print(current_spacing);
      }
      expand <- current_spacing$expand;
   }

   # optionally expand
   if (any(expand != 0) && length(useG) > 1) {
      if (verbose) {
         jamba::printDebug("adjust_cnet_nodeset(): ",
            "expand x,y coordinates.");
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
      if (verbose) {
         jamba::printDebug("adjust_cnet_nodeset(): ",
            "expand x,y coordinates - complete.");
      }
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
#' When a nodeset contains only one node, there is no calculated
#' distance so it is reported as `NA` when `verbose=TRUE` or
#' `debug=TRUE`. However for internal purposes the nodeset is
#' treated as if its nodeset spacing is exactly `percent_spacing`
#' so that no subsequent adjustments are performed on the nodeset.
#'
#' @family jam cnet igraph functions
#'
#' @return `igraph` Cnet object, whose nodes are adjusted
#'    according to the arguments given.
#'
#' @param cnet `igraph` object that contains attribute
#'    `"nodeType"` and values `nodeType=c("Gene", "Set")`,
#'    consistent with use of `get_cnet_nodeset()` and
#'    `adjust_cnet_nodeset()`.
#'    The node layout is expected to be stored in
#'    `igraph::graph_attr(cnet, "layout")`.
#' @param percent_spacing `numeric` indicating the target spacing in relation
#'    to the maximum x-axis and y-axis span for the layout coordinates.
#'    For example `percent_spacing=4` will require that the median spacing
#'    between nodes in each nodeset is at least 4 percent of the larger of
#'    the x-axis and y-axis ranges for the layout of input `cnet`.
#' @param apply_negative `logical` indicating whether to shrink node spacing
#'    when the spacing is larger than the target `percent spacing`. By default
#'    `apply_negative=FALSE` which allows nodes to exceed the spacing, and
#'    only adjusts node spacing when a nodeset is below this threshold.
#' @param metric `character` string indicating the metric used to summarize
#'    nodeset spacing:
#'    * `"median"` - takes the median spacing between nodes in the nodeset,
#'    so the end result has at least half nodes with `percent_spacing`
#'    distance, "typical distance between nodes". This option is useful
#'    when there may be one node pair close together, but most nodes have
#'    consistent distance from their nearest node neighbor.
#'    * `"minimum"` - takes the minimum spacing between nodes in the nodeset,
#'    so the end result requires all nodes to have at least `percent_spacing`
#'    distance.
#'    * `function` - a `function` can be provided which takes a `numeric`
#'    vector of node distances, and produces a single `numeric` summary value,
#'    which will be scaled relative to the maximum xy axis range. This option
#'    would allow using `mean()`, or a quantile function such as
#'    `function(x){quantile(x, 0.3)}`.
#' @param nodesets `character` with optional subset of nodesets to adjust,
#'    or the default `nodesets=NULL` will adjust all nodesets. Values can be
#'    specified either in form or nodeset names, or any one or more members
#'    of a nodeset.
#' @param tolerance `numeric` value used to apply a small tolerance to the
#'    distance threshold, for example the default `tolerance=1.05` requires
#'    nodes to be more than 5% lower than the actual threshold before
#'    adjusting nodes. This additional activation threshold is experimental
#'    and may not be necessary. It is mainly intended to reduce calculations
#'    for nodesets that are only a small rounding error from being acceptable,
#'    which can happen when the overall plot range is adjusted during rounds
#'    of small node coordinate adjustments like with this function,
#'    `nudge_igraph_node()`, `apply_nodeset_spacing()`, or
#'    `adjust_cnet_set_relayout_gene()`.
#' @param cnet_nodesets `list` optional to provide the known nodesets
#'    upfront, in order to avoid re-determining the nodesets.
#' @param debug `logical` indicating whether to print debug information
#'    without adjusting the nodeset spacing. This option is useful to see
#'    summary information about nodeset spacing before choosing the
#'    `percent_spacing` threshold.
#' @param verbose `logical` indicating whether to print verbose output.
#'    The output includes text printed with `debug=TRUE` except that nodes
#'    are adjusted when `debug=FALSE`.
#' @param ... additional arguments are ignored.
#'
#' @export
apply_nodeset_spacing <- function
(cnet,
 percent_spacing=3,
 apply_negative=FALSE,
 metric=c("median", "minimum"),
 nodesets=NULL,
 tolerance=1.05,
 cnet_nodesets=NULL,
 debug=FALSE,
 verbose=FALSE,
 ...)
{
   # validate metric argument
   if (length(metric) == 0) {
      metric <- "median"
   }
   if ("character" %in% class(metric)) {
      metric <- match.arg(metric);
      if ("median" %in% metric) {
         metric <- median;
      } else if ("min" %in% metric) {
         metric <- min;
      }
   } else if (!"function" %in% class(metric)) {
      stop("metric must either be class 'character' or 'function'");
   }

   # get node layout
   layout_xy <- igraph::graph_attr(cnet, "layout")
   if (length(layout_xy) == 0) {
      stop("apply_nodeset_spacing() requires layout stored in graph_attr(cnet, 'layout')");
   }
   rownames(layout_xy) <- igraph::V(cnet)$name;

   # get nodesets
   if (length(cnet_nodesets) == 0) {
      cnet_nodesets <- get_cnet_nodeset(cnet)
   }
   if (length(nodesets) > 0) {
      ikeep_nodeset <- sapply(jamba::nameVectorN(cnet_nodesets), function(iname){
         if (iname %in% nodesets) {
            return(TRUE)
         }
         inodeset <- cnet_nodesets[[iname]];
         if (any(inodeset %in% nodesets)) {
            return(TRUE)
         }
         FALSE
      })
      cnet_nodesets <- cnet_nodesets[ikeep_nodeset];
      keep_nodesets <- names(cnet_nodesets);
      if (length(keep_nodesets) == 0) {
         stop("After filtering by nodesets argument, no data remained.")
      }
      if (verbose) {
         jamba::printDebug("apply_nodeset_spacing(): ",
            "Subset nodesets to keep: ",
            sep="; ",
            paste0("'", keep_nodesets, "'"))
      }
   }

   # get x- and y-range for overall scale of the layout
   max_xy <- max(c(
      diff(range(layout_xy[,1])),
      diff(range(layout_xy[,2]))));
   if (verbose) {
      jamba::printDebug("apply_nodeset_spacing(): ",
         "Plot range max_xy: ", max_xy);
   }

   # determine nodeset spacing relative to layout size
   cnet_spacing <- sapply(jamba::nameVectorN(cnet_nodesets), function(iname){
      i <- cnet_nodesets[[iname]];
      if (length(i) <= 1) {
         imed <- NA;
         imin <- NA;
         imetric <- NA;
      } else {
         use_xy <- layout_xy[i,,drop=FALSE]
         # calculate distance then scale as percentage of the max xy range
         i1 <- dist(use_xy / max_xy * 100)
         i2 <- as.matrix(i1);
         diag(i2) <- NA;
         nearest_distances <- jamba::rmNA(matrixStats::colMins(i2, na.rm=TRUE));
         imed <- median(nearest_distances);
         imin <- min(nearest_distances);
         imetric <- metric(nearest_distances);
      }
      # extensive debug output logic here
      if (verbose) {
         jamba::printDebug("apply_nodeset_spacing(): ",
            c("Median spacing for nodeset '", iname, "' = "),
            imed,
            ", (minimum spacing = ", imin, ")",
            sep="");
         if (!all(imetric == imin | is.na(imetric)) && !all(imetric == imed | is.na(imetric))) {
            jamba::printDebug("apply_nodeset_spacing(): ",
               c("spacing by metric for nodeset '", iname, "' = "),
               imetric,
               sep="");
         }
      }
      if (length(i) <= 1) {
         imed <- percent_spacing;
         imin <- percent_spacing;
         imetric <- percent_spacing;
      }
      unname(imetric);
   })
   cnet_spacing_df <- data.frame(check.names=FALSE,
      nodeset=names(cnet_spacing),
      cnet_spacing=cnet_spacing,
      expand=0)

   # determine nodesets to expand
   # use tolerance to allow some small wiggle room
   adjust_nodesets <- names(which(cnet_spacing < (percent_spacing / tolerance)));

   adjust_nodesets_expand <- percent_spacing / cnet_spacing;
   names(adjust_nodesets_expand) <- names(cnet_spacing);
   adjust_nodesets_expand <- ifelse(adjust_nodesets_expand >= 1,
      adjust_nodesets_expand - 1,
      -1 * (1/adjust_nodesets_expand - 1));

   # optionally ignore nodesets that are already above the threshold
   if (apply_negative) {
      adjust_nodesets_neg <- names(which(!is.infinite(cnet_spacing) &
            cnet_spacing > (percent_spacing * 1.05)))
      adjust_nodesets <- c(adjust_nodesets,
         adjust_nodesets_neg);
   }
   if (length(adjust_nodesets) > 0) {
      imatch <- match(adjust_nodesets, rownames(cnet_spacing_df));
      cnet_spacing_df[imatch, "expand"] <- adjust_nodesets_expand[adjust_nodesets];
   }
   if (debug) {
      return(cnet_spacing_df)
   }
   if (length(adjust_nodesets) == 0) {
      if (verbose) {
         jamba::printDebug("apply_nodeset_spacing(): ",
            "No nodesets required adjustment.");
      }
      return(cnet)
   }


   # apply adjustments
   if (verbose) {
      jamba::printDebug("apply_nodeset_spacing(): ",
         "Calling adjust_cnet_nodeset()");
      jamba::printDebug("apply_nodeset_spacing(): ",
         "expand: ",
         adjust_nodesets_expand[adjust_nodesets]);
   }
   cnet_new <- adjust_cnet_nodeset(cnet,
      set_nodes=strsplit(adjust_nodesets, ","),
      x=0,
      y=0,
      expand=adjust_nodesets_expand[adjust_nodesets],
      cnet_nodesets=cnet_nodesets,
      verbose=verbose)
   return(cnet_new)
}

