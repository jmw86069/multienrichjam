
#' Jam igraph vectorized plot function (internal)
#'
#' Jam igraph vectorized plot internal function called by `jam_igraph()`
#'
#' Note that this function is intended to be called by `jam_igraph()`,
#' and is an internal function not intended to be called directly.
#'
#' The `jam_igraph()` handles the overall plot equivalent of
#' `igraph::plot.igraph()`, however it calculates layout coordinates,
#' and defines more useful x- and y-axis ranges, and then
#' adjusts node and label sizes relevant to the layout data range.
#' Specifically `vertex.size=15` is only useful when the layout range
#' is rescaled between -1 and 1; however when using `jam_igraph()`
#' the vertex is scaled relative to the actual layout ranges.
#'
#' The steps here are a reproduction of `igraph:::plot.igraph()` with
#' four changes:
#'
#' 1. Default `rescale=FALSE`, and `asp=1` which means igraph layout is
#' drawn true to the layout coordinates without distortion. To use
#' default `igraph::plot.igraph()` behavior, use `rescale=TRUE`.
#' The new default may not be appropriate for bipartite layout
#' algorithms that generate two columns, and seems most useful
#' with organic layouts where aspect ratio 1 helps convey important
#' meaning in the graph structure, namely by enforcing consistent
#' x- and y-axis visual distance between nodes.
#'
#'    * Related: the `xlim` and `ylim` values are automatically adjusted
#'    to include the layout coordinate range. The default
#'    `igraph::plot.igraph(..., rescale=FALSE)` does not adjust the
#'    `xlim` and `ylim` ranges, which can be problematic when supplying
#'    layout as a function, and therefore the output node coordinates
#'    are not known until the plot rendering step.
#'
#' When `vectorized_node_shapes=TRUE` by default:
#'
#' 2. When there are multiple different vertex `"shape"` attributes, the
#' nodes are rendered vectorized one shape at a time. The original
#' `igraph::plot.igraph()` draws each individual vertex one by one,
#' which is substantially slower (minutes compared to 1-2 seconds)
#' for large `igraph` objects.
#' 3. When there are multiple font families, the default plot function
#' draws each label one by one. The `jam_plot_igraph()` draws
#' labels in groups of font family, in order to comply with limitations
#' in `graphics::text()`. This situation is fairly rare, however
#' the speed improvement is substantial, again roughly minutes down
#' to 1-2 seconds.
#'
#' The fourth difference involves edge bundling:
#'
#' 4. When `edge_bundling` is used, it renders edges differently
#' than the approach in `igraph::plot.igraph()`, by drawing curved
#' splines for each bundle of edges.
#'
#' Some other distinctive features include:
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
#' rendering time.
#'
#' Default colors for marked node groups `mark.col` and `mark.border`
#' when not defined upfront, will call `colorjam::rainbowJam()`
#' and not `grDevices::rainbow(). The `colorjam::rainbowJam()`
#' produces more visually distinct categorical colors.
#' This behavior can be controlled by supplying a `character`
#' vector with specific colors for `mark.col` and `mark.border`. Note
#' that the border should match the colors, or it can be set to `"grey45"`
#' for a generally visible border.
#'
#' Optional argument `nodegroups` can be supplied, which is a `list`
#' of vectors, where each vector represents a group of nodes. The
#' `nodegroups` can be used with `edge_bundling="nodegroups"` to
#' define custom edge bundling.
#'
#' Finally, individual plot components can be individually disabled:
#'
#' * `render_nodes=FALSE`
#' * `render_edges=FALSE`
#' * `render_groups=FALSE`
#' * `render_nodelabels=FALSE`
#'
#'
#' Note that this function is not called by default, and is only called
#' by `multienrichjam::jam_igraph()`.
#'
#'
#' All other arguments are documented in `igraph::plot.igraph()`.
#'
#' @family jam igraph internal functions
#'
#' @inheritParams igraph::plot.igraph
#' @param xlim,ylim default x and y axis limits. When either value is `NULL`
#'    the range is defined by the layout coordinate ranges, respectively,
#'    then expanded by adding `expand` to each side of the range.
#' @param mark.alpha `numeric` value between 0 (transparent) and 1 (opaque)
#'    indicating the transparency of `mark.col` color fill values,
#'    used only when `mark.groups` is defined, and `mark.col` is not defined.
#' @param mark.lwd,mark.lty line with and line type parameters for each
#'    `mark.groups` polygon.
#' @param mark.cex `numeric` adjustment for mark label font size, used
#'    when `mark.groups` is supplied and has `names(mark.groups)`.
#' @param mark.x.nudge,mark.y.nudge `numeric` values in units of the
#'    maximum x-axis or y-axis range for the layout coordinates,
#'    used to adjust each label displayed when `names(mark.groups)`
#'    is defined. These arguments are passed to `make_point_hull()`
#'    as `label.x.nudge`, `label.y.nudge`, respectively.
#' @param pie_to_jampie `logical` indicating whether to convert
#'    vertex shape `"pie"` to `"jampie"` in order to use vectorized
#'    plotting.
#' @param use_shadowText `logical` indicating whether to use
#'    `jamba::shadowText()` instead of `graphics::text()`, in order
#'    to render text labels with a subtle shadow-like outline around
#'    each label. This change improves legibility of labels at
#'    the expense of slightly longer plot rendering time.
#' @param vectorized_node_shapes `logical` indicating whether to plot
#'    vertex node shapes using vectorized operations. It is substantially
#'    faster, however the one drawback is that nodes are plotted in
#'    order of their shape, which affects the positioning of nodes
#'    when there are node overlaps. This tradeoff is relatively minor,
#'    and it is recommended either to reposition nodes to reduce or
#'    prevent overlaps, or adjust node sizes to reduce overlaps.
#' @param edge_bundling `character` string or `function`, where:
#'    * `"default"` will try to detect an appropriate method: when
#'    `nodegroups` or `mark.groups` are defined, it chooses the matching
#'    option (see below); otherwise it chooses `"connections"`.
#'    * `"none"` will perform no edge bundling. This method is best when
#'    rendering straight edges, or for rendering multiple identical edges
#'    with curvature as defined by `igraph::igraph.plotting()`.
#'    * `"connections"` will perform graph edge bundling by
#'    shared connections by calling `edge_bundle_bipartite()` then
#'    `edge_bundle_nodegroups()`. This option is particularly good
#'    for bipartite graphs such as concept networks (cnet plots).
#'    * `"mark.groups"` will perform graph edge bundling
#'    using `mark.groups` by calling `edge_bundle_nodegroups()`.
#'    This option is equivalent to `"nodegroups"` except that it
#'    uses `mark.groups` to define node groupings.
#'    * `"nodegroups"` will perform graph edge bundling
#'    using `nodegroups` by calling `edge_bundle_nodegroups()`.
#'    This option is equivalent to `"mark.groups"` except that it
#'    uses `nodegroups` to define node groupings.
#'    * `function` will call a custom edge bundling function using
#'    the `igraph` object `x` and the igraph parameters `param`
#'    as input. This output is currently untested, and is intended
#'    to enable alternative edge bundling functions which may exist
#'    outside this package. The custom function should be able
#'    to use the node layout coordinates in `graph_attr(x, "layout")`,
#'    and render edges between nodes.
#' @param nodegroups `list` object as output by `edge_bundle_bipartite()`
#'    where each list element is a `character` vector of vertex node
#'    names present in `igraph::V(x)$name`. If no `"name"` vertex node
#'    attribute exists, then integer index values are used as names.
#'    Note that all vertex nodes must be represented in `nodegroup`
#'    in order for the corresponding edges to be plotted.
#' @param render_nodes,render_edges `logical` indicating whether to
#'    render vertex nodes, or edges, respectively. Sometimes it can
#'    be useful to call this function for other byproduct outputs,
#'    for example, `jam_plot_igraph(graph, add=FALSE, render_nodes=FALSE, render_edges=FALSE)`
#'    will create a new plot device with appropriate axis ranges,
#'    and can be used to render edge bundling results for example.
#' @param render_groups `logical` indicating whether to render groups
#'    when `mark.groups` is supplied. Groups are rendered with a
#'    shaded polygon and border.
#' @param render_nodelabels `logical` indicating whether to draw node
#'    labels, which is typically the last operation in the plot sequence.
#'    Note that node labels can be rendered without also rendering
#'    the nodes or edges.
#' @param plot_grid `logical` indicating whether to plot a background grid
#'    indicating units of 2% across the layout of the network graph. The
#'    units are calculated consistent with `nudge_igraph_nodes()`,
#'    `adjust_cnet_nodeset()` and other functions, scaled relative to the
#'    maximum x- or y-coordinate range of layout of the graph. Layout
#'    is obtained by `get_igraph_layout()` which by default uses
#'    supplied `layout`, or graph attribute `igraph::graph_attr(x, "layout")`.
#'    Note that by default, `jam_igraph()` represents the layout with
#'    aspect ratio = 1, so x-coordinates and y-coordiantes are represented
#'    with the same spacing per unit.
#'    This function calls `plot_layout_scale()` to render the grid lines.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param debug `logical` or `character` vector of attributes for
#'    which debug output will be plotted onscreen.
#'
jam_plot_igraph <- function
(x,
 ...,
 axes=FALSE,
 add=FALSE,
 xlim=NULL,
 ylim=NULL,
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
 pie_to_jampie=TRUE,
 use_shadowText=FALSE,
 vectorized_node_shapes=TRUE,
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
 render_groups=TRUE,
 render_nodelabels=TRUE,
 params=NULL,
 plot_grid=FALSE,
 verbose=FALSE,
 debug=NULL)
{
   graph <- x;
   if (!igraph::is_igraph(graph)) {
      stop("Not an igraph object")
   }

   # use mark.groups==FALSE as hard sign not to include mark.groups
   if (isFALSE(mark.groups)) {
      mark.groups <- NULL;
      mark.col <- NULL;
   }
   # validate edge_bundling input
   if (length(edge_bundling) == 0 || "default" %in% edge_bundling) {
      # default will try to detect an appropriate method
      # maybe some vcount ceiling should disable edge bundling?
      if (length(nodegroups) > 0) {
         edge_bundling <- "nodegroups";
      } else if (length(mark.groups) > 0) {
         edge_bundling <- "mark.groups";
      } else {
         edge_bundling <- "connections";
      }
      # edge_bundling <- "none";
   } else if (is.atomic(edge_bundling)) {
      edge_bundling <- head(intersect(edge_bundling,
         eval(formals(jam_plot_igraph)$edge_bundling)), 1);
      if (length(edge_bundling) == 0) {
         stop("edge_bundling method not recognized")
      }
   } else if (is.function(edge_bundling)) {
      edge_function <- edge_bundling;
      edge_bundling <- "function";
   }

   # create mark.col if needed
   if (length(mark.groups) > 0) {
      if (length(mark.col) == 0) {
         mark.col <- colorjam::rainbowJam(length(mark.groups),
            alpha=mark.alpha);
      }
      if (length(mark.border) == 0) {
         mark.border <- jamba::alpha2col(mark.col,
            alpha=1);
      }
   }

   ## avoid this technique, but all text() must occur here
   # if (use_shadowText) {
   #    text <- jamba::shadowText;
   # }

   if (length(params) == 0) {
      if (verbose) {
         jamba::printDebug("jam_plot_igraph(): ",
            "Parsing igraph params.");
      }
      params <- parse_igraph_plot_params(graph, list(...));
   }

   vertex.size <- 1/200 * params("vertex", "size");
   # jamba::printDebug("jam_plot_igraph(): ", "vertex.size: ", vertex.size);# debug
   label.family <- params("vertex", "label.family")
   label.font <- params("vertex", "label.font")
   label.fontsize <- params("vertex", "label.fontsize")
   if (verbose &&
         length(label.fontsize) > 0 &&
         any(!is.na(label.fontsize))) {
      jamba::printDebug("jam_plot_igraph(): ",
         "Applying fixed fontsize: ", "label.fontsize");
   }
   label.cex <- params("vertex", "label.cex")
   label.degree <- params("vertex", "label.degree")
   label.color <- params("vertex", "label.color")
   label.dist <- params("vertex", "label.dist")
   labels <- params("vertex", "label")
   if (length(debug) > 0 && any(c("labels","label.dist") %in% debug)) {
      jamba::printDebug("jam_plot_igraph(): ",
         "head(label.dist, 20):",
         head(label.dist, 20));
      jamba::printDebug("jam_plot_igraph(): ",
         "head(vertex.size, 20):",
         head(vertex.size, 20));
      jamba::printDebug("jam_plot_igraph(): ",
         "xlim:",
         xlim);
   }
   if (length(debug) > 0 && any(c("labels","label.degree") %in% debug)) {
      jamba::printDebug("jam_plot_igraph(): ",
         "head(label.degree, 20):",
         head(label.degree, 20));
   }

   valid_shapes <- igraph::shapes();
   shape <- params("vertex", "shape");
   ## Optionally convert shape "pie" to "jampie" for vectorized plotting
   if (pie_to_jampie && "jampie" %in% valid_shapes) {
      shape <- ifelse(shape %in% "pie",
         "jampie",
         shape);
   }
   # validate shapes exist in igraph framework
   if (!all(shape %in% valid_shapes)) {
      stop(paste0("Bad vertex shape(s): ",
         jamba::cPasteSU(setdiff(shape, valid_shapes), sep=", ")));
   }

   edge.color <- params("edge", "color")
   edge.width <- params("edge", "width")
   edge.lty <- params("edge", "lty")
   arrow.mode <- params("edge", "arrow.mode")
   edge.labels <- params("edge", "label")
   loop.angle <- params("edge", "loop.angle")
   edge.label.font <- params("edge", "label.font")
   edge.label.family <- params("edge", "label.family")
   edge.label.cex <- params("edge", "label.cex")
   edge.label.fontsize <- params("edge", "label.fontsize")
   if (verbose &&
         length(edge.label.fontsize) > 0 &&
         any(!is.na(edge.label.fontsize))) {
      jamba::printDebug("jam_plot_igraph(): ",
         "Applying fixed fontsize: ", "edge.label.fontsize");
   }
   edge.label.color <- params("edge", "label.color")
   elab.x <- params("edge", "label.x")
   elab.y <- params("edge", "label.y")
   arrow.size <- params("edge", "arrow.size")[1]
   arrow.width <- params("edge", "arrow.width")[1]
   curved <- params("edge", "curved")
   if (is.function(curved)) {
      curved <- curved(graph)
   }
   layout <- params("plot", "layout")

   margin <- params("plot", "margin")
   margin <- rep(margin,
      length.out=4);
   rescale <- params("plot", "rescale")
   asp <- params("plot", "asp")
   frame <- params("plot", "frame")
   main <- params("plot", "main")
   sub <- params("plot", "sub")
   xlab <- params("plot", "xlab")
   ylab <- params("plot", "ylab")

   ## color palette, is this chunk needed? palette is not called again
   palette <- params("plot", "palette")
   if (!is.null(palette)) {
      old_palette <- palette(palette)
      on.exit(palette(old_palette), add=TRUE)
   }

   arrow.mode <- get_igraph_arrow_mode(graph, arrow.mode)
   maxv <- max(vertex.size, na.rm=TRUE)
   medv <- median(vertex.size, na.rm=TRUE)

   ## Optional axis range scaling
   if (TRUE %in% rescale) {
      layout <- igraph::norm_coords(layout, -1, 1, -1, 1)
      xlim <- c(xlim[1] - margin[2] - maxv,
         xlim[2] + margin[4] + maxv);
      ylim <- c(ylim[1] - margin[1] - maxv,
         ylim[2] + margin[3] + maxv)
   }
   if (verbose) {
      jamba::printDebug("jam_plot_igraph(): ",
         "Parsed igraph params.");
   }

   ## initial plot device
   if (!add) {
      if (length(debug) > 0 && any(c("labels","label.dist") %in% debug)) {
         jamba::printDebug("jam_plot_igraph(): ",
            "xlim applied:",
            xlim);
         jamba::printDebug("jam_plot_igraph(): ",
            "ylim applied:",
            ylim);
         jamba::printDebug("jam_plot_igraph(): ",
            "asp applied:",
            asp);
      }
      # create blank plot space
      plot(0,
         0,
         type="n",
         xlab=xlab,
         ylab=ylab,
         xlim=xlim,
         ylim=ylim,
         axes=axes,
         frame=frame,
         asp=asp,
         main=main,
         sub=sub)
      # optionally display background grid
      if (length(plot_grid) > 0 && TRUE %in% plot_grid) {
         plot_layout_scale(layout=layout,
            ...)
      }
   }

   ## Optional mark groups
   if (render_groups && length(mark.groups) > 0) {
      if (!is.list(mark.groups) && is.numeric(mark.groups)) {
         mark.groups <- list(mark.groups)
      }
      if (length(mark.groups) > 0) {
         mark.shape <- rep(mark.shape,
            length.out=length(mark.groups))
         mark.border <- rep(mark.border,
            length.out=length(mark.groups))
         mark.col <- rep(mark.col,
            length.out=length(mark.groups))
         mark.lwd <- rep(mark.lwd,
            length.out=length(mark.groups));
         mark.lty <- rep(mark.lty,
            length.out=length(mark.groups));
         mark.cex <- rep(mark.cex,
            length.out=length(mark.groups));

         # mark.expand is scaled relative to the maximum x-/y-axis range
         max_xy_range <- max(c(
            abs(diff(range(xlim))),
            abs(diff(range(ylim)))));

         # check for mark.expand=NULL
         if (length(mark.expand) == 0) {
            ## use half the pre-adjusted vertex.size
            # medv * 200 * (1.8 / max_xy_range) / 2
            ## using 1.9 since the xy limits expand a little along the way
            # mark.expand <- (medv * 200 * (1.8 / max_xy_range) / 2);
            #
            ## Todo:
            ## 1. Consider adjustment based upon number of nodes,
            ##    fewer nodes should slightly increase the mark.expand,
            ##    Some type of exponent factor perhaps?
            ## 2. Consider calculating mark.expand per group?
            med_adj_size <- medv * 200 * (1.9 / max_xy_range);
            mark.expand <- (med_adj_size + 0) / 2 + 0.0;
            if (verbose) {
               jamba::printDebug("Auto mark.expand: ", mark.expand);
               jamba::printDebug("medv (postadjusted): ", medv);
               jamba::printDebug("median vertex.size:  ", med_adj_size);
            }
         }

         mark.expand <- rep(mark.expand,
            length.out=length(mark.groups));
         expand.by <- (mark.expand / 200) * max_xy_range;

         # scale relative to max x-/y-axis range
         mark.x.nudge <- rep(mark.x.nudge,
            length.out=length(mark.groups)) * max_xy_range;
         mark.y.nudge <- rep(mark.y.nudge,
            length.out=length(mark.groups)) * max_xy_range;

         # optional text labels
         mark.group.cluster.names <- NULL;
         if (length(names(mark.groups)) > 0) {
            if (!"membership" %in% names(mark.groups)) {
               mark.group.cluster.names <- names(mark.groups)
            } else if ("cluster_names" %in% names(mark.groups)) {
               mark.group.cluster.names <- mark.groups$cluster_names;
            }
         }
         for (g in seq_along(mark.groups)) {
            v <- igraph::V(graph)[mark.groups[[g]]]
            if (length(vertex.size) == 1) {
               vs <- vertex.size
            } else {
               vs <- rep(vertex.size,
                  length=igraph::vcount(graph))[v]
            }
            # use new method make_point_hull()
            hull_label <- mark.group.cluster.names[g];
            phxy <- make_point_hull(
               x=layout[v, , drop=FALSE],
               buffer=expand.by[g] * 2,
               do_plot=TRUE,
               hull_method="default",
               add=TRUE,
               xpd=TRUE,
               col=mark.col[g],
               border=mark.border[g],
               lwd=mark.lwd[g],
               lty=mark.lty[g],
               smooth=mark.smooth,
               shape=mark.shape[g],
               label.cex=mark.cex[g],
               label.x.nudge=mark.x.nudge[g],
               label.y.nudge=mark.y.nudge[g],
               label=hull_label);
         }
      }
   }

   ## Render edges
   if (render_edges) {
      if ("none" %in% edge_bundling) {
         ## process edges with default igraph methods
         if (verbose) {
            jamba::printDebug("jam_plot_igraph(): ",
               "Processing igraph edges");
         }
         el <- igraph::as_edgelist(graph, names=FALSE);
         loops.e <- which(el[, 1] == el[, 2])
         nonloops.e <- which(el[, 1] != el[, 2])
         loops.v <- el[, 1][loops.e]
         loop.labels <- edge.labels[loops.e]
         loop.labx <- if (is.null(elab.x)) {
            rep(NA, length(loops.e))
         } else {
            elab.x[loops.e]
         }
         loop.laby <- if (is.null(elab.y)) {
            rep(NA, length(loops.e))
         } else {
            elab.y[loops.e]
         }
         edge.labels <- edge.labels[nonloops.e]
         elab.x <- if (is.null(elab.x)) {
            NULL
         } else {
            elab.x[nonloops.e]
         }
         elab.y <- if (is.null(elab.y)) {
            NULL
         } else {
            elab.y[nonloops.e]
         }
         el <- el[nonloops.e, , drop=FALSE]
         edge.coords <- matrix(0,
            nrow=nrow(el),
            ncol=4)
         edge.coords[, 1] <- layout[, 1][el[, 1]]
         edge.coords[, 2] <- layout[, 2][el[, 1]]
         edge.coords[, 3] <- layout[, 1][el[, 2]]
         edge.coords[, 4] <- layout[, 2][el[, 2]]

         ## clip edge coordinates using node shape functions
         if (length(unique(shape)) == 1) {
            ec <- igraph::shapes(shape[1])$clip(
               coords=edge.coords,
               el=el,
               params=params,
               end="both")
         } else {
            shape <- rep(shape, length=igraph::vcount(graph))
            ec <- edge.coords;
            ec[, 1:2] <- t(sapply(seq_len(nrow(el)), function(x) {
               igraph::shapes(shape[el[x, 1]])$clip(
                  edge.coords[x, , drop=FALSE],
                  el[x, , drop=FALSE],
                  params=params,
                  end="from")
            }))
            ec[, 3:4] <- t(sapply(seq_len(nrow(el)), function(x) {
               igraph::shapes(shape[el[x, 2]])$clip(
                  edge.coords[x, , drop=FALSE],
                  el[x, , drop=FALSE],
                  params=params,
                  end="to")
            }))
         }
         x0 <- ec[, 1]
         y0 <- ec[, 2]
         x1 <- ec[, 3]
         y1 <- ec[, 4]

         ## Plot edges that are self-loops
         if (length(loops.e) > 0) {
            if (verbose) {
               jamba::printDebug("jam_plot_igraph(): ",
                  "Processing igraph edge self-loops");
            }
            ec <- edge.color;
            if (length(ec) > 1) {
               ec <- ec[loops.e]
            }
            point.on.cubic.bezier <- function(cp, t) {
               c <- 3 * (cp[2, ] - cp[1, ])
               b <- 3 * (cp[3, ] - cp[2, ]) - c
               a <- cp[4, ] - cp[1, ] - c - b
               t2 <- t * t
               t3 <- t * t * t
               a * t3 + b * t2 + c * t + cp[1, ]
            }
            compute.bezier <- function(cp, points) {
               dt <- seq(0, 1, by=1/(points - 1))
               sapply(dt, function(t) {
                  point.on.cubic.bezier(cp, t)
               })
            }
            plot.bezier <- function
            (cp,
               points,
               color,
               width,
               arr,
               lty,
               arrow.size,
               arr.w)
            {
               p <- compute.bezier(cp, points)
               polygon(p[1, ],
                  p[2, ],
                  border=color,
                  lwd=width,
                  lty=lty)
               if (arr == 1 || arr == 3) {
                  jam_igraph_arrows(p[1, ncol(p) - 1],
                     p[2, ncol(p) - 1],
                     p[1, ncol(p)],
                     p[2, ncol(p)],
                     sh.col=color,
                     sh.lwd=width,
                     sh.lty=lty,
                     h.col=color,
                     h.lwd=width,
                     h.lty=lty,
                     size=arrow.size,
                     open=FALSE,
                     code=2,
                     width=arr.w)
               }
               if (arr == 2 || arr == 3) {
                  jam_igraph_arrows(p[1, 2],
                     p[2, 2],
                     p[1, 1], p[2, 1],
                     sh.col=color,
                     h.col=color,
                     size=arrow.size,
                     sh.lwd=width,
                     h.lwd=width,
                     open=FALSE,
                     code=2,
                     width=arr.w)
               }
            }
            loop <- function
            (x0,
               y0,
               cx=x0,
               cy=y0,
               color,
               angle=0,
               label=NA,
               width=1,
               arr=2,
               lty=1,
               arrow.size=arrow.size,
               arr.w=arr.w,
               lab.x,
               lab.y)
            {
               rad <- angle
               center <- c(cx, cy)
               cp <- matrix(
                  c(x0,
                     y0,
                     x0 + 0.4,
                     y0 + 0.2,
                     x0 + 0.4,
                     y0 - 0.2,
                     x0,
                     y0),
                  ncol=2,
                  byrow=TRUE)
               phi <- atan2(
                  cp[, 2] - center[2],
                  cp[, 1] - center[1])
               r <- sqrt(
                  (cp[, 1] - center[1])^2 +
                     (cp[, 2] - center[2])^2)
               phi <- phi + rad
               cp[, 1] <- cx + r * cos(phi)
               cp[, 2] <- cy + r * sin(phi)
               plot.bezier(cp,
                  50,
                  color,
                  width,
                  arr=arr,
                  lty=lty,
                  arrow.size=arrow.size,
                  arr.w=arr.w)
               if (is.language(label) || !is.na(label)) {
                  lx <- x0 + 0.3
                  ly <- y0
                  phi <- atan2(ly - center[2],
                     lx - center[1])
                  r <- sqrt((lx - center[1])^2 + (ly - center[2])^2)
                  phi <- phi + rad
                  lx <- cx + r * cos(phi)
                  ly <- cy + r * sin(phi)
                  if (!is.na(lab.x)) {
                     lx <- lab.x
                  }
                  if (!is.na(lab.y)) {
                     ly <- lab.y
                  }
                  text_subsets <- rep(edge.label.family,
                     length.out=length(lx));
                  if (length(edge.label.fontsize) > 0) {
                     new.edge.label.cex <- edge.label.fontsize /
                        (par("ps") * par("cex"));
                     edge.label.cex <- ifelse(
                        is.na(new.edge.label.cex),
                        edge.label.cex,
                        new.edge.label.cex)
                  }
                  for (k in split(seq_along(text_subsets), text_subsets)) {
                     if (use_shadowText) {
                        jamba::shadowText(lx[k],
                           ly[k],
                           label[k],
                           col=edge.label.color[k],
                           font=edge.label.font[k],
                           family=head(edge.label.family[k], 1),
                           cex=edge.label.cex[k])
                     } else {
                        text(lx[k],
                           ly[k],
                           label[k],
                           col=edge.label.color[k],
                           font=edge.label.font[k],
                           family=head(edge.label.family[k], 1),
                           cex=edge.label.cex[k])
                     }
                  }
               }
            }
            ec <- edge.color
            if (length(ec) > 1) {
               ec <- ec[loops.e]
            }
            vs <- vertex.size
            if (length(vertex.size) > 1) {
               vs <- vs[loops.v]
            }
            ew <- edge.width
            if (length(edge.width) > 1) {
               ew <- ew[loops.e]
            }
            la <- loop.angle
            if (length(loop.angle) > 1) {
               la <- la[loops.e]
            }
            lty <- edge.lty
            if (length(edge.lty) > 1) {
               lty <- lty[loops.e]
            }
            arr <- arrow.mode
            if (length(arrow.mode) > 1) {
               arr <- arrow.mode[loops.e]
            }
            asize <- arrow.size
            if (length(arrow.size) > 1) {
               asize <- arrow.size[loops.e]
            }
            xx0 <- layout[loops.v, 1] + cos(la) * vs
            yy0 <- layout[loops.v, 2] - sin(la) * vs
            mapply(loop,
               xx0,
               yy0,
               color=ec,
               angle=-la,
               label=loop.labels,
               lty=lty,
               width=ew,
               arr=arr,
               arrow.size=asize,
               arr.w=arrow.width,
               lab.x=loop.labx,
               lab.y=loop.laby)
         }

         ## Plot edges that are not self-loops
         if (length(x0) != 0) {
            if (verbose) {
               jamba::printDebug("jam_plot_igraph(): ",
                  "Processing igraph edge non-loops");
            }
            if (length(edge.color) > 1) {
               edge.color <- edge.color[nonloops.e]
            }
            if (length(edge.width) > 1) {
               edge.width <- edge.width[nonloops.e]
            }
            if (length(edge.lty) > 1) {
               edge.lty <- edge.lty[nonloops.e]
            }
            if (length(arrow.mode) > 1) {
               arrow.mode <- arrow.mode[nonloops.e]
            }
            if (length(arrow.size) > 1) {
               arrow.size <- arrow.size[nonloops.e]
            }
            if (length(curved) > 1) {
               curved <- curved[nonloops.e]
            }
            if (length(unique(arrow.mode)) == 1) {
               lc <- jam_igraph_arrows(x0,
                  y0,
                  x1,
                  y1,
                  h.col=edge.color,
                  sh.col=edge.color,
                  sh.lwd=edge.width,
                  h.lwd=1,
                  open=FALSE,
                  code=head(arrow.mode, 1),
                  sh.lty=edge.lty,
                  h.lty=1,
                  size=arrow.size,
                  width=arrow.width,
                  curved=curved)
               lc.x <- lc$lab.x
               lc.y <- lc$lab.y
            } else {
               curved <- rep(curved,
                  length=igraph::ecount(graph))[nonloops.e]
               lc.x <- lc.y <- numeric(length(curved))
               for (code in 0:3) {
                  valid <- (arrow.mode == code)
                  if (!any(valid)) {
                     next;
                  }
                  ec <- edge.color
                  if (length(ec) > 1) {
                     ec <- ec[valid]
                  }
                  ew <- edge.width
                  if (length(ew) > 1) {
                     ew <- ew[valid]
                  }
                  el <- edge.lty
                  if (length(el) > 1) {
                     el <- el[valid]
                  }
                  # lc <- igraph:::igraph.Arrows(
                  lc <- jam_igraph_arrows(
                     x0[valid],
                     y0[valid],
                     x1[valid],
                     y1[valid],
                     code=code,
                     sh.col=ec,
                     h.col=ec,
                     sh.lwd=ew,
                     h.lwd=1,
                     h.lty=1,
                     sh.lty=el,
                     open=FALSE,
                     size=arrow.size,
                     width=arrow.width,
                     curved=curved[valid])
                  lc.x[valid] <- lc$lab.x
                  lc.y[valid] <- lc$lab.y
               }
            }
            if (!is.null(elab.x)) {
               lc.x <- ifelse(is.na(elab.x),
                  lc.x,
                  elab.x)
            }
            if (!is.null(elab.y)) {
               lc.y <- ifelse(is.na(elab.y),
                  lc.y,
                  elab.y)
            }
            ## edge labels
            if (length(lc.x) > 0 &&
                  length(lc.y) > 0 &&
                  length(edge.labels) > 0) {
               text_subsets <- rep(edge.label.family,
                  length.out=length(lc.x));
               if (length(edge.label.fontsize) > 0) {
                  new.edge.label.cex <- edge.label.fontsize /
                     (par("ps") * par("cex"));
                  edge.label.cex <- ifelse(
                     is.na(new.edge.label.cex),
                     edge.label.cex,
                     new.edge.label.cex)
               }
               for (k in split(seq_along(text_subsets), text_subsets)) {
                  if (use_shadowText) {
                     jamba::shadowText(lc.x[k],
                        lc.y[k],
                        labels=edge.labels[k],
                        col=edge.label.color[k],
                        family=head(edge.label.family[k], 1),
                        font=edge.label.font[k],
                        cex=edge.label.cex[k])
                  } else {
                     text(lc.x[k],
                        lc.y[k],
                        labels=edge.labels[k],
                        col=edge.label.color[k],
                        family=head(edge.label.family[k], 1),
                        font=edge.label.font[k],
                        cex=edge.label.cex[k])
                  }
               }
            }
         }
         ## remove x0, x1, y0, y1 (why?)
         rm(x0, y0, x1, y1)
      } else if ("connections" %in% edge_bundling) {
         if (verbose) {
            jamba::printDebug("jam_plot_igraph(): ",
               "Edge bundling with edge_bundle_bipartite()");
         }
         # bipartite nodegroups
         igraph::graph_attr(graph, "layout") <- layout;
         nodegroups <- edge_bundle_bipartite(graph,
            verbose=FALSE,
            ...);
         if (verbose) {
            jamba::printDebug("jam_plot_igraph(): ",
               "Rendering edge bundling with edge_bundle_nodegroups()");
         }
         igraph::E(graph)$color <- edge.color;
         igraph::E(graph)$width <- edge.width;
         edge_bundle_nodegroups(graph,
            nodegroups=nodegroups,
            shape=shape,
            params=params,
            bundle_self=bundle_self,
            ...);

      } else if (any(c("mark.groups", "nodegroups") %in% edge_bundling)) {
         if (verbose) {
            jamba::printDebug("jam_plot_igraph(): ",
               "Rendering edge bundling with edge_bundle_nodegroups()");
         }
         igraph::graph_attr(graph, "layout") <- layout;
         igraph::E(graph)$color <- edge.color;
         igraph::E(graph)$width <- edge.width;
         if ("nodegroups" %in% edge_bundling) {
            edge_bundle_nodegroups(graph,
               nodegroups=nodegroups,
               params=params,
               bundle_self=bundle_self,
               ...);
         } else {
            edge_bundle_nodegroups(graph,
               nodegroups=mark.groups,
               params=params,
               bundle_self=bundle_self,
               ...);
         }
      } else if ("function" %in% edge_bundling) {
         if (verbose) {
            jamba::printDebug("jam_plot_igraph(): ",
               "Edge bundling with custom function.");
         }
         edge_function(x,
            param,
            verbose=verbose>1,
            ...)
      }
   }
   ## End edge plotting

   ## plot vertex/node shapes
   if (render_nodes) {
      if (length(unique(shape)) == 1) {
         # igraph:::.igraph.shapes[[shape[1]]]$plot(
         igraph::shapes(shape[1])$plot(
            layout,
            params=params)
      } else {
         ## Loop each unique shape here instead of individual nodes
         unique_shapes <- unique(shape);
         nodes_by_shape <- split(seq_len(igraph::vcount(graph)), shape);
         if (vectorized_node_shapes) {
            if (verbose) {
               jamba::printDebug("jam_plot_igraph(): ",
                  "Using vectorized node shape plotting.");
            }
            sapply(nodes_by_shape, function(x){
               shape1 <- shape[x[1]];
               # igraph:::.igraph.shapes[[shape1]]$plot(
               igraph::shapes(shape1)$plot(
                  layout[x, , drop=FALSE],
                  v=x,
                  params=params)
            })
         } else {
            if (verbose) {
               jamba::printDebug("jam_plot_igraph(): ",
                  "Using linear non-vectorized node shape plotting.");
            }
            sapply(seq_len(igraph::vcount(graph)), function(x) {
               # igraph:::.igraph.shapes[[shape[x]]]$plot(
               igraph::shapes(shape[x])$plot(
                  layout[x, , drop=FALSE],
                  v=x,
                  params=params)
            })
         }
      }
   }

   # forxce xpd=TRUE?

   ## Define node label positions
   if (render_nodelabels) {
      opar <- par(xpd=TRUE);
      x <- layout[, 1] +
         (label.dist *
               cos(-label.degree) *
               (vertex.size + 6 * 8 * log10(2)) / 200 * 200/180);
      y <- layout[, 2] +
         (label.dist *
               sin(-label.degree) *
               (vertex.size + 6 * 8 * log10(2)) / 200 * 200/180);
      if ("labels" %in% debug) {
         jamba::printDebug("jam_plot_igraph(): ",
            "labels x,y coords:");
         print(head(
            data.frame(
               labels=labels,
               x=layout[, 1],
               y=layout[, 2],
               label_x=x,
               label_y=y,
               vertex.size=vertex.size,
               label.degree=label.degree,
               label.dist=label.dist), Inf));
      }
      if (length(label.family) == 1) {
         if ("labels" %in% debug) {
            jamba::printDebug("jam_plot_igraph(): ",
               "labels all have one font family, one call to text()");
         }
         # Note the actual font size in points is derived by
         # multiplying par("ps") * par("cex") * vertex.label.cex.
         # So it is possible to define an absolute font size in points
         # but it would need to be an alternative to using the ps*cex method.
         if (length(label.fontsize) > 0) {
            # Logic to apply fontsize:
            # - NA values use original label.cex
            # - any non-NA values override label.cex
            vertex.label.cex <- label.fontsize / (par("ps") * par("cex"))
            label.cex <- ifelse(is.na(vertex.label.cex),
               label.cex,
               vertex.label.cex)
         }
         if (use_shadowText) {
            jamba::shadowText(x,
               y,
               labels=labels,
               col=label.color,
               family=label.family,
               font=label.font,
               cex=label.cex)
         } else {
            text(x,
               y,
               labels=labels,
               col=label.color,
               family=label.family,
               font=label.font,
               cex=label.cex)
         }
      } else {
         if ("labels" %in% debug) {
            jamba::printDebug("jam_plot_igraph(): ",
               "labels all have different font families, multiple calls to text()");
         }
         if1 <- function(vect, idx) if (length(vect) == 1) {
            vect
         } else {
            vect[idx]
         }
         ## Also loop through each family here instead of each node
         unique_family <- unique(label.family);
         label.family <- rep(label.family,
            length.out=igraph::vcount(graph));
         nodes_by_family <- split(seq_len(igraph::vcount(graph)), label.family);
         if (length(label.fontsize) > 0) {
            # Logic to apply fontsize:
            # - NA values use original label.cex
            # - any non-NA values override label.cex
            vertex.label.cex <- label.fontsize / (par("ps") * par("cex"))
            label.cex <- ifelse(is.na(vertex.label.cex),
               label.cex,
               vertex.label.cex)
         }
         # iterate each font family
         for (v in nodes_by_family) {
            family1 <- label.family[head(v, 1)];
            if (use_shadowText) {
               jamba::shadowText(
                  x[v],
                  y[v],
                  labels=if1(labels, v),
                  col=if1(label.color, v),
                  family=family1,
                  font=if1(label.font, v),
                  cex=if1(label.cex, v))
            } else {
               text(
                  x[v],
                  y[v],
                  labels=if1(labels, v),
                  col=if1(label.color, v),
                  family=family1,
                  font=if1(label.font, v),
                  cex=if1(label.cex, v))
            }
         }
      }
      par(opar)
   }
   ## avoid this process
   # if (use_shadowText) {
   #    rm(text)
   # }
   ## remnant from igraph:::plot.igraph() not needed here
   # rm(x, y)
   invisible(NULL)
}
