
# jamenrich-color-nodes-edges.R

#' Colorize igraph edges by nodes
#'
#' Colorize igraph edges using node colors
#'
#' This function colorizes edges by blending colors for the
#' nodes involved, by calling `colorjam::blend_colors()`.
#'
#' The color for each node depends upon the node shape, so
#' the color or colors used to render each node shape will
#' be used for the edge. For example:
#'
#' * `shape="pie"` uses the average color from `V(g)$pie.color`
#' * `shape="coloredrectangle"` uses the avereage color from
#' `V(g)$coloredrect.color`
#' * everything else uses `V(g)$color`
#'
#' @family jam igraph utilities
#'
#' @return `igraph` object with edge color attribute updated to
#'    represent the result of blending node colors, seen by
#'    `igraph::edge_attr(g)$color`.
#'
#' @param g `igraph` object that contains vertex node attribute `"color"`
#'    as seen with `igraph::vertex_attr(g, "color")`.
#' @param edge_alpha `numeric` or `NULL`, where numeric value sets
#'    the edge alpha transparency, where `edge_alpha=0` is completely
#'    transparent, `edge_alpha=0.5` is 50% transparent, and `edge_alpha=1`
#'    is completely not transparent, and is opaque. When `edge_alpha=NULL`
#'    the alpha values are supplied by `colorjam::blend_colors()`
#'    which blends the two values.
#' @param ... additional arguments are passed to `colorjam::blend_colors()`.
#'
#' @export
color_edges_by_nodes <- function
(g,
 edge_alpha=NULL,
 Crange=c(0, 100),
 Lrange=c(0, 65),
 ...)
{
   # confirm color is present
   if (!"igraph" %in% class(g)) {
      stop("Input must be class 'igraph'");
   }
   if (!"color" %in% igraph::vertex_attr_names(g)) {
      stop("Node attributes must contain 'color' in V(g)$color");
   }
   # get edge data.frame
   if (!"name" %in% igraph::vertex_attr_names(g)) {
      vname <- as.character(seq_len(igraph::vcount(g)));
   } else {
      vname <- igraph::V(g)$name;
   }

   edge_df <- data.frame(check.names=FALSE,
      stringsAsFactors=FALSE,
      igraph::as_edgelist(g));

   # fix bug when shape is not defined
   if (length(igraph::V(g)$shape) == 0) {
      igraph::V(g)$shape <- "circle";
   }
   node_colors_list <- ifelse(
      igraph::V(g)$shape %in% "pie",
      igraph::V(g)$pie.color,
      ifelse(
         igraph::V(g)$shape %in% "coloredrectangle",
         igraph::V(g)$coloredrect.color,
         igraph::V(g)$color))
   if (!is.list(node_colors_list)) {
      node_colors <- node_colors_list;
      node_alphas <- jamba::col2alpha(
         jamba::rmNA(naValue=0.01, node_colors));
   } else {
      node_colors <- colorjam::blend_colors(node_colors_list);
      node_alphas <- sapply(node_colors_list, function(i){
         mean(jamba::col2alpha(jamba::rmNA(naValue=0.01, i)))
      })
   }
   if (any(node_alphas < 1)) {
      newalpha <- (node_alphas < 1);
      node_colors[newalpha] <- jamba::alpha2col(node_colors[newalpha],
         alpha=node_alphas[newalpha]);
   }

   edge_df$color1 <- node_colors[match(edge_df[,1], vname)];
   edge_df$color2 <- node_colors[match(edge_df[,2], vname)];

   color1_2 <- split(
      c(edge_df$color1, edge_df$color2),
      rep(seq_len(nrow(edge_df)), 2));
   blended1_2 <- colorjam::blend_colors(color1_2,
      ...);

   # apply Crange and Lrange to restrict output chroma and luminance
   if (max(Crange) < 200 || max(Lrange) < 100) {
      blended1_2 <- jamba::applyCLrange(blended1_2,
         Crange=Crange,
         Lrange=Lrange,
         CLmethod="floor");
   }

   # apply edge alpha transparency
   if (length(edge_alpha) > 0) {
      igraph::E(g)$color <- jamba::alpha2col(blended1_2,
         alpha=rep(edge_alpha,
            length.out=length(blended1_2)));
   } else {
      igraph::E(g)$color <- blended1_2;
   }
   return(g);
}

#' Color edges by nodegroups
#'
#' Color edges by nodegroups
#'
#' @family jam igraph utilities
#'
#' @return `igraph` object with edge color attribute updated to
#'    represent the result of blending node colors, seen by
#'    `igraph::edge_attr(g)$color`.
#'
#' @param g `igraph` object that contains vertex node attribute `"color"`
#'    as seen with `igraph::vertex_attr(g, "color")`.
#' @param nodegroups `list` or `communities` object that references
#'    nodes in `g` and assigns one or more to nodegroups.
#' @param edge_alpha `numeric` or `NULL`, where numeric value sets
#'    the edge alpha transparency, where `edge_alpha=0` is completely
#'    transparent, `edge_alpha=0.5` is 50% transparent, and `edge_alpha=1`
#'    is completely not transparent, and is opaque. When `edge_alpha=NULL`
#'    the alpha values are supplied by `colorjam::blend_colors()`
#'    which blends the two values.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to `colorjam::blend_colors()`.
#'
#' @export
color_edges_by_nodegroups <- function
(g,
 nodegroups,
 nodegroup_colors=NULL,
 edge_alpha=NULL,
 Crange=c(60, 100),
 Lrange=c(45, 85),
 verbose=TRUE,
 ...)
{
   #
   # confirm color is present
   if (!"igraph" %in% class(g)) {
      stop("Input must be class 'igraph'");
   }
   if (!"color" %in% igraph::vertex_attr_names(g)) {
      stop("Node attributes must contain 'color' in V(g)$color");
   }
   # get edge data.frame
   if (!"name" %in% igraph::vertex_attr_names(g)) {
      vname <- as.character(seq_len(igraph::vcount(g)));
   } else {
      vname <- igraph::V(g)$name;
   }

   if ("communities" %in% class(nodegroups)) {
      if (verbose) {
         jamba::printDebug("color_edges_by_nodegroups(): ",
            "Converted communities to nodegroups format.")
      }
      nodegroups <- communities2nodegroups(nodegroups)
   }
   if (length(names(nodegroups)) == 0) {
      names(nodegroups) <- seq_along(nodegroups)
   }
   nodegroup_df <- data.frame(name=unlist(unname(nodegroups)),
      nodegroup=rep(names(nodegroups), lengths(nodegroups)))

   if (length(nodegroup_colors) == 0) {
      nodegroup_colors <- colorjam::rainbowJam(
         n=length(nodegroups),
         Crange=Crange,
         Lrange=Lrange,
         ...)
      names(nodegroup_colors) <- names(nodegroups);
      if (verbose) {
         jamba::printDebug("color_edges_by_nodegroups(): ",
            "Defined nodegroup_colors:")
         jamba::printDebugI(jamba::nameVector(
            nodegroup_colors,
            paste0(names(nodegroup_colors), ":", nodegroup_colors)))
      }
   }

   edge_df <- data.frame(check.names=FALSE,
      stringsAsFactors=FALSE,
      igraph::as_edgelist(g));
   ematch1 <- match(edge_df[,1], nodegroup_df$name);
   ematch2 <- match(edge_df[,2], nodegroup_df$name);
   edge_df$nodegroup1 <- nodegroup_df$nodegroup[ematch1]
   edge_df$nodegroup2 <- nodegroup_df$nodegroup[ematch2]
   edge_df$nodegroup1_color <- nodegroup_colors[edge_df$nodegroup1];
   edge_df$nodegroup2_color <- nodegroup_colors[edge_df$nodegroup2];
   nodecolor_df <- unique(
      edge_df[,c("nodegroup1_color", "nodegroup2_color"), drop=FALSE])

   edgecolors <- sapply(seq_len(nrow(nodecolor_df)), function(i){
      colorjam::blend_colors(
         x=unname(unlist(nodecolor_df[i,])),
         ...)
   })

   # optionally apply edge_alpha if defined
   if (length(edge_alpha) > 0 && is.numeric(edge_alpha)) {
      edgecolors <- jamba::alpha2col(edgecolors,
         alpha=edge_alpha);
      if (verbose) {
         jamba::printDebug("color_edges_by_nodegroups(): ",
            "Applied edge_alpha:",
            edge_alpha)
      }
   }

   nodecolor_match <- match(
      jamba::pasteByRow(
         edge_df[,c("nodegroup1_color", "nodegroup2_color"), drop=FALSE]),
      jamba::pasteByRow(nodecolor_df))
   edge_df$edgecolor <- edgecolors[nodecolor_match];
   igraph::edge_attr(g, "color") <- edge_df$edgecolor;
   return(g);
}


#' Color edges by nodegroups
#'
#' Color edges by nodegroups
#'
#' @family jam igraph utilities
#'
#' @return `igraph` object with edge color attribute updated to
#'    represent the result of blending node colors, seen by
#'    `igraph::edge_attr(g)$color`.
#'
#' @param g `igraph` object that contains vertex node attribute `"color"`
#'    as seen with `igraph::vertex_attr(g, "color")`.
#' @param nodegroups `list` or `communities` object that references
#'    nodes in `g` and assigns one or more to nodegroups.
#' @param edge_alpha `numeric` or `NULL`, where numeric value sets
#'    the edge alpha transparency, where `edge_alpha=0` is completely
#'    transparent, `edge_alpha=0.5` is 50% transparent, and `edge_alpha=1`
#'    is completely not transparent, and is opaque. When `edge_alpha=NULL`
#'    the alpha values are supplied by `colorjam::blend_colors()`
#'    which blends the two values.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to `colorjam::blend_colors()`.
#'
#' @export
color_nodes_by_nodegroups <- function
(g,
 nodegroups,
 nodegroup_colors=NULL,
 node_alpha=NULL,
 Crange=c(60, 100),
 Lrange=c(45, 85),
 color_attributes=c("color"),
 verbose=TRUE,
 ...)
{
   #
   # confirm color is present
   if (!"igraph" %in% class(g)) {
      stop("Input must be class 'igraph'");
   }
   if (!"color" %in% igraph::vertex_attr_names(g)) {
      stop("Node attributes must contain 'color' in V(g)$color");
   }
   # get edge data.frame
   if (!"name" %in% igraph::vertex_attr_names(g)) {
      vname <- as.character(seq_len(igraph::vcount(g)));
   } else {
      vname <- igraph::V(g)$name;
   }

   if ("communities" %in% class(nodegroups)) {
      if (verbose) {
         jamba::printDebug("color_edges_by_nodegroups(): ",
            "Converted communities to nodegroups format.")
      }
      nodegroups <- communities2nodegroups(nodegroups)
   }
   if (length(names(nodegroups)) == 0) {
      names(nodegroups) <- seq_along(nodegroups)
   }
   nodegroup_df <- data.frame(name=unlist(unname(nodegroups)),
      nodegroup=rep(names(nodegroups), lengths(nodegroups)))

   if (length(nodegroup_colors) == 0) {
      nodegroup_colors <- colorjam::rainbowJam(
         n=length(nodegroups),
         Crange=Crange,
         Lrange=Lrange,
         ...)
      names(nodegroup_colors) <- names(nodegroups);
      if (verbose) {
         jamba::printDebug("color_edges_by_nodegroups(): ",
            "Defined nodegroup_colors:")
         jamba::printDebugI(jamba::nameVector(
            nodegroup_colors,
            paste0(names(nodegroup_colors), ":", nodegroup_colors)))
      }
   }
   nodegroup_df$new_color <- nodegroup_colors[nodegroup_df$nodegroup];

   for (color_attribute in color_attributes) {
      igraph::vertex_attr(g, name=color_attribute) <- nodegroup_df$new_color;
   }

   return(g);
}

