

# highlight_igraph_node <- function
# (g,
#  node=NULL,
#  highlight.border="orange",
#  ...)
# {
#    #
# }

#' Highlight edges connected to a node or nodes
#'
#' Highlight edges connected to a node or nodes
#'
#' This function highlights edges connected to one of more nodes
#' by applying a color, and adjusting the edge width.
#'
#' Note that when edge attributes 'width' and 'color' are not yet
#' defined in the `igraph` object, they will be defined and stored
#' into the object by using `default_igraph_values()`.
#' For example: `default_igraph_values()$edge$width` and
#' `default_igraph_values()$edge$color`.
#'
#' @returns `igraph` object after adjusting edge attributes.
#'
#' @family jam igraph functions
#'
#' @param g `igraph` object
#' @param node `character` node name, or `integer` node index.
#' @param highlight_color `character` color used for highlighting,
#'    default 'royalblue'.
#' @param highlight_width `numeric` default NULL, set a fixed edge width.
#' @param highlight_cex `numeric` default 2, multiplies the edge width
#'    by this expansion factor.
#' @param nonhighlight_alpha `numeric` default NULL, applies an optional
#'    alpha transparency to non-highlighted edges. Values must be
#'    between 0 (transparent) and 1 (opaque), where a value 0.5 is
#'    recommended when used.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' cnet1 <- make_cnet_test(seed=1);
#' # highlight edges by node
#' cnet1h <- highlight_edges_by_node(cnet1, "CA", nonhighlight_alpha=0.5)
#' # plot the original
#' jam_igraph(cnet1, node_factor=2, label_factor=2)
#' # plot the highlighted variant
#' jam_igraph(cnet1h, node_factor=2, label_factor=2, label_dist_factor=3)
#'
#' # highlight edges using nodeset
#' nodesets <- get_cnet_nodeset(cnet1)
#' cnet1hns <- highlight_edges_by_node(cnet1,
#'    nodesets[[grep("CA", nodesets)]],
#'    nonhighlight_alpha=0.5)
#' jam_igraph(cnet1hns, node_factor=2, label_factor=2, label_dist_factor=3,
#'    mark.groups=unname(nodesets[grep("CA", nodesets)]))
#' 
#' @export
highlight_edges_by_node <- function
(g,
 node=NULL,
 highlight_color="royalblue3",
 highlight_width=NULL,
 highlight_cex=2,
 nonhighlight_alpha=NULL,
 ...)
{
   #
   # get neighbors
   if (length(node) == 0) {
      return(g);
   }
   use_edge_ids_list <- lapply(igraph::incident_edges(g, node),
      as.numeric);
   use_edge_ids <- unique(unlist(use_edge_ids_list))

   # add "width" if not present
   if (!"width" %in% igraph::edge_attr_names(g)) {
      igraph::E(g)$width <- default_igraph_values()$edge$width;
   }
   # when provided, use highlight_width
   if (length(highlight_width) > 0) {
      igraph::E(g)[use_edge_ids]$width <- rep(highlight_width,
         length.out=length(use_edge_ids))
   }
   # when provided, apply highlight_cex
   if (length(highlight_cex) > 0) {
      igraph::E(g)[use_edge_ids]$width <- (rep(highlight_cex,
         length.out=length(use_edge_ids)) *
            igraph::E(g)[use_edge_ids]$width);
   }

   # add "color" if not present
   if (!"color" %in% igraph::edge_attr_names(g)) {
      igraph::E(g)$color <- default_igraph_values()$edge$color;
   }
   # when provided, use highlight_width
   if (length(highlight_color) > 0) {
      igraph::E(g)[use_edge_ids]$color <- rep(highlight_color,
         length.out=length(use_edge_ids))
   }

   # optional alpha for non-highlighted edges
   if (length(nonhighlight_alpha) > 0 && is.numeric(nonhighlight_alpha)) {
      non_edge_ids <- setdiff(seq_len(igraph::ecount(g)), use_edge_ids);
      igraph::E(g)[non_edge_ids]$color <- jamba::alpha2col(
         x=igraph::E(g)[non_edge_ids]$color,
         alpha=rep(nonhighlight_alpha, length.out=length(non_edge_ids)));
   }

   return(g)
}
