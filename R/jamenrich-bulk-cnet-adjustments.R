

#' Bulk Cnet plot adjustments
#'
#' Bulk Cnet plot adjustments
#'
#' This function applies a series of adjustments to nodesets and nodes in
#' a Cnet `igraph` object. Nodesets are adjusted first, then nodes
#' are adjusted afterward. This order is important because the x,y
#' coordinates may be altered by nodeset first, before the node adjustment
#' is applied.
#'
#' All x,y coordinate adjustments are performed using scaled units,
#' using the max x- or y-range for the plot layout.
#' All Cnet `igraph` plots are assumed to use `asp=1` to enforce
#' aspect ratio 1, in order to maintain the intended visual distances
#' used by the layout algorithms.
#'
#' Therefore it is possible to expand a plot layout such that the
#' overall scale is wider than before, thereby affecting the scaled unit
#' space.
#'
#' In any event, this function is intended to replicated the effect done
#' by `launch_shinycat()`.
#'
#' @param g `igraph` object, assumed to be Cnet plot with vertex attribute
#'    'nodeType' with nodes assigned as 'Set' or 'Gene'.
#' @param node_adj `data.frame` with columns:
#'    * 'node' - `character` name of each node to adjust
#'    * 'x', 'y' - `numeric` values to add to existing node coordinates.
#'    * 'label.degree' - `numeric` value in degrees (not radians) to add
#'    to the existing node label angle 'label.dist'.
#'    It is stored in the igraph object using radians.
#'    * 'label.dist' - `numeric` distance added to 'label.dist' for each node.
#' @param nodeset_adj `data.frame` with columns:
#'    * 'nodeset' - `character` name of the nodeset
#'    * 'x', 'y' - `numeric` values to add to existing node coordinates for
#'    all nodes in each nodeset.
#'    * 'rotate_degrees' - `numeric` rotation to apply to all nodes in
#'    a given nodeset, units are degrees and not radians.
#'    * 'percent_spacing' - `numeric` spacing to apply to each nodeset.
#'    Note the argument `apply_negative` which controls whether node spacing
#'    would be contracted/shrunk by this process.
#' @param nodesets `list` named by nodeset, with `character` vectors
#'    containing node names, as output from `get_cnet_nodesets()`.
#'    Any node grouping will work.
#' @param apply_negative `logical` default TRUE, whether the `percent_spacing`
#'    is able to force node spacing to be reduced.
#' @param ... additional arguments are ignored.
#'
#' @returns `igraph` object after applying adjustments.
#'
#' @family jam cnet utilities
#'
#' @param g `igraph` object expected to be a Cnet object with vertex attribute
#'    'nodeType' containing values 'Gene' and 'Set'.
#' @param nodeset_adj `data.frame` recognizing five colnames:
#'    * nodeset: character name of each nodeset
#'    * x,y: numeric x,y coordinate adjustments to add to existing coordinates
#'    * rotate_degrees: numeric angle in degrees to rotate existing nodes
#'    * percent_spacing: absolute spacing to enforce for nodes in the nodeset.
#' @param node_adj `data.frame` recognizing five colnames:
#'    * node: `character` name of each node
#'    * x,y: `numeric` x,y coordinate adjustments to add to existing coordinates
#'    * label.dist: `numeric` value to add to existing label distance.
#'    * label.degree: `numeric` angle in degrees to add to the existing angle.
#' @param apply_negative `logical` whether to allow shrinking nodeset spacing.
#' @param nodesets `list` as output from `get_cnet_nodesets()`, or any
#'    list named by nodeset, containing `character` vector of node names.
#' @param ... additional arguments are ignored.
#'
#' @export
bulk_cnet_adjustments <- function
(g,
 nodeset_adj,
 node_adj,
 apply_negative=TRUE,
 nodesets=NULL,
 ...)
{
   #
   # operate on a copy of the data
   adj_cnet <- g;

   # nodeset adjustments
   if (inherits(nodeset_adj, "data.frame") && nrow(nodeset_adj) > 0) {

      if (length(nodesets) == 0) {
         nodesets <- get_cnet_nodeset(g)
      }

      ## test the difference from starting point
      # get nodeset spacing
      default_spacings <- summarize_node_spacing(g,
         scaled=TRUE,
         each_group=TRUE,
         node_groups=nodesets,
         dist_type="nearest_node")
      # default values
      default_spacing <- jamba::nameVector(
         round(
            jamba::rmInfinite(default_spacings$nearest_within[, 'Median'],
               infiniteValue=0) * 10) / 10,
         rownames(default_spacings$nearest_within))
      percent_spacing_diff <- (
         nodeset_adj$percent_spacing -
            default_spacing[rownames(nodeset_adj)]);
      # change only nodes that require it
      nodeset_which <- which(nodeset_adj$x != 0 |
            nodeset_adj$y != 0 |
            (percent_spacing_diff != 0 & nodeset_adj$percent_spacing != 0) |
            nodeset_adj$rotate_degrees != 0);
      # adjust each
      for (i in nodeset_which) {
         adj_cnet <- adjust_cnet_nodeset(adj_cnet,
            set_nodes=strsplit(nodeset_adj$nodeset[i], ",")[[1]],
            x=nodeset_adj$x[i],
            y=nodeset_adj$y[i],
            apply_negative=apply_negative,
            percent_spacing=nodeset_adj$percent_spacing[i],
            rotate_degrees=nodeset_adj$rotate_degrees[i])
      }
   }

   # node adjustments
   if (inherits(node_adj, "data.frame") && nrow(node_adj) > 0) {
      # node x,y coordinates
      node_which <- which(node_adj$x != 0 |
            node_adj$y != 0);
      if (length(node_which) > 0) {
         adj_cnet <- nudge_igraph_node(adj_cnet,
            nodes=node_adj$node[node_which],
            x=node_adj$x[node_which],
            y=node_adj$y[node_which])
      }

      # label.dist
      node_which_l <- which(node_adj$label.degree != 0 |
            node_adj$label.dist != 0);
      if (length(node_which_l) > 0) {
         matchv <- match(node_adj$node[node_which_l],
            igraph::V(adj_cnet)$name);
         # distance
         igraph::V(adj_cnet)[matchv]$label.dist <- (
            igraph::V(adj_cnet)[matchv]$label.dist +
               node_adj$label.dist[node_which_l])
         # angle
         igraph::V(adj_cnet)[matchv]$label.degree <- (
            igraph::V(adj_cnet)[matchv]$label.degree +
               jamba::deg2rad(
                  node_adj$label.degree[node_which_l])) %% (pi*2);
      }
   }

   return(adj_cnet)
}
