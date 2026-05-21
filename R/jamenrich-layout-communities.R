
#' Layout igraph communities
#' 
#' Layout igraph communities, applies layout to each node community
#' 
#' This function mimics `igraph::layout_components()`, in that it
#' splits the input `igraph` into a `list`, applies layout to
#' each graph, then merges the layouts together.
#' 
#' @returns `igraph`
#' 
#' @param g `igraph` with graph attributes that include either
#'    'mark.groups' or 'nodegroups'.
#' @param layout `function` to apply to each igraph community.
#' @param ... additional arguments are ignored.
#' 
#' @export
layout_communities <- function
(g,
 layout=layout_with_qfrf(repulse=3),
 ...)
 {
   #
   split_igraph_by_nodegroups <- function
   (g,
    nodegroups)
   {
      if (is.list(nodegroups)) {
         nodegroups <- jamba::nameVector(
            rep(names(nodegroups), lengths(nodegroups)),
            unname(unlist(nodegroups)))
      }
      nodegroups <- nodegroups[igraph::V(g)$name];
      lapply(split(igraph::V(g), nodegroups), function(vids) {
         igraph::subgraph(g, vids)
      })
   }
   mark.groups <- igraph::graph_attr(g, "mark.groups");
   nodegroups <- igraph::graph_attr(g, "nodegroups")
   if (length(mark.groups) > 0) {
      nodegroups <- multienrichjam::communities2nodegroups(mark.groups)
   }
   # split igraph into list of igraph
   g_list <- split_igraph_by_nodegroups(g=g,
      nodegroups=nodegroups)
      g_list_layout <- igraph::merge_coords(g_list,
      layouts=lapply(g_list, layout))
   rownames(g_list_layout) <- unlist(
      lapply(unname(g_list), function(i){igraph::V(i)$name}))

   # match original igraph
   v_match <- match(igraph::V(g)$name, rownames(g_list_layout))
   g_list_layout_new <- g_list_layout[v_match, , drop=FALSE];
   # assign to original g (the o.g. haha) for consistency
   igraph::graph_attr(g, "layout") <- g_list_layout_new
   return(g);
   
   # recombine igraph list
   g_new <- igraph::disjoint_union(g_list)
   # assign layout
   igraph::graph_attr(g_new, "layout") <- g_list_layout
   # re-assign mark.groups
   if (length(mark.groups) > 0) {
      igraph::graph_attr(g_new, "mark.groups") <- mark.groups
   }
   if (length(nodegroups) > 0) {
      igraph::graph_attr(g_new, "nodegroups") <- nodegroups
   }
}
