
#' Summarize spacing between igraph nodes given a layout
#'
#' Summarize spacing between igraph nodes given a layout
#'
#' This function is a simple wrapper to calculate typical distances
#' between nodes in a given network layout. It is experimental,
#' and intended to provide helpful information when determining
#' an appropriate value for `percent_spacing` to use for example
#' with `apply_nodeset_spacing()`. The optimal value depends upon
#' the number of nodes overall, also the number of nodes in
#' each nodeset, and the relative position of each nodeset
#' in layout coordinates.
#'
#' The `node_groups` argument is intended to provide summary data
#' for each node group (for example Cnet `nodesets`) so that
#' individual node groups can be adjusted accordingly.
#'
#' @family jam utility functions
#'
#' @return `data.frame` with summary information about node distances
#'    for connected nodes (only where edges connect any two nodes),
#'    and unconnected nodes (only where two nodes are not connected
#'    by an edge). When `node_groups` is defined, the summary
#'    also includes each individual node group.
#'
#' @param g `igraph` object
#' @param layout passed together with `g` to `get_igraph_layout()`
#' @param nodes `character` with optional node names, or `integer`
#'    index of nodes in `g` to define a subset of nodes for which
#'    statistics are calculated. Useful to focus on a specific subset
#'    of nodes, for example one or two Cnet nodesets.
#' @param node_groups `list` implemented to use nodesets. The intent is to
#'    define node groups, then calculate statistics of node spacing across
#'    and within node groups.
#' @param each_group `logical` indicating whether to include each node group
#'    when `node_groups` is also supplied.
#' @param scaled `logical` indicating whether to report spacing
#'    relative to the max x-axis/y-axis range, similar to the `min_percent`
#'    and `percent_spacing` argument units in other node spacing functions.
#' @param dist_type `character` string indicating the type of distance to
#'    summarize:
#'    * `"nearest_node"` only uses the nearest node to each node, which is
#'    helpful when trying to ensure all nodes have a minimum distance
#'    from other nodes.
#'    * `"all_nodes"` uses all node distances from each node, which is
#'    helpful when assessing the overall spacing between nodes.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to internal functions.
#'
#' @export
summarize_node_spacing <- function
(g,
 layout=NULL,
 nodes=NULL,
 node_groups=NULL,
 each_group=FALSE,
 scaled=FALSE,
 dist_type=c("nearest_node",
    "all_nodes"),
 verbose=FALSE,
 ...)
{
   #
   dist_type <- match.arg(dist_type);
   xy_new <- get_igraph_layout(g=g,
      layout=layout,
      verbose=verbose,
      ...)

   # distance between all nodes
   d <- as.matrix(dist(xy_new));

   # optionally scale by the overall layout range
   max_xy <- max(apply(xy_new, 2, function(i){
      diff(range(i, na.rm=TRUE))
   }))
   if (scaled) {
      if (verbose) {
         jamba::printDebug("summarize_node_spacing(): ",
            "scaling distances by max range (",
            format(max_xy, digits=3), ") to produce units as percent layout.");
      }
      d <- d / max_xy * 100;
   }

   # incidence matrix of graph edges
   im <- as.matrix(as_adjacency_matrix(g)) * 1;

   # custom function to take nearest node distance
   # then summarize
   row_min_summary <- function(d) {
      nearest_d <- apply(d, 1, function(id){
         min(id[id > 0], na.rm=TRUE)
      })
      summary(nearest_d, na.rm=TRUE);
   }

   # optional subset of nodes
   if (length(nodes) > 0) {
      nmatch <- match(nodes, rownames(d));
      if (any(is.na(nmatch))) {
         stop("Input nodes were not all present in V(g)$name")
      }
      d <- d[nmatch, nmatch, drop=FALSE]
      im <- im[nmatch, nmatch, drop=FALSE]
   }

   # distance only of edges
   de_new <- d * im;

   # summarize distance between connected nodes
   edge_summary <- summary(de_new[de_new > 0])
   edge_nearest_summary <- row_min_summary(de_new);

   # summarize distance between all nodes
   all_summary <- summary(d[d > 0])
   all_nearest_summary <- row_min_summary(d);

   # summarize distance between unconnected nodes
   dn_new <- d * (1 - im);
   nonedge_summary <- summary(dn_new[dn_new > 0])
   nonedge_nearest_summary <- row_min_summary(dn_new);

   retvals <- list();
   if ("nearest_node" %in% dist_type) {
      retvals$edge_summary <- edge_nearest_summary;
      retvals$nonedge_summary <- nonedge_nearest_summary;
      retvals$all_summary <- all_nearest_summary;
   } else {
      retvals$edge_summary <- edge_summary;
      retvals$nonedge_summary <- nonedge_summary;
      retvals$all_summary <- all_summary;
   }

   # iterate optional node_groups
   for (iname in names(node_groups)) {
      node_group <- node_groups[[iname]]
      # summary within group
      nmatch <- match(node_group, rownames(d));
      d_ng <- d[nmatch, nmatch, drop=FALSE];
      ng_all_summary <- summary(d_ng[d_ng > 0]);
      # nearest node distance within group
      nearest_summary_within <- row_min_summary(d_ng);
      # summary across groups
      d_nong <- d[nmatch, -nmatch, drop=FALSE];
      nong_all_summary <- summary(d_nong[d_nong > 0])
      # nearest node distance across group
      nearest_summary_across <- row_min_summary(d_nong);
      # encode into the retvals list
      if ("nearest_node" %in% dist_type) {
         nearest_ng_name <- paste0(iname, "_nearest_within");
         nearest_nong_name <- paste0(iname, "_nearest_across");
         retvals[[nearest_ng_name]] <- nearest_summary_within;
         retvals[[nearest_nong_name]] <- nearest_summary_across;
      } else {
         ng_name <- paste0(iname, "_within_summary");
         nong_name <- paste0(iname, "_across_summary");
         retvals[[ng_name]] <- ng_all_summary;
         retvals[[nong_name]] <- nong_all_summary;
      }
   }

   retvals <- tryCatch({
      rdf <- jamba::rbindList(retvals);
      within_rows <- jamba::vigrep("within", rownames(rdf));
      if (length(within_rows) > 0) {
         within_mean <- apply(rdf[within_rows, , drop=FALSE], 2, mean, na.rm=TRUE);
         within_median <- apply(rdf[within_rows, , drop=FALSE], 2, median, na.rm=TRUE);
         retvals <- c(retvals, list(within_mean=within_mean,
            within_median=within_median));
         rdf <- jamba::rbindList(retvals);
      }
      across_rows <- jamba::vigrep("across", rownames(rdf));
      if (length(across_rows) > 0) {
         across_mean <- apply(rdf[across_rows, , drop=FALSE], 2, mean, na.rm=TRUE);
         across_median <- apply(rdf[across_rows, , drop=FALSE], 2, median, na.rm=TRUE);
         retvals <- c(retvals, list(across_mean=across_mean,
            across_median=across_median));
         rdf <- jamba::rbindList(retvals);
      }
      head_rows <- setdiff(rownames(rdf), c(within_rows, across_rows));
      if (each_group) {
         rdf_rows <- c(head_rows,
            within_rows,
            across_rows);
      } else {
         rdf_rows <- head_rows;
      }
      rdf <- rdf[rdf_rows, , drop=FALSE];
      names(dimnames(rdf)) <- c("node_group",
         paste0(dist_type,
            ifelse(TRUE %in% scaled, "_percent", "_dist")));
      rdf
   }, error=function(e){
      print(e)
      retvals
   });

   return(retvals)
}
