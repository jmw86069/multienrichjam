
#' Summarize spacing between igraph nodes or node groups
#'
#' Summarize spacing between igraph nodes or node groups
#'
#' This function is a simple wrapper to calculate typical distances
#' between nodes or node groups in a given network layout.
#'
#' It is intended to help determine values to use for `percent_spacing`
#' with `apply_nodeset_spacing()`.
#' The optimal value is subjective, and depends upon
#' the total number of nodes, and the relative spacing of nodes in
#' subgroups such as Cnet nodesets, or igraph community clusters.
#'
#' Specific to Cnet plots and bipartite graphs, there are two
#' important metrics:
#'
#' 1. Within-group spacing
#'
#'    * 'nearest_within' matrix, column 'Median' is useful
#'    to assess the spacing within each group of nodes.
#'
#' 2. Across-group spacing
#'
#'    * 'nearest_across' matrix, column 'Min' is useful to assess
#'    the minimum distance between two node groups.
#'    Column 'Q1' may be useful to determine whether a node group has
#'    substantial overlaps with another node group.
#'
#' Even still, sometimes the node groups overlap each other, which
#' is not summarized here.
#'
#' # Todo
#'
#' Re-evaluate how to represent summary values that cannot be calculated.
#' For example, a nodeset with only one node cannot calculate within-nodeset
#' distances.
#'
#' Background:
#'
#' * Typically, `min(NULL)` returns `Inf` (infinite).
#' If there were no node-node distances, the calculation would return `Inf`.
#' * Currently, this function returns `NA` for such cases.
#' The main motivation is for downstream summary statistics, for example
#' taking column median or column mean with argument `na.rm=TRUE`
#' for convenience.
#'
#'
#' @family jam utility functions
#'
#' @return `list` of `numeric` matrix data, named by the summary used.
#'    * Columns: 'Min', 'Q1', 'Median', 'Mean', 'Q3', 'Max'.
#'    * Rows when `each_group=TRUE` include each `names(node_groups)`.
#'    * Rows with `each_group=FALSE` include 'all', 'edge', 'nonedge'.
#'    The 'edge' represents only node-node connected by an edge;
#'    'nonedge' represents only node-node not connected by an edge;
#'    and 'all' includes all node-node pairings.
#'    * The list names when `each_group=FALSE`:
#'
#'       * 'summary' with overall summary values
#'
#'    * The list names when `each_group=TRUE`:
#'
#'       * 'nearest_within': nearest-node, within each node group
#'       * 'nearest_across': nearest-node, across each node group
#'       * 'all_within': all-nodes, within each node group
#'       * 'all_across': all-nodes, across each node group
#'
#' @param g `igraph` object
#' @param layout `numeric` matrix with x,y coordinates, or NULL (default)
#'    to call `get_igraph_layout()` which uses either:
#'    `igraph::graph_attr(g, 'layout')` (preferred) or
#'    `igraph::V(g)$x` and `igraph::V(g)$y` (deprecated in igraph).
#' @param nodes `character` with optional node names, or `integer`
#'    index of nodes in `g` to define a subset of nodes for which
#'    statistics are calculated. Useful to focus on a specific subset
#'    of nodes, for example one or two Cnet nodesets.
#' @param scaled `logical`, default TRUE, whether to report percent spacing
#'    relative to the max x-axis/y-axis range, similar to the `min_percent`
#'    and `percent_spacing` argument units in other node spacing functions.
#'    Note that percent spacing is scaled from 0 to 100.
#' @param each_group `logical`, default TRUE, whether to summarize each
#'    node group. When `node_groups` is NULL, by default it will use
#'    `get_cnet_nodesets()`.
#' @param node_groups `list` of vectors named by node group, where each
#'    vector contains vertex names `igraph::V(g)$name`.
#'    * The purpose is to help summarize spacing within and across node groups.
#'    * For Cnet plot data, when `each_group=TRUE` and `node_groups=NULL`,
#'    it will call `get_cnet_nodesets()` to use Cnet nodesets as node groups.
#'    * Cnet nodesets were the primary motivation for this function,
#'    however it also works well using network communities.
#' @param dist_type `character` string, default 'nearest_node' with
#'    the distance summary to provide:
#'    * `"nearest_node"` only the nearest node distance,
#'    helpful to assess whether all nodes have a minimum distance
#'    from other nodes.
#'    * `"all_nodes"` all node distances from each node,
#'    helpful to assess the overall spacing between nodes.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to internal functions.
#'
#' @examples
#' cnet <- make_cnet_test(num_sets=2)
#' igraph::V(cnet)$size <- igraph::V(cnet)$size * 2;
#' igraph::V(cnet)$label.cex <- igraph::V(cnet)$label.cex * 2;
#' jam_igraph(cnet)
#' summarize_node_spacing(cnet)
#' summarize_node_spacing(cnet, dist_type="all_nodes")
#'
#' @export
summarize_node_spacing <- function
(g,
 layout=NULL,
 nodes=NULL,
 scaled=TRUE,
 each_group=TRUE,
 node_groups=NULL,
 dist_type=c("all",
    "nearest_node",
    "all_nodes"),
 verbose=FALSE,
   debug=FALSE,
 ...)
{
   # validate arguments
   dist_type <- match.arg(dist_type);

   # iterate "all"
   if ("all" %in% dist_type) {
      sns_list <- lapply(c("nearest_node", "all_nodes"), function(idist_type){
         sns <- summarize_node_spacing(
            g=g,
            layout=layout,
            nodes=nodes,
            scaled=scaled,
            each_group=each_group,
            node_groups=node_groups,
            dist_type=idist_type,
            verbose=verbose,
            ...)
      })
      return(unlist(recursive=FALSE, sns_list));
   }

   xy_new <- get_igraph_layout(g=g,
      layout=layout,
      verbose=verbose,
      ...)

   # distance between all nodes
   d <- as.matrix(dist(xy_new));

   # try to populate node_groups when each_group==TRUE
   if (TRUE %in% each_group && length(node_groups) == 0) {
      if ("nodeType" %in% igraph::vertex_attr_names(g)) {
         node_groups <- get_cnet_nodeset(g)
      } else {
         each_group <- FALSE;
      }
   }

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
   im <- as.matrix(igraph::as_adjacency_matrix(g)) * 1;

   # custom function to take nearest node distance
   # then summarize
   row_min_summary <- function(d) {
      nearest_d <- apply(d, 1, function(id){
         if (length(id[id > 0]) > 0) {
            min(id[id > 0], na.rm=TRUE)
         } else {
            NA_real_;
         }
      })
      head(summary(nearest_d), 6)
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
   edge_summary <- head(summary(de_new[de_new > 0]), 6)
   edge_nearest_summary <- row_min_summary(de_new);

   # summarize distance between all nodes
   all_summary <- head(summary(d[d > 0]), 6)
   all_nearest_summary <- row_min_summary(d);

   # summarize distance between unconnected nodes
   dn_new <- d * (1 - im);
   nonedge_summary <- head(summary(dn_new[dn_new > 0]), 6)
   nonedge_nearest_summary <- row_min_summary(dn_new);

   retvals <- list();

   ## summary values - not useful for each_group=TRUE
   if (!each_group) {
      if ("nearest_node" %in% dist_type) {
         retvals$all_summary <- all_nearest_summary;
         retvals$edge_summary <- edge_nearest_summary;
         retvals$nonedge_summary <- nonedge_nearest_summary;
      } else {
         retvals$all_summary <- all_summary;
         retvals$edge_summary <- edge_summary;
         retvals$nonedge_summary <- nonedge_summary;
      }
   }

   # iterate optional node_groups
   for (iname in names(node_groups)) {
      node_group <- node_groups[[iname]]

      ## within group
      nmatch <- match(node_group, rownames(d));
      d_ng <- d[nmatch, nmatch, drop=FALSE];
      ng_all_summary <- head(summary(d_ng[d_ng > 0]), 6)
      # nearest node distance within group
      nearest_summary_within <- row_min_summary(d_ng);

      ## across groups
      d_nong <- d[nmatch, -nmatch, drop=FALSE];
      nong_all_summary <- head(summary(d_nong[d_nong > 0]), 6)
      # nearest node distance across group
      nearest_summary_across <- row_min_summary(d_nong);
      # encode into the retvals list
      if ("nearest_node" %in% dist_type) {
         nearest_ng_name <- paste0(iname, "_nearest_within");
         nearest_nong_name <- paste0(iname, "_nearest_across");
         retvals[[nearest_ng_name]] <- nearest_summary_within;
         retvals[[nearest_nong_name]] <- nearest_summary_across;
      } else {
         ng_name <- paste0(iname, "_all_within");
         nong_name <- paste0(iname, "_all_across");
         retvals[[ng_name]] <- ng_all_summary;
         retvals[[nong_name]] <- nong_all_summary;
      }
   }

   if (debug) {
      return(retvals);
   }
   retvals <- tryCatch({
      rdf <- jamba::rbindList(retvals);
      rdf <- jamba::renameColumn(rdf,
         from=c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."),
         to=c("Min", "Q1", "Median", "Mean", "Q3", "Max"))

      # split
      if (any(grepl("_within|_across|_summary", rownames(rdf)))) {
         rdf_split <- gsub("^(.+)_((nearest|all)_(within|across)|summary)$", "\\2",
            rownames(rdf));
         rdf_split <- factor(rdf_split,
            levels=jamba::provigrep(c("within", "."), unique(rdf_split)));
         rdf <- lapply(split(rownames(rdf), rdf_split), function(i){
            im <- rdf[i, , drop=FALSE];
            rownames(im) <- gsub("^(.+)_((nearest|all)_(within|across)|summary)$", "\\1",
               rownames(im))
            im
         })
      }

      # Previous summary values "within" rows mean,median
      if (FALSE) {
         # "within" rows mean,median
         within_rows <- jamba::vigrep("within", rownames(rdf));
         if (length(within_rows) > 0) {
            within_mean <- apply(rdf[within_rows, , drop=FALSE], 2, mean, na.rm=TRUE);
            within_median <- apply(rdf[within_rows, , drop=FALSE], 2, median, na.rm=TRUE);
            retvals <- c(retvals, list(within_mean=within_mean,
               within_median=within_median));
            rdf <- jamba::rbindList(retvals);
         }

         # "across" rows mean,median
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
      }

      rdf
   }, error=function(e){
      print(e)
      retvals
   });
   # Replace NaN with NA
   retvals <- lapply(retvals, jamba::rmNA, naValue=NA)

   return(retvals)
}
