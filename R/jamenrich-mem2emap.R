
#' Convert multiEnrichMap mem output to EnrichmentMap emap
#'
#' Convert multiEnrichMap mem output to EnrichmentMap emap
#'
#' This function is currently In development.
#'
#' This function is intended to convert `mem` output from
#' `multiEnrichMap()` into an EnrichmentMap `igraph` format
#' which represents the statistical enrichment support from each
#' pathway enrichment.
#'
#' This function can apply P-value thresholds using the input `mem`,
#' or using a custom value.
#'
#' A node clustering function is applied by default, which may help
#' define suitable subgroups of nodes. When defined, the clusters
#' are used to define nodegroups for edge bundling.
#'
#' @family jam igraph functions
#'
#' @param mem `Mem` or legacy `list` mem output from `multiEnrichMap()`
#' @param overlap `numeric` value between 0 and 1 indicating the Jaccard
#'    overlap coefficient required between any two pathways in order to
#'    create a network edge connecting these two pathways. Typically,
#'    `overlap=0.2` is default, which specifies roughly 20% overlap
#'    in genes shared between two pathway nodes. Note these genes must be
#'    involved in enrichment, and therefore does not use all possible
#'    genes annotated to a pathway. Therefore connections are only
#'    created with enriched genes are shared between pathways.
#' @param p_cutoff `numeric` threshold used for significant enrichment
#'    P-value, usually defined in the `mem` object.
#' @param min_count `integer` threshold for minimum genes involved in
#'    enrichment in order for a pathway to be considered significant
#'    during this analysis.
#' @param colorV `character` vector of R colors used for each enrichment.
#' @param cluster_function `function` used to cluster nodes in the resulting
#'    `igraph` object, used to help generate a visual summary.
#' @param num_keep_terms `integer` number of terms to keep from each
#'    pathway cluster, when `cluster_function` is supplied above. Common
#'    terms are removed from each pathway cluster, then remaining terms
#'    are sorted by decreasing occurrence, and used as a straightforward
#'    summary of pathways in each cluster.
#' @param keep_terms_sep `character` string used to separated multiple
#'    pathway terms defined by `num_keep_terms` above.
#' @param repulse `numeric` value passed to `layout_with_qfr()`.
#' @param remove_singlets `logical` indicating whether to remove pathway
#'    singlets, which have no connections to any other pathways. It can
#'    help simplify busy figures, however removing a singlet pathway
#'    is not recommended because it may imply the pathways were not
#'    statistically significant, and in fact they were.
#' @param color_by_nodes `logical` indicating whether to colorize pathway
#'    clusters based upon blending the node colors within each cluster.
#'    Note that a mix of colors often turns brown, so this feature has
#'    unpredictable benefit.
#' @param do_plot `logical` indicating whether to render the resulting
#'    plot.
#' @param ... additional arguments are passed to `jam_igraph()` to customize
#'    the network plot.
#'
#' @export
mem2emap <- function
(mem,
 overlap=0.2,
 p_cutoff=NULL,
 min_count=4,
 colorV=NULL,
 cluster_function=igraph::cluster_walktrap,
 cluster_list=NULL,
 num_keep_terms=3,
 keep_terms_sep=",\n",
 repulse=3.3,
 remove_singlets=TRUE,
 color_by_nodes=FALSE,
 apply_direction=TRUE,
 direction_max=2,
 direction_floor=0.5,
 do_plot=TRUE,
 ...)
{
   #
   # determine nodes to include
   if (is.list(mem) && "enrichIM" %in% names(mem)) {
      mem <- list_to_Mem(mem)
   }
   if (length(p_cutoff) == 0) {
      if ("p_cutoff" %in% names(thresholds(mem))) {
         p_cutoff <- thresholds(mem)$p_cutoff;
      } else if ("cutoffRowMinP" %in% names(thresholds(mem))) {
         p_cutoff <- thresholds(mem)[["cutoffRowMinP"]];
      } else {
         p_cutoff <- 1;
      }
   }
   if (length(colorV) == 0) {
      colorV <- mem@colorV
   }
   
   if (length(cluster_list) > 0) {
      enrich_rows_use <- unname(unlist(cluster_list));
   } else {
      enrich_use <- (enrichIM(mem) <= p_cutoff) * 1;

      # optionally apply min_count when defined, and when present in the mem data
      if (min_count > 1) {
         enrich_use <- (enrich_use > 0 &
               mem@enrichIMgeneCount >= min_count) * 1;
      }
      rownames(enrich_use) <- rownames(enrichIM(mem));
      enrich_rows_use <- rownames(enrich_use)[rowSums(enrich_use) > 0]
   }
   if (length(enrich_rows_use) < nrow(enrichIM(mem))) {
      jamba::printDebug("mem2emap(): ",
         "Note that ",
         (nrow(enrichIM(mem)) - length(enrich_rows_use)),
         " pathways were removed due to filtering thresholds.");
   }
   if (length(enrich_rows_use) == 0) {
      stop(paste0(
         "No enrichment results met the thresholds given: p_cutoff<=",
         format(digits=3, p_cutoff),
         ", min_count>=", min_count));
   }

   # convert memIM to Jaccard overlap matrix
   col_match <- match(enrich_rows_use,
      colnames(memIM(mem)));
   if (any(is.na(col_match))) {
      stop("There is a mismatch in rownames(enrichIM(mem)) and colnames(memIM(mem)).")
   }
   jacc_overlap <- 1 - as.matrix(
      dist(t(memIM(mem)[, col_match, drop=FALSE]),
         method="binary"));

   # convert to graph using Jaccard overlap threshold
   jacc_overlap_filtered <- jacc_overlap * (jacc_overlap >= overlap)
   jacc_g <- igraph::graph_from_adjacency_matrix(
      jacc_overlap_filtered,
      mode="undirected",
      diag=FALSE,
      weighted=TRUE,
      add.colnames=NULL)
   igraph::V(jacc_g)$size <- 5;

   # optionally remove singlets
   if (TRUE %in% remove_singlets && length(cluster_list) == 0) {
      jacc_g <- removeIgraphSinglets(jacc_g)
   }

   # define layout
   igraph::graph_attr(jacc_g, "layout") <- layout_with_qfr(
      jacc_g,
      repulse=repulse)


   imatch <- match(igraph::V(jacc_g)$name, rownames(enrich_use))
   igraph::V(jacc_g)$pie.color <- lapply(imatch, function(i){
      colorV[colnames(enrich_use)[enrich_use[i, ] != 0]]
   })
   igraph::V(jacc_g)$pie <- lapply(igraph::V(jacc_g)$pie.color, function(i){
      rep(1, length(i))
   })
   igraph::V(jacc_g)$shape <- "jampie";

   # optionally apply direction to pie frame color
   has_negative <- any(
      jamba::rmNA(naValue=0, enrichIMdirection(mem)) < 0);
   # Todo: decide if these criteria are flexible enough
   if (TRUE %in% apply_direction && has_negative) {
      # define color gradient for border color
      dir_col_fn <- colorjam::col_div_xf(
         direction_max,
         floor=direction_floor,
         colramp=jamba::getColorRamp(
            "RdBu_r",
            divergent=TRUE,
            trimRamp=c(1, 1)))
      # apply color to borders
      pie_borders <- lapply(seq_len(igraph::vcount(jacc_g)), function(i){
         j <- match(igraph::V(jacc_g)$name[i],
            rownames(enrichIMdirection(mem)));
         k <- igraph::V(jacc_g)$pie.color[[i]];
         knames <- names(k);
         dir_col_fn(enrichIMdirection(mem)[j, knames])
      })
      igraph::V(jacc_g)$pie.border <- pie_borders;
      igraph::V(jacc_g)$pie.lwd <- 3;
      igraph::V(jacc_g)$frame.lwd <- 0.5;
   } else {
      apply_direction <- FALSE
   }

   # cluster nodes
   if (length(cluster_list) > 0) {
      wc <- cluster_list;
      nodegroups_wc <- cluster_list;
      mark.colors <- jamba::alpha2col(alpha=0.3,
         colorjam::rainbowJam(n=length(nodegroups_wc),
            Crange=c(60, 90),
            Lrange=c(50, 85)))

   } else if (is.function(cluster_function)) {
      wc <- cluster_function(jacc_g)
      # define list
      nodegroups_wc <- split(igraph::V(jacc_g)$name, wc$membership)
      nodegroups_wc <- communities2nodegroups(wc);
      # assign most common terms as a cluster label
      wc <- label_communities(wc,
         keep_terms_sep=keep_terms_sep,
         num_keep_terms=num_keep_terms);
      nodegroups_wc_labels <- wc$cluster_names;
      names(nodegroups_wc) <- nodegroups_wc_labels;

      # bonus points: define mark.group colors by node colors
      if (TRUE %in% color_by_nodes) {
         mark.colors <- sapply(nodegroups_wc, function(i){
            ic1 <- igraph::V(jacc_g)[i]$pie.color;
            ic2 <- jamba::alpha2col(alpha=0.3,
               colorjam::blend_colors(unname(unlist(ic1))))
            ic2
         })
      } else {
         mark.colors <- jamba::alpha2col(alpha=0.3,
            colorjam::rainbowJam(n=length(nodegroups_wc),
               Crange=c(60, 90),
               Lrange=c(50, 85)))
      }

      # bonus points: color edges by mark.group colors
   } else {
      wc <- NULL;
      nodegroups_wc <- NULL;
      mark.colors <- NULL;
   }

   if (TRUE %in% do_plot) {
      jam_igraph(jacc_g,
         mark.groups=wc,
         mark.col=mark.colors,
         nodegroups=nodegroups_wc,
         # render_nodelabels=FALSE,
         ...)
      mem_legend("bottomleft",
         mem=mem,
         do_direction=apply_direction)
   }
   return(invisible(jacc_g));
}

