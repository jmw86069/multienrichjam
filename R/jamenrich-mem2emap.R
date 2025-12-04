
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
#' @family jam Mem utilities
#' @family jam igraph functions
#'
#' @param mem `Mem` or legacy `list` mem output from `multiEnrichMap()`
#' @param overlap `numeric`, default 0.2, value between 0 and 1 with Jaccard
#'    overlap coefficient required between any two pathways in order to
#'    create a network edge connecting these two pathways.
#'    * Jaccard overlap coefficient is reciprocal of Jaccard distance,
#'    using `(1 - Jaccard_distance)`.
#'    * `overlap=0.2` is default, which specifies roughly 20% overlap
#'    in genes shared between two pathway nodes.
#'    * Note these genes must be involved in enrichment, and
#'    therefore does not use all possible genes annotated to a pathway.
#'    Therefore connections are only created with enriched genes
#'    are shared between pathways.
#' @param p_cutoff `numeric` threshold used for significant enrichment
#'    P-value. The default NULL uses the value in the `mem` object.
#' @param min_count `integer` threshold for minimum genes involved in
#'    enrichment in order for a pathway to be considered significant
#'    during this analysis. When NULL it uses the threshold in `mem`.
#' @param colorV `character` vector of R colors used for each enrichment.
#'    Default NULL uses colors defined in `mem`.
#' @param cluster_function `function`, default `igraph::cluster_walktrap`,
#'    used to cluster `igraph` nodes in the resulting network graph.
#'    It is used to enhance the visual summary.
#'    * This function is not used when `cluster_list` is supplied.
#' @param cluster_list `list` default NULL, optional list of pathway
#'    clusters, containing `character` vectors of `igraph` node/vertex
#'    names. It will be used instead of `cluster_function` when supplied.
#' @param num_keep_terms `integer`, default 3, number of text terms to
#'    keep for each pathway cluster, used together with `cluster_function`.
#'    Common words are removed, remaining terms are sorted by decreasing
#'    occurrence, then used to summarize each cluster.
#' @param keep_terms_sep `character` string, default is comma-newline,
#'    used to separate terms, used together with `num_keep_terms`.
#' @param repulse `numeric` value passed to `layout_with_qfr()`,
#'    default 3.3.
#'    * Use repulse 'FALSE', 'NULL', or '0' will skip the layout, which
#'    is used by `mem_find_overlap()` to iterate numerous overlaps
#'    without spending time on layout each iteration.
#' @param remove_singlets `logical`, default TRUE, whether to remove pathway
#'    singlets which have no connections to any other pathways.
#'    * Using TRUE will help simplify busy figures, at the expense of
#'    being less complete.
#' @param color_by_nodes `logical`, default TRUE, whether to colorize pathway
#'    clusters using node colors in each cluster.
#'    Note that a mix of colors often turns brown, so this feature has
#'    unpredictable benefit.
#' @param apply_direction `logical` default TRUE, whether to apply direction
#'    to the node border color, when supplied. Used with `direction_max`
#'    and `direction_floor`.
#' @param direction_max `numeric`, default 2, indicating the directional
#'    score at which the maximum color is applied. Typically when using
#'    some derivative of z-score to represent directionality, values at or
#'    above 2 are effectively "maximum".
#' @param direction_floor `numeric` default 0.5, the minimum directional score
#'    for any color to be applied, below which values are considered
#'    "non-directional". Typically when using some form of z-score, a value
#'    at or below 0.5 may be considered "not directional."
#' @param seed `numeric` passed to `set.seed()` via the `layout_with_qfr()`
#'    layout algorithm, default 123.
#' @param do_express `logical` default FALSE, when TRUE it skips a number
#'    of aesthetic steps. This option is intended mainly for
#'    `mem_find_overlap()` to perform more rapid iterative evaluation
#'    of overlap thresholds.
#' @param do_plot `logical`, default FALSE, whether to render the resulting
#'    plot using `jam_igraph()`.
#' @param verbose `logical` whether to print verbose output.
#' @param ... additional arguments are passed to `jam_igraph()` to customize
#'    the network plot, used when `do_plot=TRUE`.
#' 
#' @examples
#' emap <- mem2emap(Memtest, min_count=2,)
#' jam_igraph(emap)
#' 
#' emap2 <- relayout_nodegroups(emap, final_repulse=4, y_bias=10, label_min_dist=2.5)
#' jam_igraph(emap2)
#' 
#' emapB <- mem2emap(fixSetLabels(Memtest), min_count=2, overlap=0.35)
#' jam_igraph(emapB)
#' 
#' @export
mem2emap <- function
(mem,
 overlap=0.2,
 p_cutoff=NULL,
 min_count=NULL,
 colorV=NULL,
 cluster_function=igraph::cluster_walktrap,
 cluster_list=NULL,
 num_keep_terms=3,
 keep_terms_sep=",\n",
 repulse=3.3,
 remove_singlets=FALSE,
 color_by_nodes=FALSE,
 size_by_genes=TRUE,
 median_size=5,
 apply_direction=TRUE,
 direction_max=2,
 direction_floor=0.5,
 seed=123,
 spread_labels=TRUE,
 label_min_dist=2.5,
 y_bias=10,
 vertex.label.font=2,
 use_shadowText=TRUE,
 do_express=FALSE,
 do_plot=FALSE,
 verbose=FALSE,
 ...)
{
   #
   # determine nodes to include
   if (is.list(mem) && "enrichIM" %in% names(mem)) {
      mem <- list_to_Mem(mem)
   }
   if (!inherits(mem, "Mem")) {
      stop("Input 'mem' must be class 'Mem' or coercible to 'Mem'.");
   }
   # apply default values from mem when not supplied
   if (length(p_cutoff) == 0) {
      if ("p_cutoff" %in% names(thresholds(mem))) {
         p_cutoff <- thresholds(mem)$p_cutoff;
      } else if ("cutoffRowMinP" %in% names(thresholds(mem))) {
         p_cutoff <- thresholds(mem)[["cutoffRowMinP"]];
      } else {
         p_cutoff <- 1;
      }
   }
   if (verbose && p_cutoff < 1) {
      jamba::printDebug("mem2emap(): ",
         "p_cutoff: ", p_cutoff);
   }
   if (length(min_count) == 0) {
      if ("min_count" %in% names(thresholds(mem))) {
         min_count <- thresholds(mem)$min_count;
      } else {
         min_count <- 1;
      }
   }
   if (verbose && min_count > 1) {
      jamba::printDebug("mem2emap(): ",
         "min_count: ", min_count);
   }
   if (length(colorV) == 0) {
      colorV <- colorV(mem)
   }
   
	# define enrich_use while applying filtering
   enrich_use <- (enrichIM(mem) <= p_cutoff) * 1;
   # optionally apply min_count when defined, and when present in the mem data
   if (min_count > 1) {
      enrich_use <- (enrich_use > 0 &
            enrichIMgeneCount(mem) >= min_count) * 1;
   }
   rownames(enrich_use) <- rownames(enrichIM(mem));
   
   # determine which rows to use
   if (length(cluster_list) > 0) {
      enrich_rows_use <- unname(unlist(cluster_list));
   } else {
      enrich_rows_use <- rownames(enrich_use)[rowSums(enrich_use) > 0]
   }
   
   # indicate whether the filtering reduced the overall pathways
   if (TRUE %in% verbose && length(enrich_rows_use) < nrow(enrichIM(mem))) {
      jamba::printDebug("mem2emap(): ",
         "Note that ",
         (nrow(enrichIM(mem)) - length(enrich_rows_use)),
         " pathways were removed due to filtering thresholds.");
   }
   if (length(enrich_rows_use) == 0) {
      stop_msg <- paste0(
         "No enrichment results met the thresholds given: p_cutoff<=",
         format(digits=3, p_cutoff),
         ", min_count>=", min_count);
      stop(stop_msg);
   }

   # convert memIM to Jaccard overlap matrix
   col_match <- match(enrich_rows_use,
      colnames(memIM(mem)));
   if (any(is.na(col_match))) {
      stop_msg <- paste0("There is a mismatch in ",
      	"rownames(enrichIM(mem)) and colnames(memIM(mem)).");
      stop(stop_msg);
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
   
   if (isTRUE(do_express)) {
   	size_by_genes <- FALSE;
   }
   if (!isTRUE(do_express)) {
	   # calculate gene count per overlap
	   # - note convert non-zero, non-NA to 1, all else 0
	   genecount_matrix <- crossprod((!memIM(mem) == 0 & !is.na(memIM(mem))) * 1);
	   # - multiply by non-zero Jaccard overlap so that it uses identical edges
	   genecount_matrix <- genecount_matrix * (jacc_overlap_filtered > 0);
	   # Create the igraph using genecount weights
	   jacc_g_ct <- igraph::graph_from_adjacency_matrix(
	   	genecount_matrix,
	   	mode="undirected",
	   	diag=FALSE,
	   	weighted=TRUE,
	   	add.colnames=NULL)
	   igraph::E(jacc_g)$gene_count <- igraph::E(jacc_g_ct)$weight;
	
	   # add gene count to vertex attributes
	   genecount_set <- colSums((!memIM(mem) == 0 & !is.na(memIM(mem))) * 1);
	   names(genecount_set) <- colnames(memIM(mem));
	   igraph::V(jacc_g)$gene_count <- genecount_set[igraph::V(jacc_g)$name];
   }

   # Define size by vcount
   if (isTRUE(size_by_genes)) {
   	gc_sizes <- sqrt(igraph::V(jacc_g)$gene_count / 
   			median(igraph::V(jacc_g)$gene_count)) * median_size;
   	igraph::V(jacc_g)$size <- gc_sizes;
   } else {
	   igraph::V(jacc_g)$size <- median_size;
   }

   # define some default aesthetics
   if (length(vertex.label.font) == 1) {
   	igraph::V(jacc_g)$label.font <- vertex.label.font;
   }
   if (isTRUE(use_shadowText)) {
   	igraph::graph_attr(jacc_g, "use_shadowText") <- use_shadowText;
   }
   
   # optionally remove singlets
   if (isTRUE(remove_singlets) && length(cluster_list) == 0) {
      if (verbose) {
         jamba::printDebug("mem2emap(): ",
            "Removing singlet nodes.");
      }
      jacc_g <- removeIgraphSinglets(jacc_g)
   }
   # 0.0.108.900: check for empty graph
   if (igraph::vcount(jacc_g) == 0) {
   	return(NULL)
   }

   # define layout
   # Todo: consider using igraph::layout_components()
   new_layout <- NULL;
   if (!isTRUE(do_express)) {
	   jacc_g_comp <- igraph::components(jacc_g);
	   if (jacc_g_comp$no > 1) {
	      if (verbose) {
	         jamba::printDebug("mem2emap(): ",
	            "Applying layout_components() with layout_with_qfrf().");
	      }
	   	# skip layout when repulse == 0, used by mem_find_overlap()
	      if (length(repulse) == 0 || isFALSE(repulse) || repulse == 0) {
	         new_layout <- NULL;
	      } else {
	         new_layout <- igraph::layout_components(jacc_g,
	            layout=layout_with_qfrf(repulse=repulse,
	               seed=seed,
	               ...));
	         if (!all(rownames(new_layout) %in% igraph::V(jacc_g)$name)) {
	            rownames(new_layout) <- igraph::V(jacc_g)$name;
	         }
	      }
	   } else {
	      if (verbose) {
	         jamba::printDebug("mem2emap(): ",
	            "Applying layout_with_qfr().");
	      }
	      if (length(repulse) == 0 || isFALSE(repulse) || repulse == 0) {
	         new_layout <- NULL;
	      } else {
	         new_layout <- layout_with_qfr(
	            jacc_g,
	            seed=seed,
	            repulse=repulse,
	            ...)
	      }
	   }
   }
   if (!is.null(new_layout)) {
	   igraph::graph_attr(jacc_g, "layout") <- new_layout;
   }

   # assign node colors
   wc <- NULL;
   mark.colors <- NULL;
   nodegroups_wc <- NULL;
   if (!isTRUE(do_express)) {
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
	   wc <- NULL;
	   mark.colors <- NULL;
	   nodegroups_wc <- NULL;
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
	      # nodegroups_wc <- split(igraph::V(jacc_g)$name, wc$membership)
	      nodegroups_wc <- communities2nodegroups(wc);
	      
	      # assign most common terms as a cluster label
	      nodegroups_wc <- tryCatch({
	      	wc <- label_communities(wc,
		         keep_terms_sep=keep_terms_sep,
		         num_keep_terms=num_keep_terms);
		      nodegroups_wc_labels <- wc$cluster_names;
		      names(nodegroups_wc) <- jamba::cPaste(nodegroups_wc_labels,
		         keep_terms_sep)
		      nodegroups_wc
	      }, error=function(e){
	      	# jamba::printDebug("vcount:", igraph::vcount(jacc_g));# debug
	      	# jamba::printDebug("Using nodegroups_wc fallback:");print(nodegroups_wc);# debug
	      	nodegroups_wc;
	      })
	
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
	      mark.colors <- NULL;
	      nodegroups_wc <- NULL;
	   }
   }
   
   # Check for NULL layout
   if ("layout" %in% igraph::graph_attr_names(jacc_g)) {
   	use_layout <- igraph::graph_attr(jacc_g, "layout");
   	if (is.null(use_layout) || length(use_layout) == 0) {
   		jacc_g <- igraph::delete_graph_attr(jacc_g, "layout")
   	}
   }
   # optionally spread labels
   if ("layout" %in% igraph::graph_attr_names(jacc_g) &&
   		isTRUE(spread_labels)) {
   	jacc_g <- spread_igraph_labels(jacc_g,
   		label_min_dist=label_min_dist,
   		y_bias=y_bias,
   		...)
   }

   if (!isTRUE(do_express) && isTRUE(do_plot)) {
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
   if (length(wc) > 0) {
      igraph::graph_attr(jacc_g, "mark.groups") <- wc;
      igraph::graph_attr(jacc_g, "mark.colors") <- mark.colors;
      igraph::graph_attr(jacc_g, "nodegroups") <- nodegroups_wc;
   }
   
   return(invisible(jacc_g));
}

