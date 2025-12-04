
#' Average geometric angles
#'
#' Average geometric angles
#'
#' @family jam utility functions
#'
#' @examples
#' x <- c(90, 180);
#' avg_angles(x);
#'
#' @export
avg_angles <- function
(x,
 w=1,
 maxAngle=360,
 na.rm=TRUE,
 ...)
{
   ## Purpose is to average angular values, which otherwise don't average properly,
   ## e.g. 359 degrees and 1 degrees should have average angle of 0, not 180.
   if (length(x) == 0) {
      return(NULL);
   }
   w <- rep(w, length.out=length(x));
   x_values <- cos(x * 2*pi / maxAngle);
   y_values <- sin(x * 2*pi / maxAngle);
   xComp <- weighted.mean(x_values,
      w=w,
      na.rm=na.rm,
      ...);
   yComp <- weighted.mean(y_values,
      w=w,
      na.rm=na.rm,
      ...);
   newAngle <- atan2(y=yComp, x=xComp) * (maxAngle / (2*pi));
   tryCatch({
      if (newAngle < 0) {
         newAngle <- newAngle + maxAngle;
      }
   }, error=function(e){
      jamba::printDebug("avg_angles(): ",
         "Error:", e);
      jamba::printDebug("xComp:", xComp);
      jamba::printDebug("yComp:", yComp);
   });
   newAngle;
}

#' Average colors by list
#'
#' Average colors by list
#'
#' This is a simple wrapper function intended to provide a rapid
#' average color, when supplied a list of color vectors in hex or
#' R color name format.
#'
#' This function simply converts each color to HCL, determines the
#' color hue angle (from 0 to 360) then calculates the average angular
#' color hue using `avg_angles()`, then applies that to the maximum
#' C and L values to determine the new color. It is deliberately intended
#' to ignore muddiness when averaging multiple colors.
#'
#' Colors are only modified for elements with 2 or more entries.
#'
#' This method also only operates on the unique set of colors,
#' so it should be substantially more efficient on large lists that
#' contain only a few unique subsets of colors.
#'
#' @family jam utility functions
#'
#' @return vector of R colors.
#'
#' @param x `list` of character vectors.
#' @param useWeightedHue logical indicating whether to weight the
#'    hue wheel using `colorjam::h2hw()` and `colorjam::hw2h()`
#'    which effectively converts the RGB angles to RYB (red-yellow-blue),
#'    and therefore makes additive color blending more sensible.
#'    Specifically, "yellow and blue makes green".
#' @param ... additional arguments are ignored.
#'
#' @examples
#' x <- list(input1=c(red="red", blue="blue"),
#'    input2=c(blue="blue", gold="gold"),
#'    input3=c(red="red", yellow="yellow"));
#' x_avg <- avg_colors_by_list(x, useWeightedHue=TRUE);
#' jamba::showColors(list(
#'    input1=c(x[[1]], x_avg[1]),
#'    input2=c(x[[2]], x_avg[2]),
#'    input2=c(x[[3]], x_avg[3])),
#'    main="With weighted hue")
#'
#' x_avg <- avg_colors_by_list(x, useWeightedHue=FALSE);
#' jamba::showColors(list(
#'    input1=c(x[[1]], x_avg[1]),
#'    input2=c(x[[2]], x_avg[2]),
#'    input2=c(x[[3]], x_avg[3])),
#'    main="Without weighted hue")
#'
#' @export
avg_colors_by_list <- function
(x,
 useWeightedHue=TRUE,
 Cmethod=c("mean", "max", "min"),
 Lmethod=c("mean", "max", "min"),
 c_min=4,
 grey_hue=359,
 ...)
{
   if (length(x) == 0) {
      return(x);
   }
   Cmethod <- match.arg(Cmethod);
   Lmethod <- match.arg(Lmethod);
   x_full <- x;
   x <- unique(x_full);
   x_match <- match(x_full, x);
   xdo <- which(lengths(x) > 1);
   k <- 1;
   if (length(xdo) > 0) {
      xh <- lapply(unname(x[xdo]), function(i){
         #ihcl <- jamba::col2hcl(i);
         ihcl <- farver::decode_colour(i, to="hcl");
         if (any(ihcl[,"c"] <= c_min)) {
            ihcl[ihcl[,"c"] <= c_min,"h"] <- grey_hue;
            ihcl[ihcl[,"c"] <= c_min,"c"] <- 1;
         }
         w <- jamba::noiseFloor(ihcl[,"c"], minimum=0.01, ceiling=100);
         if (useWeightedHue) {
            H <- colorjam::hw2h(
               avg_angles(colorjam::h2hw(ihcl[,"h"]),
                  w=w,
                  maxAngle=360));
         } else {
            H <- avg_angles(ihcl[,"h"],
               w=w,
               maxAngle=360);
         }
         if ("mean" %in% Lmethod) {
            L <- mean((ihcl[,"l"])^k)^(1/k);
         } else if ("min" %in% Lmethod) {
            L <- min((ihcl[,"l"])^k)^(1/k);
         } else {
            L <- max(ihcl[,"l"]);
         }
         if ("mean" %in% Cmethod) {
            C <- mean((ihcl[,"c"])^k)^(1/k);
            #C <- weighted.mean((ihcl[,"c"])^k, w=sqrt(w))^(1/k);
         } else if ("min" %in% Cmethod) {
            C <- min((ihcl[,"c"])^k)^(1/k);
         } else {
            C <- max(ihcl[,"c"]);
         }
         ahcl <- as.matrix(data.frame(h=H, c=C, l=L));
         #jamba::nameVector(jamba::hcl2col(H=H,
         #   C=C,
         #   L=L),
         #   jamba::cPaste(names(i)));
         jamba::nameVector(
            farver::encode_colour(ahcl, from="hcl"),
            jamba::cPaste(names(i)));
      });
      x[xdo] <- xh;
   }
   x_new <- unlist(unname(x));
   x_extended <- x_new[x_match];
   names(x_extended) <- names(x_full);
   return(x_extended);
}

#' Find recommended overlap threshold for EnrichMap, experimental
#'
#' Find recommended overlap threshold for EnrichMap, experimental
#'
#' It implements a straightforward approach to determine
#' a reasonable Jaccard overlap threshold for Enrichment Map data,
#' and is still very much open to improvement after more
#' experience using it on varied datasets.
#' 
#' The premise is that two pathways that have Jaccard overlap above
#' a threshold are connected by a network "edge".
#' * With extremely low threshold, most pathways would be connected,
#' even if they have only one gene in common.
#' * With an extremely high threshold, pathways would only be
#' connected if nearly all genes were in common.
#' * A moderate threshold is intended to balance the two extremes.
#' * The aesthetic and biological interesting threshold appears
#' to be dependent upon the type and number of pathways returned
#' from enrichment analysis. For example, immunology pathways
#' may favor a different threshold than metabolic pathways.
#' (Purely hypothetical.)
#' * As a result, this function is intended to find a middle ground
#' based upon the pathway data used for analysis at the time,
#' where some but not all pathways are connected.
#' 
#' The method finds the overlap threshold at which the first connected
#' component is no more than `max_cutoff` fraction of the whole
#' network. This fraction is defined by the number of nodes in the
#' largest connected component, divided by the total number of
#' non-singlet nodes.
#'
#' We found that `max_cutoff=0.4`, the point at which the
#' largest connected component contains no more than 40% of all nodes,
#' seems to be a reasonably good threshold.
#'
#' @family jam utility functions
#' 
#' @returns `numeric` value with recommended Jaccard overlap coefficient.
#'
#' @param mem `list` output from `multiEnrichMap()`
#' @param overlap_range `numeric` range of Jaccard overlap values,
#'    default `0.1, 0.99` using step `0.01`.
#' @param max_cutoff `numeric` value between 0 and 1, to define the
#'    maximum fraction of nodes in the largest connected component,
#'    compared to the total number of non-singlet nodes.
#' @param adjust `numeric` used to adjust the final overlap, default
#'    `-0.01` will use the overlap one step before the max O score.
#' @param debug `logical` indicating whether to return full debug
#'    data, which is used internally to determine the best overlap
#'    cutoff to use.
#' @param ... additional arguments are passed to `mem2emap()`.
#'
#' @export
mem_find_overlap <- function
(mem,
 overlap_range=c(0.1, 0.99),
 max_cutoff=0.4,
 adjust=-0.01,
 debug=FALSE,
 ...)
{
   # define the range of values to attempt
   overlap_range <- range(overlap_range);
   oseq <- seq(from=overlap_range[1],
      to=overlap_range[2],
      by=0.01);
   odata <- lapply(jamba::nameVector(oseq), function(i){
      g <- mem2emap(mem,
         overlap=i,
         repulse=0,
      	remove_singlets=TRUE,
      	do_express=TRUE,
         do_plot=FALSE,
      	spread_labels=FALSE,
         ...)
      if (!inherits(g, "igraph")) {
      	return(NULL);
      }
      if (igraph::vcount(g) < 1) {
         return(NULL);
      }
      # Todo: consider cluster_walktrap() or something similar
      gc <- igraph::components(g);
      k <- rev(sort(gc$csize));
      names(k) <- jamba::colNum2excelName(seq_along(k))
      c(k, frac_max=unname(k[1] / sum(k)));
   });

   ## Remove empty entries
   odata <- jamba::rmNULL(odata);
   if (length(odata) == 0) {
      return(0);
   }

   ## Calculate total nodes remaining in each case
   oct <- sapply(odata, function(i){sum(head(i, -1))});
   o_fraction <- oct / max(oct);
   omax <- sapply(odata, function(i){
      unname(i["frac_max"])
   });
   o_score <- (1 - omax) * o_fraction;
   o_choose <- as.numeric(names(which.max(o_score)));
   ## quick plot to review scores
   # plot(x=as.numeric(names(o_score)), y=o_score,
   #    xlab="overlap", ylab="o_score");# debug

   # old method
   iseq <- jamba::nameVector(seq(from=max_cutoff, to=1, by=0.01));
   old_scores <- jamba::rmNULL(lapply(iseq, function(max_cutoff_i) {
   # for (max_cutoff_i in seq(from=max_cutoff, to=1, by=0.01)) {
   	omet <- sapply(odata, function(i){
   		unname(i["frac_max"] <= max_cutoff_i);
   	});
   	if (any(omet)) {
   		use_score <- as.numeric(names(head(omet[omet], 1)));
   		adjusted_score <- use_score + adjust;
   		# jamba::printDebug("max_cutoff_i: ", max_cutoff_i, " found adjusted_score: ", adjusted_score);# debug
   		# return(use_score + adjust);
   		return(adjusted_score)
   	}
   	return(NULL)
   }))
   
   if (length(debug) > 0 && debug) {
      return(list(odata=odata,
         remaining_nodes=oct,
         fraction_max=o_fraction,
         fraction_in_largest_group=omax,
         composite_score=o_score,
      	old_scores=unlist(old_scores)));
   }
   if (length(o_choose) > 0) {
      return(o_choose + adjust);
   } else {
      return(min(igraph::E(g)$overlap, na.rm=TRUE));
   }
   return(NULL)
}

#' Extract gene hit list from list of enrichResult
#'
#' Extract gene hit list from list of enrichResult
#'
#' This function is mainly for internal use in multienrichjam,
#' it takes a list of `enrichResult` objects, and determines
#' the full set of genes involved in each `enrichResult`.
#' 
#' Note that genes are sorted using `jamba::mixedSort()` for
#' alpha-numeric sorting, based upon version sorting.
#' 
#' This function also works with `ComplexHeatmap::HeatmapList`
#' objects.
#'
#' @family jam utility functions
#'
#' @returns `list` of character vectors, containing the unique
#' set of genes involved in each enrichment.
#' 
#' @param enrichList `list` of `enrichResult` objects
#' @param geneColname `character` string with the column name containing
#'    delimited gene identifiers.
#' @param geneDelim `character` regular expression used to split delimited
#'    gene values. By default, `enrichResult` uses '/' forward slash
#'    as delimiter, however the default here will split any space or
#'    comma as well.
#' @param make_unique `logical` default TRUE, whether to return only unique
#'    genes per set, or potentially multiple genes per set. Typically there
#'    should only ever be one instance of a gene per set, but through a
#'    variety of other mechanisms they may exist, for example if two gene
#'    identifiers are resolved into the same gene symbol.
#' @param verbose `logical` whether to print verbose output.
#' @param ... additional arguments are ignored.
#' 
#' @export
enrichList2geneHitList <- function
(enrichList,
 geneColname,
 geneDelim="[,/ ]",
 make_unique=TRUE,
 verbose=FALSE,
 ...)
{
   if (verbose) {
      jamba::printDebug("enrichList2geneHitList(): ",
         "geneColname: '",
         geneColname,
         "'")
      jamba::printDebug("enrichList2geneHitList(): ",
         "geneDelim: '",
         geneDelim,
         "'")
   };
   geneHitList <- lapply(enrichList, function(iDF){
      ## Split text field of delimited genes into proper vector
      if (!inherits(iDF, "data.frame")) {
         if (verbose) {
            jamba::printDebug("enrichList2geneHitList(): ",
               "head(iDF, 3):");
            print(head(iDF, 3));
         }
         k <- jamba::mixedSort(unlist(
            strsplit(
               as.character(
                  as.data.frame(iDF)[[geneColname]]),
               geneDelim)));
      } else {
         k <- jamba::mixedSort(unlist(
            strsplit(
               as.character(
                  iDF[[geneColname]]),
               geneDelim)));
      }
      if (make_unique) {
         k <- unique(k);
      }
      k;
   });
   return(geneHitList);
}


#' Order colors
#'
#' @family jam utility functions
#' 
#' @param x `character` vector of colors
#' @param sort_by `character` vector with color colnames used for sorting,
#'    passed as argument 'byCols' to `colorjam::sort_colors()`.
#' @param ... additional arguments are passed to `colorjam::sort_colors()`.
#' 
#' @examples
#' set.seed(123);
#' x <- sample(grep("royal|^golden", colors(), value=TRUE))
#' x_sorted <- x[order_colors(x)];
#' jamba::showColors(list(input=x, sorted=x_sorted))
#' 
#' @export
order_colors <- function
(x,
 sort_by=c("H", "-L", "-C"),
 ...)
{
   ## hue_offset=-12 makes red the first color,
   ## and moves pink to the end with purple
   x_sorted <- colorjam::sort_colors(unique(x),
      byCols=sort_by,
      ...);
   x_factor <- factor(x, levels=x_sorted);
   order(x_factor);
}


#' Rank Multienrichment clusters
#'
#' Rank Multienrichment clusters
#'
#' This function takes `list` output from `multiEnrichMap()`, and
#' a `list` of clusters, and returns a `data.frame` that contains
#' several rank order metrics. It is intended to be used with
#' column clusters following `mem_gene_path_heatmap()`,
#' see examples.
#'
#' The argument `per_cluster` is intended to make it convenient
#' to pick the top exemplar pathways, especially when argument
#' `byCols` is defined so that it sorts by the rank columns.
#'
#' The argument `choose` is intended to make it easy to retrieve
#' pathways from specific clusters.
#'
#' @return `data.frame` sorted by the criteria defined by `byCols`,
#'    with colname `"set"` to indicate the pathway/set name.
#'
#' @family jam utility functions
#'
#' @param mem `Mem` or legacy `list` mem output from `multiEnrichMap()`
#' @param clusters `list` containing set names, that must match
#'    `colnames(mem$memIM)` and `rownames(mem$enrichIM)`.
#' @param choose optional vector indicating which clusters to
#'    return. If an integer vector, it refers to the elements
#'    in the input `clusters`. If a character vector, it must
#'    contain values in `names(clusters)`. When `choose` is `NULL`,
#'    all clusters are returned.
#' @param per_cluster integer vector with the number of entries
#'    to return for each cluster. Values will be recycled to the
#'    length of the clusters to be returned, defined by `choose`
#'    or by `length(clusters)` when `choose` is `NULL`.
#' @param byCols character vector used to sort the resulting
#'    `data.frame` within each cluster. This argument is passed
#'    directly to `jamba::mixedSortDF()`.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' ## Start with mem
#' # mem <- multiEnrichMap(...);
#' # gp_hm <- mem_gene_path_heatmap(mem, column_split=4);
#' ## Retrieve clusters from the Heatmap output, there should be 4 clusters
#' # clusters <- heatmap_column_order(gp_hm)
#' # clusters_df <- rank_mem_clusters(mem, clusters)
#'
#' @export
rank_mem_clusters <- function
(mem,
 clusters,
 choose=NULL,
 per_cluster=Inf,
 byCols=c("composite_rank", "minp_rank", "gene_count_rank"),
 verbose=FALSE,
 ...)
{
   Mem <- NULL;
   if (inherits(mem, "Mem")) {
      Mem <- mem;
      mem <- Mem_to_list(Mem);
   } else {
      Mem <- list_to_Mem(mem);
   }
   
   if (!all(c("memIM") %in% names(mem))) {
      stop("mem must be a list with components: 'memIM'");
   }
   if (length(names(clusters)) == 0) {
      if (verbose) {
         jamba::printDebug("rank_mem_clusters(): ",
            "Assigning names(clusters) using ",
            "jamba::colNum2excelName");
      }
      names(clusters) <- jamba::colNum2excelName(seq_along(clusters));
   }
   if (length(choose) == 0) {
      choose <- seq_along(clusters);
   }
   clusters <- jamba::rmNULL(clusters[choose]);
   if (length(names(clusters)) == 0) {
      stop("No clusters remained after subsetting with argument 'choose'.");
   }
   if (length(per_cluster) == 0) {
      per_cluster <- Inf;
   }
   per_cluster <- rep(per_cluster, length.out=length(clusters))
   names(per_cluster) <- names(clusters);
   clusters_dfs <- lapply(names(clusters), function(iname){
      i <- clusters[[iname]];
      j <- mem$memIM[,i,drop=FALSE];
      pval_m <- mem$enrichIM[i,,drop=FALSE];
      colnames(pval_m) <- paste0("P_", colnames(pval_m));
      minp <- jamba::rmNA(naValue=1,
         apply(mem$enrichIM[i,,drop=FALSE], 1, min, na.rm=TRUE));
      names(minp) <- i;
      gene_count <- colSums(j);
      minp_rank <- order(do.call(order, list(minp, -gene_count)));
      gene_count_rank <- order(do.call(order, list(-gene_count, minp)));
      composite_rank <- order(do.call(order,
         list(
            rank(floor(log10(minp))),
            gene_count_rank)));
      ijdf <- data.frame(
         stringsAsFactors=FALSE,
         check.names=FALSE,
         cluster=iname,
         set=i,
         gene_count=gene_count,
         pval_m,
         minp=minp,
         gene_count_rank=gene_count_rank,
         minp_rank=minp_rank,
         composite_rank=composite_rank);
      if (length(byCols) > 0) {
         ijdf <- jamba::mixedSortDF(ijdf,
            byCols=byCols);
      }
      ijdf$cluster_rank <- paste0(iname, "_", seq_len(nrow(ijdf)));
      rownames(ijdf) <- ijdf$cluster_rank;
      head(ijdf, per_cluster[[iname]]);
   });
   clusters_df <- jamba::rbindList(clusters_dfs);
   return(clusters_df);
}


#' Collapse Multienrichment clusters
#'
#' Collapse Multienrichment clusters
#'
#' This function is similar to `rank_mem_clusters()` in that it
#' starts with `mem` results from `multiEnrichMap()` and a list
#' of `clusters` of pathways/sets. Instead of ranking and choosing
#' exemplars from each clusters, it simply collapses each cluster
#' into one super-set that contains union of all genes.
#'
#' It does also run `rank_mem_clusters()` in the event one would
#' want to collapse only the top `per_cluster` entries for each
#' cluster, but the default is to include them all.
#'
#' @returns By default `return_type="cnet"` and this function returns
#'    an augmented Cnet plot. The labels of each cluster are defined
#'    by the input `names(clusters)`, however an `igraph` attribute
#'    `"set_labels"` includes an abbreviated label of the top ranked
#'    sets for each cluster. This label is probably the closest thing
#'    to summarizing the composition of each cluster.
#'
#' @family jam utility functions
#' 
#' @param mem `Mem` or legacy `list` mem object
#' @param clusters `list` named by cluster name, containing `character`
#'    vectors all of whose values should be present in `sets(mem)`.
#' @param choose `character` default NULL, passed to `rank_mem_clusters()`
#'    to choose a subset of clusters to return. Values must be found in
#'    `names(clusters)`.
#' @param per_cluster `integer` default 'Inf' with the maximum sets
#'    to include in each cluster. Use an integer like `per_cluster=1`
#'    to choose one "exemplar" per cluster, ranked according
#'    to `byCols` which is described in `rank_mem_clusters()`.
#' @param byCols `character` vector of summary values to use when
#'    sorting sets within each cluster, passed to and described
#'    in `rank_mem_clusters()`.
#' @param return_type `character` default "cnet" with the output format:
#'    * "cnet" returns an `igraph` Cnet plot.
#'    * "Mem" returns a `Mem` object after collapsing by cluster.
#'    * "mem" returns legacy `list` mem output.
#' @param max_labels `integer` default 4, the number of pathway labels
#'    to include in each cluster, ordered according to 'byCols'.
#' @param max_nchar_labels `integer` default 25, max character length
#'    for each pathway label.
#' @param include_cluster_title `logical` default TRUE, whether to include
#'    the cluster title as a prefix.
#'    For example if the cluster name is 'A' the new title would be:
#'    'A: pathway 1; pathway 2; pathway 3; pathway 4'.
#' @param cluster_color_min_fraction `numeric` value between 0 and 1
#'    used to define the cnet colors for each cluster. The number of
#'    significant pathways is calculated per cluster for each enrichment
#'    using `mem$enrichIMcolors`, and the fraction of significant
#'    pathways versus the max number per enrichment is used to filter.
#'    For example a cluster with 10 pathways, might have 8 significant
#'    pathways in one enrichment result, and 3 significant pathways in
#'    another enrichment result. The fraction for each enrichment is
#'    `8/8 == 1` for the first enrichment, and `3/8 = 0.375` for the
#'    second enrichment. When `cluster_color_min_fraction=0.5` (default)
#'    the first enrichment color would be included, but not the second
#'    enrichment color. The intent is to represent enrichment colors
#'    that have at least half (`0.5`) the pathways of the most
#'    representative (max) enrichment color. Therefore, a cluster with
#'    only one significant pathway from a given enrichment would
#'    typically not be representive of that enrichment, and its
#'    enrichment color would not be included.
#' @param cluster_color_min_fraction `numeric` value default 0.5,
#'    the fraction of pathways in a cluster which has "non blank"
#'    assigned color in `enrichIMcolors()` in order for the color
#'    to be assigned to the resulting node, in the case of Cnet
#'    output with `return_type="cnet"`.
#' @param verbose `logical` whether to print verbose output.
#' @param ... additional arguments are passed to `cnet2mem()` which
#'    includes thresholds used by `apply_cnet_direction()` such as
#'    `direction_col_fn` (color function applied to enrichIMdirection),
#'    `gene_direction_col_fn` (color function applied to geneIMdirection),
#'    and others.
#'    To avoid applying directional color to Set nodes (nodeType='Set')
#'    use `direction_cutoff=100` and `direction_max=101`, which will
#'    consider any enrichIMdirection value at or below 100 to be
#'    non-directional.
#' 
#' @export
collapse_mem_clusters <- function
(mem,
 clusters,
 choose=NULL,
 per_cluster=Inf,
 byCols=c("composite_rank", "minp_rank", "gene_count_rank"),
 return_type=c("cnet", "Mem", "mem"),
 max_labels=4,
 max_nchar_labels=25,
 include_cluster_title=TRUE,
 cluster_color_min_fraction=0.5,
 verbose=FALSE,
 ...)
{
   #
   return_type <- match.arg(return_type);
   ##
   if (inherits(mem, "Mem")) {
      Mem <- mem;
      mem <- Mem_to_list(Mem);
   }
   clusters_df <- rank_mem_clusters(mem=mem,
      clusters=clusters,
      choose=choose,
      per_cluster=per_cluster,
      byCols=byCols,
      verbose=verbose);
   cluster_sets_l <- split(clusters_df$set,
      clusters_df$cluster);

   max_labels <- rep(max_labels,
      length.out=length(cluster_sets_l));
   names(max_labels) <- names(cluster_sets_l);
   max_nchar_labels <- rep(max_nchar_labels,
      length.out=length(cluster_sets_l));
   names(max_nchar_labels) <- names(cluster_sets_l);

   cluster_sets <- jamba::cPaste(cluster_sets_l, sep="; ");

   cluster_labels <- jamba::cPaste(sep=";\n",
      lapply(jamba::nameVectorN(cluster_sets_l), function(j){
         i <- cluster_sets_l[[j]];
         if (length(i) > max_labels[j]) {
            more <- paste0("(", length(i) - max_labels[j], " more)");
         } else {
            more <- NULL;
         }
         c(head(ifelse(nchar(i) <= max_nchar_labels[j],
            i,
            paste0(substr(i, 1, max_nchar_labels[j] - 3), "...")),
            max_labels[j]), more);
      }));
   ## Optionally prepend the cluster title to the cluster_labels
   if (include_cluster_title) {
      cluster_labels <- paste0(
         "Cluster ",
         names(cluster_labels),
         ":\n",
         cluster_labels);
   }
   #for (j in cluster_labels){cat("\n");cat(j, "\n\n");}

   clusters_df$cluster <- factor(clusters_df$cluster,
      levels=unique(clusters_df$cluster));
   cluster_memIM <- do.call(cbind,
      lapply(split(clusters_df, clusters_df$cluster), function(idf){
         # im1 <- subset(mem$memIM[, idf$set, drop=FALSE]);
         im1 <- mem$memIM[, idf$set, drop=FALSE];
         # iv <- rowMaxs(im1);
         iv <- apply(im1, 1, max, na.rm=TRUE);
         names(iv) <- rownames(im1);
         iv;
      }));

   cluster_enrichIMgeneCount <- do.call(cbind,
      lapply(jamba::nameVector(colnames(mem$geneIM)), function(i){
         xgenes <- intersect(rownames(mem$geneIM)[mem$geneIM[,i] > 0],
            rownames(cluster_memIM));
         if (verbose) {
            jamba::printDebug("length(xgenes):", length(xgenes));
            jamba::printDebug("head(xgenes, 30):");
            print(head(xgenes, 30));
         }
         colSums(cluster_memIM[xgenes,,drop=FALSE] > 0);
      }));
   
   ## aggregate enrichment P-values
   ## - take mean of log10 P-values, then exponentiate
   cluster_enrichIM <- jamba::rbindList(
      lapply(split(clusters_df, clusters_df$cluster), function(idf){
         # im1 <- subset(mem$enrichIM[idf$set,,drop=FALSE]);
         im1 <- mem$enrichIM[idf$set, , drop=FALSE];
         10^colMeans(log10(im1))
      }));

   ## aggregate enrichment directionality
   ## - take weighted mean?
   ##   weighted_mean <- sum(z * abs(z)) / sum(abs(z))
   ## - directional coherence (similar to concordance) but with thresholding
   ##
   cluster_enrichIMdirection <- jamba::rbindList(
      lapply(split(clusters_df, clusters_df$cluster), function(idf){
         # im1 <- subset(mem$enrichIMdirection[idf$set, , drop=FALSE]);
         im1 <- mem$enrichIMdirection[idf$set, , drop=FALSE];
         
         ## directional coherence with thresholding (above 1)
         ## - then scale by 4 to make it approximately similar to z-score
         # abs(colSums(sign(im1 * (abs(im1) >= 1)))) / nrow(im1) * 4
         
         ## weighted mean
         jamba::rmNA(naValue=0, colSums(im1 * abs(im1)) / colSums(abs(im1)))
      }));
   
   ## Method to determine grouped colors per cluster
   ## - start with enrichIMcolors
   ## - count non-blank colors per column
   ## - determine fraction versus max column count
   cluster_enrichIMcolors <- jamba::rbindList(lapply(cluster_sets_l, function(iset1){
      isetm1 <- mem$enrichIMcolors[iset1, , drop=FALSE]
      # fix issue when only one row or one column
      isetm1_blank <- matrix(apply(isetm1, 2, isColorBlank), ncol=ncol(isetm1));
      color_counts <- colSums(!isetm1_blank);
      names(color_counts) <- colnames(isetm1);
      sapply(colnames(isetm1), function(k1){
         if (max(color_counts) == 0) {
            # should not occur, in theory
            return("#FFFFFFFF");
         }
         if ((color_counts[k1] / max(color_counts)) >= cluster_color_min_fraction) {
            unname(mem$colorV[k1])
         } else {
            "#FFFFFFFF";
         }
      })
   }))
   #cluster_enrichIMcolors <- do.call(cbind, lapply(jamba::nameVector(colnames(mem$enrichIMcolors)), function(i){
   #   avg_colors_by_list(split(mem$enrichIMcolors[clusters_df$set,i], clusters_df$cluster))
   #}))
   ## avg_colors_by_list
   cluster_mem <- mem;
   cluster_mem$memIM <- cluster_memIM;
   cluster_mem$enrichIM <- cluster_enrichIM;
   
   # Todo: make this directional?
   cluster_mem$enrichIMdirection <- cluster_enrichIMdirection;
   
   cluster_mem$enrichIMcolors <- cluster_enrichIMcolors;
   cluster_mem$enrichIMgeneCount <- cluster_enrichIMgeneCount;

   ## 0.0.101.900 - do not edit nor remove these objects bc
   ## they are difficult to recreate in this context but
   ## are required in the Mem object
   # cluster_mem$multiEnrichDF <- data.frame()
   # cluster_mem$multiEnrichResult <- list();

   ## Make Cnet plot
   if (verbose) {
      jamba::printDebug("collapse_mem_clusters(): ",
         "Calling mem2cnet()");
   }
   cnet <- mem2cnet(cluster_mem,
      verbose=verbose,
      ...);
   
   set_match <- match(names(cluster_sets),
      igraph::V(cnet)$name)
   igraph::V(cnet)$set_names <- "";
   igraph::V(cnet)$set_names[set_match] <- cluster_sets;
   igraph::V(cnet)$set_labels <- "";
   igraph::V(cnet)$set_labels[set_match] <- cluster_labels;
   cluster_mem$multiCnetPlot <- cnet;

   if ("cnet" %in% return_type) {
      return(cnet);
   }
   if ("Mem" %in% return_type) {
      return(list_to_Mem(cluster_mem));
   }
   return(cluster_mem);
}
