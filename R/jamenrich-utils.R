
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

#' Find best overlap threshold for EnrichMap
#'
#' Find best overlap threshold for EnrichMap
#'
#' This function implements a straightforward approach to determine
#' a reasonable Jaccard overlap threshold for EnrichMap data.
#' It finds the overlap threshold at which the first connected
#' component is no more than `max_cutoff` fraction of the whole
#' network. This fraction is defined as the number of nodes in the
#' largest connected component, divided by the total number of
#' non-singlet nodes. When all nodes are connected, this fraction == 1.
#'
#' We found empirically that a `max_cutoff=0.4`, the point at which the
#' largest connected component contains no more than 40% of all nodes,
#' seems to be a reasonably good place to start.
#'
#' @family jam utility functions
#'
#' @param mem `list` output from `multiEnrichMap()`
#' @param overlap_range numeric range of Jaccard overlap values
#' @param overlap_count numeric value passed to `mem_multienrichplot()`
#'    which is used to filter the multienrichmap by Jaccard overlap
#'    and by overlap_count.
#' @param max_cutoff numeric value between 0 and 1, to define the
#'    maximum fraction of nodes in the largest connected component,
#'    compared to the total number of non-singlet nodes.
#' @param debug logical indicating whether to return full debug
#'    data, which is used internally to determine the best overlap
#'    cutoff to use.
#'
#' @export
mem_find_overlap <- function
(mem,
 overlap_range=c(0.1,0.99),
 overlap_count=2,
 node_fraction=0.5,
 max_cutoff=0.4,
 debug=FALSE,
 ...)
{
   overlap_range <- range(overlap_range);
   oseq <- seq(from=overlap_range[1],
      to=overlap_range[2],
      by=0.01);
   odata <- lapply(jamba::nameVector(oseq), function(i){
      g <- mem_multienrichplot(mem,
         overlap=i,
         remove_blanks=FALSE,
         remove_singlets=TRUE,
         spread_labels=FALSE,
         do_plot=FALSE);
      if (igraph::vcount(g) <= 1) {
         return(NULL);
      }
      gc <- igraph::components(g);
      k <- rev(sort(gc$csize));
      c(k, frac_max=k[1] / sum(k));
   });

   ## Remove empty entries
   odata <- jamba::rmNULL(odata);
   if (length(odata) == 0) {
      return(0);
   }

   ## Calculate total nodes remaining in each case
   oct <- sapply(odata, function(i){sum(head(i, -1))});
   o_fraction <- oct / max(oct);
   omax <- sapply(odata, function(i){unname(i["frac_max"])});
   o_score <- (1 - omax) * o_fraction;
   o_choose <- as.numeric(names(which.max(o_score)));

   if (length(debug) > 0 && debug) {
      return(list(odata=odata,
         oct=oct,
         omax=omax,
         o_fraction=o_fraction,
         o_score=o_score));
   }
   if (length(o_choose) > 0) {
      return(o_choose);
   } else {
      return(min(igraph::E(g)$overlap, na.rm=TRUE));
   }
   for (max_cutoff_i in seq(from=max_cutoff, to=1, by=0.01)) {
      omet <- sapply(odata, function(i){
         unname(i["frac_max"] <= max_cutoff_i);
      });
      if (any(omet)) {
         return(as.numeric(names(head(omet[omet], 1))));
      }
   }
}

#' Extract gene hit list from list of enrichResult
#'
#' Extract gene hit list from list of enrichResult
#'
#' This function is mainly for internal use in multienrichjam,
#' it takes a list of `enrichResult` objects, and determines
#' the full set of genes involved in each `enrichResult`.
#'
#' This function also works with `ComplexHeatmap::HeatmapList`
#' objects.
#'
#' @family jam utility functions
#'
#' @return `list` of character vectors, containing the unique
#' set of genes involved in each enrichment.
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
      if (!jamba::igrepHas("data.frame", class(iDF))) {
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

#' Return Heatmap row order
#'
#' Return Heatmap row order as a list of character vectors
#'
#' This function is a helpful utility to return the fully
#' qualified list of rownames in a `ComplexHeatmap::Heatmap`
#' object.
#'
#' This function also works with `ComplexHeatmap::HeatmapList`
#' objects.
#'
#' @family jam utility functions
#'
#' @export
heatmap_row_order <- function
(hm)
{
   ##
   if ("HeatmapList" %in% class(hm)) {
      hm <- hm@ht_list[[1]];
   }
   x <- lapply(ComplexHeatmap::row_order(hm), function(i){
      hm@row_names_param$labels[i]
   });
   if (length(hm@row_title) == length(x)) {
      names(x) <- hm@row_title;
   }
   x;
}

#' Return Heatmap column order
#'
#' Return Heatmap column order as a list of character vectors
#'
#' This function is a helpful utility to return the fully
#' qualified list of colnames in a `ComplexHeatmap::Heatmap`
#' object.
#'
#' @family jam utility functions
#'
#' @export
heatmap_column_order <- function
(hm)
{
   ##
   if ("HeatmapList" %in% class(hm)) {
      hm <- hm@ht_list[[1]];
   }
   x <- lapply(ComplexHeatmap::column_order(hm), function(i){
      hm@column_names_param$labels[i]
   });
   if (length(hm@column_title) == length(x)) {
      names(x) <- hm@column_title;
   }
   x;
}

#' Sort colors (deprecated)
#'
#' Sort colors (deprecated)
#'
#' This function is deprecated, please use `colorjam::sort_colors()`.
#'
#' This function is intended to be a very rapid method to sort
#' colors, based upon hue, then chroma descending, then luminance
#' descending.
#'
#' @family jam utility functions
#'
#' @examples
#' x <- jamba::nameVector(colors());
#'
#' ## Basic color sort
#' c2 <- sort_colors_deprecated(x);
#' jamba::showColors(c2, main="cmin=4");
#'
#' ## Increase filtering of unsaturated colors
#' c3 <- sort_colors_deprecated(x, c_min=20);
#' jamba::showColors(c3, main="cmin=20");
#'
#' ## Increase filtering of unsaturated colors
#' c4 <- sort_colors_deprecated(x, c_min=50);
#' jamba::showColors(c4, main="cmin=50");
#'
#' @export
sort_colors_deprecated <- function
(x,
 sort_by=c("h", "-c", "-l"),
 c_min=4,
 grey_hue=359,
 hue_offset=0,
 ...)
{
   ## hue_offset=-12 makes red the first color,
   ## and moves pink to the end with purple
   if (length(x) == 0) {
      return(x);
   }
   x_sort <- data.frame(input=x);
   if (any(c("h","c","l","-c", "-l") %in% sort_by)) {
      x_hcl <- farver::decode_colour(x, to="hcl");
      x_hcl[,"h"] <- jamba::rmNA(naValue=0, x_hcl[,"h"]);
      x_hcl[,"c"] <- jamba::rmNA(naValue=c_min, x_hcl[,"c"]);
      x_hcl[,"l"] <- jamba::rmNA(naValue=0, x_hcl[,"l"]);
      x_hcl[,"h"] <- (x_hcl[,"h"] + hue_offset) %% 360;
      if (any(x_hcl[,"c"] <= c_min)) {
         c_is_min <- which(x_hcl[,"c"] <= c_min);
         x_hcl[c_is_min,"h"] <- grey_hue;
         x_hcl[c_is_min,"c"] <- c_min;
      }
      x_hcl <- round(x_hcl/4)*4;
      x_sort <- cbind(x_sort, x_hcl);
   }
   if (any(c("s","v","-s","-v") %in% sort_by)) {
      x_hcl[,"s"] <- jamba::rmNA(naValue=0, x_hcl[,"s"]);
      x_hcl[,"v"] <- jamba::rmNA(naValue=0, x_hcl[,"v"]);
      x_sv <- 100*farver::decode_colour(x, to="hsv")[,c("s","v"),drop=FALSE];
      x_sort <- cbind(x_sort, x_sv);
   }
   if (any(c("a","b","-a","-b") %in% sort_by)) {
      x_hcl[,"a"] <- jamba::rmNA(naValue=0, x_hcl[,"a"]);
      x_hcl[,"b"] <- jamba::rmNA(naValue=0, x_hcl[,"b"]);
      x_ab <- round(farver::decode_colour(x, to="lab")[,c("a","b"),drop=FALSE]);
      x_sort <- cbind(x_sort, x_ab);
   }
   x_sorted <- jamba::mixedSortDF(data.frame(x_sort),
      byCols=sort_by);
   #x[x_sorted$i];
   x_sorted$input;
}

#' Order colors
#'
#' @family jam utility functions
#'
#' @export
order_colors <- function
(x,
 sort_by=c("H", "-C", "-L"),
 c_min=4,
 grey_hue=359,
 hue_offset=0,
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

#' Apply color channel numeric range cap
#'
#' Apply color channel numeric range cap
#'
#' This function extends `farver::set_channel()` by enforcing
#' a numeric range for any numeric channel accessible by
#' the `farver` package. Common use is to restrict the luminance
#' `"l"` channel of an `"hcl"` color from `c(0,100)` to `c(70,100)`
#' in order to make all colors brighter.
#'
#' @family jam utility functions
#'
#' @param x vector colors, or list of color vectors. When input
#'    is a `list`, the list is flattened, operations are performed,
#'    then the result is split back into the original structure.
#' @param channel character channel recognized by `farver::set_channel()`.
#' @param range numeric range allowed for values returned by
#'    `farver::get_channel()`.
#' @param space character name of a color space recognized by
#'    `farver::set_channel()`.
#' @param ... additional arguments are passed to `farver::get_channel()`
#'    and `farver::set_channel()`.
#'
#' @export
apply_color_cap <- function
(x,
 channel="l",
 range=c(0, 100),
 space="hcl",
 ...)
{
   if (is.list(x)) {
      x_ext <- unlist(x);
      x_len <- lengths(x);
      x_fac <- rep(factor(seq_along(x)), x_len);
      x_new <- apply_color_cap(x_ext,
         channel=channel,
         range=range,
         space=space,
         ...);
      x_list <- split(x_new, x_fac);
      names(x_list) <- names(x);
      return(x_list);
   }
   channel_range <- range(range, na.rm=TRUE);
   channel_values <- farver::get_channel(x,
      channel=channel,
      space=space,
      ...);
   channel_values <- jamba::noiseFloor(channel_values,
      minimum=min(range),
      ceiling=max(range));
   x_new <- farver::set_channel(x,
      channel=channel,
      space=space,
      value=channel_values,
      ...);
   names(x_new) <- names(x);
   return(x_new);
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
#' @param mem `list` object output from `multiEnrichMap()`
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
      ijdf <- data.frame(cluster=iname,
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
#' @return By default `return_type="cnet"` and this function returns
#'    an augmented Cnet plot. The labels of each cluster are defined
#'    by the input `names(clusters)`, however an `igraph` attribute
#'    `"set_labels"` includes an abbreviated label of the top ranked
#'    sets for each cluster. This label is probably the closest thing
#'    to summarizing the composition of each cluster.
#'
#' @family jam utility functions
#'
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
#'
#' @export
collapse_mem_clusters <- function
(mem,
 clusters,
 choose=NULL,
 per_cluster=Inf,
 byCols=c("composite_rank", "minp_rank", "gene_count_rank"),
 return_type=c("cnet", "mem"),
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
   cluster_memIM <- do.call(cbind, lapply(split(clusters_df, clusters_df$cluster), function(idf){
      im1 <- subset(mem$memIM[,idf$set,drop=FALSE]);
      iv <- rowMaxs(im1);
      names(iv) <- rownames(im1);
      iv;
   }));

   cluster_enrichIMgeneCount <- do.call(cbind, lapply(jamba::nameVector(colnames(mem$geneIM)), function(i){
      xgenes <- intersect(rownames(mem$geneIM)[mem$geneIM[,i] > 0],
         rownames(cluster_memIM));
      if (verbose) {
         jamba::printDebug("length(xgenes):", length(xgenes));
         jamba::printDebug("head(xgenes, 30):");
         print(head(xgenes, 30));
      }
      colSums(cluster_memIM[xgenes,,drop=FALSE] > 0);
   }));
   cluster_enrichIM <- jamba::rbindList(lapply(split(clusters_df, clusters_df$cluster), function(idf){
      im1 <- subset(mem$enrichIM[idf$set,,drop=FALSE]);
      10^colMeans(log10(im1))
   }));

   ## Method to determine grouped colors per cluster
   ## - start with enrichIMcolors
   ## - count non-blank colors per column
   ## - determine fraction versus max column count
   cluster_enrichIMcolors <- jamba::rbindList(lapply(cluster_sets_l, function(iset1){
      isetm1 <- mem$enrichIMcolors[iset1,,drop=FALSE]
      # fix issue when only one row or one column
      isetm1_blank <- matrix(apply(isetm1, 2, isColorBlank), ncol=ncol(isetm1));
      color_counts <- colSums(!isetm1_blank);
      names(color_counts) <- colnames(isetm1);
      sapply(colnames(isetm1), function(k1){
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
   cluster_mem$enrichIMcolors <- cluster_enrichIMcolors;
   cluster_mem$enrichIMgeneCount <- cluster_enrichIMgeneCount;

   ## For now remove some difficult-to-update objects
   cluster_mem$multiEnrichDF <- NULL;
   cluster_mem$multiEnrichResult <- NULL;
   cluster_mem$multiEnrichMap <- NULL;
   cluster_mem$multiEnrichMap2 <- NULL;
   cluster_mem$multiCnetPlot <- NULL;
   cluster_mem$multiCnetPlot1 <- NULL;
   cluster_mem$multiCnetPlot1b <- NULL;
   cluster_mem$multiCnetPlot2 <- NULL;

   ## Make Cnet plot
   if (verbose) {
      jamba::printDebug("collapse_mem_clusters(): ",
         "Calling memIM2cnet()");
   }
   cnet <- memIM2cnet(cluster_mem,
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
   return(cluster_mem);
}
