
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
   xComp <- mean(cos(x * 2*pi / maxAngle),
      na.rm=na.rm,
      ...);
   yComp <- mean(sin(x * 2*pi / maxAngle),
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
 ...)
{
   if (length(x) == 0) {
      return(x);
   }
   Cmethod <- match.arg(Cmethod);
   Lmethod <- match.arg(Lmethod);
   xdo <- which(lengths(x) > 1);
   k <- 1;
   if (length(xdo) > 0) {
      xh <- lapply(unname(x[xdo]), function(i){
         ihcl <- jamba::col2hcl(i);
         if (useWeightedHue) {
            H <- colorjam::hw2h(
               avg_angles(colorjam::h2hw(ihcl["H",]),
                  maxAngle=360));
         } else {
            H <- avg_angles(ihcl["H",],
                  maxAngle=360);
         }
         if ("mean" %in% Lmethod) {
            L <- mean((ihcl["L",])^k)^(1/k);
         } else if ("min" %in% Lmethod) {
            L <- min((ihcl["L",])^k)^(1/k);
         } else {
            L <- max(ihcl["L",]);
         }
         if ("mean" %in% Cmethod) {
            C <- mean((ihcl["C",])^k)^(1/k);
         } else if ("min" %in% Cmethod) {
            C <- min((ihcl["C",])^k)^(1/k);
         } else {
            C <- max(ihcl["C",]);
         }
         jamba::nameVector(jamba::hcl2col(H=H,
            C=C,
            L=L),
            jamba::cPaste(names(i)));
      });
      x[xdo] <- xh;
   }
   x_new <- unlist(unname(x));
   return(x_new);
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
 max_cutoff=0.4,
 debug=FALSE,
 ...)
{
   overlap_range <- range(overlap_range);
   oseq <- seq(from=overlap_range[1],
      to=overlap_range[2],
      by=0.01);
   odata <- lapply(nameVector(oseq), function(i){
      g <- mem_multienrichplot(mem,
         overlap=i,
         remove_blanks=FALSE,
         remove_singlets=TRUE,
         spread_labels=FALSE,
         do_plot=FALSE);
      if (vcount(g) <= 1) {
         return(NULL);
      }
      gc <- components(g);
      k <- rev(sort(gc$csize));
      c(k, frac_max=k[1] / sum(k));
   });
   if (length(debug) > 0 && debug) {
      return(odata);
   }
   ## Remove empty entries
   odata <- jamba::rmNULL(odata);
   if (length(odata) == 0) {
      return(0);
   }
   omet <- sapply(odata, function(i){
      unname(i["frac_max"] <= max_cutoff);
   });
   return(as.numeric(names(head(omet[omet], 1))));
}

#' Extract gene hit list from list of enrichResult
#'
#' Extract gene hit list from list of enrichResult
#'
#' This function is mainly for internal use in multienrichjam,
#' it takes a list of `enrichResult` objects, and determines
#' the full set of genes involved in each `enrichResult`.
#'
#' @return `list` of character vectors, containing the unique
#' set of genes involved in each enrichment.
#'
#' @export
enrichList2geneHitList <- function
(enrichList,
 geneColname,
 geneDelim,
 make_unique=TRUE,
 ...)
{
   geneHitList <- lapply(enrichList, function(iDF){
      ## Split text field of delimited genes into proper vector
      if (!jamba::igrepHas("data.frame", class(iDF))) {
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
#' @jam utility functions
#'
#' @export
heatmap_row_order <- function
(hm)
{
   ##
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
#' @jam utility functions
#'
#' @export
heatmap_column_order <- function
(hm)
{
   ##
   x <- lapply(ComplexHeatmap::column_order(hm), function(i){
      hm@column_names_param$labels[i]
   });
   if (length(hm@column_title) == length(x)) {
      names(x) <- hm@column_title;
   }
   x;
}
