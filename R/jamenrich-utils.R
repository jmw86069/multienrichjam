
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
      return(list(odata=odata, oct=oct, omax=omax, o_fraction=o_fraction, o_score=o_score));
   }
   if (length(o_choose) > 0) {
      return(o_choose);
   } else {
      return(min(E(g)$overlap, na.rm=TRUE));
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
#' @family jam utility functions
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
#' @family jam utility functions
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
#' @family jam utility functions
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

#' Sort colors
#'
#' @family jam utility functions
#'
#' @examples
#' x <- jamba::nameVector(colors());
#'
#' ## Basic color sort
#' c2 <- sort_colors(x);
#' jamba::showColors(c2, main="cmin=4");
#'
#' ## Increase filtering of unsaturated colors
#' c3 <- sort_colors(x, c_min=20);
#' jamba::showColors(c3, main="cmin=20");
#'
#' ## Increase filtering of unsaturated colors
#' c4 <- sort_colors(x, c_min=50);
#' jamba::showColors(c4, main="cmin=50");
#'
#' @export
sort_colors <- function
(x,
 sort_by=c("h", "-c", "-l"),
 c_min=4,
 grey_hue=359,
 hue_offset=0,
 ...)
{
   ## hue_offset=-12 makes red the first color,
   ## and moves pink to the end with purple
   x_sort <- data.frame(input=x);
   if (any(c("h","c","l") %in% sort_by)) {
      x_hcl <- farver::decode_colour(x, to="hcl");
      x_hcl[,"h"] <- (x_hcl[,"h"] + hue_offset) %% 360;
      if (any(x_hcl[,"c"] <= c_min)) {
         c_is_min <- which(x_hcl[,"c"] <= c_min);
         x_hcl[c_is_min,"h"] <- grey_hue;
         x_hcl[c_is_min,"c"] <- c_min;
      }
      x_hcl <- round(x_hcl/4)*4;
      x_sort <- cbind(x_sort, x_hcl);
   }
   if (any(c("s","v") %in% sort_by)) {
      x_sv <- 100*farver::decode_colour(x, to="hsv")[,c("s","v"),drop=FALSE];
      x_sort <- cbind(x_sort, x_sv);
   }
   if (any(c("a","b") %in% sort_by)) {
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
 sort_by=c("h", "-c", "-l"),
 c_min=4,
 grey_hue=359,
 hue_offset=0,
 ...)
{
   ## hue_offset=-12 makes red the first color,
   ## and moves pink to the end with purple
   x_sorted <- sort_colors(unique(x),
      sort_by=sort_by,
      c_min=c_min,
      grey_hue=grey_hue,
      hue_offset=hue_offset,
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
 range=c(0,100),
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
