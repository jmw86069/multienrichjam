
#' Define gene score color function with fixed midpoint
#' 
#' Define gene score color function with fixed midpoint and
#' scalable maximum value
#' 
#' The gene score color function tries to accomplish two goals:
#' 1. Color gene hits with values `-1` and `1`, blue and red,
#' respectively.
#' 2. Color gene scores, with values that may range from `-8`
#' to `8`, with a color gradient to indicate the magnitude of
#' the score.
#' 
#' This color function is applied to ORA enrichment data,
#' using `-1` and `1` to indicate 'down' and 'up', respectively.
#' This color function is also applied to GSEA enrichment data,
#' using gene scores that are typically derived from t-statistics
#' or signed-P-values: `-log10(P-value) * sign(logFC)`.
#' 
#' The argument `x` defines the range beyond `1` for colors
#' available to the new range '1 to `ceiling`'. It therefore
#' also defines the color at `1`.
#' 
#' The argument `ceiling` defines how far to stretch the colors
#' originally used for '1 to `x`', so that they span '1 to `ceiling`'.
#' 
#' In general, higher `x` should consider using higher `lens`
#' so the colors at `-1` and `1` are still visibly blue and
#' red, respectively.
#' 
#' ## Todo
#' 
#' * Consider how to handle direction values being `-1` and `1`
#' while the signed-P used for GSEA may range much higher,
#' and therefore `-1` and `1` for signed-P scores might be
#' presented as too strong compared to the direction hits.
#' For example, maybe the signed-P scores should be scaled
#' somehow?
#' 
#' @param x `numeric` step 1 maximum value, must be greater than 1.
#' @param ceiling `numeric` step 2 maximum value
#' @param lens `numeric` passed to `colorjam::col_div_xf()`
#' @param colramp `character` string passed to `colorjam::col_div_xf()`
#'    with default color ramp.
#' @param floor `numeric` optional floor, below which the middle
#'    color is assigned.  
#'    **Note: `floor` is not well-implemented yet.**
#'    It imposes a very pale color at the interface from middle
#'    color to the first non-middle color.
#' @param ... additional arguments are ignored.
#' 
#' @keywords internal
#' @examples
#' get_gene_color_fun(do_plot=TRUE, x=2, ceiling=5)
#' 
#' @noRd
get_gene_color_fun <- function(
   x=2,
   ceiling=8,
   lens=2,
   colramp="RdBu_r",
   floor=0,
   do_plot=FALSE,
   ...
) {
   #
   col1 <- colorjam::col_div_xf(x=x, floor=floor, lens=lens, colramp=colramp)
   col1at <- jamba::nameVector(
      attr(col1, "legend_at"),
      attr(col1, "legend_labels")
   )
   # cl <- jamba::showColors(list(
   #    col1=col1(col1at)
   # ))

   col1breaks <- c(
      -x,
      round(digits=5, seq(from=-1, to=1, by=0.2)),
      x
   )
   col1values <- col1(col1breaks)
   col2breaks <- c(
      -ceiling,
      round(digits=5, seq(from=-1, to=1, by=0.2)),
      ceiling
   )
   col2 <- circlize::colorRamp2(
      breaks=col2breaks,
      colors=col1values
   )
   
   ceil1 <- pretty(c(-ceiling, col1at[1]), n=3)
   ceil1 <- ceil1[ceil1 > -ceiling & ceil1 < col1at[1]];
   col2at <- c(
      -ceiling, ceil1, col1at, rev(-ceil1), ceiling);
   names(col2at) <- c(
      -ceiling, ceil1, names(col1at), rev(-ceil1), ceiling);
   if (isTRUE(do_plot)) {
      demoat <- seq(-ceiling, ceiling, by=0.02);
      labelat <- ifelse(
         round(digits=4, demoat)*5 == round(demoat*5),
         round(digits=4, demoat),
         ""
      )
      names(demoat) <- labelat;
      # discrete colors
      n <- floor(length(demoat) / 3);
      dn <- rep_len(col2(-1), n);
      up <- rep_len(col2(1), n);
      mid <- rep_len(col2(0), length(demoat) - n * 2);
      dnn <- rep_len("", n);
      upn <- rep_len("", n);
      dnn[ceiling(n/2)] <- "-1";
      upn[ceiling(n/2)] <- "1";
      midn <- rep_len("", length(mid));
      midn[ceiling(length(mid)/2)] <- "0";
      dnup <- jamba::nameVector(
         c(dn, mid, up),
         c(dnn, midn, upn),
         makeNamesFunc=c
      )
      jamba::showColors(
         labelCells=TRUE,
         groupCellnotes=FALSE,
         list(
            down_mid_up=dnup,
            col1=col1(demoat),
            col2=col2(demoat)
         )
      )
   }

   col2
}
