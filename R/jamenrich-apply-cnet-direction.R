
#' Apply Cnet border color by directionality
#'
#' Apply Cnet border color by directionality
#'
#' This function specifically requires the `"pie.color"`
#' node attribute of the `cnet` object `igraph` is populated
#' with a `list` of colors whose name is the enrichment,
#' and matches `colnames(hitim)` of the supplied hit
#' incidence matrix.
#'
#' The node names `V(cnet)$name` also must match `rownames(hitim)`,
#' otherwise the `"pie.border"` and `"coloredrect.border"` is assigned
#' its previous value.
#' In that case `"frame.color"` retains its original value.
#'
#' This function may be run multiple times with different `hitim`,
#' for example using `geneIMdirection` and `enrichIMdirection` in two steps.
#'
#' When there are multiple colors for a given node, they are populated in
#' node attributes `"pie.border"` in `list` form. For these nodes,
#' the attribute `"frame.color"` is set to `NULL` so it is not displayed
#' on top of the pie wedge colors.
#'
#' When there are no matching `rownames(hitim)`, or all colors
#' are identical, the `"pie.border"` is populated with `NULL`,
#' so the pie wedges do not each show a color. Also for single-wedge
#' pie nodes, this process avoids drawing a small line at the top
#' of each node.
#' Instead the attribute `"frame.color"` is populated with the
#' one unique border color, so that only the outer border is colorized.
#'
#' @family jam cnet igraph functions
#'
#' @param cnet `igraph` object with node attribute `"pie.color"`
#'    populated as a `list` of `character` vectors, named by
#'    enrichment. The enrichment names should match `colnames(hitim)`.
#' @param hitim `numeric` matrix with values centered at zero for
#'    no change, positive values `+1` for up-regulation, and
#'    negative values `-1` for down-regulation. The values are
#'    converted to color using color function `col`.
#' @param col `function` that takes `numeric` input and assigns
#'    a color. The default assigns red for positive values, blue
#'    for negative values, and white for zero, using
#'    `colorjam::col_div_xf(1.2)`.
#' @param col_l_max `numeric` maximum HCL Lightness permitted for
#'    output colors, for example because the middle color in `"RdBu_r"`
#'    is nearly white, it may be preferable to apply a darker grey
#'    color.
#' @param hide_solo_pie `logical` indicating whether a single-color
#'    border should only be applied to `"frame.color"` and not to
#'    individual `"pie.border"` entries.
#'    * When `hide_solo_pie=TRUE` and all wedges (one or more) will
#'    have the same `pie.border` color, then `pie.border` is defined
#'    as `NA`, and instead the `frame.color` is assigned to this color.
#'    Both `pie.lwd` and `frame.lwd` will be assigned `border_lwd`,
#'    since `pie.border` is `NA` it will not be rendered, only the
#'    `frame.color` will be rendered.
#'    * When `hide_solo_pie=FALSE` each pie wedge border color is assigned
#'    to `pie.border`, `frame.color` will be assigned `frame_blank`,
#'    and `frame.lwd` will be assigned `frame_lwd_blank` which is useful
#'    for displaying a small outer frame for each node.
#' @param frame_blank `character` string to define the color used
#'    for `"frame.color"` when there are multiple colors in `"pie.border"`
#'    and therefore there should be no visible `"frame.color"`.
#' @param do_reorder `logical` indicating whether to call
#'    `reorder_igraph_nodes()` on the resulting `igraph`, so that
#'    the border color can be used in the sort conditions.
#' @param ... additional arguments are ignored.
#'
#' @export
apply_cnet_direction <- function
(cnet,
 hitim,
 col=colorjam::col_div_xf(1.2),
 col_l_max=65,
 hide_solo_pie=TRUE,
 frame_blank=default_igraph_values()$vertex$frame.color,
 frame_lwd_blank=1,
 border_lwd=2,
 do_reorder=FALSE,
 ...)
{
   #
   frame_blank <- head(frame_blank, 1);
   # pie.border, coloredrect.border
   attr_names <- c("pie.border",
      "pie.lwd",
      "coloredrect.border",
      "coloredrect.lwd");
   cap_color_l <- c;
   if (col_l_max < 100) {
      cap_color_l <- function(x, l){
         l <- jamba::noiseFloor(l, minimum=0, ceiling=100);
         xhcl <- jamba::col2hcl(x);
         xhcl["L",] <- jamba::noiseFloor(xhcl["L",],
            minimum=0,
            ceiling=l)
         jamba::hcl2col(xhcl)
      }
   }
   for (attr_name in attr_names) {
      attr_values <- lapply(seq_along(igraph::V(cnet)), function(i){
         iname <- igraph::vertex_attr(cnet, "name")[[i]];
         attr_name2 <- gsub("[.](lwd|border)$", ".color", attr_name);
         attr_name1 <- gsub("[.]lwd$", ".border", attr_name);
         icolor <- igraph::vertex_attr(cnet, attr_name2)[[i]];
         if (length(names(icolor)) == 0) {
            if (length(icolor) == ncol(hitim)) {
               ienrich <- colnames(hitim);
            } else {
               ienrich <- NULL;
            }
         } else {
            ienrich <- names(icolor);
         }
         iborder <- igraph::vertex_attr(cnet, attr_name1)[[i]];
         # only apply direction when the node name matches rownames(hitim)
         if (length(ienrich) > 0 &&
               iname %in% rownames(hitim)) {
            ipieborder <- jamba::nameVector(
               cap_color_l(col(hitim[iname, ienrich]), col_l_max),
               ienrich)
            if (TRUE %in% hide_solo_pie &&
                  length(unique(ipieborder)) == 1) {
               # hide_solo_pie is active
               if (grepl("border", attr_name)) {
                  return(rep(NA,
                     length.out=length(ipieborder)))
               } else {
                  return(rep(border_lwd,
                     length.out=length(ipieborder)))
               }
            }
            # hide_solo_pie is not active
            if (grepl("border", attr_name)) {
               return(ipieborder)
            } else {
               return(rep(border_lwd,
                  length.out=length(ipieborder)))
            }
         }
         # no direction is applied
         if (grepl("border", attr_name)) {
            return(iborder)
         }
         # assign the same pie.lwd as already present
         ipielwd <- jamba::rmNA(naValue=1,
            igraph::vertex_attr(cnet, attr_name)[[i]]);
         if (length(ipielwd) == 0) {
            ipielwd <- 1;
         }
         if (length(iborder) == 0 || all(is.na(iborder))) {
            return(border_lwd)
         }
         return(rep(ifelse(ipielwd == 0, 1, ipielwd),
            length.out=length(iborder)))
      })
      igraph::vertex_attr(cnet, attr_name) <- attr_values;
   }

   # frame.color
   cnet_framecolor <- sapply(seq_along(igraph::V(cnet)), function(i){
      iname <- igraph::vertex_attr(cnet, "name")[[i]];
      icolor <- igraph::vertex_attr(cnet, "pie.color")[[i]];
      if (length(names(icolor)) == 0) {
         if (length(icolor) == ncol(hitim)) {
            ienrich <- colnames(hitim);
         } else {
            ienrich <- NULL;
         }
      } else {
         ienrich <- names(icolor);
      }
      iframecolor <- igraph::vertex_attr(cnet, "frame.color")[[i]];
      # only apply direction when the node name matches rownames(hitim)
      if (length(ienrich) > 0 &&
            iname %in% rownames(hitim)) {
         ipieborder <- jamba::nameVector(
            cap_color_l(col(hitim[iname, ienrich]), col_l_max),
            ienrich)
         if (TRUE %in% hide_solo_pie &&
               length(unique(ipieborder)) == 1) {
            return(unique(ipieborder))
         }
         return(frame_blank)
      }
      # keep the original frame.color if no direction was applied
      iframecolor
   })
   igraph::vertex_attr(cnet, "frame.color") <- cnet_framecolor

   # frame.lwd
   cnet_framelwd <- sapply(seq_along(igraph::V(cnet)), function(i){
      iname <- igraph::vertex_attr(cnet, "name")[[i]];
      icolor <- igraph::vertex_attr(cnet, "pie.color")[[i]];
      if (length(names(icolor)) == 0) {
         if (length(icolor) == ncol(hitim)) {
            ienrich <- colnames(hitim);
         } else {
            ienrich <- NULL;
         }
      } else {
         ienrich <- names(icolor);
      }
      iframecolor <- igraph::vertex_attr(cnet, "frame.color")[[i]];
      iframelwd <- igraph::vertex_attr(cnet, "frame.lwd")[[i]];
      if (length(iframelwd) == 0) {
         iframelwd <- 1;
      }
      if (length(ienrich) > 0 &&
            iname %in% rownames(hitim)) {
         ipieborder <- jamba::nameVector(
            cap_color_l(col(hitim[iname, ienrich]), col_l_max),
            ienrich)
         if (TRUE %in% hide_solo_pie &&
               length(unique(ipieborder)) == 1) {
            return(border_lwd)
         }
         return(frame_lwd_blank)
      }
      # if no directionality, then keep the existing frame.lwd
      return(iframelwd)
   })
   igraph::vertex_attr(cnet, "frame.lwd") <- cnet_framelwd

   # optionally reorder nodes using new border color
   if (TRUE %in% do_reorder) {
      # apply node re-ordering step
      cnet <- reorderIgraphNodes(cnet,
         ...);
   }

   return(cnet);
}
