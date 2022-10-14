
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
#'    `colorjam::col_div_xf(1.5)`.
#' @param hide_solo_pie `logical` indicating whether a single-color
#'    border should only be applied to `"frame.color"` and not to
#'    individual `"pie.border"` entries. When `hide_solo_pie=FALSE`
#'    the border color is applied to each pie wedge, and if all
#'    pie wedges are identical, it is also applied as `"frame.color"`.
#' @param frame_blank `character` string to define the color used
#'    for `"frame.color"` when there are multiple colors in `"pie.border"`
#'    and therefore there should be no visible `"frame.color"`.
#' @param border_blank
#' @param ... additional arguments are ignored.
#'
#' @export
apply_cnet_direction <- function
(cnet,
 hitim,
 col=colorjam::col_div_xf(1.2),
 hide_solo_pie=TRUE,
 frame_blank="#FFFFFF00",
 ...)
{
   #
   frame_blank <- head(frame_blank, 1);
   # pie.border, coloredrect.border
   for (attr_name in c("pie.border", "coloredrect.border")) {
      igraph::vertex_attr(cnet, attr_name) <- lapply(seq_along(igraph::V(cnet)), function(i){
         iname <- igraph::vertex_attr(cnet, "name")[i];
         ienrich <- names(igraph::vertex_attr(cnet, attr_name)[[i]]);
         iborder <- igraph::vertex_attr(cnet, attr_name)[i];
         if (iname %in% rownames(hitim)) {
            ipieborder <- jamba::nameVector(
               col(hitim[iname, ienrich]),
               ienrich)
            if (!TRUE %in% hide_solo_pie || length(unique(ipieborder)) > 1) {
               return(ipieborder)
            }
         }
         iborder
      })
   }
   # frame.color
   cnet_framecolor <- sapply(seq_along(igraph::V(cnet)), function(i){
      iframecolor <- igraph::vertex_attr(cnet, "frame.color")[[i]];
      iname <- igraph::vertex_attr(cnet, "name")[i];
      ienrich <- names(igraph::vertex_attr(cnet, attr_name)[[i]]);
      if (iname %in% rownames(hitim)) {
         ipieborder <- jamba::nameVector(
            col(hitim[iname, ienrich]),
            ienrich)
         if (length(unique(ipieborder)) == 1) {
            return(unique(ipieborder))
         } else if (length(unique(ipieborder)) > 1) {
            return(frame_blank)
         }
      }
      iframecolor
   })
   vertex_attr(cnet, "frame.color") <- cnet_framecolor
   return(cnet);
}
