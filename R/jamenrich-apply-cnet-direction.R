
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
#' ## Todo
#' 
#' * Consider permitting input to be a simple igraph, then creating
#' pie.color, pie.border, pie.border.lwd attributes even if they
#' previously did not exist. It would probably optionally also use
#' `colorV` to colorize each component of multi-component nodes.
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
#' @param col_l_max `numeric` maximum HCL Lightness, default 80, for
#'    output colors. For example, the middle color in `"RdBu_r"`
#'    is nearly white, the `col_l_max` can be used to apply a darker grey.
#' @param hide_solo_pie `logical` default TRUE, whether a single-color
#'    border for a multi-part pie node should only apply the color
#'    to the overall node with 'frame.color', and not apply the color
#'    to each pie wedge using 'pie.border'.
#'    * Default `hide_solo_pie=TRUE`: when all wedges (one or more)
#'    have the same `pie.border` color, the 'pie.border' is defined
#'    as `NA`, and 'frame.color' is assigned to this color.
#'    The effect is to display the outline color and not each wedge.
#'    Both 'pie.lwd' and 'frame.lwd' will be assigned `border_lwd`,
#'    and since 'pie.border' is `NA` it will not be rendered. Only the
#'    'frame.color' will be rendered.
#'    * When `hide_solo_pie=FALSE` each pie wedge border color is assigned
#'    to 'pie.border', `frame.color` will be assigned 'frame_blank',
#'    and 'frame.lwd' will be assigned 'frame_lwd_blank' which is useful
#'    for displaying a small outer frame for each node.
#' @param frame_blank `character` string to define the color used
#'    for 'frame.color' when colors are defined in 'pie.border'.
#'    The default uses the igraph defaults, currently 'black'.
#'    In this case, the frame is drawn around the inner `pie.border`
#'    colors, and only serves to add visual clarity. The frame border
#'    can be blank `frame_blank="transparent"` or can be a thinner line,
#'    controlled with `frame_lwd_blank=0.2`.
#' @param frame_lwd_blank `numeric` line width, default 0.2,
#'    for nodes that have "blank" frame, which also means the
#'    'pie.border' colors must also defined. In this case
#'    the frame border can be invisible (`frame_lwd_blank=0`) or
#'    a very thin line (default `frame_lwd_blank=0.2`) to surround the
#'    inner borders drawn with 'pie.border'.
#' @param border_lwd `numeric` line width, default 2, used when a node
#'    matches `rownames(hitim)`.
#'    When the colors are applied to 'pie.color',
#'    the border is defined with 'pie.lwd'.
#'    When colors are applied to 'frame.color', the border is defined
#'    with 'frame.lwd'.
#' @param do_reorder `logical` default FALSE, whether to reorder nodes
#'    by node attributes such as color and border, by calling
#'    `reorder_igraph_nodes()`.
#'    When `do_reorder=TRUE`, other relevant arguments are passed
#'    through `...` to `reorder_igraph_nodes()` in particular:
#'    * `colorV`: to control the expected order of colors.
#'    It should be supplied if known upfront.
#'    * `sortAttributes`: to customize the default attribute sort order.
#'    * `nodeSortBy`: to customize the x-/y- axis arrangement.
#'    * `orderByAspect`: to enable x-/y- sorting by the aspect ratio of
#'    nodes in the group. For example, tall-skinny node groups should sort by
#'    y-axis first, short-wide node groups should sort by x-axis first.
#' @param ... additional arguments are passed to `reorder_igraph_nodes()`
#'    when `do_reorder=TRUE`.
#'
#' @export
apply_cnet_direction <- function
(cnet,
 hitim=NULL,
 col=circlize::colorRamp2(breaks=c(-1, 0, 1),
    colors=c("blue", "grey80", "firebrick3")),
 col_l_max=80,
 hide_solo_pie=TRUE,
 frame_blank=default_igraph_values()$vertex$frame.color,
 frame_lwd_blank=0.2,
 border_lwd=2,
 do_reorder=FALSE,
 verbose=FALSE,
 ...)
{
   # if input hit incidence matrix is empty, return input cnet
   if (length(hitim) == 0) {
      return(cnet)
   }
   # if input cnet is not "igraph" then return unchanged
   if (!"igraph" %in% class(cnet)) {
      return(cnet)
   }
   frame_blank <- head(frame_blank, 1);
   if (length(frame_blank) == 0) {
      # when supplied as NULL, interpret as NA for no frame.color
      frame_blank <- NA
   }
   # pie.border, coloredrect.border
   attr_names <- c("pie.border",
      "coloredrect.border");
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
   
   # 0.0.104.900: vectorized approach - the rest will be removed soon
   use_approach <- "vectorized";
   
   for (attr_name in attr_names) {
      # attr_name1 is the border color(s)
      attr_name1 <- attr_name;
      # attr_name2 is the fill color(s)
      attr_name2 <- gsub("[.](border)$", ".color", attr_name);
      if (verbose) {
         jamba::printDebug("apply_cnet_direction(): ",
            "applying attribute: ", attr_name,
            ", attr_name2:", attr_name2, ", attr_name1:", attr_name1);
      }
      
      # 0.0.104.900: vectorized logic
      # use only rows represented in hitim
      whichv <- which(igraph::V(cnet)$name %in% rownames(hitim));
      
      if ("vectorized" %in% use_approach) {
         # node names
         inames <- igraph::vertex_attr(cnet, "name")[whichv]
         
         # node pie.color or coloredrect.color (list of color vectors)
         icolor <- igraph::vertex_attr(cnet, attr_name2)[whichv]
         if (is.atomic(icolor)) {
            # it can only have one color per node, convert to list
            icolorl <- rep(1, length(icolor))
            # but it cannot have names, required to match colnames(hitim)
            icolornames <- NULL;
         } else {
            icolorl <- lengths(icolor);
            icolornames <- lapply(icolor, names)
         }
         icolordf <- data.frame(check.names=FALSE,
            inames=rep(inames, icolorl),
            icolor=unlist(unname(icolor)),
            whichv=factor(rep(whichv, icolorl), levels=whichv))
         
         # icolornames must not contain NULL
         if (length(icolornames) > 0 && all(lengths(icolornames) > 0)) {
            icolordf$icolornames <- unlist(unname(icolornames))
         } else {
            # skip the rest
            if (verbose) {
               jamba::printDebug("apply_cnet_direction(): ",
                  "no names('", attr_name,
                  "') to assign direction, skipping.");
            }
            next;
         }
         
         # node pie.border or coloredrect.border (list of color vectors)
         # (Is this section necessary?)
         iborder <- igraph::vertex_attr(cnet, attr_name1)[whichv]
         if (is.atomic(iborder)) {
            iborderl <- rep(1, length(iborder))
         } else {
            iborderl <- lengths(iborder);
         }
         if (length(iborder) == 0) iborder <- NA;
         
         if (length(unlist(iborder)) == nrow(icolordf)) {
            icolordf$iborder <- unlist(unname(iborder));
         } else {
            # when in doubt use NA to apply no color
            icolordf$iborder <- NA;
         }
         
         ## Add numeric value
         # matrix index to hitim using inames,icolornames
         midx <- as.matrix(icolordf[, c("inames", "icolornames"), drop=FALSE]);
         ipieborder <- hitim[midx];
         icolordf$ipieborder <- ipieborder;
         
         ## Calculate expected border color upfront so we can recognize dupes
         icolordf$ipiebordercolor <- cap_color_l(
            col(icolordf$ipieborder), col_l_max)
         zerocolor <- jamba::unalpha(
            cap_color_l(col(0), col_l_max));
         
         # find which pie have multiple directions
         if (TRUE %in% hide_solo_pie) {
            # count unique ipieborder values per node
            # Note: use expected border color and not the numeric value
            # - otherwise use "ipieborder" for just the numeric value
            ipietc <- jamba::tcount(
               unique(icolordf[, c("ipiebordercolor", "inames"),
                     drop=FALSE])$inames);

            # apply to the data.frame summary
            icolordf$ipieborderct <- ipietc[icolordf$inames];

            # determine which nodes to apply border to the node overall
            apply_frame_only <- (
               icolordf$ipieborderct == 1);
            
            # determine section color
            # Note: color that matches zero is not shown
            section_color <- ifelse(
               icolordf$ipieborder %in% c(0, NA) |
                  icolordf$ipiebordercolor %in% zerocolor,
               NA,
               cap_color_l(col(icolordf$ipieborder), col_l_max));
            icolordf$use_border <- ifelse(
               apply_frame_only,
               NA,
               section_color)
            # frame.lwd is tied to this value during rendering
            # so all pie.lwd must be at least the frame border lwd
            icolordf$use_border_lwd <- border_lwd;

            # frame color and lwd
            icolordf$use_frame <- ifelse(
               apply_frame_only,
               section_color,
               frame_blank)
            icolordf$use_frame_lwd <- ifelse(
               apply_frame_only,
               border_lwd,
               frame_lwd_blank) #
         } else {
            section_color <- ifelse(
               icolordf$ipieborder %in% c(0, NA),
               NA,
               cap_color_l(col(icolordf$ipieborder), col_l_max));
            icolordf$use_border <- section_color;
            icolordf$use_border_lwd <- ifelse(
               is.na(section_color),
               frame_lwd_blank,
               border_lwd)
            icolordf$use_frame <- ifelse(
               icolordf$ipieborderct > 1,
               frame_blank,
               NA);
            icolordf$use_frame_lwd <- ifelse(
               icolordf$ipieborderct > 1,
               frame_lwd_blank,
               border_lwd);
         }
         
         # update numeric direction for convenience
         if (!"direction.values" %in% igraph::vertex_attr_names(cnet)) {
            igraph::vertex_attr(cnet, name="direction.values") <- rep(
               list(), length.out=igraph::vcount(cnet));
         }
         use_dir_values <- split(
            # icolordf$ipieborder,
            ## consider assigning names
            jamba::nameVector(icolordf$ipieborder,
               icolordf$icolornames,
               makeNamesFunc=c),
            icolordf$whichv);
         igraph::vertex_attr(cnet,
            index=as.integer(names(use_dir_values)),
            name="direction.values") <- use_dir_values;
         
         # update sectioned border color values
         use_border <- split(icolordf$use_border, icolordf$whichv);
         use_border_lwd <- split(icolordf$use_border_lwd, icolordf$whichv);
         igraph::vertex_attr(cnet,
            index=as.integer(names(use_border)),
            name=attr_name1) <- use_border;
         
         # use pie.lwd and not pie.border.lwd
         attr_name1_border <- paste0(
            gsub("pie.border", "pie", attr_name1),
            ".lwd");
         igraph::vertex_attr(cnet,
            index=as.integer(names(use_border_lwd)),
            name=attr_name1_border) <- use_border_lwd;
         
         icolordf_frame <- subset(icolordf, !duplicated(inames));
         igraph::vertex_attr(cnet,
            index=icolordf_frame$whichv,
            name="frame.color") <- icolordf_frame$use_frame;
         # igraph::vertex_attr(cnet, "frame.color",
         #    index=icolordf_frame$whichv) <- icolordf_frame$use_frame;
         igraph::vertex_attr(cnet,
            index=icolordf_frame$whichv,
            name="frame.lwd") <- icolordf_frame$use_frame_lwd;
         # jamba::printDebug("frame.lwd:");print(igraph::vertex_attr(cnet,
         #    index=icolordf_frame$whichv,
         #    name="frame.lwd"));# debug
         
         # igraph::vertex_attr(cnet, "frame.lwd",
         #    index=icolordf_frame$whichv) <- icolordf_frame$use_frame_lwd;
         igraph::vertex_attr(cnet, "frame.width") <- (
            igraph::vertex_attr(cnet, "frame.lwd"));
      }
      
      # looping each vector (below) seems inefficient, and is quite slow
      if (!"vectorized" %in% use_approach) {
         attr_values <- lapply(seq_along(igraph::V(cnet)), function(i){
            iname <- igraph::vertex_attr(cnet, "name")[[i]];
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
   }

   # frame.color
   # 0.0.104.900: no longer necessary
   if (!"vectorized" %in% use_approach) {
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
   }

   # optionally reorder nodes using new border color
   if (TRUE %in% do_reorder) {
      # apply node re-ordering step
      cnet <- reorder_igraph_nodes(cnet,
         ...);
   }

   return(cnet);
}
