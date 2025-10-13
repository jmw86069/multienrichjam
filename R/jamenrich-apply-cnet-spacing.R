
#' Apply minimum node spacing for each Cnet node set
#'
#' Apply minimum node spacing for each Cnet node set
#'
#' This function automates the process of spacing nodesets
#' in a Cnet plot by calling `adjust_cnet_nodeset()` on
#' all nodesets whose median nearest node distance is
#' above or below a threshold.
#'
#' By default nodes must be at least `percent_spacing`
#' distance apart, relative to the maximum x-axis and
#' y-axis layout coordinate range. Roughly speaking,
#' nodes are clearly visible with distance between nodes
#' of at least 3 percent of the coordinate range.
#'
#' Nodes are not compressed down to this distance unless
#' using the argument `apply_negative=TRUE`.
#'
#' Note that this function does not detect when two
#' clusters overlap. For other specific layout adjustments
#' use `adjust_cnet_nodeset()`. In future, this function
#' may attempt to re-position nodesets relative to
#' each other, but that process tends to be iterative
#' and is not guaranteed to have a solution.
#'
#' When a nodeset contains only one node, there is no calculated
#' distance so it is reported as `NA` when `verbose=TRUE` or
#' `debug=TRUE`. However for internal purposes the nodeset is
#' treated as if its nodeset spacing is exactly `percent_spacing`
#' so that no subsequent adjustments are performed on the nodeset.
#'
#' @family jam cnet igraph functions
#'
#' @return `igraph` Cnet object, whose nodes are adjusted
#'    according to the arguments given.
#'
#' @param cnet `igraph` object that contains attribute
#'    `"nodeType"` and values `nodeType=c("Gene", "Set")`,
#'    consistent with use of `get_cnet_nodeset()` and
#'    `adjust_cnet_nodeset()`.
#'    The node layout is expected to be stored in
#'    `igraph::graph_attr(cnet, "layout")`.
#' @param percent_spacing `numeric` indicating the target spacing in relation
#'    to the maximum x-axis and y-axis span for the layout coordinates.
#'    For example `percent_spacing=4` will require that the median spacing
#'    between nodes in each nodeset is at least 4 percent of the larger of
#'    the x-axis and y-axis ranges for the layout of input `cnet`.
#' @param apply_negative `logical` indicating whether to shrink node spacing
#'    when the spacing is larger than the target `percent spacing`. By default
#'    `apply_negative=FALSE` which allows nodes to exceed the spacing, and
#'    only adjusts node spacing when a nodeset is below this threshold.
#' @param metric `character` string indicating the metric used to summarize
#'    nodeset spacing:
#'    * `"median"` - takes the median spacing between nodes in the nodeset,
#'    so the end result has at least half nodes with `percent_spacing`
#'    distance, "typical distance between nodes". This option is useful
#'    when there may be one node pair close together, but most nodes have
#'    consistent distance from their nearest node neighbor.
#'    * `"minimum"` - takes the minimum spacing between nodes in the nodeset,
#'    so the end result requires all nodes to have at least `percent_spacing`
#'    distance.
#'    * `function` - a `function` can be provided which takes a `numeric`
#'    vector of node distances, and produces a single `numeric` summary value,
#'    which will be scaled relative to the maximum xy axis range. This option
#'    would allow using `mean()`, or a quantile function such as
#'    `function(x){quantile(x, 0.3)}`.
#' @param nodesets `character` with optional subset of nodesets to adjust,
#'    or the default `nodesets=NULL` will adjust all nodesets. Values can be
#'    specified either in form or nodeset names, or any one or more members
#'    of a nodeset.
#' @param tolerance `numeric` value used to apply a small tolerance to the
#'    distance threshold, for example the default `tolerance=1.05` requires
#'    nodes to be more than 5% lower than the actual threshold before
#'    adjusting nodes. This additional activation threshold is experimental
#'    and may not be necessary. It is mainly intended to reduce calculations
#'    for nodesets that are only a small rounding error from being acceptable,
#'    which can happen when the overall plot range is adjusted during rounds
#'    of small node coordinate adjustments like with this function,
#'    `nudge_igraph_node()`, `apply_nodeset_spacing()`, or
#'    `adjust_cnet_set_relayout_gene()`.
#' @param cnet_nodesets `list` optional to provide the known nodesets
#'    upfront, in order to avoid re-determining the nodesets.
#' @param debug `logical` indicating whether to print debug information
#'    without adjusting the nodeset spacing. This option is useful to see
#'    summary information about nodeset spacing before choosing the
#'    `percent_spacing` threshold.
#' @param verbose `logical` indicating whether to print verbose output.
#'    The output includes text printed with `debug=TRUE` except that nodes
#'    are adjusted when `debug=FALSE`.
#' @param ... additional arguments are ignored.
#'
#' @export
apply_nodeset_spacing <- function
(cnet,
 percent_spacing=3,
 apply_negative=FALSE,
 metric=c("median", "minimum"),
 nodesets=NULL,
 tolerance=1.05,
 cnet_nodesets=NULL,
 debug=FALSE,
 verbose=FALSE,
 ...)
{
   # validate metric argument
   if (length(metric) == 0) {
      metric <- "median"
   }
   use_metric <- "median";
   if ("character" %in% class(metric)) {
      metric <- match.arg(metric);
      use_metric <- metric;
      if ("median" %in% metric) {
         metric <- median;
      } else if ("min" %in% metric) {
         metric <- min;
      }
   } else if (!"function" %in% class(metric)) {
      stop("metric must either be class 'character' or 'function'");
   }

   # get node layout
   layout_xy <- igraph::graph_attr(cnet, "layout")
   if (length(layout_xy) == 0) {
      stop("apply_nodeset_spacing() requires layout stored in graph_attr(cnet, 'layout')");
   }
   rownames(layout_xy) <- igraph::V(cnet)$name;

   # get nodesets
   if (length(cnet_nodesets) == 0) {
      cnet_nodesets <- get_cnet_nodeset(cnet)
   }
   if (length(nodesets) > 0) {
      ikeep_nodeset <- sapply(jamba::nameVectorN(cnet_nodesets), function(iname){
         if (iname %in% nodesets) {
            return(TRUE)
         }
         inodeset <- cnet_nodesets[[iname]];
         if (any(inodeset %in% nodesets)) {
            return(TRUE)
         }
         FALSE
      })
      cnet_nodesets <- cnet_nodesets[ikeep_nodeset];
      keep_nodesets <- names(cnet_nodesets);
      if (length(keep_nodesets) == 0) {
         stop("After filtering by nodesets argument, no data remained.")
      }
      if (verbose) {
         jamba::printDebug("apply_nodeset_spacing(): ",
            "Subset nodesets to keep: ",
            sep="; ",
            paste0("'", keep_nodesets, "'"))
      }
   }

   # determine nodeset spacing relative to layout size
   sns <- summarize_node_spacing(cnet,
      layout=layout_xy,
      node_groups=cnet_nodesets,
      dist_type="nearest_node",
      scaled=TRUE);
   cnet_spacing_df <- data.frame(check.names=FALSE,
      nodeset=rownames(sns$nearest_within),
      cnet_spacing=sns$nearest_within[, 'Median'],
      expand=0)
   rownames(cnet_spacing_df) <- cnet_spacing_df[, 1];
   if ("minimum" %in% use_metric) {
      cnet_spacing_df$cnet_spacing <- sns$nearest_within[, 'Min'];
   }

   # jamba::printDebug("cnet_spacing_df:");print(cnet_spacing_df);# debug
   cnet_spacing <- cnet_spacing_df$cnet_spacing;
   names(cnet_spacing) <- cnet_spacing_df$nodeset;

   # determine nodesets to expand
   # use tolerance to allow some small wiggle room
   adjust_nodesets <- names(which(cnet_spacing < (percent_spacing / tolerance)));

   adjust_nodesets_expand <- percent_spacing / cnet_spacing;
   names(adjust_nodesets_expand) <- names(cnet_spacing);
   adjust_nodesets_expand <- ifelse(adjust_nodesets_expand >= 1,
      adjust_nodesets_expand - 1,
      -1 * (1/adjust_nodesets_expand - 1));

   # optionally ignore nodesets that are already above the threshold
   if (apply_negative) {
      adjust_nodesets_neg <- names(which(!is.infinite(cnet_spacing) &
            cnet_spacing > (percent_spacing * tolerance)))
      adjust_nodesets <- c(adjust_nodesets,
         adjust_nodesets_neg);
   }
   if (length(adjust_nodesets) > 0) {
      imatch <- match(adjust_nodesets, rownames(cnet_spacing_df));
      cnet_spacing_df[imatch, "expand"] <- adjust_nodesets_expand[adjust_nodesets];
   }
   if (debug) {
      return(cnet_spacing_df)
   }
   if (length(adjust_nodesets) == 0) {
      if (verbose) {
         jamba::printDebug("apply_nodeset_spacing(): ",
            "No nodesets required adjustment.");
      }
      return(cnet)
   }


   # apply adjustments
   if (verbose) {
      jamba::printDebug("apply_nodeset_spacing(): ",
         "Calling adjust_cnet_nodeset()");
      jamba::printDebug("apply_nodeset_spacing(): ",
         "expand: ",
         adjust_nodesets_expand[adjust_nodesets]);
   }
   cnet_new <- adjust_cnet_nodeset(cnet,
      set_nodes=strsplit(adjust_nodesets, ","),
      x=0,
      y=0,
      expand=adjust_nodesets_expand[adjust_nodesets],
      cnet_nodesets=cnet_nodesets,
      verbose=verbose)
   return(cnet_new)
}

