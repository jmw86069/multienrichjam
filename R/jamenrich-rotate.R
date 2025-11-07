
# jam coordinate rotation functions

#' Rotate igraph layout coordinates
#'
#' Rotate igraph layout coordinates, optionally after reflecting coordinates
#' along one or more coordinate axes.
#'
#' This function rotates igraph layout coordinates by calling
#' the function `rotate_coordinates()`. The input can either be
#' `g` as `igraph` object, or `layout` as a numeric `matrix`.
#'
#' Note that the `reflect` is applied before `degrees`. To change
#' the order, call this function multiple times.
#'
#' When both `g` and `layout` are supplied, the coordinates are
#' used from `layout`, rotated, then stored in the `g` `igraph` object
#' as a graph attribute, using `igraph::graph_attr(g, "layout")`.
#'
#' When only `g` is supplied, it is expected to contain
#' layout coordinates in graph attributes, obtained with
#' `igraph::graph_attr(g, "layout")`.
#'
#' When only `layout` is supplied, and no `g` `igraph` object
#' is supplied, this function serves only as a wrapper to
#' `rotate_coordinates()`.
#'
#' Rotation code kindly contributed by Don MacQueen to the `maptools`
#' package, and is reproduced here to avoid a dependency on `maptools`
#' and therefore the `sp` package.
#'
#' This function also calls other useful helper functions,
#' when `spread_labels=TRUE` it calls `spread_igraph_labels()` to
#' position labels around each node based upon the angle of
#' incoming edges, which has the effect of reducing label overlaps.
#' When `do_reorder=TRUE` it calls `reorderIgraphNodes()` which
#' sorts nodes within a nodeset by color then by label, to help
#' visually group similar nodes together in equivalent positions
#' in the layout.
#'
#' @family jam cnet igraph functions
#' @family jam igraph layouts
#'
#' @return `igraph` when input `g` is supplied, otherwise `numeric matrix`.
#'
#' @param g `igraph` object that contains layout coordinates in
#'    graph attributes, for example `igraph::graph_attr(g, "layout")`.
#' @param degrees numeric value indicating the degrees to rotate
#'    layout coordinates, where 360 degrees is one complete rotation.
#' @param reflect `character` string indicating one or more axes
#'    to reflect coordinates, or `"none"` to reflect no axis.
#' @param center `numeric` coordinates to use as the center, or
#'    `center=NULL` to calculate the center using `center_rule`.
#' @param center_rule `character` string indicating which rule to
#'    apply to determine the center coordinates, when `center=NULL`:
#'    `"origin"` uses c(0, 0); `"mean"` uses the mean of each axis;
#'    `"median"` uses the median of each axis; `"min"` uses the minimum
#'    of each axis; `"max"` uses the max of each axis.
#' @param rotation_axes `integer` vector indicating which axis
#'    coordinates to rotate, by default `c(1, 2)` uses the first
#'    two axes.
#' @param spread_labels,do_reorder `logical` indicating whether to
#'    call `spread_igraph_labels()`, and subsequently
#'    whether to call `reorderIgraphNodes()`.
#' @param layout `matrix` with 2 or more columns, when defined this
#'    layout is used and not the layout from the `g` `igraph` object.
#' @param ... additional arguments are passed to
#'    `spread_igraph_labels()` which calls `reorderIgraphNodes()` when
#'    `spread_labels=TRUE` and `do_reorder=TRUE`, or to
#'    `reorderIgraphNodes()` when `spread_labels=FALSE` and
#'    `do_reorder=TRUE`. Notably, the optional argument `sortAttributes`
#'    can be passed through those functions to affect the node sort
#'    priority.
#'
#' @examples
#' layout <- cbind(0:10, 0:10);
#' layout_rot50 <- rotate_igraph_layout(layout=layout, degrees=50);
#' layout_rot40_ctrmean <- rotate_igraph_layout(layout=layout, degrees=40, center_rule="mean");
#' plot(rbind(layout, layout_rot50, layout_rot40_ctrmean),
#'    col=rep(c("darkorchid", "darkorange1", "dodgerblue"), each=11),
#'    pch=rep(c(17, 20, 18), each=11),
#'    cex=2);
#'
#' if (require(igraph)) {
#'    g <- igraph::make_graph( ~ A-B-C-D-A, E-A:B:C:D,
#'       F-G-H-I-F, J-F:G:H:I,
#'       K-L-M-N-K, O-K:L:M:N,
#'       P-Q-R-S-P, T-P:Q:R:S,
#'       B-F, E-J, C-I, L-T, O-T, M-S,
#'       C-P, C-L, I-L, I-P);
#'    g <- relayout_with_qfr(g, repulse=8);
#'    g2 <- rotate_igraph_layout(g, degrees=45);
#'    opar <- par("mfrow"=c(1,2));
#'    on.exit(par(opar));
#'    jam_igraph(g,
#'       main="original layout",
#'       node_factor=0.6,
#'       label_dist_factor=7);
#'    jam_igraph(g2,
#'       main="rotated 45 degrees",
#'       node_factor=0.6,
#'       label_dist_factor=7);
#' }
#'
#' @export
rotate_igraph_layout <- function
(g=NULL,
 degrees=0,
 reflect=c("none", "x", "y", "z"),
 center=NULL,
 center_rule=c("median",
    "origin",
    "mean",
    "min",
    "max"),
 rotation_axes=c(1, 2),
 spread_labels=TRUE,
 do_reorder=FALSE,
 layout=NULL,
 verbose=FALSE,
...)
{
   if (length(layout) == 0) {
      if (length(g) == 0) {
         stop("Input must contain either igraph object g, or layout");
      }
      if ("igraph" %in% class(g)) {
         layout <- igraph::graph_attr(g, "layout");
         # ensure rownames match node names
         rownames(layout) <- igraph::V(g)$name;
      } else {
         stop("Input g must be an igraph object, no layout was supplied.");
      }
   }
   # when g and layout are supplied,
   # make sure V(g)$name matches the order of rownames(layout)
   if (length(g) > 0 &&
         "igraph" %in% class(g) &&
         all(igraph::V(g)$name %in% rownames(layout)) &&
         !all(igraph::V(g)$name == rownames(layout))) {
      # confirm order matches node order
      if (verbose) {
         jamba::printDebug("rotate_igraph_layout(): ",
            "Re-ordered layout so rownames(layout) match V(g)$name.")
      }
      layout <- layout[match(igraph::V(g)$name, rownames(layout)), , drop=FALSE];
   }
   if (length(dim(layout)) < 2) {
      stop("Layout must contain at least 2 columns.");
   }
   center_rule <- match.arg(center_rule);
   layout <- rotate_coordinates(layout,
      degrees=degrees,
      reflect=reflect,
      center=center,
      center_rule=center_rule,
      rotation_axes=rotation_axes,
      ...);

   if (length(g) > 0) {
      g <- igraph::set_graph_attr(g,
         name="layout",
         value=layout);
      if (spread_labels) {
         g <- multienrichjam::spread_igraph_labels(g,
            do_reorder=do_reorder,
            update_g_coords=TRUE,
            ...);
      } else if (do_reorder) {
         g <- reorderIgraphNodes(g,
            ...);
      }
      return(g);
   }
   return(layout);
}

#' Rotate numeric coordinates
#'
#' Rotate numeric coordinates, optionally after reflecting coordinates
#' along one or more coordinate axes.
#'
#' This function rotates coordinates in two axes, by the angle
#' defined in `degrees`. It optionally reflects coordinates in
#' one or more axes, which occurs before rotation.
#'
#' Note that the `reflect` is applied before `degrees`.
#'
#' Rotation code kindly contributed by Don MacQueen to the `maptools`
#' package, and is reproduced here to avoid a dependency on `maptools`
#' and therefore the `sp` package.
#'
#' @return `numeric matrix` with the same number of columns as the
#'    input `x`.
#'
#' @family jam utility functions
#'
#' @param x `matrix` with 2 or more columns.
#' @param degrees numeric value indicating the degrees to rotate
#'    layout coordinates, where 360 degrees is one complete rotation.
#' @param reflect `character` string indicating one or more axes
#'    to reflect coordinates, which flips the position of coordinates
#'    along that axis. It is usually called to flip x-axis or y-axis
#'    coordinates, for example with `reflect="x"` or `reflect=1`.
#'    Input is handled as follows:
#'    * if `reflect` contains `"none"`, then reflect is applied to
#'    none of the coordinate axes, therefore the default
#'    `reflect=c("none", "x", "y", "z")` will apply no reflection.
#'    * `character` input: `reflect` values are matched to `colnames(x)`.
#'    When there are no `colnames(x)`, then `reflect` values of
#'    `c("x", "y", "z")` are automatically recognized as columns
#'    `c(1, 2, 3)` respectively.
#'    * `integer` input is treated as a vector of column index positions,
#'    for example `reflect=c(2)` will reflect values on the second
#'    coordinate column.
#' @param center `numeric` coordinates to use as the center, or
#'    `center=NULL` to calculate the center using `center_rule`.
#' @param center_rule `character` string indicating which rule to
#'    apply to determine the center coordinates when `center=NULL`.
#'    Note that it has little effect on most downstream plotting
#'    assuming the plot function adjusts x- and y-axis ranges to
#'    the data range, but may modify the axis ranges as a result.
#'    * `"origin"` uses c(0, 0);
#'    * `"mean"` uses the mean of each axis;
#'    * `"median"` uses the median of each axis;
#'    * `"min"` uses the minimum of each axis;
#'    * `"max"` uses the max of each axis.
#' @param rotation_axes `integer` vector indicating which axis
#'    coordinates to rotate, by default `c(1, 2)` uses the first
#'    two axes in `x`. Note that `rotation_axes` must represent
#'    columns present in x.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' layout <- cbind(0:10, 0:10);
#' layout_rot50 <- rotate_coordinates(x=layout, degrees=50);
#' layout_rot40_ctrmean <- rotate_coordinates(x=layout, degrees=40, center_rule="mean");
#' layout_reflectx_ctrmean <- rotate_coordinates(x=layout, reflect="x", center_rule="mean");
#' plot(rbind(layout, layout_rot50, layout_rot40_ctrmean, layout_reflectx_ctrmean),
#'    col=rep(c("darkorchid", "darkorange1", "dodgerblue", "red4"), each=11),
#'    pch=rep(c(17, 20, 18, 17), each=11),
#'    cex=2);
#'
#' @export
rotate_coordinates <- function
(x,
 degrees=0,
 reflect=c("none",
    "x",
    "y",
    "z"),
 center=NULL,
 center_rule=c("median",
    "origin",
    "mean",
    "min",
    "max"),
 rotation_axes=c(1, 2),
 verbose=FALSE,
 ...)
{
   center_rule <- match.arg(center_rule);
   degrees <- head(degrees, 1);
   if (length(degrees) == 0) {
      degrees <- 0;
   }
   if (length(dim(x)) < 2) {
      stop("Input x must contain at least two columns.");
   }
   if (!is.matrix(x)) {
      x <- tryCatch({
         as.matrix(x);
      }, error=function(e){
         as(x, "matrix");
      });
      if (!is.matrix(x)) {
         stop("Input x must be a matrix or coercible to matrix with as.matrix(x) or as(x, 'matrix')");
      }
   }
   if ("character" %in% class(rotation_axes)) {
      rotation_axes <- intersect(jamba::rmNA(rotation_axes),
         colnames(x));
   } else {
      rotation_axes <- intersect(jamba::rmNA(rotation_axes),
         seq_len(ncol(x)));
   }
   if (length(rotation_axes) != 2) {
      stop("rotation_axes must contain two unique non-NA columns in input x");
   }
   if (verbose) {
      jamba::printDebug("rotate_coordinates(): ",
         "center_rule: ", center_rule,
         ", rotation_axes: ", rotation_axes)
   }

   ## Calculate center coordinates if not supplied
   if (length(center) == 0) {
      if ("origin" %in% center_rule) {
         center <- rep(0, ncol(x));
      } else if ("mean" %in% center_rule) {
         center <- colMeans(x, na.rm=TRUE);
      } else if ("median" %in% center_rule) {
         center <- apply(x, 2, median, na.rm=TRUE);
      } else if ("min" %in% center_rule) {
         center <- apply(x, 2, min, na.rm=TRUE);
      } else if ("max" %in% center_rule) {
         center <- apply(x, 2, max, na.rm=TRUE);
      }
   } else {
      center <- rep(center,
         length.out=ncol(x));
   }
   center_m <- matrix(rep(center, nrow(x)),
      ncol=ncol(x),
      byrow=TRUE);
   if (length(colnames(x)) > 0) {
      colnames(center_m) <- colnames(x);
   }

   ## Optionally reflect about each axis
   if (length(reflect) > 0 && !"none" %in% reflect) {
      if (any(c("character","factor") %in% class(reflect))) {
         reflect <- as.character(reflect);
         if (length(colnames(x)) == 0) {
            reflect <- unique(ifelse(reflect %in% "x", 1,
               ifelse(reflect %in% "y", 2,
                  ifelse(reflect %in% "z", 3, NA))));
            if (any(is.na(reflect))) {
               stop("reflect must contain only c('x','y','z') when there are no colnames(x).");
            }
            for (i in reflect) {
               x[,i] <- (x[,i] - center_m[,i]) * -1 + center_m[,i];
            }
         } else {
            if (!all(reflect %in% colnames(x))) {
               stop("reflect was supplied as character vector, but is not fully contained in colnames(x).");
            }
            reflect <- intersect(reflect, colnames(x));
            for (i in reflect) {
               x[,i] <- (x[,i] - center_m[,i]) * -1 + center_m[,i];
            }
         }
      } else {
         if (!any(reflect %in% seq_len(ncol(x)))) {
            stop("reflect supplied as numeric must only contain column index values in seq_len(ncol(x)).");
         }
         for (i in unique(reflect)) {
            x[,i] <- (x[,i] - center_m[,i]) * -1 + center_m[,i];
         }
      }
   }
   if (degrees == 0) {
      return(x);
   }

   ## Prepare rotation matrix math
   co <- cos(-jamba::deg2rad(degrees));
   si <- sin(-jamba::deg2rad(degrees));
   ## Apply rotation matrix math
   x <- x - center_m;
   axis1 <- rotation_axes[1];
   axis2 <- rotation_axes[2];
   new_x <- cbind(
      co * x[,axis1] - si * x[,axis2],
      si * x[,axis1] + co * x[,axis2]) + center_m;
   return(new_x);
}

