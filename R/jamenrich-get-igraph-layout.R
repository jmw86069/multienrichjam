
#' Obtain or create layout for igraph object
#'
#' Obtain or create layout for igraph object
#'
#' This function is a simple helper function intended to retrieve
#' the node layout for an `igraph` object. When `layout` is supplied
#' as a `function` it is used to define a specific layout matrix.
#'
#' @family jam utility functions
#' 
#'
#' @return `numeric` matrix of x,y coordinates with nrow equal
#'    to the number of nodes in the input, from `igraph::vcount(g)`.
#'    Note that `rownames()` are defined to match node names
#'    `igraph::V(g)$name`, unlike default `igraph` layouts.
#'
#' @param g `igraph` object
#' @param layout one of:
#'    * `numeric` matrix of layout coordinates, with `nrow(layout)`
#'    equal to the number of nodes `igraph::vcount(g)`.
#'    * `function` that takes `igraph` input, and returns `numeric` matrix
#'    of layout coordinates.
#'    * `NULL`: to use the layout embedded in `g` using
#'    `igraph::graph_attr(g, "layout")` if it exists. If it does not exist,
#'    then `make_circular=TRUE` causes a new layout to be created with
#'    arbitrary circular coordinates; or `make_circular=FALSE` will return
#'    `NULL`.
#' @param make_circular `logical` indicating whether to create a makeshift
#'    layout when the input data does not already contain layout coordinates,
#'    and when `layout` is not supplied.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
get_igraph_layout <- function
(g,
 layout=NULL,
 make_circular=TRUE,
 verbose=FALSE,
 ...)
{
   #
   # obtain layout
   xy <- NULL;
   if (length(layout) == 0) {
      # igraph layout
      if ("layout" %in% igraph::graph_attr_names(g) &&
            length(igraph::graph_attr(g, "layout")) > 0) {
         if (verbose) {
            jamba::printDebug("get_igraph_layout(): ",
               "taking graph_attr layout.");
         }
         xy <- igraph::graph_attr(g, "layout");
      } else {
         # makeshift circular layout
         if (make_circular) {
            if (verbose) {
               jamba::printDebug("get_igraph_layout(): ",
                  "creating initial circular layout.");
            }
            n <- igraph::vcount(g);
            xy <- matrix(0, nrow=n, ncol=2);
            xy[,1] <- head(sin(seq(0, 2*pi, length.out=n + 1)) * n * 0.8, n);
            xy[,2] <- head(cos(seq(0, 2*pi, length.out=n + 1)) * n * 0.8, n);
         } else {
            if (verbose) {
               jamba::printDebug("get_igraph_layout(): ",
                  "no layout is available, returning NULL.");
            }
            return(NULL);
         }
      }
      rownames(xy) <- igraph::V(g)$name;
   } else if ("function" %in% class(layout)) {
      #
      if (verbose) {
         jamba::printDebug("get_igraph_layout(): ",
            "creating layout using supplied ",
            "layout()");
      }
      xy <- layout(g, ...)
   } else {
      xy <- layout;
      if (!"matrix" %in% class(xy)) {
         # try two common methods to coerce to matrix
         xy <- tryCatch({
            as.matrix(xy);
         }, error=function(e){
            as(xy, "matrix");
         });
      }
   }

   # validate coordinates versus input g
   if (nrow(xy) != igraph::vcount(g)) {
      stop(paste0(
         "layout nrow (", nrow(xy),
         ") does not equal vcount (", igraph::vcount(g), ")."))
   }
   if (length(rownames(xy)) > 0) {
      if (!all(rownames(xy) %in% igraph::V(g)$name)) {
         stop("rownames(layout) do not match V(g)$name.")
      }
      xy <- xy[match(igraph::V(g)$name, rownames(xy)), , drop=FALSE];
   } else {
      rownames(xy) <- igraph::V(g)$name;
   }

   if (verbose) {
      jamba::printDebug("get_igraph_layout(): ",
         "head(xy, 3):");
      print(head(xy, 3));
   }
   return(xy);
}
