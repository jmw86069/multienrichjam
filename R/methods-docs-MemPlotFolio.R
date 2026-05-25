# methods-docs-MemPlotFolio.R

#' Cnet-methods for MemPlotFolio
#' 
#' @description
#' Methods to create Concept Network (Cnet) plots
#' from a `MemPlotFolio` object.
#' 
#' @details
#' The following methods are defined for `MemPlotFolio`:
#' 
#' \describe{
#'    \item{`CnetCollapsed(x, type, do_plot, legend_x, legend_y)`}{
#'       Returns `igraph` Cnet object with pathways collapsed
#'       using `Clusters()`.  
#'       Argument `type` controls the pathway cluster label:
#'       'title' uses the cluster title; 'set' uses the abbreviated
#'       set names; 'cluster' uses `ClusterLabels()` if
#'       defined.  
#'       Append **`'2'`** suffix to hide gene labels, for example
#'       `type='title2'` will display pathway cluster title, and
#'       hide gene labels.
#'    }
#'    \item{`CnetExemplar(x, num, do_plot, legend_x, legend_y, main)`}{
#'       Returns `igraph` Cnet object with exemplar pathways
#'       for each cluster in `Clusters()`.  
#'       Argument 'num' defines the number of exemplar pathways
#'       used for each cluster. It can be a single value, or
#'       vector of values applied in order of `Clusters()`.
#'    }
#'    \item{`CnetCluster(x, cluster, do_plot, legend_x, legend_y, main)`}{
#'       Returns `igraph` Cnet object foran individual cluster
#'       defined in `Clusters()`.  
#'       Argument 'cluster' defines the cluster, and must match
#'       an entry in `Clusters()` or be an `integer` index to
#'       the cluster.
#'    }
#' }
#' 
#' @param x `MemPlotFolio` object
#' @param type `character` used for `CnetCollapsed()`
#' @param num `integer` number of exemplar pathways, used with
#'    `CnetExemplar()`
#' @param cluster `character` name of pathway cluster which
#'    matches an entry in `Clusters(x)`.
#' @param do_plot `logical` whether to render the igraph using
#'    `jam_igraph()`.
#' @param legend_x,legend_y passed to `mem_legend()`. The default
#'    `legend_x='bottomleft'`. Use `legend_x=FALSE` to hide the
#'    legend.
#' @param main `character` with plot title.
#'    * For `CnetCluster()` if `ClusterLabels(x)` is defined, it
#'    will use the corresponding label as plot title, otherwise
#'    if it uses `names(Clusters(x))`.
#'    To hide the title, use main=FALSE.
#' @param ... additional arguments may be recognized.
#' 
#' @returns `igraph` Cnet object as defined for each method above.
#' 
#' @name Cnet-methods
#' @docType methods
NULL
