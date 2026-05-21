# methods-MemPlotFolio.R
#
# Todo:
# heatmap_row_order()
# heatmap_column_order()
# geneOrder()
# setOrder()

#' @returns `list_to_MemPlotFolio()` returns a `MemPlotFolio` S4 object,
#'    from 'list' or 'MemPlotFolio' input.
#'
#' @describeIn MemPlotFolio-class Convert legacy `list` to S4 `MemPlotFolio`
#' @param mpf `list` output from `mem_plot_folio()`
#' 
#' @examples
#' # list_to_MemPlotFolio examples
#' data(Memtest)
#' mpf <- mem_plot_folio(Memtest, do_plot=FALSE,
#'    do_which=c(1, 2, 3, 4), returnType="list")
#' Mpf <- list_to_MemPlotFolio(mpf)
#'
#' @export
list_to_MemPlotFolio <- function
(mpf)
{
   #
   if (length(mpf$gp_hm) > 0 &&
         inherits(mpf$gp_hm, "Heatmap") &&
         "draw_caption" %in% names(attributes(mpf$gp_hm))) {
      mpf$caption_function <- attributes(mpf$gp_hm)[["draw_caption"]];
   }
   Mpf <- new("MemPlotFolio",
      enrichment_hm=mpf$enrichment_hm,
      gp_hm=mpf$gp_hm,
      clusters=mpf$clusters_mem,
      gene_clusters=mpf$gene_clusters_mem,
      caption=list(
         caption=mpf$gp_hm_caption,
         caption_legendlist=mpf$gp_hm_caption_legendlist,
         caption_fn=mpf$caption_fn
      ),
      cnet_collapsed=list(
         collapsed=mpf$cnet_collapsed,
         collapsed_set=mpf$cnet_collapsed_set,
         collapsed_set2=mpf$cnet_collapsed_set2
      ),
      cnet_exemplars=mpf$cnet_exemplars,
      cnet_clusters=mpf$cnet_clusters,
      thresholds=mpf$thresholds,
      metadata=mpf$metadata
   )
   return(Mpf);
}


#' @describeIn MemPlotFolio-class Coerce S4 `MemPlotFolio` to `list` format
#' @param x `MemPlotFolio` object
#' @param ... additional arguments are ignored
#'
#' @returns `MemPlotFolio_to_list()` returns a `list`
#'
#' @export
MemPlotFolio_to_list <- function
(x,
   ...)
{
   #
   if (!inherits(x, "MemPlotFolio")) {
      stop("Input must be 'MemPlotFolio'")
   }
   
   mpf <- list(
      enrichment_hm=x@enrichment_hm,
      gp_hm=x@gp_hm,
      clusters_mem=x@clusters,
      gene_clusters_mem=x@gene_clusters,
      
      gp_hm_caption=x@caption$caption,
      gp_hm_caption_legendlist=x@caption$caption_legendlist,
      
      cnet_collapsed=x@cnet_collapsed$collapsed,
      cnet_collapsed_set=x@cnet_collapsed$collapsed_set,
      cnet_collapsed_set2=x@cnet_collapsed$collapsed_set2,
      
      cnet_exemplars=x@cnet_exemplars,
      cnet_clusters=x@cnet_clusters,
      thresholds=x@thresholds,
      metadata=x@metadata
   )
   
   # et voila
   mpf
}


#' @describeIn MemPlotFolio-class Show summary of a MemPlotFolio object,
#'    dimensions defined by genes, sets, enrichments.
setMethod("show", "MemPlotFolio",
   function(object) {
      has_plots <- character(0)
      if (inherits(object@enrichment_hm, "Heatmap"))
         has_plots <- c(has_plots, "Enrichment Heatmap");
      if (inherits(object@gp_hm, "Heatmap"))
         has_plots <- c(has_plots, "Gene-Pathway Heatmap");
      if (length(object@cnet_collapsed) > 0 &&
            inherits(object@cnet_collapsed[[1]], "igraph"))
         has_plots <- c(has_plots, "Cnet Collapsed");
      if (length(object@cnet_exemplars) > 0) {
         cnet_ex_list <- object@cnet_exemplars;
         any_igraph <- any(unlist(lapply(cnet_ex_list, function(i){
            length(i) > 0 && inherits(i, "igraph")
         })))
         if (TRUE %in% any_igraph)
            has_plots <- c(has_plots,
               paste0("Cnet Exemplars: ",
                  jamba::cPaste(names(object@cnet_exemplars), sep=", ")));
      }
      if (length(object@cnet_clusters) > 0) {
         has_plots <- c(has_plots,
            paste0("Cnet Clusters: ",
               jamba::cPaste(names(object@cnet_clusters), sep=", ")));
      }
      use_plots <- paste0("Contains Plots:\n",
         paste0(collapse="",
            paste0("  ", has_plots, "\n")))

      add_text <- character(0);
      if (is.list(object@clusters) && length(object@clusters) > 0) {
         add_text <- c(add_text,
            paste0("Clusters:\n  ",
               jamba::cPaste(names(object@clusters), sep=", "),
               "\n"));
      }
      use_add_text <- paste0(add_text, collapse="\n");
      
      if (length(object@caption$caption) > 0) {
         caption_lines <- strsplit(object@caption$caption, "\n")[[1]];
         caption_lines <- gsub("^(.+[^:])$", "  \\1", caption_lines);
         use_captions <- paste0(caption_lines,
            collapse="\n");
      } else {
         use_captions <- character(0);
      }

      all_text <- paste0(use_plots,
         use_add_text,
         use_captions);
      
      cat(all_text, sep="");
   }
)


#' @param x `MemPlotFolio` object
#' @docType methods
#' @describeIn MemPlotFolio-class Returns the order of gene sets from
#'    `MemPlotFolio` results as a `list` named by pathway cluster,
#'    containing `character` vectors of pathway gene sets.
#' @returns `Clusters(MemPlotFolio)` returns a `list` of `character` vectors,
#'    named by cluster.
#' @export
setMethod("Clusters", "MemPlotFolio", function(x) {
   x@clusters;
})


#' @param x `MemPlotFolio` object
#' @docType methods
#' @describeIn MemPlotFolio-class Returns the order of genes from
#'    `MemPlotFolio` results as a `list` named by gene cluster,
#'    containing `character` vectors of genes.
#' @returns `GeneClusters(MemPlotFolio)` returns a `list` of
#'    `character` vectors, named by cluster.
#' @export
setMethod("GeneClusters", "MemPlotFolio", function(x) {
   return(x@gene_clusters);
})


#' @param x `MemPlotFolio` object
#' @docType methods
#' @describeIn MemPlotFolio-class Returns cluster labels stored in the
#'    `metadata` slot of `MemPlotFolio`.
#' @aliases ClusterLabels
#' @returns `ClusterLabels(MemPlotFolio)` returns a named `character` vector
#'    of cluster labels, or `NULL` if none are defined.
#' @export
setMethod("ClusterLabels", "MemPlotFolio", function(x) {
   x@metadata$cluster_labels
})


#' @param x `MemPlotFolio` object
#' @param value named `character` vector of cluster labels, where names
#'    correspond to cluster names as returned by `Clusters(x)`.
#' @docType methods
#' @describeIn MemPlotFolio-class Sets cluster labels in the `metadata` slot
#'    of `MemPlotFolio`.
#' @aliases ClusterLabels
#' @returns `ClusterLabels<-(MemPlotFolio)` returns the updated `MemPlotFolio`
#'    object with new cluster labels stored in `metadata$cluster_labels`.
#' @export
setReplaceMethod("ClusterLabels", "MemPlotFolio", function(x, value) {
   # validate length
   xClusters <- Clusters(x)
   if (length(value) == 0) {
      # blank cluster_labels
      x@metadata$cluster_labels <- NULL;
   } else {
      if (!length(value) == length(xClusters)) {
         stop(paste0("length(ClusterLabels(x)) must equal ",
            "length(Clusters(x))."))
      }
      # validate names(value) is equal to names(Clusters(x))
      # it empty, assign consistent names
      if (length(names(value)) == 0) {
         names(value) <- names(xClusters);
      }
      if (!all(names(value) == names(xClusters))) {
         stop(paste0("names(ClusterLabels(x)) do not match ",
            "names(Clusters(x))."))
      }
      x@metadata$cluster_labels <- value;
   }
   validObject(x)
   x
})


#' @param x `MemPlotFolio` object
#' @docType methods
#' @describeIn MemPlotFolio-class Returns cluster data stored in the
#'    `metadata` slot of `MemPlotFolio`, intended to store further
#'    information to describe each cluster.
#' @aliases ClusterData
#' @returns `ClusterData(MemPlotFolio)` returns a named `list`
#'    of cluster data, or `NULL` if none are defined.
#' @export
setMethod("ClusterData", "MemPlotFolio", function(x) {
   x@metadata$cluster_data
})


#' @param x `MemPlotFolio` object
#' @param value named `character` vector of cluster labels, where names
#'    correspond to cluster names as returned by `Clusters(x)`.
#' @docType methods
#' @describeIn MemPlotFolio-class Sets cluster data in the `metadata` slot
#'    of `MemPlotFolio`, intended to store further information to
#'    describe each cluster. It should be the same length as `Clusters(x)`
#'    and have no names, or have names equal to `names(Clusters(x))`.
#' @aliases ClusterData
#' @returns `ClusterData<-(MemPlotFolio)` returns the updated `MemPlotFolio`
#'    object with new cluster data stored in `metadata$cluster_data`.
#' @export
setReplaceMethod("ClusterData", "MemPlotFolio", function(x, value) {
   # validate length
   xClusters <- Clusters(x)
   if (!length(value) == length(xClusters)) {
      stop(paste0("length(ClusterData(x)) must equal ",
      "length(Clusters(x))."))
   }
   # validate names(value) is equal to names(Clusters(x))
   # it empty, assign consistent names
   if (length(names(value)) == 0) {
      names(value) <- names(xClusters);
   }
   if (!all(names(value) == names(xClusters))) {
      stop(paste0("names(ClusterData(x)) do not match ",
      "names(Clusters(x))."))
   }
   x@metadata$cluster_data <- value;
   validObject(x)
   x
})


#' @param x `MemPlotFolio` object
#' @docType methods
#' @describeIn MemPlotFolio-class Returns the thresholds used with
#'    `MemPlotFolio`.
#' @aliases thresholds
#' @returns `thresholds(MemPlotFolio)` returns a `list` of
#'    thresholds used with `mem_plot_folio()`.
#' @export
setMethod("thresholds", "MemPlotFolio", function(x) {
   return(x@thresholds);
})


#' @param x `MemPlotFolio` object
#' @docType methods
#' @describeIn MemPlotFolio-class Returns the metadata used with
#'    `MemPlotFolio`.
#' @aliases metadata
#' @returns `metadata(MemPlotFolio)` returns a `list` of
#'    metadata used with `mem_plot_folio()`.
#' @export
setMethod("metadata", "MemPlotFolio", function(x) {
   return(x@metadata);
})


#' @param x `MemPlotFolio` object
#' @docType methods
#' @describeIn MemPlotFolio-class Returns the caption summary for
#'    `MemPlotFolio`.
#' @aliases Caption
#' @returns `metadata(MemPlotFolio)` returns a `character` string
#'    with caption summary used with `mem_plot_folio()`.
#'    Multiple lines are delimited by newline characters.
#' @export
setMethod("Caption", "MemPlotFolio", function(x, ...) {
   return(x@caption$caption);
})


#' @param x `MemPlotFolio` object
#' @docType methods
#' @describeIn MemPlotFolio-class Returns the caption summary for
#'    `MemPlotFolio` as `ComplexHeatmap::Legends`.
#' @aliases CaptionLegendList
#' @returns `metadata(MemPlotFolio)` returns the caption summary
#'    in the form of `ComplexHeatmap::Legends` suitable to `draw()`
#'    as R grid graphics.
#' @export
setMethod("CaptionLegendList", "MemPlotFolio", function(x, ...) {
   return(x@caption$caption_legendlist);
})


#' @param x `MemPlotFolio` object
#' @docType methods
#' @describeIn MemPlotFolio-class Draws the enrichment heatmap
#'    from `MemPlotFolio` results.  
#'    Optional arguments '...' include:
#'    * 'main' or 'column_title' for overall plot title;
#'    * 'column_title_gp' as `grid::gpar()` to adjust title font;
#'    * 'use_cluster_labels' `logical` whether to use `ClusterLabels()`
#'    as row title entries.
#' @aliases EnrichmentHeatmap
#' @returns `EnrichmentHeatmap(MemPlotFolio)` returns a
#'    `ComplexHeatmap::HeatmapList` when do_plot is TRUE (default),
#'    in addition to rendering the heatmap. It returns
#'    `ComplexHeatmap::Heatmap` when do_plot is FALSE, containing
#'    enrichment P-values by enrichment, and pathway rows in clusters.  
#'    Optional arguments '...' include:
#'    * 'main' or 'column_title' for overall plot title;
#'    * 'column_title_gp' as `grid::gpar()` to adjust title font;
#'    * 'use_cluster_labels' `logical` whether to use `ClusterLabels()`
#'    as row title entries.
#' 
#' @examples
#' data(Memtest)
#' mpf <- mem_plot_folio(Memtest, do_plot=FALSE, returnType="list")
#' Mpf <- list_to_MemPlotFolio(mpf)
#' 
#' # enrichment heatmap
#' EnrichmentHeatmap(Mpf, column_title="Enrichment Heatmap")
#' 
#' # Gene-path heatmap
#' GenePathHeatmap(Mpf, column_title="Gene-Path Heatmap")
#' 
#' # Cnet collapsed sets
#' CnetCollapsed(Mpf, type="set", use_shadowText=TRUE, main="Cnet Collapsed Sets")
#' 
#' # Cnet exemplar plot
#' CnetExemplar(Mpf, num=2, use_shadowText=TRUE, main="Cnet Exemplars, num=2")
#' 
#' # Cnet cluster plot
#' CnetCluster(Mpf, cluster="B", use_shadowText=TRUE, main="Cnet Cluster 'B'")
#' 
#' @export
setMethod("EnrichmentHeatmap", "MemPlotFolio", function(x, do_plot, ...) {
   if (missing(do_plot)) {
      do_plot <- TRUE;
   }
   arglist <- list(...);
   if ("width" %in% names(arglist)) {
      suppressWarnings(width <- jamba::rmNA(naValue=30,
         as.numeric(arglist$width)));
      arglist[["width"]] <- NULL;
   } else {
      width <- 30;
   }

   if ("use_cluster_labels" %in% names(arglist)) {
      use_cluster_labels <- arglist$use_cluster_labels;
   } else {
      use_cluster_labels <- TRUE;
   }
   if (isTRUE(use_cluster_labels) &&
      length(ClusterLabels(x)) > 0) {
      hm_row_title <- fixSetLabels(
         paste0(names(Clusters(x)), ": ",
            ClusterLabels(x)),
         width=width,
         lowercaseAll=FALSE,
         adjustCase=FALSE,
         removeGrep=NULL,
         do_abbreviations=FALSE)
      x@enrichment_hm@row_title <- hm_row_title;
   }

   if (isTRUE(do_plot)) {
      # check for off-book arguments in '...'
      column_title <- NULL;
      column_title_gp <- grid::gpar(fontsize=18);
      if ("main" %in% names(arglist)) {
         column_title <- arglist$main;
      } else if ("column_title" %in% names(arglist)) {
         column_title <- arglist$column_title;
      }
      if ("column_title_gp" %in% names(arglist)) {
         column_title_gp <- arglist$column_title_gp;
      }

      annotation_legend_list <- c(
         attr(x@enrichment_hm, "annotation_legend_list"),
         x@caption$caption_legendlist);
      ComplexHeatmap::draw(x@enrichment_hm,
         newpage=mem_do_newpage(),
         annotation_legend_list=annotation_legend_list,
         column_title=column_title,
         column_title_gp=column_title_gp,
         merge_legends=TRUE)
      # return invisibly?
      return(invisible(x@enrichment_hm));
   } else {
      x@enrichment_hm;
   }
})


#' @param x `MemPlotFolio` object
#' @docType methods
#' @describeIn MemPlotFolio-class Draws the gene-pathway set heatmap
#'    from `MemPlotFolio` results.
#' @aliases GenePathHeatmap
#' @returns `GenePathHeatmap(MemPlotFolio)` returns a
#'    `ComplexHeatmap::HeatmapList` when do_plot is TRUE (default),
#'    `ComplexHeatmap::Heatmap` when do_plot is FALSE, containing
#'    the genes-pathways incidence matrix and associated caption.
#' @export
setMethod("GenePathHeatmap", "MemPlotFolio", function(x, do_plot, ...) {
   if (missing(do_plot)) {
      do_plot <- TRUE;
   }
   if (isTRUE(do_plot)) {
      # check for off-book arguments in '...'
      arglist <- list(...);
      use_thresholds <- thresholds(x);
      use_ht_opts <- list();
      if ("column_anno_padding" %in% names(arglist)) {
         use_ht_opts$COLUMN_ANNO_PADDING <- arglist$column_anno_padding;
      } else if ("column_anno_padding" %in% names(use_thresholds)) {
         use_ht_opts$COLUMN_ANNO_PADDING <- use_thresholds$column_anno_padding;
      }
      if ("row_anno_padding" %in% names(arglist)) {
         use_ht_opts$ROW_ANNO_PADDING <- arglist$row_anno_padding;
      } else if ("row_anno_padding" %in% names(use_thresholds)) {
         use_ht_opts$ROW_ANNO_PADDING <- use_thresholds$row_anno_padding;
      }
      if (length(use_ht_opts) > 0) {
         local_ht_opts(use_ht_opts)
      }
      # check for known attributes
      title_list <- attr(x@gp_hm, "title_list");
      
      column_title <- title_list$column_title;
      column_title_gp <- title_list$column_title_gp;
      if ("main" %in% names(arglist)) {
         column_title <- arglist$main;
      } else if ("column_title" %in% names(arglist)) {
         column_title <- arglist$column_title;
      }
      if ("column_title_gp" %in% names(arglist)) {
         column_title_gp <- arglist$column_title_gp;
      }
      
      caption_legendlist <- x@caption$caption_legendlist;
      ComplexHeatmap::draw(x@gp_hm,
         newpage=mem_do_newpage(),
         annotation_legend_list=caption_legendlist,
         column_title=column_title,
         column_title_gp=column_title_gp,
         merge_legends=TRUE)
   } else {
      x@gp_hm;
   }
})


# Todo: Use only cnet_collapsed and change V(cnet)$label by 'type'

#' @param x `MemPlotFolio` object
#' @docType methods
#' @describeIn MemPlotFolio-class Draws the Cnet collapsed
#'    network from `MemPlotFolio` results. Note that '...' arguments are
#'    passed to `jam_igraph()` and `mem_legend()` when `do_plot=TRUE`.
#'    When 'type' is not defined: it uses  type='cluster' if
#'    `ClusterLabels()` is defined, otherwise type='set';
#'    then it appends '2' to the end if there are more than 1000 nodes.
#'    Argument 'type' can be:
#'    * `type='title'` or `type=''` to use cluster title
#'    * `type='set'` to use abbreviated pathway names
#'    * `type='cluster'` to use cluster labels
#'    * `type='attribute'` to use any existing vertex attribute.
#'    * Append **'2'** to the end to hide gene labels. For example `type='set2'`.
#' 
#'    Recognized arguments in '...':
#'    * `'width'` to apply word-wrap to cluster labels, or vertex attribute.
#'    * `'maxNchar'` to set maximum string character length, passed to
#'    `fixSetLabels()`.
#'    * `'layout'` passed to `set_igraph_layout(cnet, layout)` to update
#'    the node layout.
#'    * `'rotate_degrees'` with `numeric` value. When set, it calls
#'    `rotate_igraph_layout()` with defaults.
#'    * All other '...' arguments are passed to `jam_igraph()`.
#'    
#'    The legend includes direction if encoded in the 'MemPlotFolio', but
#'    can be forced with `do_directional=TRUE` or FALSE.
#' @aliases CnetCollapsed
#' @returns `CnetCollapsed(MemPlotFolio)` returns an `igraph` object invisibly,
#'    with Gene and Set nodes representing the collapsed pathway clusters.
#'    The legend position can be adjusted using 'legend_x', 'legend_y' which
#'    are passed to `mem_legend()` as 'x' and 'y'.
#' @export
setMethod("CnetCollapsed", "MemPlotFolio",
   function(x, type, do_plot, legend_x="bottomleft", legend_y=NULL, ...) {
   if (missing(do_plot)) {
      do_plot <- TRUE;
   }
   # type indicates which cnet data to use
   cnet <- x@cnet_collapsed$collapsed;
   if (!inherits(cnet, "igraph")) {
      stop_msg <- paste0("No Cnet data were prepared by collapsed sets, ",
         "use:\nmem_plot_folio(Mem, do_which=3)");
      stop(stop_msg);
   }
   arglist <- list(...);
   if ("width" %in% names(arglist)) {
      suppressWarnings(width <- jamba::rmNA(naValue=30,
         as.numeric(arglist$width)));
      arglist[["width"]] <- NULL;
   } else {
      width <- 30;
   }
   maxNchar <- Inf;
   if ("maxNchar" %in% names(arglist)) {
      maxNchar <- arglist$maxNchar;
      arglist[["maxNchar"]] <- NULL;
   }

   if ("layout" %in% names(arglist)) {
      cnet <- set_igraph_layout(g=cnet,
         layout=arglist$layout)
      arglist[["layout"]] <- NULL;
   }
   if ("rotate_degrees" %in% names(arglist)) {
      cnet <- rotate_igraph_layout(cnet,
         degrees=arglist$rotate_degrees)
      arglist[["rotate_degrees"]] <- NULL;
   }

   if (missing(type)) {
      type <- "title";
      if (length(ClusterLabels(Mpf)) > 0) {
         type <- "cluster";
      }
      if (igraph::vcount(cnet) > 1000) {
         type <- paste0(type, "2");
      }
   }
   if (length(type) == 0) {
      type <- "";
   }
   type <- head(type, 1);
   if (type %in% c(NA, "") || grepl("^title[ \t\n2]*$", type)) {
      if (!"label" %in% igraph::vertex_attr_names(cnet)) {
         igraph::V(cnet)$label <- igraph::V(cnet)$name;
      }
   } else if (any(grepl("^(set|cluster)[ \t\n2]*$", type))) {
      # add set_names to labels
      if (any(grepl("^cluster", type))) {
         # add cluster_labels to labels
         cluster_labels <- ClusterLabels(x);
         if (length(cluster_labels) == 0) {
            ## no cluster_labels, what to do? Error, inaction, warning?
            # return(invisible(NULL))
            ## for testing use "Cluster A", "Cluster B", etc.
            ## In future use "A", "B", etc.
            cluster_labels <- jamba::nameVector(
               # names(x@clusters),
               paste("Cluster", names(x@clusters)),
               names(x@clusters))
         }
         # match igraph Set nodes to names(cluster_labels)
         isset <- which(igraph::V(cnet)$nodeType %in% "Set");
         clmatch <- match(igraph::V(cnet)$name[isset], names(cluster_labels))
         if (!"label" %in% igraph::vertex_attr_names(cnet)) {
            igraph::V(cnet)$label <- igraph::V(cnet)$name;
         }
         # any non-NA match gets updated
         # NA match re-uses the original node name, there is no label
         igraph::V(cnet)$label[isset] <- ifelse(!is.na(clmatch),
            cluster_labels[clmatch],
            igraph::V(cnet)$label[isset])
      } else if ("set_labels" %in% igraph::vertex_attr_names(cnet)) {
         # add set_names to labels
         use_labels <- ifelse(
            nchar(jamba::rmNA(naValue="", igraph::V(cnet)$set_labels)) > 0,
            igraph::V(cnet)$set_labels,
            igraph::V(cnet)$name);
         igraph::V(cnet)$label <- use_labels;
      }
   } else {
      # check if type matches vertex attribute names
      test_type <- gsub("[ \n]*2$", "", type);
      if (!test_type %in% igraph::vertex_attr_names(cnet)) {
         stop_msg <- paste0(
            "The 'type' did not contain 'title', 'set', 'cluster', ",
            "nor any vertex attribute names.")
         stop(stop_msg);   
      }
      use_type <- type;
      if (!type %in% igraph::vertex_attr_names(cnet)) {
         use_type <- test_type;
      }
      use_labels <- ifelse(
         nchar(jamba::rmNA(naValue="",
            igraph::vertex_attr(cnet, name=use_type))) > 0,
         igraph::vertex_attr(cnet, name=use_type),
         igraph::V(cnet)$name);
      # apply word wrap
      if (length(width) > 0) {
         use_labels <- fixSetLabels(
            x=use_labels,
            width=width,
            lowercaseAll=FALSE,
            adjustCase=FALSE,
            removeGrep=NULL,
            do_abbreviations=FALSE,
            maxNchar=maxNchar)
      }
      igraph::V(cnet)$label <- use_labels;
   }
   # hide gene labels when type ends with '2'
   if (grepl("2$", type)) {
      isgene <- jamba::igrep("^gene$",
         igraph::vertex_attr(cnet, "nodeType"));
      if (length(isgene) > 0) {
         igraph::vertex_attr(cnet, index=isgene, name="label") <- "";
      }
   }
   
   if (isTRUE(do_plot)) {
      do.call(jam_igraph,
         c(
            alist(x=cnet),
            arglist))
   	
   	# determine whether to include direction in the legend
      arglist$x <- legend_x;
      arglist$y <- legend_y;
   	if (length(arglist) > 0 && "do_directional" %in% names(arglist)) {
         do_directional <- arglist[["do_directional"]];
   		# arglist <- arglist[-match("do_directional", names(arglsit))];
      } else {
         hasDirection <- ifelse(isTRUE(metadata(x)[["hasDirection"]]),
            TRUE, FALSE)
   		arglist[["do_directional"]] <- hasDirection;
   	}
   	# draw the legend
   	do.call(mem_legend,
   		c(
   			alist(mem=x),
   			arglist))
   }
   invisible(cnet);
})


#' @param x `MemPlotFolio` object
#' @docType methods
#' @describeIn MemPlotFolio-class Draws the Cnet exemplar network
#'    from `MemPlotFolio` results for 'num' exemplar per cluster.
#'    Note that '...' arguments are
#'    passed to `jam_igraph()` and `mem_legend()` when `do_plot=TRUE`.
#'    
#'    The legend includes direction if encoded in the 'MemPlotFolio', but
#'    can be forced with `do_directional=TRUE` or FALSE. Legend position
#'    is controlled by legend_x,legend_y, default legend_x="bottomleft".
#' @aliases CnetExemplar
#' @returns `CnetExemplar(MemPlotFolio)` returns an `igraph` object
#'    with Gene and Set nodes for the 'num' number of exemplars per cluster.
#' @export
setMethod("CnetExemplar", "MemPlotFolio",
   function(x, num, do_plot, legend_x="bottomleft", legend_y=NULL, main=NULL, ...) {
   if (missing(do_plot)) {
      do_plot <- TRUE;
   }
   # type indicates which cnet data to use
   if (missing(num)) {
      cnet <- x@cnet_exemplars[[1]];
      if (!inherits(cnet, "igraph")) {
         # Todo: prepare dynamically
         stop("No Cnet exemplar networks were prepared.")
      }
   } else {
      num <- as.character(head(num, 1));
      if (!num %in% names(x@cnet_exemplars)) {
         stop_msg <- paste0("No Cnet exemplar network with ",
            "exemplar number ", num, ".");
         stop(stop_msg);
      }
      cnet <- x@cnet_exemplars[[num]];
   }

   if (isTRUE(do_plot)) {
      jam_igraph(cnet,
         ...)
   	
   	# determine whether to include direction in the legend
   	hasDirection <- ifelse(isTRUE(metadata(x)[["hasDirection"]]),
   		TRUE, FALSE)
   	arglist <- list(...);
      arglist$x <- legend_x;
      arglist$y <- legend_y;
   	if (length(arglist) > 0 && "do_directional" %in% names(arglist)) {
   		do_directional <- arglist[["do_directional"]];
   		# arglist <- arglist[-match("do_directional", names(arglsit))];
   	} else {
   		arglist[["do_directional"]] <- hasDirection;
   	}
   	# draw the legend
   	do.call(mem_legend,
   		c(
   			alist(mem=metadata(x)),
   			arglist))
   }
   invisible(cnet);
   })


#' @param x `MemPlotFolio` object
#' @docType methods
#' @describeIn MemPlotFolio-class Draws the Cnet network
#'    from `MemPlotFolio` results for a specific cluster.
#'    Note that '...' arguments are
#'    passed to `jam_igraph()` and `mem_legend()` when `do_plot=TRUE`.
#'    
#'    The legend includes direction if encoded in the 'MemPlotFolio', but
#'    can be forced with `do_directional=TRUE` or FALSE.
#'    Plot title 'main', default NULL uses the cluster label if present,
#'    then cluster name. Use main=FALSE or main="" to hide the title.
#' @aliases CnetCluster
#' @returns `CnetCluster(MemPlotFolio)` returns an `igraph` object invisibly,
#'    using all pathways in the cluster defined with argument `cluster`.
#' @export
setMethod("CnetCluster", "MemPlotFolio",
   function(x, cluster, do_plot, legend_x="bottomleft", legend_y=NULL, main=NULL, ...) {
   if (missing(do_plot)) {
      do_plot <- TRUE;
   }
   if (length(main) == 0) {
      if (length(x@metadata$cluster_labels) > 0) {
         main <- x@metadata$cluster_labels[cluster];
      } else {
         main <- names(Clusters(x)[cluster]);
      }
   }
      if (isFALSE(main)) {
         main <- NULL;
      }

   # type indicates which cnet data to use
   if (missing(cluster) || length(cluster) == 0) {
      cluster <- head(Clusters(x), 1);
      if (length(cluster) == 0) {
         stop("No clusters are defined in the MemPlotFolio provided.")
      }
   }
   
   # Todo: Consider calling mem_plot_folio() to prepare when not present
   # - mem_plot_folio(Mpf) would require the MemPlotFolio to store Mem also
   #   therefore it seems inappropriate.
   if (is.numeric(cluster)) {
      cluster <- head(names(Clusters(x))[as.integer(cluster)], 1);
      if (length(cluster) == 0) {
         stop("Argument 'cluster' as integer did not match Clusters(x).");
      }
   }
   if (!cluster %in% names(x@cnet_clusters)) {
      # cluster Cnet does not exist, describe steps to create
      stop_msg <- paste0("Cluster '", cluster, "' is not present. ",
         "Create:\n",
         "Mpf <- mem_plot_folio(Mem, do_which=7:(7 + n))\n",
         "where 'n' is total number of clusters. Then plot:\n",
         "CnetCluster(Mpf, cluster=cluster)");
      stop(stop_msg);
   }
   
   # obtain the Cnet data
   cluster <- as.character(head(cluster, 1));
   if (!cluster %in% names(x@cnet_clusters)) {
      stop_msg <- paste0("No Cnet cluster network prepared for ",
         cluster, ".");
      stop(stop_msg);
   }
   cnet <- x@cnet_clusters[[cluster]];

   if (isTRUE(do_plot)) {
   	# draw the igraph
      jam_igraph(cnet,
         main=main,
         ...)
   	
   	# determine whether to include direction in the legend
   	hasDirection <- ifelse(isTRUE(metadata(x)[["hasDirection"]]),
   		TRUE, FALSE)
   	arglist <- list(...);
      arglist$x <- legend_x;
      arglist$y <- legend_y;
   	if (length(arglist) > 0 && "do_directional" %in% names(arglist)) {
   		do_directional <- arglist[["do_directional"]];
   		# arglist <- arglist[-match("do_directional", names(arglsit))];
   	} else {
   		arglist[["do_directional"]] <- hasDirection;
   	}
   	# draw the legend
   	do.call(mem_legend,
   		c(
   			alist(mem=metadata(x)),
	      	arglist))
   }
   invisible(cnet)
   })

#' @describeIn MemPlotFolio-class Plot a `MemPlotFolio` object calling
#'    `plot_mpf()`
#' @param x `MemPlotFolio` object
#' @param y ignored
#' @param ... additional arguments passed to `plot_mpf()`
#' @export
setMethod("plot", signature(x = "MemPlotFolio"), function(x, y, ...) {
   plot_mpf(x, y, ...)
})
