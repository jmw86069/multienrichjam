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
#' @describeIn MemPlotFolio-class Returns the thresholds used with
#'    `MemPlotFolio`.
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
#' @returns `EnrichmentHeatmap(MemPlotFolio)` returns a
#'    `ComplexHeatmap::HeatmapList` when do_plot is TRUE (default),
#'    `ComplexHeatmap::Heatmap` when do_plot is FALSE, containing
#'    enrichment P-values by enrichment, and pathway rows in clusters.
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
   if (isTRUE(do_plot)) {
      # check for off-book arguments in '...'
      arglist <- list(...);
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
         annotation_legend_list=annotation_legend_list,
         column_title=column_title,
         column_title_gp=column_title_gp,
         merge_legends=TRUE)
   } else {
      x@enrichment_hm;
   }
})


#' @param x `MemPlotFolio` object
#' @docType methods
#' @describeIn MemPlotFolio-class Draws the gene-pathway set heatmap
#'    from `MemPlotFolio` results.
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
#'    Argument `type` can be:
#'    * `type=''` (default) to use cluster title
#'    * `type='set'` to use abbreviated pathway names
#'    * `type='set2'` to use abbreviated pathway names, with hidden gene labels
#' @returns `CnetCollapsed(MemPlotFolio)` returns an `igraph` object invisibly,
#'    with Gene and Set nodes representing the collapsed pathway clusters.
#' @export
setMethod("CnetCollapsed", "MemPlotFolio", function(x, type, do_plot, ...) {
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
   if (missing(type)) {
      type <- "";
   }
   type <- head(type, 1);
   if (type %in% c(NA, "")) {
      # no changes required
   } else if (type %in% c("set", "set2")) {
      # add set_names to labels
      if ("set_labels" %in% igraph::vertex_attr_names(cnet)) {
         use_labels <- ifelse(
            nchar(jamba::rmNA(naValue="", igraph::V(cnet)$set_labels)) > 0,
            igraph::V(cnet)$set_labels,
            igraph::V(cnet)$name);
         igraph::V(cnet)$label <- use_labels;
      }
      if ("set2" %in% type) {
         # remove gene labels
         isgene <- which(tolower(igraph::vertex_attr(cnet, "nodeType")) %in%
               "gene");
         igraph::vertex_attr(cnet, index=isgene, name="label") <- "";
      }
   } else if (type %in% igraph::vertex_attr_names(cnet)) {
      use_labels <- ifelse(
         nchar(jamba::rmNA(naValue="",
            igraph::vertex_attr(cnet, name=type))) > 0,
         igraph::vertex_attr(cnet, name=type),
         igraph::V(cnet)$name);
   } else {
      stop_msg <- paste0("The 'type' did not match '', 'set', 'set2', nor ",
         "any vertex attribute names.")
      stop(stop_msg);
   }
   
   if (isTRUE(do_plot)) {
      jam_igraph(cnet,
         ...)
      mem_legend(metadata(x),
         ...)
   }
   invisible(cnet);
})


#' @param x `MemPlotFolio` object
#' @docType methods
#' @describeIn MemPlotFolio-class Draws the Cnet exemplar network
#'    from `MemPlotFolio` results for 'num' exemplar per cluster.
#'    Note that '...' arguments are
#'    passed to `jam_igraph()` and `mem_legend()` when `do_plot=TRUE`.
#' @returns `CnetExemplar(MemPlotFolio)` returns an `igraph` object
#'    with Gene and Set nodes for the 'num' number of exemplars per cluster.
#' @export
setMethod("CnetExemplar", "MemPlotFolio", function(x, num, do_plot, ...) {
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
      mem_legend(metadata(x),
         ...)
   }
   invisible(cnet);
})


#' @param x `MemPlotFolio` object
#' @docType methods
#' @describeIn MemPlotFolio-class Draws the Cnet network
#'    from `MemPlotFolio` results for a specific cluster.
#'    Note that '...' arguments are
#'    passed to `jam_igraph()` and `mem_legend()` when `do_plot=TRUE`.
#' @returns `CnetCluster(MemPlotFolio)` returns an `igraph` object invisibly,
#'    using all pathways in the cluster defined with argument `cluster`.
#' @export
setMethod("CnetCluster", "MemPlotFolio", function(x, cluster, do_plot, ...) {
   if (missing(do_plot)) {
      do_plot <- TRUE;
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
      jam_igraph(cnet,
         ...)
      mem_legend(metadata(x),
         ...)
      
   }
   invisible(cnet)
})
