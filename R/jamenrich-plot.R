
#' MultiEnrichment Heatmap of Genes and Pathways
#'
#' MultiEnrichment Heatmap of Genes and Pathways
#'
#' @param mem `list` object created by `multiEnrichMap()`. Specifically
#'    the object is expected to contain `colorV`, `enrichIM`,
#'    `memIM`, `geneIM`.
#' @param genes character vector of genes to include in the heatmap,
#'    all other genes will be excluded.
#' @param sets character vector of sets (pathways) to include in the heatmap,
#'    all other sets will be excluded.
#' @param min_gene_ct minimum number of occurrences of each gene
#'    across the pathways, all other genes are excluded.
#' @param min_set_ct minimum number of genes required for each set,
#'    all other sets are excluded.
#' @param column_fontsize,row_fontsize arguments passed to
#'    `ComplexHeatmap::Heatmap()` to size column and row labels.
#' @param row_method,column_method character string of the distance method
#'    to use for row and column clustering. The clustering is performed
#'    by `amap::hcluster()`.
#' @param name character value passed to `ComplexHeatmap::Heatmap()`,
#'    used as a label above the heatmap color legend.
#' @param p_cutoff numeric value of the enrichment P-value cutoff,
#'    above which P-values are not colored, and are therefore white.
#'    The enrichment P-values are displayed as an annotated heatmap
#'    at the top of the main heatmap. Any cell that has a color meets
#'    at least the minimum P-value threshold.
#' @param column_title optional character string with title to display
#'    above the heatmap.
#' @param ... additional arguments are passed to `ComplexHeatmap::Heatmap()`
#'    for customization.
#'
#' @export
mem_gene_path_heatmap <- function
(mem,
 genes=NULL,
 sets=NULL,
 min_gene_ct=1,
 min_set_ct=1,
 column_fontsize=6,
 row_fontsize=8,
 row_method="spearman",
 column_method="spearman",
 cluster_columns=NULL,
 cluster_rows=NULL,
 name="gene_ct",
 p_cutoff=0.05,
 column_title=NULL,
 ...)
{
   #
   if (length(colnames(mem$memIM)) == 0) {
      colnames(mem$memIM) <- rownames(mem$enrichIM);
   }
   memIM <- mem$memIM;
   if (length(genes) > 0) {
      memIM <- memIM[intersect(rownames(memIM), genes),,drop=FALSE];
   }
   if (length(sets) > 0) {
      memIM <- memIM[,intersect(colnames(memIM), sets),drop=FALSE];
   }
   if (any(dim(memIM) == 0)) {
      stop("No remaining data after filtering.");
   }
   memIMgenect <- rowSums(memIM > 0);
   genes <- rownames(memIM)[memIMgenect >= min_gene_ct];
   memIMsetct <- colSums(memIM > 0);
   sets <- colnames(memIM)[memIMsetct >= min_set_ct];
   memIM <- memIM[genes,sets,drop=FALSE];
   if (any(dim(memIM) == 0)) {
      stop("No remaining data after filtering.");
   }
   col_iml1 <- lapply(nameVectorN(mem$colorV), function(i){
      j <- colorV[[i]];
      circlize::colorRamp2(breaks=c(0,1),
         colors=jamba::getColorRamp(j, n=3)[1:2])
   });
   col_iml4 <- lapply(nameVectorN(mem$colorV), function(i){
      j <- colorV[[i]];
      circlize::colorRamp2(
         #breaks=c(0,4),
         breaks=c(-log10(p_cutoff+1e-5), -log10(p_cutoff), 4, 6),
         colors=c("white", jamba::getColorRamp(j, trimRamp=c(4,2), n=3, lens=3)))
   });
   ## Cluster columns and rows
   if (length(cluster_columns) == 0) {
      cluster_columns <- amap::hcluster(
         link="ward",
         cbind(
            noiseFloor(
               -log10(mem$enrichIM[sets,,drop=FALSE]),
               minimum=-log10(p_cutoff+1e-5),
               newValue=0,
               ceiling=3),
            t(mem$memIM[genes,sets,drop=FALSE])),
         method=column_method);
   }
   if (length(cluster_rows) == 0) {
      cluster_rows <- amap::hcluster(
         link="ward",
         cbind(
            (mem$geneIM[genes,,drop=FALSE]),
            (mem$memIM[genes,sets,drop=FALSE])),
         method=row_method);
   }
   ComplexHeatmap::Heatmap(memIM[genes,sets,drop=FALSE],
      border=TRUE,
      name=name,
      cluster_columns=cluster_columns,
      cluster_rows=cluster_rows,
      clustering_distance_columns=column_method,
      clustering_distance_rows=row_method,
      top_annotation=ComplexHeatmap::HeatmapAnnotation(
         which="column",
         border=TRUE,
         #gp=gpar(col="black"),
         col=col_iml4,
         df=-log10(mem$enrichIM[sets,,drop=FALSE])),
      col=getColorRamp("Reds", lens=2),
      left_annotation=ComplexHeatmap::HeatmapAnnotation(which="row",
         border=TRUE,
         col=col_iml1,
         #gp=gpar(col="black"),
         df=mem$geneIM[genes,,drop=FALSE]),
      row_names_gp=gpar(fontsize=row_fontsize),
      column_names_gp=gpar(fontsize=column_fontsize),
      column_names_rot=90,
      column_title=column_title,
      ...)
}

#' MultiEnrichment Heatmap of enrichment P-values
#'
#' MultiEnrichment Heatmap of enrichment P-values
#'
#' This function is a lightweight wrapper to `ComplexHeatmap::Heatmap()`
#' intended to visualize the enrichment P-values from multiple
#' enrichment results. The P-value threshold is used to colorize
#' every cell whose P-value meets the threshold, while all other
#' cells are therefore white.
#'
#' @family jam plot functions
#'
#' @param mem `list` object created by `multiEnrichMap()`. Specifically
#'    the object is expected to contain `enrichIM`.
#' @param p_cutoff numeric value of the enrichment P-value cutoff,
#'    above which P-values are not colored, and are therefore white.
#'    The enrichment P-values are displayed as an annotated heatmap
#'    at the top of the main heatmap. Any cell that has a color meets
#'    at least the minimum P-value threshold.
#' @param p_floor numeric minimum P-value used for the color gradient.
#'    P-values below this floor are colored with the maximum color gradient.
#' @param row_method character string of the distance method
#'    to use for row and column clustering. The clustering is performed
#'    by `amap::hcluster()`.
#' @param name character value passed to `ComplexHeatmap::Heatmap()`,
#'    used as a label above the heatmap color legend.
#' @param row_dend_reorder logical indicating whether to reorder the
#'    row dendrogram using the method described in
#'    `ComplexHeatmap::Heatmap()`. The end result is minor reshuffling of
#'    leaf nodes on the dendrogram based upon mean signal in each row,
#'    which can sometimes look cleaner.
#' @param row_fontsize,column_fontsize arguments passed to
#'    `ComplexHeatmap::Heatmap()` to size row and column labels.
#' @param cluster_columns logical indicating whether to cluster heatmap
#'    columns, by default columns are not clustered.
#' @param sets character vector of sets (pathways) to include in the heatmap,
#'    all other sets will be excluded.
#' @param color_by_column logical indicating whether to colorize the
#'    heatmap using `mem$colorV` as defined for each comparison. This
#'    option is currently experimental, and produces a base R heatmap
#'    using `jamba::imageByColors()`.
#' @param cex.axis sent to `jamba::imageByColors()` only when
#'    `color_by_column=TRUE`.
#' @param column_title optional character string with title to display
#'    above the heatmap.
#' @param ... additional arguments are passed to `ComplexHeatmap::Heatmap()`
#'    for customization.
#'
#' @export
mem_enrichment_heatmap <- function
(mem,
 p_cutoff=0.05,
 p_floor=1e-6,
 row_method="euclidean",
 name="-log10P",
 row_dend_reorder=TRUE,
 row_fontsize=8,
 column_fontsize=12,
 cluster_columns=FALSE,
 sets=NULL,
 color_by_column=FALSE,
 cex.axis=1,
 column_title=NULL,
 ...)
{
   #
   col_logp <- circlize::colorRamp2(
      breaks=c(-log10(p_cutoff+1e-5),
         seq(from=-log10(p_cutoff), to=-log10(p_floor), length.out=25)),
      colors=c("white",
         getColorRamp("Reds", lens=4, n=25, trimRamp=c(2,0)))
   )
   if (length(sets) > 0) {
      mem$enrichIM <- mem$enrichIM[intersect(rownames(mem$enrichIM), sets),,drop=FALSE];
   } else {
      sets <- rownames(mem$enrichIM);
   }
   if (any(dim(mem$enrichIM) == 0)) {
      stop("No remaining data after filtering.");
   }
   er_hc2 <- amap::hcluster(
      link="ward",
      noiseFloor(
         -log10(mem$enrichIM[sets,,drop=FALSE]),
         minimum=-log10(p_cutoff+1e-5),
         newValue=0,
         ceiling=3),
      method=row_method);
   hm <- ComplexHeatmap::Heatmap(
      -log10(mem$enrichIM),
      name=name,
      col=col_logp,
      cluster_rows=er_hc2,
      row_dend_reorder=row_dend_reorder,
      border=TRUE,
      row_names_gp=gpar(fontsize=row_fontsize),
      column_names_gp=gpar(fontsize=column_fontsize),
      cluster_columns=cluster_columns,
      row_dend_width=grid::unit(30, "mm"),
      row_names_max_width=grid::unit(8, "cm"),
      column_title=column_title,
      ...);
   if (color_by_column) {
      hm_sets <- rownames(mem$enrichIM)[row_order(hm)];
      jamba::imageByColors(mem$enrichIMcolors[hm_sets,],
         cellnote=sapply(mem$enrichIM[hm_sets,],
            format.pval,
            eps=1e-50,
            digits=2),
         adjustMargins=TRUE,
         flip="y",
         cexCellnote=0.7,
         cex.axis=cex.axis,
         main=column_title,
         groupCellnotes=FALSE,
         ...);
   } else {
      return(hm);
   }
}

#' MultiEnrichMap plot
#'
#' MultiEnrichMap plot
#'
#' This function takes output from `multiEnrichMap()` and produces
#' customized "multiple enrichMap" plots using an igraph network.
#' It differs from the data provided from `multiEnrichMap()` mostly
#' by enabling different overlap filters, and by automating
#' several steps that help with network layout, and node label
#' placement.
#'
#' For the most flexible exploration of data, run `multiEnrichMap()`
#' using a lenient `overlapThreshold`, for example `overlapThreshold=0.1`.
#' Then call this function with increasing `overlap` until the
#' visualization has insightful structure.
#'
#' @return invisibly returns the `igraph` object used for plotting,
#'    a by-product of this function when `do_plot=TRUE` is that
#'    the igraph object is also visualized. All custom plot elements
#'    are updated in the `igraph` object, so in principle a
#'    simple call to `plot(...)` should suffice.
#'
#' @family jam igraph functions
#' @family jam plot functions
#'
#' @param mem `list` object output from `multiEnrichMap()`, specifically
#'    containing elements `"multiEnrichMap","multiEnrichMap2"` which
#'    are expected to be `igraph` objects.
#' @param overlap numeric value between 0 and 1, indicating the Jaccard
#'    overlap filter to use for edge nodes. The value is used to
#'    delete edges whose values `E(g)$overlap` are below this threshold.
#' @param overlap_count numeric value indicating the minimum overlap count
#'    below which edges are removed. The value `E(g)$overlap_count` is
#'    used for this filter.
#' @param do_plot logical indicating whether to plot the final result.
#' @param remove_blanks logical indicating whether to call
#'    `removeIgraphBlanks()` which removes blank/empty colors in
#'    igraph nodes.
#' @param remove_singlets logical indicating whether to remove singlet
#'    nodes, which are nodes that have no connections to other nodes.
#'    By default, singlets are removed, in order to help visualize the
#'    node connections that remain after filtering by `overlap`.
#' @param spread_labels logical indicating whether to call
#'    `spread_igraph_labels()`, which places node labels at an angle offset
#'    from the node, in order to improve default label positions.
#' @param y_bias numeric value passed to `spread_igraph_labels()` when
#'    `spread_labels=TRUE`.
#' @param repulse numeric value used for network layout when
#'    `layout_with_qfrf()` is used.
#' @param sets optional character vector of enriched sets to include,
#'    all other sets will be excluded. These values are matched with
#'    `V(g)$name`.
#' @param rescale logical indicating whether the `igraph` layout
#'    coordinates are scaled to range `c(-1, 1)` before plotting.
#'    In practice, when `rescale=FALSE` the function `jam_igraph()`
#'    is called because it does much better at properly handling
#'    other settings during the change. The effect is mainly to keep
#'    layout aspect intact, in cases where the x- and y-axis ranges
#'    are not roughly the same size, for example a short-wide
#'    layout.
#' @param main character string used as the title when `do_plot=TRUE`.
#'    This character string is passed to `glue::glue()` in order to
#'    include certain argument values such as `overlap`.
#' @param ... additional arguments are passed to `removeIgraphBlanks()`,
#'    `removeIgraphSinglets()`, and `spread_igraph_labels()` as needed.
#'
#' @export
mem_multienrichplot <- function
(mem,
 overlap=0.1,
 overlap_count=2,
 do_plot=TRUE,
 remove_blanks=TRUE,
 remove_singlets=TRUE,
 spread_labels=TRUE,
 y_bias=1,
 label_edges=c("overlap_count","count","overlap","label","none"),
 edge_cex=1,
 node_cex=0.8,
 node_size=5,
 edge_color="#55555588",
 frame_color="#55555500",
 shape="pie",
 repulse=3.7,
 sets=NULL,
 rescale=TRUE,
 main="MultiEnrichMap\noverlap >= {overlap}, overlap_count >= {overlap_count}",
 ...)
{
   ##
   if (!is.list(mem) || !"multiEnrichMap2" %in% names(mem)) {
      stop("Input mem must be a list with element 'multiEnrichMap2'.");
   }
   g <- mem$multiEnrichMap2;

   ## Optionally filter for specific sets
   if (length(sets) > 0) {
      keep_nodes <- which(V(g)$name %in% sets);
      if (length(keep_nodes) == 0) {
         stop("No sets remain after filtering V(g)$name for sets.");
      }
      g <- subgraph(g, v=keep_nodes);
   }

   ## Filter to remove nodes as needed
   if (length(overlap) > 0 && overlap > 0) {
      delete_edges <- which(E(g)$overlap < overlap);
      if (length(delete_edges) > 0) {
         g <- igraph::delete.edges(g, delete_edges);
      }
   }
   if (length(overlap_count) > 0 && "overlap_count" %in% list.edge.attributes(g)) {
      delete_edges <- which(E(g)$overlap_count < overlap_count);
      if (length(delete_edges) > 0) {
         g <- igraph::delete.edges(g, delete_edges);
      }
   }
   if (remove_blanks) {
      g <- removeIgraphBlanks(g, ...);
   }
   if (remove_singlets) {
      g <- removeIgraphSinglets(g, ...);
   }
   if (spread_labels) {
      g <- spread_igraph_labels(g,
         do_reorder=FALSE,
         y_bias=y_bias,
         repulse=repulse,
         ...);
   }
   ## Optionally label edges
   label_edges <- head(intersect(label_edges, list.edge.attributes(g)), 1);
   if (length(label_edges) > 0) {
      if (!"label" %in% label_edges) {
         edge_text <- edge_attr(g, label_edges);
         if (is.numeric(edge_text)) {
            edge_text <- format(edge_text,
               trim=TRUE,
               digits=2);
         }
         E(g)$label <- edge_text;
      }
   }
   if (length(edge_cex) > 0) {
      E(g)$label.cex <- edge_cex;
   }
   if (length(edge_color) > 0) {
      E(g)$color <- edge_color;
   }
   if (length(shape) > 0) {
      V(g)$shape <- shape;
   }
   if (length(node_size) > 0) {
      V(g)$size <- node_size;
   }
   if (length(node_cex) > 0) {
      V(g)$label.cex <- node_cex;
   }
   if (length(frame_color) > 0) {
      V(g)$frame.color <- frame_color;
   }
   if (length(main) > 0) {
      main <- glue::glue(main);
   }
   if (do_plot) {
      if (rescale) {
         plot(g,
            main=main,
            rescale=rescale,
            ...);
      } else {
         jam_igraph(g,
            main=main,
            rescale=rescale,
            ...);
      }
   }
   invisible(g);
}

#' MultiEnrichMap color legend
#'
#' MultiEnrichMap color legend
#'
#' This function is a simple wrapper around `legend()` to add
#' a color key to a plot, typically for `igraph` plots.
#'
#' @family jam plot functions
#'
#' @param mem `list` object output from `multiEnrichMap()`, specifically
#'    expected to contain element `"colorV"`.
#' @param x,y,bg,box.col,title,cex,ncol,pch,pt.cex arguments passed
#'    to `legend()`.
#' @param ... additional arguments are passed to `legend()`.
#'
#' @export
mem_legend <- function
(mem,
 x="bottomleft",
 y=NULL,
 bg="#FFFFFF99",
 box.col="transparent",
 title="Color Key",
 cex=0.8,
 ncol=1,
 pch=21,
 pt.cex=2,
 ...)
{
   ##
   if (!is.list(mem) || !"colorV" %in% names(mem)) {
      stop("Input mem must be a list with element 'colorV'");
   }
   colorV <- mem[["colorV"]];
   colorVb <- makeColorDarker(colorV, darkFactor=1.5);

   legend(x=x,
      y=y,
      title=title,
      ncol=ncol,
      legend=names(colorV),
      pch=pch,
      pt.cex=pt.cex,
      pt.bg=colorV,
      col=colorVb,
      bg=bg,
      box.col=box.col,
      cex=cex,
      ...);
}

#' Jam wrapper to plot igraph
#'
#' Jam wrapper to plot igraph
#'
#' This function is a lightweight wrapper around `igraph::plot.igraph()`
#' intended to handle `rescale=FALSE` properly, which is not done
#' in the former function (as of January 2020). The desired outcome
#' is for the `xlim` and `ylim` defaults to be scaled according to the
#' `igraph` layout. Similarly, the `vertex.size` and `vertex.label.dist`
#' parameters should also be scaled proportionally.
#'
#' @family jam igraph functions
#' @family jam plot functions
#'
#' @param x `igraph` object to be plotted
#' @param ... additional arguments are passed to `igraph::plot.igraph()`
#' @param xlim,ylim default x and y axis limits
#' @param expand numeric value used to expand the x and y axis ranges,
#'    where `0.03` expands each size `3%`.
#' @param rescale logical indicating whether to rescale the layout
#'    coordinates to `c(-1, 1)`. When `rescale=FALSE` the original
#'    layout coordinates are used as-is without change.
#'
#' @export
jam_igraph <- function
(x,
 ...,
 xlim=c(-1,1),
 ylim=c(-1,1),
 expand=0.03,
 rescale=FALSE)
{
   ##
   params <- igraph:::i.parse.plot.params(x, list(...));
   layout <- params("plot", "layout");
   vertex.size <- params("vertex", "size");
   label.dist <- params("vertex", "label.dist");
   if (!rescale) {
      xlim <- range(layout[,1]);
      vertex.size <- vertex.size * diff(xlim) / 2;
      label.dist <- label.dist * diff(xlim) / 2;
      xlim <- xlim + diff(xlim) * c(-1,1) * expand;
      ylim <- range(layout[,2]);
      ylim <- ylim + diff(ylim) * c(-1,1) * expand;
   }
   plot(x=x,
      ...,
      rescale=rescale,
      vertex.size=vertex.size,
      vertex.label.dist=label.dist,
      xlim=xlim,
      ylim=ylim);
}

