
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
