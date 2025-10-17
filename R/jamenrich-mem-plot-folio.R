

#' Multienrichment folio of summary plots
#'
#' Multienrichment folio of summary plots
#'
#' This function is intended to create multiple summary plots
#' using the output data from `multiEnrichMap()`. By default
#' it creates all plots one by one, sufficient for including
#' in a multi-page PDF document with `cairo_pdf(..., onefile=TRUE)`
#' or `pdf(..., onefile=TRUE)`.
#'
#' The data for each plot object can be created and visualized later
#' with argument `do_plot=FALSE`.
#'
#' Note: Since version `0.0.76.900` the first step in the workflow is
#' to cluster the underlying gene-pathway incidence matrix.
#' This step defines a consistent dendrogram driven by underlying
#' gene content in each pathway.
#' The dendrogram is used by each subsequent plot
#' including the enrichment heatmap.
#'
#' There are two recommended strategies for visualizing multienrichment
#' results:
#'
#' 1. Pathway clusters viewed as a concept network (Cnet) plot.
#'
#'    * Given numerous statistically enriched pathways,
#'    this process defines pathway clusters using the underlying gene-pathway
#'    incidence matrix.
#'    * Within each pathway cluster, the pathways typically share a
#'    high proportion of the same genes, and therefore are expected
#'    to represent very similar functions. Ideally, each cluster represents
#'    some distinct biological function, or a functional theme.
#'    * Benefit: Reducing a large number of pathways to a small number
#'    of clusters greatly improves the options for visualization,
#'    while retaining a comprehensive view of all genes and pathways
#'    involved.
#'    * Benefit: This option is recommended when there are numerous pathways,
#'    and when including more pathways is beneficial to understanding
#'    the overall functional effects of the experimental study.
#'    * Limitation: The downside with this approach is that sometimes this
#'    comprehensive content can be too much detail to interpret in one
#'    figure, overshadowing individual pathways in each cluster.
#'    * Limitation: It may be difficult to recognize a functional theme for
#'    each pathway cluster, unfortunately that process is not (yet) automated
#'    and requires some domain expertise of the pathways and
#'    functions involved.
#'    * Limitation: It may not be possible for one Cnet plot to represent all
#'    functional effects of an experimental study.
#'
#' 2. Exemplar pathways are viewed as a Cnet plot.
#'
#'    * As described above, given numerous statistically enriched pathways,
#'    pathways are clustered using the gene-pathway incidence matrix. One
#'    "exemplar" pathway is selected from each cluster to represent the typical
#'    pathway content in each cluster, usually the most significant pathway in
#'    the cluster, but optionally the pathway containing the most total genes.
#'    * Benefit: This process can produce a cleaner figure than Option 1
#'    PathwayClusters, because fewer pathways and their associated genes are
#'    included in the figure.
#'    * Limitation: This cleaner figure is understandably somewhat
#'    less comprehensive, and may be subject to bias when selecting
#'    exemplar pathways. However the selection of relevant pathways may
#'    be very effective within the context of the experimental study.
#'    * Benefit: The resulting Cnet plot can often improve focus on specific
#'    genes and pathways, which can be advantageous when including numerous
#'    "synonyms" for the same or similar pathways is not beneficial.
#'    * Benefit: This strategy also works particularly well when there are
#'    relatively few enriched pathways, or when argument `topEnrichN` used
#'    with `multiEnrichMap()` was relatively small.
#'
#'
#' The folio of plots includes:
#'
#' 1. **Enrichment Heatmap** (Plot #1), enrichment P-values created using
#' `mem_enrichment_heatmap()`. Note that by default, the Gene-Pathway
#' incidence matrix is also created (invisibly) in order to define
#' consistent pathway clusters.
#' Output list name: `"enrichment_hm"`
#' 2. **Gene-Pathway Incidence Matrix Heatmap** (Plot #2) is created
#' using `mem_gene_path_heatmap()`.
#' This step defines and visualizes the pathway clustering used by all
#' plots in the folio.
#' Output list name: `"gp_hm"`
#' 3. **Cnet Cluster Plot** (Plots #3,#4,#5) creates a collapsed
#' Concept network (Cnet) of Genes with Pathway clusters,
#' using `collapse_mem_clusters()`, then plotted with `jam_igraph()`.
#'
#'    * Plot #3 labels the pathway clusters with the first N pathways.
#'    Output list name: `"cnet_collapsed"`
#'    * Plot #4 labels the pathway clusters with LETTERS.
#'    **This file is typically used for other plots.**
#'    Output list name: `"cnet_collapsed_set"`
#'    * Plot #5 hides all gene labels.
#'    Output list name: `"cnet_collapsed_set2"`
#'
#' 4. **Cnet Exemplar Plots** (Plots #6,#7,#8) creates smaller pathway
#' Cnet plots, as opposed to pathway-cluster Cnets in #3,#4,#5 above,
#' using exemplar pathways from each gene-pathway cluster.
#' Output list name: `"cnet_exemplars"` with a `list` of `igraph` objects:
#'
#'    * Plots #6 includes one exemplar pathway per pathway cluster.
#'    * Plots #7 includes two exemplar pathways per pathway cluster.
#'    * Plots #8 includes three exemplar pathways per pathway cluster.
#'
#' 5. **Cnet Individual Cluster Plots** (Plots #9,#10,#11,etc.) create one
#' pathway Cnet plot per individual pathway cluster, showing only
#' the pathways in that cluster. The number of plots are defined by
#' the number of pathway cluters, usually `pathway_column_split`.
#' These plots may be useful to explore pathways in detail within each
#' pathway cluster, for example when there are many pathways which are
#' not well-defined for a particular pathway cluster in the Gene-Pathway
#' heatmap.
#' Output list name `"cnet_clusters"`
#'
#' The specific plots to be created are controlled with `do_which`:
#'
#' * `do_which=1` will create the enrichment heatmap.
#' * `do_which=2` will create the gene-pathway heatmap.
#' * `do_which=3` will create the Cnet Cluster Plot using
#' pathway cluster labels for each pathway node, by default it uses `LETTERS`:
#' `"A", "B", "C", "D"`, etc.
#' * `do_which=4` will create the Cnet Cluster Plot using abbreviated
#' pathway labels for each pathway cluster node.
#' * `do_which=5` will create the Cnet Cluster Plot with no node labels.
#' * `do_which=6` begins the series of Cnet Exemplar Plots for each value
#' in argument `exemplar_range`, whose default is `c(1, 2, 3)`.
#' * `do_which=9` (by default) begins the series of Cnet individual
#' cluster plots, which includes all pathways from each cluster.
#'
#' The most frequently used plots are `do_which=2` for the
#' gene-pathway heatmap, and `do_which=4` for the collapsed Cnet
#' plot, where Cnet clusters are based upon the gene-pathway heatmap.
#'
#' Arguments `p_cutoff` and `min_set_ct_each` can be used to
#' apply more stringent thresholds than the original `mem` data.
#' For example, applying `p_cutoff=0.05` during `multiEnrichMap()`
#' will colorize pathways in `mem$enrichIMcolors`, however when
#' calling `mem_plot_folio()` with `p_cutoff=0.001` will use blank
#' color in the color gradient for pathways that do not
#' have `mem$enrichIM` value at or below `0.001`.
#'
#' Our experience is that the pathway clustering does not need to
#' be perfect to be useful and valid. The pathway clusters
#' are valid based upon the parameters used for clustering,
#' and provide insight into the genes that help define each
#' cluster distinct from other clusters.
#' Sometimes the clustering results are more or less effective
#' based upon the type of pattern observed in the data, so it
#' can be helpful to adjust parameters to drill down to
#' the most effective patterns.
#'
#' # Gene-Pathway clustering
#'
#' The clustering is performed by combining the gene-pathway incidence
#' matrix `mem$memIM` with the `-log10(mem$enrichIM)` enrichment P-values.
#' The relative weight of each matrix is controlled by
#' `enrich_im_weight`, where `enrich_im_weight=0` assigns weight=0
#' to the enrichment P-values, and thus clusters only using the
#' gene-pathway matrix. Similarly, `enrich_im_weight=1` will assign
#' full weight to the enrichment P-value matrix, and will ignore
#' the gene-pathway matrix data.
#'
#' The corresponding weight for gene (rows) is controlled by
#' `gene_im_weight`, which balances row clustering with the
#' `mem$geneIM` matrix, and the gene-pathway matrix `mem$memIM`.
#'
#' The argument `column_method` defines the distance method,
#' for example `"euclidean"` and `"binary"` are two immediate choices.
#' The method also adds `"correlation"` from `amap::hcluster()` which
#' can be very useful especially with large datasets.
#'
#' The number of pathway clusters is controlled by
#' `pathway_column_split`, by default when `pathway_column_split=NULL`
#' and `auto_cluster=TRUE` the number of clusters is defined based
#' upon the total number of pathways. In practice, `pathway_column_split=4`
#' or `pathway_column_split=3` is recommended, as this number of
#' clusters is most convenient to visualize as a Cnet plot.
#'
#' To define your own pathway cluster labels, define `pathway_column_title`
#' as a vector with length equal to `pathway_column_split`. These labels
#' become network node labels in subsequent plots, and in the
#' resulting `igraph` object.
#'
#' The pathway clusters are dependent upon the genes and pathways
#' used during clustering, which are also controlled by
#' `min_set_ct` and `min_gene_ct`.
#' * `min_set_ct` filters the matrix by the number of times a Set is
#' represented in the matrix,
#' which can be helpful when there are pathways with large number of
#' genes, with some pathways with very low number of genes.
#' * `min_gene_ct` filters the matrix by the number of times a gene is
#' represented in the matrix. It can be helpful for requiring a gene
#' be represented in more than one enriched pathway.
#' * `min_set_ct_each` filters the matrix to require each Set to
#' contain at least this many entries from one enrichment result,
#' rather than using the combined incidence matrix. It is mostly
#' helpful to increase the value used in `multiEnrichMap()` argument
#' `min_count`, which already filters pathways for minimum number
#' of genes involved.
#' * Note: These filters are only recommended when the gene-pathway
#' matrix is very large, perhaps 100 pathways, or 500 genes.
#'
#' # Cnet pathway clusters
#'
#' The resulting Cnet pathway clusters are single nodes in the
#' network, and these nodes are colorized based upon the enrichment
#' tests involved. The threshold for including the color for
#' each enrichment test is defined by `cluster_color_min_fraction`,
#' which requires at least this fraction of pathways in a
#' pathway cluster meets the significance criteria for that
#' enrichment test.
#'
#' To adjust the coloration filter to include any enrichment
#' test with at least one significant result, use
#' `cluster_color_min_fraction=0.01`.
#' In the gene-pathway heatmap,
#' these colors are shown across the top of the heatmap.
#' The default `cluster_color_min_fraction=0.4` requires 40%
#' of pathways in a cluster for each enrichment test.
#'
#' Note: Prior to version `0.0.76.900`
#' the enrichment heatmap was clustered only using enrichment
#' P-values, transformed with `log10(Pvalue)`. The clustering was
#' inconsistent with other plots in the folio, and was not effective
#' at clustering pathways based upon similar content, which is the
#' primary goal of the `multienrichjam` R package.
#'
#' @family jam plot functions
#'
#' @returns `list` returned using `invisible()`, containing each
#'    plot object enabled by the argument `do_which`:
#'    * `enrichment_hm` is a Heatmap object from `ComplexHeatmap`
#'    that contains the enrichment P-value heatmap. Note that this
#'    data is not used directly in subsequent plots, the pathway
#'    clusters shown here are based upon `-log10(Pvalue)` and not
#'    the underlying gene content of each pathway. This plot is
#'    a useful overview that answers the question "How many
#'    pathways are significantly enriched across the different
#'    enrichment tests?"
#'    * `gp_hm` is a Heatmap object from `ComplexHeatmap` with
#'    the gene-pathway incidence matrix heatmap. This heatmap and
#'    the column/pathway clusters are the subject of subsequent
#'    Cnet plots.
#'    * `gp_hm_caption` is a text caption that describes the gene
#'    and set filter criteria, and the row and column distance methods
#'    used for clustering. Because the filtering and clustering
#'    options have substantial impact on clustering, and the
#'    pathway clusters are the key for all subsequent plots,
#'    these values are important to keep associated with the
#'    output of this function.
#'    * `clusters_mem` is a `list` with the pathways contained
#'    in each pathway cluster shown by the gene-pathway heatmap,
#'    obtained by `heatmap_column_order(gp_hm)`. The pathway names
#'    should also be present in `colnames(mem$memIM)` and
#'    `rownames(mem$enrichIM)`, for follow-up inspection.
#'    * `cnet_collapsed` is an `igraph` object with Cnet plot data,
#'    where the pathways have been collapsed by cluster, using the
#'    gene-pathway heatmap clusters defined in `clusters_mem`. Each
#'    pathway cluster is labeled by cluster name, and the first few
#'    pathway names.
#'    This data can be plotted using `jam_igraph(cnet_collapsed)`.
#'    * `cnet_collapsed_set` is the same as `cnet_collapsed` except the
#'    pathways are labeled by the cluster name only, for example
#'    `c("A", "B", "C", "D")`.
#'    This data can be plotted using `jam_igraph(cnet_collapsed_set)`.
#'    * `cnet_collapsed_set2` is the same as `cnet_collapsed_set` except the
#'    gene labels are hidden, useful when there are too many genes to label
#'    clearly. The gene symbols are still stored in `V(g)$name` but the labels
#'    in `V(g)$label` are updated to hide the genes.
#'    This data can be plotted using `jam_igraph(cnet_collapsed_set2)`.
#'    * `cnet_exemplars` is a `list` of `igraph` Cnet objects, each
#'    one contains only the number of exemplar pathways from each cluster
#'    defined by argument `exemplar_range`. By default it uses `1` exemplar
#'    per cluster, then `2` exemplars per cluster, then `3` exemplars
#'    per cluster. A number of published figures use `1` exemplar per
#'    pathway cluster.
#'    This data can be plotted using `jam_igraph(cnet_exemplars[[1]])`,
#'    which will plot only the first `igraph` object from the list.
#'    * `cnet_clusters` is a `list` of `igraph` Cnet objects, each one
#'    contains all the pathways in one pathway cluster.
#'    This data can be plotted using `jam_igraph(cnet_clusters[[1]])`,
#'    or by calling a specific cluster `jam_igraph(cnet_clusters[["A"]])`.
#'
#'
#' @param mem `list` object created by `multiEnrichMap()`. Specifically
#'    the object is expected to contain `colorV`, `enrichIM`,
#'    `memIM`, `geneIM`.
#' @param do_which `integer` vector of plots to produce. When `do_which`
#'    is `NULL`, then all plots are produced. This argument is intended
#'    to help produce one plot from a folio, therefore each plot is referred
#'    by the number of the plot, in order.
#' @param p_cutoff `numeric` value indicating the enrichment P-value threshold
#'    used for `multiEnrichMap()`, but when `NULL` this value is taken
#'    from the `mem` input, or `0.05` is used by default.
#' @param p_floor `numeric` value indicating the lowest enrichment P-value
#'    used in the color gradient on the Enrichment Heatmap.
#' @param main `character` string used as a title on Cnet plots.
#' @param use_raster `logical` default FALSE, deprecated,
#'    whether to use raster heatmaps, passed to `ComplexHeatmap::Heatmap()`.
#'    * Note that `use_raster=TRUE` may produce visual artifacts, and
#'    changing this argument is no longer supported
#' @param min_gene_ct,min_set_ct `integer` values passed to
#'    `mem_gene_path_heatmap()`. The `min_gene_ct` requires each set
#'    to contain `min_gene_ct` genes, and `min_set_ct` requires each gene
#'    to be present in at least `min_set_ct` sets.
#' @param min_set_ct_each `integer` minimum genes required for each set,
#'    required for at least one enrichment test.
#' @param column_method,row_method `character` arguments passed to
#'    `ComplexHeatmap::Heatmap()` which indicate the distance method used
#'    to cluster columns and rows, respectively.
#' @param exemplar_range `integer` vector (or `NULL`) used to create Cnet
#'    exemplar plots, using this many exemplars per cluster.
#' @param pathway_column_split,gene_row_split `integer` value passed
#'    as `column_split` and `row_split`, respectively, to
#'    `mem_gene_path_heatmap()`, indicating the number of pathway
#'    clusters, and gene clusters, to create in the gene-pathway heatmap.
#'    When either value is `NULL` then auto-split logic is used.
#' @param pathway_column_title,gene_row_title `character` vectors
#'    passed to `mem_gene_path_heatmap()` as `column_title` and
#'    `row_title`, respectively. When one value is supplied, it is
#'    displayed and centered across all the respective splits. When
#'    multiple values are supplied, values are used to the number
#'    of splits, and recycled as needed. In that case, repeated
#'    values are made unique by `jamba::makeNames()`.
#' @param cex.main,cex.sub `numeric` values passed to `title()` which
#'    size the default title and sub-title in Cnet plots.
#' @param row_cex,column_cex `numeric` character expansion factor, used
#'    to adjust the relative size of row and column labels,
#'    respectively. A value of `1.1` will make row font size 10%
#'    larger.
#' @param color_by_column `logical` indicating whether to colorize
#'    the enrichment heatmap columns using `colorV` in the input `mem`.
#'    This argument is only relevant when `do_which` include `1`.
#' @param enrich_im_weight,gene_im_weight `numeric` value between 0 and
#'    1, passed to `mem_gene_path_heatmap()`, used to apply relative
#'    weight to clustering columns and rows, respectively, when
#'    combining the gene-pathway incidence matrix with either column
#'    enrichment P-values, or row gene incidence matrix data.
#' @param colorize_by_gene `logical` passed to `mem_gene_path_heatmap()`
#'    indicating whether the heatmap body for the gene-pathway heatmap
#'    will be colorized using the enrichment colors for each gene.
#' @param cluster_color_min_fraction `numeric` value passed to
#'    `collapse_mem_clusters()` used to determine which enrichment
#'    colors to associate with each Cnet cluster.
#' @param byCols `character` vector describing how to sort the
#'    pathways within Cnet clusters. This argument is passed
#'    to `rank_mem_clusters()`.
#' @param edge_bundling `character` string passed to `jam_igraph()`
#'    to control edge bundling. The default `edge_bundling="connections"`
#'    will bundle Cnet plot edges for genes that share the same pathway
#'    connections.
#' @param apply_direction `logical` or `NULL` indicating whether to
#'    indicate directionality in the `mem_enrichment_heatmap()` which is
#'    the first plot in the series. The default `apply_direction=NULL`
#'    will auto-detect whether there is directionality present in the
#'    data, and will set `apply_direction=TRUE` only when there are non-NA
#'    values that differ from zero.
#' @param rotate_heatmap `logical` passed to `mem_gene_path_heatmap()()`
#'    and only this function, default `FALSE`. It indicates whether to
#'    rotate the heatmap to have gene columns and pathway rows.
#'    If you find most people tilt their head to read the pathways,
#'    it might be preferable.
#' @param row_anno_padding,column_anno_padding `grid::unit` or `numeric`
#'    which will be converted to "mm" units. These values control the
#'    space between the heatmap body and row/column annotations,
#'    respectively, only relevant for `mem_gene_path_heatmap()()`.
#'    The value is only applied during `draw()` and cannot be
#'    defined in the `Heatmap` object itself, which is why it is
#'    included here and not `mem_gene_path_heatmap()()`.
#' @param do_plot `logical` indicating whether to render each plot.
#'    When `do_plot=FALSE` the plot objects will be created and returned,
#'    but the plot itself will not be rendered. This option may be
#'    useful to generate the full set of figures in one set, then
#'    review each figure one by one in an interactive session.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to downstream functions.
#'    Some useful examples:
#'    * `sets` is passed to `mem_gene_path_heatmap()` which
#'    allows one to define a specific subset of sets to use in the
#'    gene-pathway heatmap.
#'    * `cell_size` is passed to `mem_enrichment_heatmap()` with the
#'    option to define square cell size in the heatmap dotplot. However,
#'    the resulting heatmap will be at least `ncol * cell_height` width,
#'    and `nrow * cell_size[2]` height, in addition to the heights of
#'    the title and column labels, and widths of the color key and
#'    dendrogram.
#'
#' @export
mem_plot_folio <- function
(mem,
 do_which=NULL,
 p_cutoff=NULL,
 p_floor=1e-10,
 main="",
 use_raster=FALSE,
 min_gene_ct=1,
 min_set_ct=1,
 min_set_ct_each=4,
 column_method="euclidean",
 row_method="euclidean",
 exemplar_range=c(1, 2, 3),
 pathway_column_split=NULL,
 pathway_column_title=LETTERS,
 gene_row_split=NULL,
 gene_row_title=letters,
 edge_color=NULL,
 cex.main=2,
 cex.sub=1.5,
 row_cex=1,
 column_cex=1,
 max_labels=4,
 max_nchar_labels=25,
 include_cluster_title=TRUE,
 repulse=4,
 use_shadowText=FALSE,
 color_by_column=FALSE,
 style="dotplot_inverted",
 enrich_im_weight=0.3,
 gene_im_weight=0.5,
 colorize_by_gene=TRUE,
 cluster_color_min_fraction=0.4,
 byCols=c("composite_rank",
    "minp_rank",
    "gene_count_rank"),
 edge_bundling="connections",
 apply_direction=NULL,
 rotate_heatmap=FALSE,
 row_anno_padding=grid::unit(3, "mm"),
 column_anno_padding=grid::unit(3, "mm"),
 do_plot=TRUE,
 verbose=FALSE,
 ...)
{
   # accept Mem input, convert to list for now
   if (inherits(mem, "Mem")) {
      mem <- Mem_to_list(Mem)
   }
   # validate arguments
   if (length(p_cutoff) == 0) {
      if ("p_cutoff" %in% names(mem)) {
         p_cutoff <- mem$p_cutoff;
      } else {
         p_cutoff <- 0.05;
      }
   }
   ret_vals <- list();
   plot_num <- 0;
   if (length(do_which) == 0) {
      do_which <- seq_len(50);
   }

   #############################################################
   # Apply optional heatmap padding parameters
   rowpad <- NULL;
   columnpad <- NULL;
   if (length(row_anno_padding) > 0) {
      if (!grid::is.unit(row_anno_padding)) {
         row_anno_padding <- grid::unit(row_anno_padding, "mm");
      }
      rowpad <- ComplexHeatmap::ht_opt("ROW_ANNO_PADDING")
      ComplexHeatmap::ht_opt(ROW_ANNO_PADDING=row_anno_padding)
   }
   if (length(column_anno_padding) > 0) {
      if (!grid::is.unit(column_anno_padding)) {
         column_anno_padding <- grid::unit(column_anno_padding, "mm");
      }
      columnpad <- ComplexHeatmap::ht_opt("COLUMN_ANNO_PADDING")
      ComplexHeatmap::ht_opt(COLUMN_ANNO_PADDING=column_anno_padding)
   }
   # revert changes later
   revert_hm_padding <- function(...){
      if (length(rowpad) > 0) {
         ComplexHeatmap::ht_opt(ROW_ANNO_PADDING=rowpad)
      }
      if (length(columnpad) > 0) {
         ComplexHeatmap::ht_opt(COLUMN_ANNO_PADDING=columnpad)
      }
   }

   #############################################################
   # First step: Define the Gene-Pathway Heatmap
   #
   # All subsequent plots depend upon mem_gene_path_heatmap()
   if (verbose) {
      jamba::printDebug("mem_plot_folio(): ",
         "Gene-pathway heatmap (pre-emptive)");
   }
   use_raster <- FALSE;
   gp_hm <- mem_gene_path_heatmap(mem,
      p_cutoff=p_cutoff,
      p_floor=p_floor,
      row_method=row_method,
      column_split=pathway_column_split,
      column_title=pathway_column_title,
      row_split=gene_row_split,
      row_title=gene_row_title,
      column_method=column_method,
      row_cex=row_cex,
      column_cex=column_cex,
      use_raster=use_raster,
      min_gene_ct=min_gene_ct,
      min_set_ct=min_set_ct,
      min_set_ct_each=min_set_ct_each,
      enrich_im_weight=enrich_im_weight,
      gene_im_weight=gene_im_weight,
      colorize_by_gene=colorize_by_gene,
      rotate_heatmap=rotate_heatmap,
      ...);

   # Extract hclust object for re-use in the enrichment heatmap
   use_sets <- colnames(gp_hm@matrix);
   gphm_pw_order <- NULL;
   if (rotate_heatmap) {
      gp_hm_hclust <- gp_hm@row_dend_param$obj;
      if (length(gp_hm_hclust) == 0) {
         gp_hm_hclust <- gp_hm@row_dend_param$fun(t(gp_hm@matrix));
      }
      use_sets <- rownames(gp_hm@matrix);
      gphm_pw_order <- jamba::heatmap_row_order(gp_hm);
   } else {
      gp_hm_hclust <- gp_hm@column_dend_param$obj;
      if (length(gp_hm_hclust) == 0) {
         gp_hm_hclust <- gp_hm@column_dend_param$fun(t(gp_hm@matrix));
      }
      gphm_pw_order <- jamba::heatmap_column_order(gp_hm);
   }


   #############################################################
   ## Enrichment P-value heatmap
   plot_num <- plot_num + 1;
   if (length(do_which) == 0 || plot_num %in% do_which) {
      if (verbose) {
         jamba::printDebug("mem_plot_folio(): ",
            c("plot_num ", plot_num, ": "),
            c("Enrichment P-value Heatmap"),
            sep="");
      }
      if (length(apply_direction) == 0) {
         apply_direction <- FALSE;
         if ("enrichIMdirection" %in% names(mem) &&
               any(!is.na(mem$enrichIMdirection)) &&
               any(jamba::rmNA(mem$enrichIMdirection) < 0)) {
            apply_direction <- TRUE;
         }
      }
      # version 0.0.76.900 use fixed sets, in same order as previous gp_hm,
      # then re-use the hclust object and same pathway_column_split
      ## 0.0.89.900 - apply hclust and turn off additional row ordering
      ## so the row order will exactly match the gp_hm pathway order.
      mem_hm <- jamba::call_fn_ellipsis(
         mem_enrichment_heatmap,
         mem=mem,
         sets=use_sets,
         p_cutoff=p_cutoff,
         p_floor=p_floor,
         color_by_column=color_by_column,
         row_cex=row_cex,
         row_method=row_method,
         row_split=length(gphm_pw_order),
         row_title=names(gphm_pw_order),
         row_title_rot=0,
         cluster_rows=gp_hm_hclust,
         row_dend_reorder=FALSE,
         column_cex=column_cex,
         column_method=column_method,
         style=style,
         apply_direction=apply_direction,
         use_raster=use_raster,
         do_plot=do_plot,
         ...);
      if (do_plot) {
         if (!color_by_column) {
            if ("annotation_legend_list" %in% names(attributes(mem_hm))) {
               annotation_legend_list <- attributes(mem_hm)$annotation_legend_list;
            } else {
               annotation_legend_list <- NULL;
               if (length(main) > 0 && nchar(main) > 0) {
                  mem_hm <- ComplexHeatmap::draw(mem_hm,
                     annotation_legend_list=annotation_legend_list,
                     column_title=main);
               } else {
                  mem_hm <- ComplexHeatmap::draw(mem_hm,
                     annotation_legend_list=annotation_legend_list);
               }
            }
         }
      }
      ret_vals$enrichment_hm <- mem_hm;
   }


   #############################################################
   ## Gene-Pathway Heatmap
   if (length(do_which) > 0 && !any(do_which > plot_num)) {
      revert_hm_padding();
      return(invisible(ret_vals));
   }
   ## All subsequent plots depend upon mem_gene_path_heatmap()
   if (verbose) {
      jamba::printDebug("mem_plot_folio(): ",
         "Gene-pathway heatmap");
   }

   ## draw the heatmap
   plot_num <- plot_num + 1;

   # generate caption and include in returned results
   draw_caption <- NULL;
   caption_legendlist <- NULL;
   caption <- NULL;
   if ("caption_legendlist" %in% names(attributes(gp_hm))) {
      caption_legendlist <- attr(gp_hm, "caption_legendlist");
   }
   if ("caption" %in% names(attributes(gp_hm))) {
      caption <- attr(gp_hm, "caption");
   }
   if ("draw_caption" %in% names(attributes(gp_hm))) {
      draw_caption <- attr(gp_hm, "draw_caption");
   }
   ret_vals$gp_hm <- gp_hm;
   ret_vals$gp_hm_caption_legendlist <- caption_legendlist;
   ret_vals$gp_hm_caption <- caption;

   if (length(do_which) == 0 || plot_num %in% do_which) {
      if (verbose) {
         jamba::printDebug("mem_plot_folio(): ",
            c("plot_num ", plot_num, ": "),
            c("Gene-Pathway Heatmap"),
            sep="");
      }
      # Optionally increase padding between annotation and heatmap body
      #row_anno_padding <- ComplexHeatmap::ht_opt$ROW_ANNO_PADDING;
      #column_anno_padding <- ComplexHeatmap::ht_opt$COLUMN_ANNO_PADDING;
      if (do_plot) {
         if (length(caption_legendlist) > 0) {
            ComplexHeatmap::draw(gp_hm,
               annotation_legend_list=caption_legendlist,
               merge_legends=TRUE,
               column_title=main,
               column_title_gp=grid::gpar(fontsize=18));
         } else if (length(draw_caption) > 0) {
            ComplexHeatmap::draw(gp_hm,
               merge_legends=TRUE,
               column_title=main,
               column_title_gp=grid::gpar(fontsize=18));
            draw_caption();
         } else {
            ComplexHeatmap::draw(gp_hm,
               merge_legends=TRUE,
               column_title=main,
               column_title_gp=grid::gpar(fontsize=18));
            # grid_with_title(gp_hm,
            #    title=main,
            #    caption=caption);
         }
      }
   }
   ## Obtain heatmap pathway clusters
   if (rotate_heatmap) {
      clusters_mem <- jamba::heatmap_row_order(gp_hm);
   } else {
      clusters_mem <- jamba::heatmap_column_order(gp_hm);
   }
   ret_vals$clusters_mem <- clusters_mem;
   ## Get number of pathway clusters
   pathway_clusters_n <- length(clusters_mem);
   if (verbose) {
      jamba::printDebug("mem_plot_folio(): ",
         c("Defined ", pathway_clusters_n, " pathway clusters."),
         sep="");
   }


   #############################################################
   ## Cnet collapsed
   if (any(c(plot_num + c(1, 2, 3)) %in% do_which)) {
      if (verbose) {
         jamba::printDebug("mem_plot_folio(): ",
            "Preparing Cnet collapsed");
      }
      # first apply stat cutoffs to temporary enrichIMcolors matrix
      # to ensure colors use the defined cutoffs
      # NOTE: Consider updating mem$enrichIMcolors to add missing colors
      # where the p_cutoff, min_set_ct_each are more lenient than the
      # original data. For now it is acceptable to force this step to use
      # thresholds at least as stringent as the original.
      blank_hit_matrix <- (mem$enrichIM > p_cutoff |
            mem$enrichIMgeneCount < min_set_ct_each)
      if (any(blank_hit_matrix)) {
         mem$enrichIMcolors[blank_hit_matrix] <- "#FFFFFF";
      }

      cnet_collapsed <- tryCatch({
         collapse_mem_clusters(mem=mem,
            clusters=clusters_mem,
            verbose=verbose>1,
            max_labels=max_labels,
            byCols=byCols,
            cluster_color_min_fraction=cluster_color_min_fraction,
            max_nchar_labels=max_nchar_labels,
            include_cluster_title=include_cluster_title,
            return_type="cnet",
            ...);
      }, error=function(e){
         jamba::printDebug("Error during collapse_mem_clusters(), ",
            "returning NULL.");
         print(e);
         NULL;
      });
      if (length(cnet_collapsed) == 0) {
         revert_hm_padding();
         return(list(mem=mem,
            clusters_mem=clusters_mem,
            ret_vals=ret_vals))
      }
      igraph::V(cnet_collapsed)$pie.color <- lapply(
         igraph::V(cnet_collapsed)$pie.color, function(i){
            j <- ifelse(names(i) %in% names(mem$colorV) & !isColorBlank(i),
               mem$colorV[names(i)],
               i);
         });
      igraph::V(cnet_collapsed)$coloredrect.color <- lapply(
         igraph::V(cnet_collapsed)$coloredrect.color, function(i){
            j <- ifelse(names(i) %in% names(mem$colorV) & !isColorBlank(i),
               mem$colorV[names(i)],
               i);
         });

      if (verbose) {
         jamba::printDebug("mem_plot_folio(): ",
            "subsetCnetIgraph()");
      }
      cnet_collapsed <- tryCatch({
         cnet_collapsed <- subsetCnetIgraph(
            cnet_collapsed,
            remove_blanks=TRUE,
            repulse=repulse,
            verbose=verbose > 1);
         # cnet_collapsed %>%
         #    subsetCnetIgraph(remove_blanks=TRUE,
         #       repulse=repulse,
         #       verbose=verbose>1);
      }, error=function(e){
         if (verbose) {
            jamba::printDebug("mem_plot_folio(): ",
               "subsetCnetIgraph() error during ",
               "subsetCnetIgraph(..., remove_blanks=TRUE):");
            print(e);
         }
         cnet_collapsed <- tryCatch({
            cnet_collapsed <- subsetCnetIgraph(
               cnet_collapsed,
               remove_blanks=FALSE,
               repulse=repulse,
               verbose=verbose>1);
            # cnet_collapsed %>%
            #    subsetCnetIgraph(remove_blanks=FALSE,
            #       repulse=repulse,
            #       verbose=verbose>1);
         }, error=function(e2){
            jamba::printDebug("mem_plot_folio(): ",
               "subsetCnetIgraph() error during ",
               "subsetCnetIgraph(), skipping this operation.");
            print(e2);
            cnet_collapsed;
         })
      })
      if (length(edge_color) > 0) {
         igraph::E(cnet_collapsed)$color <- edge_color;
      }
      plot_num <- plot_num + 1;
      if (length(do_which) == 0 || plot_num %in% do_which) {
         if (verbose) {
            jamba::printDebug("mem_plot_folio(): ",
               c("plot_num ", plot_num, ": "),
               c("Cnet collapsed ", "with gene and cluster labels"),
               sep="");
         }
         cnet_title <- "Cnet plot using collapsed clusters";
         cnet_collapsed <- igraph::set_graph_attr(cnet_collapsed,
            name="title",
            value=cnet_title);
         ret_vals$cnet_collapsed <- cnet_collapsed;
         ## Draw Cnet collapsed
         if (do_plot) {
            jam_igraph(cnet_collapsed,
               use_shadowText=use_shadowText,
               edge_bundling=edge_bundling,
               ...);
            mem_legend(mem);
            title(sub="Cnet plot using collapsed clusters",
               main=main,
               cex.main=cex.main,
               cex.sub=cex.sub);
         }
      }
      ## Draw Cnet collapsed with top n labels
      #isset <- (V(cnet_collapsed)$nodeType %in% "Set");
      if ("set_labels" %in% igraph::list.vertex.attributes(cnet_collapsed)) {
         igraph::V(cnet_collapsed)$label <- ifelse(
            nchar(jamba::rmNA(naValue="", igraph::V(cnet_collapsed)$set_labels)) > 0,
            igraph::V(cnet_collapsed)$set_labels,
            igraph::V(cnet_collapsed)$name);
      }
      plot_num <- plot_num + 1;
      if (length(do_which) == 0 || plot_num %in% do_which) {
         if (verbose) {
            jamba::printDebug("mem_plot_folio(): ",
               c("plot_num ", plot_num, ": "),
               c("Cnet collapsed ", "with gene and set labels"),
               sep="");
         }
         cnet_title <- "Cnet plot using collapsed clusters\nlabeled by set";
         cnet_collapsed <- igraph::set_graph_attr(cnet_collapsed,
            name="title",
            value=cnet_title);
         ret_vals$cnet_collapsed_set <- cnet_collapsed;
         if (do_plot) {
            jam_igraph(cnet_collapsed,
               use_shadowText=use_shadowText,
               edge_bundling=edge_bundling,
               ...);
            mem_legend(mem);
            title(sub=cnet_title,
               main=main,
               cex.main=cex.main,
               cex.sub=cex.sub);
         }
      }
      ## Draw Cnet collapsed with top n labels, no gene labels
      plot_num <- plot_num + 1;
      if (length(do_which) == 0 || plot_num %in% do_which) {
         if (verbose) {
            jamba::printDebug("mem_plot_folio(): ",
               c("plot_num ", plot_num, ": "),
               c("Cnet collapsed ", "with set labels, without gene labels"),
               sep="");
         }
         if (length(igraph::V(cnet_collapsed)$label) == 0) {
            if ("set_labels" %in% igraph::list.vertex.attributes(cnet_collapsed)) {
               igraph::V(cnet_collapsed)$label <- ifelse(
                  nchar(jamba::rmNA(naValue="",
                     igraph::V(cnet_collapsed)$set_labels)) > 0,
                  igraph::V(cnet_collapsed)$set_labels,
                  igraph::V(cnet_collapsed)$name);
            } else {
               igraph::V(cnet_collapsed)$label <- igraph::V(cnet_collapsed)$name;
            }
         }
         igraph::V(cnet_collapsed)$label <- ifelse(
            igraph::V(cnet_collapsed)$nodeType %in% "Gene",
            "",
            igraph::V(cnet_collapsed)$label);
         cnet_title <- paste(sep="\n",
            "Cnet plot using collapsed clusters",
            "labeled by set",
            "gene labels hidden");
         cnet_collapsed <- igraph::set_graph_attr(cnet_collapsed,
            name="title",
            value=cnet_title);
         ret_vals$cnet_collapsed_set2 <- cnet_collapsed;
         if (do_plot) {
            jam_igraph(cnet_collapsed,
               use_shadowText=use_shadowText,
               edge_bundling=edge_bundling,
               ...);
            mem_legend(mem);
            title(sub=cnet_title,
               main=main,
               cex.main=cex.main,
               cex.sub=cex.sub);
         }
      }
   } else {
      plot_num <- plot_num + 3;
   }

   #############################################################
   ## Prepare for Cnet plots
   cnet_range <- seq_len(length(exemplar_range) + pathway_clusters_n);
   if (any(c(plot_num + cnet_range) %in% do_which)) {
      if (verbose) {
         jamba::printDebug("mem_plot_folio(): ",
            "Preparing cnet for subsetting with memIM2cnet().");
      }
      # spread_labels=FALSE because no layout is required yet
      cnet <- memIM2cnet(mem,
         remove_blanks=TRUE,
         spread_labels=FALSE,
         ...);
      # ## Freshen pie.color by using the original colorV value by name
      # igraph::V(cnet)$pie.color <- lapply(igraph::V(cnet)$pie.color, function(i){
      #    j <- ifelse(names(i) %in% names(mem$colorV) & !isColorBlank(i),
      #       mem$colorV[names(i)],
      #       i);
      # });
      # ## Freshen coloredrect.color by using the original colorV value by name
      # igraph::V(cnet)$coloredrect.color <- lapply(igraph::V(cnet)$coloredrect.color, function(i){
      #    j <- ifelse(names(i) %in% names(mem$colorV) & !isColorBlank(i),
      #       mem$colorV[names(i)],
      #       i);
      # });
      # cnet <- cnet %>%
      #    removeIgraphBlanks();
   }

   #############################################################
   ## Cnet exemplars
   if (any((plot_num + seq_along(exemplar_range)) %in% do_which)) {
      cnet_exemplars <- list();
      for (exemplar_n in exemplar_range) {
         pluralized <- "";
         if (exemplar_n > 1) {
            pluralized <- "s";
         }
         clusters_mem_n <- rank_mem_clusters(mem,
            clusters_mem,
            per_cluster=exemplar_n,
            byCols=byCols,
            ...);
         cnet_exemplar <- subsetCnetIgraph(cnet,
            includeSets=clusters_mem_n$set,
            ...);
         # cnet_exemplar <- cnet %>%
         #    subsetCnetIgraph(includeSets=clusters_mem_n$set,
         #       ...);
         plot_num <- plot_num + 1;
         if (length(do_which) == 0 || plot_num %in% do_which) {
            if (verbose) {
               jamba::printDebug("mem_plot_folio(): ",
                  c("plot_num ", plot_num, ": "),
                  c("Cnet with ", exemplar_n,
                     paste0(" exemplar", pluralized)),
                  sep="");
            }
            cnet_title <- paste0("Cnet plot using ",
               exemplar_n,
               " exemplar",
               pluralized,
               " per cluster")
            cnet_exemplar <- igraph::set_graph_attr(cnet_exemplar,
               name="title",
               value=cnet_title);
            ## Draw Cnet exemplar
            if (do_plot) {
               jam_igraph(cnet_exemplar,
                  use_shadowText=use_shadowText,
                  edge_bundling=edge_bundling,
                  ...);
               title(
                  sub=cnet_title,
                  main=main,
                  cex.main=cex.main,
                  cex.sub=cex.sub);
               mem_legend(mem);
            }
            ## Add to return list
            cnet_exemplars <- c(cnet_exemplars,
               list(cnet_exemplar));
            names(cnet_exemplars)[length(cnet_exemplars)] <- exemplar_n;
         }
      }
      ret_vals$cnet_exemplars <- cnet_exemplars;
   } else {
      plot_num <- plot_num + length(exemplar_range);
   }


   #############################################################
   ## Cnet each cluster
   if (length(do_which) == 0 || any(plot_num + seq_len(pathway_clusters_n) %in% do_which)) {
      cnet_clusters <- list();
      for (cluster_name in names(clusters_mem)) {
         plot_num <- plot_num + 1;
         if (length(do_which) == 0 || plot_num %in% do_which) {
            if (verbose) {
               jamba::printDebug("mem_plot_folio(): ",
                  c("plot_num ", plot_num, ": "),
                  c("Cnet cluster ", cluster_name),
                  sep="");
            }
            cluster_sets <- unique(unlist(clusters_mem[[cluster_name]]));
            cnet_cluster <- subsetCnetIgraph(
               cnet,
               includeSets=cluster_sets,
               repulse=repulse,
               ...);
            # cnet_cluster <- cnet %>%
            #    subsetCnetIgraph(includeSets=cluster_sets,
            #       repulse=repulse,
            #       ...);
            cnet_title <- paste0("Cnet plot for cluster ",
               cluster_name);
            cnet_cluster <- igraph::set_graph_attr(cnet_cluster,
               name="title",
               value=cnet_title);
            ## Draw Cnet cluster
            if (do_plot) {
               jam_igraph(cnet_cluster,
                  use_shadowText=use_shadowText,
                  edge_bundling=edge_bundling,
                  ...);
               title(sub=paste0("Cnet plot for cluster ",
                  cluster_name),
                  main=main,
                  cex.main=cex.main,
                  cex.sub=cex.sub);
               mem_legend(mem);
            }
            ## Add to return list
            cnet_clusters <- c(cnet_clusters,
               list(cnet_cluster));
            names(cnet_clusters)[length(cnet_clusters)] <- cluster_name;
         }
      }
      ret_vals$cnet_clusters <- cnet_clusters;
   } else {
      plot_num <- plot_num + pathway_clusters_n;
   }

   revert_hm_padding();

   invisible(ret_vals);
}
