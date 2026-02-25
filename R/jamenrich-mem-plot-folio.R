

#' Multienrichment folio of summary plots
#'
#' Multienrichment folio of summary plots
#'
#' `prepare_folio()` and `mem_plot_folio()` both prepare data
#' visualizations from `Mem` data input.
#' However, `prepare_folio()` does not render figures,
#' while `mem_plot_folio()` does render each resulting figure.
#' 
#' Multiple figures can be added to a single PDF file
#' using `pdf(file, onefile=TRUE)` or `cairo_pdf(filename, onefile=TRUE)`.
#' 
#' The data are returned as `MemPlotFolio-class` which can be used to
#' create figures.
#' 
#' Pathways are clustered using the gene-pathway incidence matrix,
#' then used to define pathway clusters.
#' This step can be customized by supplying `pathway_column_split`
#' as a `list` of `character` vectors containing
#' pathway (set) names.
#' 
#' 
#' ## Pathway Clustering
#' 
#' Pathways are hierarchical clustered by default using `amap::hcluster()`,
#' with `column_method="euclidean"`.
#' The resulting hclust/dendrogram is split using `pathway_column_split`
#' with an `integer` number of sub-clusters. The default (NULL) will
#' determine a default based upon the total Sets in the analysis, and
#' is intended to be customized by the analyst.
#' 
#' The data to be clustered include the `memIM()` matrix of Genes (rows)
#' and Sets (columns), itself an incidence matrix. It could therefore
#' be clustered using `column_method="binary"`. However, the default
#' behavior is to append `-log10()` values from the `enrichIM` matrix,
#' in order to cluster both the incidence matrix, and the enrichment
#' P-values together. The matrix data are weighted relative to one another
#' using `enrich_im_weight=0.3`. Note that the enrichment P-value matrix
#' applies the `p_cutoff` such that values that do not meet the threshold
#' are considered '1' and become '0' with `-log10()` transformation.
#' 
#' In order to cluster the gene-pathway incidence matrix with no
#' influence of the enrichment P-value matrix, use `enrich_im_weight=0`.
#' In this case it often works well to use `column_method="binary"`.
#' 
#' Caveat:
#' Clustering is itself "imperfect": The results are not definitive,
#' and there is generally no one "true" answer.
#' However, the output is mainly intended to help organize
#' information already present in the data, and not to declare
#' ground truth.
#' As such, the results are considered stable for the methods and
#' parameters used, and the interpretation is performed in that context.
#' 
#' The clustering methods and parameters are included in the legend
#' alongside the `GenePathHeatmap()` and the `EnrichmentHeatmap()`
#' for clarity.
#' 
#' Our experience is that clusters do not need to be perfect to be
#' useful, informative, and valid.
#' 
#' ### Custom Cluster Function
#' It is also possible to supply a custom function to perform clustering,
#' by using argument `cluster_columns` which is passed via '...'
#' to the underlying function `mem_gene_path_heatmap()`. This function
#' should take a `matrix` of `numeric` values as input, and return
#' either hclust or dendrogram output, suitable for use with `cutree()`.
#' 
#' An interesting custom function is `cluster::daisy(x, method='gower')`
#' whose output can be converted to `hclust`. The 'gower' method accepts
#' mixed input types such as ordinal, numeric, signed, or other values.
#' 
#' When supplying a custom function via `cluster_columns`, we recommend
#' defining a label using argument `column_method`, for example
#' `custom_method='gower'`. This label will appear in the legend of
#' the resulting visualizations.
#' 
#' The default `amap::hcluster()` was chosen because it is fast, includes
#' both the `dist()` and `hclust()` steps in one function, and
#' includes novel distance metrics such as `method='correlation'`
#' which have proven useful in other contexts.
#' 
#' 
#' ### Pathway Grouping
#' 
#' Pathways can be grouped manually by defining `pathway_column_split`
#' as a `list` of `character` vectors, where those vectors match
#' pathway names returned by `sets(Mem)` - which are also colnames
#' of the `memIM()` incidence matrix.
#' This step may also use a subset of sets, for example including only
#' pathways of relevance to the downstream analysis.
#' 
#' 
#' ## Gene Clustering
#' 
#' Genes in the Gene-Pathway Heatmap are clustered using similar logic
#' as described for Pathways. Note that gene clusters are not used
#' directly in other analysis steps (yet).
#' 
#' However, `score_gene_path_clusters()` offers some metric to identify
#' "hot spots" where a majority of genes in a gene cluster are
#' represented in the majority of pathways in a pathway cluster.
#' These "hot spots" are subject to further study, as potential
#' hubs for interpreting core components across similar biological
#' pathways.
#' 
#' 
#' Custom gene groups can be supplied with `gene_row_split` as
#' a `list` of `character` vectors whose values match `genes(Mem)`.
#' 
#' 
#' ## Enrichment Heatmap
#' 
#' The Gene-Pathway Heatmap is used to define the pathway clusters
#' (or pathway groups), and therefore the order of pathways shown
#' in the Gene-Pathway Heatmap.
#' 
#' The same pathway order is used with the Enrichment Heatmap,
#' for consistency.
#' 
#' An important point is that the Enrichment Heatmap does not
#' cluster nor order pathways using enrichment P-values,
#' and instead re-uses the dendrogram (if relevant) from the
#' Gene-Pathway Heatmap.
#' 
#' Enrichment P-values are not specific to the genes
#' contained in each pathway, and therefore clustering by enrichment
#' P-value can be misleading by grouping pathways together which
#' are not related by the underlying biological data.
#' 
#' 
#' ## Cnet Collapsed Plot
#' 
#' The pathways in each pathway cluster (or pathway group) are collapsed
#' into one virtual pathway, then used to create a Concept network (Cnet)
#' plot. This Cnet plot is intended to balance the motivation to
#' show complete pathway enrichment data, with the motivation to
#' reduce redundant information.
#' 
#' The pathway clusters (or pathway groups) are the key components
#' of the Cnet Collapsed output. Therefore adjusting the clustering
#' options `enrich_im_weight` or `column_method` are the most common
#' ways to influence and optimize the outcome.
#' 
#' 
#' ## Cnet Exemplar Plot
#' 
#' There are two general paradigms for Concept networks (Cnets).
#' 
#' 1. Pathways in clusters: `CnetCollapsed()`
#' 2. Exemplar Pathways per cluster: `CnetExemplar()`
#' 
#' The Cnet Exemplar plot includes one pathway for each pathway cluster
#' `num=1` but may include two `num=2` or three `num=3` pathways
#' per cluster if relevant. The exemplar pathway is intended as a
#' representative of the cluster, which can help simplify the number
#' of genes displayed, while also focusing the meaning of those genes
#' to the pathway shown.
#' 
#' This option may be preferred when pathway clusters are not visibly
#' distinctive, for example when one pathway cluster does not appear
#' to have a common set of genes shared with other pathways in the same
#' cluster.
#' 
#' The alternative occurs when pathways in each cluster do share
#' many of the same genes, and in this case the `CnetCollapsed()`
#' plot may be more effective.
#' 
#' Finally, a custom Cnet exemplar plot may be warranted when
#' displaying a few pathways with particular relevance to a
#' research study or experiment. The steps are described in the
#' README.Rmd and README html page for the multienrichjam package:
#' [https://jmw86069.github.io/multienrichjam](https://jmw86069.github.io/multienrichjam)
#' 
#' 
#' ## Cnet Cluster Plot
#' 
#' The last set of Concept network (Cnet) plots include all pathways
#' in each individual pathway cluster. These plots are useful for
#' pathway clusters which may have different functional sub-components,
#' or when pathways in a cluster do not appear to be very cohesive.
#' 
#' The Cnet Cluster plot helps display the relationship of pathways within
#' one pathway cluster, which may help provide the basis for interpreting
#' the results.
#' 
#' ## Order of Plots
#' 
#' When using `mem_plot_folio()` plots are produced in a specific order.
#' The argument `do_which` may also be useful to focus the output
#' only on a particular plot, thereby skipping the preparatory
#' steps used for other plot types.
#' 
#' 1. **Enrichment Heatmap** (Plot #1), `EnrichmentHeatmap()`.
#' Detailed arguments are described in `mem_enrichment_heatmap()`.
#' 2. **Gene-Pathway Heatmap** (Plot #2) `GenePathHeatmap()`.
#' Detailed arguments are described in `mem_gene_path_heatmap()`.
#' 3. **Cnet Cluster Plot** (Plots #3,#4,#5) `CnetCollapsed()`.
#' Detailed arguments are described in `collapse_mem_clusters()` and
#' `mem2cnet()`.
#'
#'    * For `mem_plot_folio()` three styles are produced,
#'    with different node labeling strategies.
#'    * Plot #3 uses pathway cluster titles.
#'    * Plot #4 uses pathway names.
#'    * Plot #5 uses pathway names, and hides gene labels.
#'
#' 4. **Cnet Exemplar Plots** (Plots #6,#7,#8) `CnetExemplar()`.
#' The number of exemplars uses `exemplar_range=c(1, 2, 3)`.
#' 5. **Cnet Cluster Plots** (Plots #9,#10,#11,etc.) `CnetCluster()`.
#' Detailed arguments are described in `mem2cnet()`.
#' 
#' ## Final Points
#' 
#' To define your own pathway cluster labels, define `pathway_column_title`
#' as a vector with length equal to `pathway_column_split`. These labels
#' become network node labels in subsequent plots, and in the
#' resulting `igraph` object.
#' 
#' Metadata are stored in `metadata(Mpf)` for MemPlotFolio output, and
#' includes:
#' 
#' * 'colorV': the color vector used
#' * 'hasDirection': whether the input data contained directionality,
#' defined by any negative value in 'geneIMdirection', 'enrichIMdirection',
#' or any non-blank column for headers(Mem) 'directionColname'.
#'
#' @family multienrichjam core functions
#'
#' @returns `MemPlotFolio` object using `invisible()`, containing each
#'    plot object enabled by the argument `do_which`.
#'    The `MemPlotFolio-class` data are accessible using common functions:
#'    * `EnrichmentHeatmap()`
#'    * `GenePathHeatmap()`
#'    * `CnetCollapsed()`
#'    * `CnetExemplar()`
#'    * `CnetCluster()`
#'    * `Clusters()`
#'    * `GeneClusters()`
#'    * `thresholds()`
#'    * `metadata()`
#'    * `Caption()`
#'    * `CaptionLegendList()`
#'
#' @param mem `Mem` or `list` object from `multiEnrichMap()`.
#' @param do_which `integer` vector of plots to produce.
#'    When `do_which` is `NULL`, then all plots are produced.
#'    This argument is intended to produce only a subset of plots.
#' @param mpf `MemPlotFolio`, default NULL, used only to re-apply the same
#'    settings as another `MemPlotFolio`.
#'    * Note: `pathway_column_split` and `gene_row_split` are also
#'    taken from `thresholds(mpf)`, unless defined in function arguments.
#'    * Note: All other values in `thresholds(mpf)` are used,
#'    and not taken from corresponding function arguments.
#' @param p_cutoff `numeric` value, default NULL is taken from `mem`,
#'    indicating the enrichment P-value threshold.
#' @param p_floor `numeric` with the lowest enrichment P-value
#'    used in the color gradient on the Enrichment Heatmap.
#'    The purpose is to prevent very low P-values from shifting the
#'    color gradient too far from the `p_cutoff` causing those colors
#'    to be pale and nearly white.
#' @param main `character` string used as a title on Cnet plots.
#' @param use_raster `logical` default FALSE, deprecated,
#'    whether to use raster heatmaps, passed to `ComplexHeatmap::Heatmap()`.
#'    * Note that `use_raster=TRUE` may produce visual artifacts especially
#'    with argument `colorize_by_gene=TRUE` in `mem_gene_path_heatmap()`.
#'    Changing this argument is no longer supported
#' @param min_gene_ct,min_set_ct `integer` values passed to
#'    `mem_gene_path_heatmap()`. The `min_gene_ct` requires each set
#'    to contain `min_gene_ct` genes, and `min_set_ct` requires each gene
#'    to be present in at least `min_set_ct` sets.
#' @param min_set_ct_each `integer`, default NULL, minimum genes per set
#'    in at least one enrichment test.
#'    * Default NULL uses `thresholds(mem)$min_count` to use the same
#'    criteria.
#'    * The distinction from `min_set_ct` is that this threshold
#'    requires this number of genes in one enrichment,
#'    while `min_set_ct` applies the threshold to the
#'    combined multi-enrichment data.
#' @param column_method,row_method `character` arguments passed to
#'    `ComplexHeatmap::Heatmap()` which indicate the distance method used
#'    to cluster columns and rows, respectively.
#' @param cluster_columns,cluster_rows `logical`, default NULL, whether to
#'    cluster columns (Sets) and rows (Genes), respectively.
#'    When NULL it uses default clustering with `amap::hcluster()` and
#'    applies `column_method` or `row_method`, respectively.
#' @param exemplar_range `integer` vector (or `NULL`) used to create Cnet
#'    exemplar plots, using this many exemplars per cluster.
#' @param pathway_column_split,gene_row_split `integer` value passed
#'    as `column_split` and `row_split`, respectively, to
#'    `mem_gene_path_heatmap()`, indicating the number of pathway
#'    clusters, and gene clusters, to create in the gene-pathway heatmap.
#'    When either value is `NULL` then auto-split logic is used.
#'    * `integer` input will define this number of clusters
#'    * `list` input should be named by pathway cluster, in order the
#'    clusters will be split; and have `character` vectors that
#'    match pathways (sets) in the 'mem' input data.
#'    * `character` input should be named by pathway (set) and have
#'    values with the name of the pathway cluster. It can be `factor`.
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
#'    Default `TRUE` for `mem_plot_folio()`, and `FALSE` for `prepare_folio()`.
#'    In either case, plot data are created and returned, but `do_plot=TRUE`
#'    will draw each plot on a unique page, suitable for use with
#'    PDF output and `onefile=TRUE` for example.
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
#' @examples
#' data(Memtest)
#' Mpf <- prepare_folio(Memtest)
#' GenePathHeatmap(Mpf)
#' 
#' # take a subset of clusters and re-Mpf
#' MpfSub <- prepare_folio(Memtest,
#'    pathway_column_split=Clusters(Mpf)[c("A", "C")])
#' GenePathHeatmap(MpfSub)
#' 
#' # adjust set names
#' Memtest2 <- fixSetLabels(Memtest)
#' Mpf2 <- prepare_folio(Memtest2,
#'    row_names_rot=0,
#'    column_names_rot=60)
#' GenePathHeatmap(Mpf2)
#' 
#' # Now add a temporary buffer for the heatmap dimnames
#' with_ht_opts(list(DIMNAME_PADDING=grid::unit(c(2), "mm")), {
#'    GenePathHeatmap(Mpf2)
#' })
#' 
#' @export
mem_plot_folio <- function
(mem,
 do_which=NULL,
 mpf=NULL,
 p_cutoff=NULL,
 p_floor=1e-10,
 main="",
 use_raster=FALSE,
 min_gene_ct=1,
 min_set_ct=1,
 min_set_ct_each=NULL,
 column_method="euclidean",
 cluster_columns=NULL,
 row_method="euclidean",
 cluster_rows=NULL,
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
 returnType=c("MemPlotFolio",
    "list"),
 do_plot=TRUE,
 verbose=FALSE,
 ...)
{
   # validate arguments
   returnType <- match.arg(returnType);
   
   # accept Mem input, convert to list for now
   Mem <- NULL;
   if (inherits(mem, "Mem")) {
      Mem <- mem;
      mem <- Mem_to_list(Mem)
   } else {
      Mem <- list_to_Mem(mem);
   }
   # validate arguments
   # fill some default arguments
   if (length(p_cutoff) == 0) {
      if ("p_cutoff" %in% names(thresholds(Mem))) {
         p_cutoff <- thresholds(Mem)$p_cutoff;
      } else if ("cutoffRowMinP" %in% names(thresholds(Mem))) {
         p_cutoff <- thresholds(Mem)[["cutoffRowMinP"]];
      } else {
         p_cutoff <- 1;
      }
   }
   if (length(min_set_ct_each) == 0) {
      if ("min_count" %in% names(thresholds(Mem))) {
         min_set_ct_each <- thresholds(Mem)$min_count;
      } else {
         min_set_ct_each <- 1;
      }
   }
   
   # do_which
   if (length(do_which) == 0) {
      if (isTRUE(do_plot)) {
         do_which <- seq_len(50)
      } else {
         # if not plotting, skip the two variant Cnet plots
         do_which <- setdiff(seq_len(50), c(4, 5));
      }
   }
   
   # prepare thresholds
   thresholds <- list(
      min_gene_ct=min_gene_ct,
      min_set_ct=min_set_ct,
      min_set_ct_each=min_set_ct_each,
      enrich_im_weight=enrich_im_weight,
      gene_im_weight=gene_im_weight,
      column_method=column_method,
      cluster_columns=cluster_columns,
      column_split=pathway_column_split,
      column_title=pathway_column_title,
      row_method=row_method,
      cluster_rows=cluster_rows,
      row_split=gene_row_split,
      row_title=gene_row_title,
      p_cutoff=p_cutoff,
      p_floor=p_floor,
      byCols=byCols,
      do_which=do_which,
      do_plot=do_plot,
      rotate_heatmap=rotate_heatmap,
      cluster_color_min_fraction=cluster_color_min_fraction,
      row_anno_padding=row_anno_padding,
      column_anno_padding=column_anno_padding
   )
   
   # Re-use mpf if supplied
   if (inherits(mpf, "MemPlotFolio")) {
      mpf_thresholds <- thresholds(mpf);
      if (verbose) {
         jamba::printDebug("Applying thresholds(mpf) named: ",
            names(mpf_thresholds));
      }
      thresholds[names(mpf_thresholds)] <- mpf_thresholds;
      # 0.0.105.900: inherit Clusters(mpf)
      if (length(pathway_column_split) == 0) {
         pathway_column_split <- Clusters(mpf);
         thresholds$column_split <- pathway_column_split;
      }
      if (length(gene_row_split) == 0) {
         gene_row_split <- GeneClusters(mpf);
         thresholds$row_split <- gene_row_split;
      }
   }
   
   # hasDirection is TRUE when any criteria are met:
   # 1. enrichIMdirection has any negative value, OR
   # 2. geneIMdirection has any negative value, OR
   # 3. directionColname is defined
   hasDirection <- FALSE;
   if (any(enrichIMdirection(Mem) < 0) ||
   		any(geneIMdirection(Mem) < 0) ||
   		!is.null(headers(Mem)[["directionColname"]])) {
   	hasDirection <- TRUE;
   }

   # metadata
   metadata <- list(
      colorV=Mem@colorV,
   	hasDirection=hasDirection
   );
   
   return_mpf <- function(ret_vals, returnType="MemPlotFolio") {
      if ("MemPlotFolio" %in% returnType) {
         return(list_to_MemPlotFolio(ret_vals))
      }
      return(ret_vals)
   }
   ret_vals <- list();
   # add thresholds and parameters used
   ret_vals$thresholds <- thresholds;
   ret_vals$metadata <- metadata;
   
   plot_num <- 0;

   #############################################################
   # Apply optional heatmap padding parameters
   rowpad <- NULL;
   columnpad <- NULL;
   use_ht_opt <- list();
   if (length(row_anno_padding) > 0) {
      if (!grid::is.unit(row_anno_padding)) {
         row_anno_padding <- grid::unit(row_anno_padding, "mm");
      }
      rowpad <- ComplexHeatmap::ht_opt("ROW_ANNO_PADDING")
      use_ht_opt$ROW_ANNO_PADDING <- row_anno_padding;
   }
   if (length(column_anno_padding) > 0) {
      if (!grid::is.unit(column_anno_padding)) {
         column_anno_padding <- grid::unit(column_anno_padding, "mm");
      }
      columnpad <- ComplexHeatmap::ht_opt("COLUMN_ANNO_PADDING")
      use_ht_opt$COLUMN_ANNO_PADDING <- column_anno_padding;
   }
   if (length(use_ht_opt) > 0) {
      local_ht_opts(use_ht_opt)
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
   gphm_gn_order <- NULL;

   if (rotate_heatmap) {
      gp_hm_hclust <- gp_hm@row_dend_param$obj;
      if (length(gp_hm_hclust) == 0) {
         gp_hm_hclust <- gp_hm@row_dend_param$fun(t(gp_hm@matrix));
      }
      use_sets <- rownames(gp_hm@matrix);
      suppressWarnings(gphm_pw_order <- jamba::heatmap_row_order(gp_hm))
      suppressWarnings(gphm_gn_order <- jamba::heatmap_column_order(gp_hm))
   } else {
      gp_hm_hclust <- gp_hm@column_dend_param$obj;
      if (length(gp_hm_hclust) == 0) {
         gp_hm_hclust <- gp_hm@column_dend_param$fun(t(gp_hm@matrix));
      }
      suppressWarnings(gphm_pw_order <- jamba::heatmap_column_order(gp_hm))
      suppressWarnings(gphm_gn_order <- jamba::heatmap_row_order(gp_hm))
   }
   ret_vals$clusters_mem <- gphm_pw_order;
   ret_vals$gene_clusters_mem <- gphm_gn_order;
   

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
   # add title to attributes
   attr(gp_hm, "title_list") <- list(
      column_title=main,
      column_title_gp=grid::gpar(fontsize=18))
   
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
   # if (rotate_heatmap) {
   #    suppressWarnings(clusters_mem <- jamba::heatmap_row_order(gp_hm))
   # } else {
   #    suppressWarnings(clusters_mem <- jamba::heatmap_column_order(gp_hm))
   # }
   # ret_vals$clusters_mem <- clusters_mem;
   clusters_mem <- ret_vals$clusters_mem;
   
   ## Get number of pathway clusters
   pathway_clusters_n <- length(clusters_mem);
   if (verbose) {
      jamba::printDebug("mem_plot_folio(): ",
         c("Defined ", pathway_clusters_n, " pathway clusters."),
         sep="");
   }
   
   ## If no other plots are necessary, exit here
   if (length(do_which) > 0 && !any(do_which > plot_num)) {
      return(invisible(return_mpf(ret_vals, returnType)));
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
            repulse=repulse,
            spread_labels=TRUE,
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
      # 0.0.104.900: this step may not be necessary
      # - if the collapse_mem_clusters() already created layout,
      #   and calling mem2cnet() should already remove blanks
      if (FALSE) {
         cnet_collapsed <- tryCatch({
            cnet_collapsed <- subsetCnetIgraph(
               cnet_collapsed,
               remove_blanks=TRUE,
               repulse=repulse,
               verbose=verbose > 1);
         }, error=function(e){
            if (verbose) {
               jamba::printDebug("mem_plot_folio(): ",
                  "subsetCnetIgraph() error during ",
                  "subsetCnetIgraph(..., remove_blanks=TRUE):");
               print(e);
            }
            # try again without removing blanks
            cnet_collapsed <- tryCatch({
               cnet_collapsed <- subsetCnetIgraph(
                  cnet_collapsed,
                  remove_blanks=FALSE,
                  repulse=repulse,
                  verbose=verbose > 1);
            }, error=function(e2){
               jamba::printDebug("mem_plot_folio(): ",
                  "subsetCnetIgraph() error during ",
                  "subsetCnetIgraph(), skipping this step.");
               print(e2);
               cnet_collapsed;
            })
         })
         }
      if (length(edge_color) > 0) {
         igraph::E(cnet_collapsed)$color <- edge_color;
      }
      plot_num <- plot_num + 1;
      cnet_title <- "Cnet plot using collapsed clusters";
      cnet_collapsed <- igraph::set_graph_attr(cnet_collapsed,
         name="title",
         value=cnet_title);
      ret_vals$cnet_collapsed <- cnet_collapsed;
      if (length(do_which) == 0 || plot_num %in% do_which) {
         if (verbose) {
            jamba::printDebug("mem_plot_folio(): ",
               c("plot_num ", plot_num, ": "),
               c("Cnet collapsed ", "with gene and cluster labels"),
               sep="");
         }
         ## Draw Cnet collapsed
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
      ## Draw Cnet collapsed with top n labels
      #isset <- (V(cnet_collapsed)$nodeType %in% "Set");
      if ("set_labels" %in% igraph::vertex_attr_names(cnet_collapsed)) {
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
         if ("list" %in% returnType) {
            # 0.0.104.900: only return when using 'list' format
            ret_vals$cnet_collapsed_set <- cnet_collapsed;
         }
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
            if ("set_labels" %in% igraph::vertex_attr_names(cnet_collapsed)) {
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
         if ("list" %in% returnType) {
            # 0.0.104.900: only return when using 'list' format
            ret_vals$cnet_collapsed_set2 <- cnet_collapsed;
         }
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
            "Preparing cnet for subsetting with mem2cnet().");
      }
      # spread_labels=FALSE because no layout is required yet
      cnet <- mem2cnet(mem,
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
   if (length(do_which) == 0 ||
         any((plot_num + seq_len(pathway_clusters_n)) %in% do_which)) {
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

   invisible(return_mpf(ret_vals, returnType));
}

#' @inheritParams mem_plot_folio
#' @rdname mem_plot_folio
#' @export
prepare_folio <- function
(mem,
 do_which=NULL,
 mpf=NULL,
 p_cutoff=NULL,
 p_floor=1e-10,
 main="",
 use_raster=FALSE,
 min_gene_ct=1,
 min_set_ct=1,
 min_set_ct_each=NULL,
 column_method="euclidean",
 cluster_columns=NULL,
 row_method="euclidean",
 cluster_rows=NULL,
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
 returnType=c("MemPlotFolio",
 	"list"),
 do_plot=FALSE,
 verbose=FALSE,
 ...)
{
   mem_plot_folio(mem=mem,
   	do_which=do_which,
   	mpf=mpf,
   	p_cutoff=p_cutoff,
   	p_floor=p_floor,
   	main=main,
   	use_raster=use_raster,
   	min_gene_ct=min_gene_ct,
   	min_set_ct=min_set_ct,
   	min_set_ct_each=min_set_ct_each,
   	column_method=column_method,
   	cluster_columns=cluster_columns,
   	row_method=row_method,
   	cluster_rows=cluster_rows,
   	exemplar_range=exemplar_range,
   	pathway_column_split=pathway_column_split,
   	pathway_column_title=pathway_column_title,
   	gene_row_split=gene_row_split,
   	gene_row_title=gene_row_title,
   	edge_color=edge_color,
   	cex.main=cex.main,
   	cex.sub=cex.sub,
   	row_cex=row_cex,
   	column_cex=column_cex,
   	max_labels=max_labels,
   	max_nchar_labels=max_nchar_labels,
   	include_cluster_title=include_cluster_title,
   	repulse=repulse,
   	use_shadowText=use_shadowText,
   	color_by_column=color_by_column,
   	style=style,
   	enrich_im_weight=enrich_im_weight,
   	gene_im_weight=gene_im_weight,
   	colorize_by_gene=colorize_by_gene,
   	cluster_color_min_fraction=cluster_color_min_fraction,
   	byCols=byCols,
   	edge_bundling=edge_bundling,
   	apply_direction=apply_direction,
   	rotate_heatmap=rotate_heatmap,
   	row_anno_padding=row_anno_padding,
   	column_anno_padding=column_anno_padding,
   	returnType=returnType,
   	do_plot=do_plot,
   	verbose=verbose,
   	...)
}
