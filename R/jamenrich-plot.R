
#' MultiEnrichment Heatmap of Genes and Pathways
#'
#' MultiEnrichment Heatmap of Genes and Pathways
#'
#' This function takes the `Mem` output from
#' `multiEnrichMap()` and creates a gene-by-pathway incidence
#' matrix heatmap, using `ComplexHeatmap::Heatmap()`.
#' The major output of this function is to define pathway
#' clusters which influences other figures produced by
#' `mem_plot_folio()`: the enrichment heatmap; Cnet cluster plots;
#' Cnet exemplars.
#' 
#' It uses three basic sources of data to annotate the heatmap:
#'
#' 1. `memIM` the gene-set incidence matrix
#' 2. `geneIM` the gene incidence matrix by dataset
#' 3. `enrichIM` the pathway enrichment P-value matrix by dataset
#' 
#' ## Pathway and gene clusters
#' 
#' When `column_split` (pathways) or `row_split` (genes) are
#' not defined, a reasonable number is assigned to split columns
#' and rows, respectively. The specific number is controlled by
#' defining an integer number. The splits are named using
#' `column_title` and `row_title`, which defaults `LETTERS` and
#' `letters`, respectively.
#' 
#' Custom pathway clusters (columns) can be defined by defining
#' `column_split` as a vector named by values in `sets(Mem)`,
#' with values used to define the name for each cluster.
#' When provided as a `factor` it will honor the factor level order.
#' 
#' When a subset of `sets(Mem)` are provided, it will also subset
#' the `Mem` object to show only the matching sets, however it does
#' not re-apply row and column filtering described above.
#' It is recommended to use argument 'sets' to subset the pathway
#' gene sets upfront, then use `column_split` with names that
#' match 'sets'.
#' 
#' Alternatively, `column_split` can be supplied as a `data.frame`
#' with rownames that match `sets(Mem)`, in which case column values
#' are combined using `jamba::pasteByRowOrdered()` which also keeps
#' factor level order.
#' 
#' Custom gene clusters (rows) can be defined with `row_split` using
#' the same mechanism as with `column_split` described above,
#' with names that match `genes(Mem)` as appropriate.
#' 
#' ### Column pathway clustering
#' 
#' Columns (pathways) are clustered using a combination of the gene-pathway
#' incidence matrix, and the -log10 P-values shown along the
#' pathway axis (usually columns). The purpose is to allow the
#' enrichment P-value to influence the clustering together with
#' the gene content. The balance is adjusted using
#' `enrich_im_weight`, default 0.3, and where 0.0 will ignore the
#' enrichment P-values during clustering. Higher values will tend
#' to create clusters that represent shared/unique *significant*
#' pathways, and less determined based solely on the gene content.
#' 
#' Note that the enrichment P-value matrix used during clustering
#' is adjusted by using `p_cutoff` and `p_floor`.
#' Values above the `p_cutoff` are converted to 1 so they do not
#' influence clustering. Values below `p_floor` are converted to
#' `p_floor` so they influence clustering only at the level of
#' other values at `p_floor`.
#' 
#' The default `cluster_columns=TRUE` will employ `amap::hcluster()`
#' as a convenient and efficient one-step distance-hierarchical clustering
#' approach, using `cluster_method` as the distance method. A custom
#' function can be supplied with `cluster_columns` as long as the
#' output is 'hclust' or can be coerced to 'hclust'.
#' 
#' The clustering distance method `column_method` uses default 'binary',
#' which treats any non-zero value as 1.
#' In this case, the relative weight is not important since all non-zero
#' values are equivalent.
#' However, `column_method='euclidean'`, current default in `mem_plot_folio()`
#' may improve the output when adjusting the relative weight of
#' enrichment P-values with incidence matrix.
#' 
#' ### Gene clustering
#' 
#' Rows (genes) are clustered using a combination of the gene-pathway
#' incidence matrix, and the 'geneIM' or 'geneIMdirection' matrix data
#' shown, as defined with `gene_annotations`. The relative weight of
#' these matrices is controlled with `gene_im_weight` with default 0.5,
#' which gives equal weight to each matrix. By default, when directional
#' gene values are shown the directional matrix is used with clustering.
#' 
#' The default `cluster_rows=TRUE` will employ `amap::hcluster()` as
#' described for Pathway clustering, or `cluster_rows` can be a custom
#' function that produces 'hclust' or can be coerced to 'hclust'.
#' 
#' When `row_method='binary'` the relative weight of gene incidence matrix
#' and the pathway-gene matrix is not important, since all non-zero
#' values are treated as 1 during clustering. The `mem_plot_folio()`
#' default uses 'euclidean' to improve the effect of the relative
#' weight of these two matrices.
#' 
#' Gene clusters are not often used for downstream analysis, for
#' example they do not form clusters in the Cnet plots, and are not
#' (yet) used for other analysis in multienrichjam.
#' However, gene clusters are quite useful when interpreting
#' pathway-gene data.
#' 
#' It is helpful in practice to refer to a gene cluster by name:
#' "The genes in cluster 'd' all appear to be cytokines."
#' (Also maybe there should be better names, but that's for
#' another time.)
#' 
#' Gene clusters may form what we call "hot spots", where
#' most of a gene cluster is colorized and associated with
#' one or more pathway clusters.
#' A hot spot indicates a set of genes shared across multiple pathways,
#' a "core" set of genes which may have serve an important functional
#' basis in several pathways.
#' 
#' An example might be mitogen-activated protein kinase (MAPK) genes,
#' which typically involve a multi-step kinase call signaling
#' cascade. Pathways which involve one MAPK would very often
#' involve each MAPK at subsequent steps in the cascade.
#' In fact, MAPK genes are often the core signaling mechanism of
#' numerous apparently unrelated pathways - we refer to it as the
#' "internal wiring" of the signaling in a cell, as an analogy to
#' a building which may use wiring to pass along any number of messages.
#' Different cell types may employ the MAPK cascade to send a message,
#' and this message may have different meaning across cell types,
#' and indeed may have meaning based upon the cell state.
#' 
#' ## Filtering
#' 
#' A subset of pathways can be defined with argument `sets`,
#' which may be useful to prepare this plot for a set of
#' "exemplar pathways" which represent key pathways of interest
#' for a study.
#' 
#' Similar for genes with argument `genes`, however this option
#' is less commonly used.
#' 
#' Pathways and genes can be subset by number of occurrences
#' of each, for example:
#' 
#' * `min_gene_ct`: minimum occurrences of a gene across pathways,
#' which also means the number of pathways in which a gene
#' occurs.
#' * `min_set_ct`: minimum occurrences of a set (pathway) across genes,
#' which means the number of genes in the set across all enrichments.
#' * `min_set_ct_each`: minimum occurrences of a set (pathway) across genes.
#' It requires at least `min_set_ct_each` genes in at least one enrichment,
#' which is consistent with `min_count` in `multiEnrichMap()`.
#' 
#' When pathways are filtered by `min_gene_ct`, `min_set_ct`,
#' and `min_set_ct_each`, the order of operations is as follows:
#'
#' 1. `min_set_ct_each`, `min_set_ct` - these filters are applied
#' before filtering genes, in order to ensure all genes are present
#' from the start.
#' 2. `min_gene_ct` - genes are filtered after pathway filtering,
#' in order to remove pathways which were not deemed "significant"
#' based upon the required number of genes. Only after those pathways
#' are removed can the number of occurrences of each gene be judged
#' appropriately.
#' 
#'
#' @family jam plot functions
#'
#' @param mem `Mem` or `list` object created by `multiEnrichMap()`.
#' @param genes `character` vector of genes to include in the heatmap,
#'    all other genes will be excluded.
#' @param sets `character` vector of sets (pathways) to include in the heatmap,
#'    all other sets will be excluded.
#' @param min_gene_ct `integer` minimum number of occurrences of each gene
#'    across the pathways, all other genes are excluded.
#' @param min_set_ct `integer` minimum number of genes required for each set,
#'    all other sets are excluded.
#' @param min_set_ct_each `integer` minimum number of genes required
#'    for each set, required for at least one enrichment test.
#' @param column_fontsize,row_fontsize `numeric`
#'    passed as `fontsize` to `ComplexHeatmap::Heatmap()`
#'    to define a specific fontsize for column and row labels. When
#'    `NULL` the nrow/ncol of the heatmap are used to infer a reasonable
#'    starting point fontsize, which can be adjusted with `column_cex`
#'    and `row_cex`.
#' @param row_method,column_method `character` string with distance method,
#'    default is 'binary' which considers all non-zero values to be 1.
#'    It is used for row and column hierarchical clustering
#'    by `amap::hcluster()`. It offers more methods than `hclust()`
#'    hence its use here.
#' @param enrich_im_weight `numeric` value between 0 and 1 (default 0.3),
#'    the relative weight of enrichment `-log10 P-value` and overall
#'    gene-pathway incidence matrix when clustering pathways.
#'    * When `enrich_im_weight=0` then only the gene-pathway incidence
#'    matrix is used for pathway clustering.
#'    * When `enrich_im_weight=1` then only the pathway significance
#'    (`-log10 P-value`) is used for pathway clustering.
#'    * The default `enrich_im_weight=0.3` balances the combination
#'    of the enrichment P-value matrix, with the gene-pathway incidence
#'    matrix.
#' @param gene_im_weight `numeric` value between 0 and 1 (default 0.5),
#'    the relative weight of the `mem$geneIM` gene incidence matrix,
#'    and overall gene-pathway incidence matrix when clustering genes.
#'    * When `gene_im_weight=0` then only the gene-pathway incidence
#'    matrix is used for gene clustering.
#'    * When `gene_im_weight=1` then only the gene incidence matrix
#'    (`mem$geneIM`) is used for gene clustering.
#'    * The default `_im_weight=0.5` balances the gene incidence matrix
#'    with the gene-pathway incidence matrix, giving each matrix equal weight
#'    (since values are typically all `(0, 1)`.
#' @param gene_annotations `character` string indicating which annotation(s)
#'    to display alongside the gene axis of the heatmap.
#'    * Default is `"im", "direction"`, and `"default"` which will hide
#'    the `"direction"` if all non-zero values are positive.
#'    * When `"im"` is present, the colored incidence matrix is displayed.
#'    * When `"direction"` is present, the directional matrix is displayed
#'    using colors defined by `colorjam::col_div_xf(1.2)`.
#'    * When both `"im"` and `"direction"` are present, they are displayed
#'    in the order defined.
#'    * When no values are given, the gene annotation is not displayed.
#' @param annotation_suffix `character` vector named by values permitted
#'    by `gene_annotations`, with optional suffix to add to the annotation
#'    labels. For example, as by default, it may be helpful to add
#'    "hit" or "dir" to distinguish the enrichment labels.
#' @param name `character` value passed to `ComplexHeatmap::Heatmap()`,
#'    used as a label above the heatmap color legend. Default `NULL`
#'    uses "Gene Hits by Enrichment".
#' @param p_cutoff `numeric` value of the enrichment P-value cutoff,
#'    above which P-values are not colored, and are therefore white.
#'    The enrichment P-values are displayed as an annotated heatmap
#'    at the top of the main heatmap. Any cell that has a color meets
#'    at least the minimum P-value threshold. The default
#'    is taken from input `mem` for
#'    consistency with the input multienrichment analysis.
#' @param column_split,row_split row and column split, default NULL,
#'    detects an appropriate number of clusters.
#'    Passed to `ComplexHeatmap::Heatmap()`.
#'    * To turn off split use `1`.
#'    * To specify a fixed number of clusters, use an `integer` value.
#'    The value may be changed if the underlying data does not support
#'    that number of clusters.
#'    * To specify fixed clusters by name, use an atomic vector
#'    named by the rownames `genes(Mem)` or colnames `sets(Mem)` with
#'    values which will become the cluster names. When supplying a
#'    `factor` the factor level order will be maintained.
#'    * Alternatively, supply a `data.frame` whose rownames match
#'    the `genes(Mem)` rows or `sets(Mem)` columns, respectively.
#'    In this case, column values are combined using
#'    `jamba::pasteByRowOrdered()` which also maintains factor level
#'    order.
#'    * When supplying either an atomic vector or `data.frame`
#'    the actual names will be used to subset the resulting `Mem` data
#'    if there are fewer names provided than exist in `Mem`.
#'    Note that the `Mem` data are not subset again by other filter
#'    criteria, for example `min_set_ct`, `min_gene_ct`, etc.
#'    * The `column_title` argument is used when
#' @param auto_split `logical` whether to determine clusters when column_split
#'    or row_split is NULL, default TRUE.
#' @param column_title optional `character` string or vector
#'    to display above the column splits. Default uses `LETTERS`
#'    with length to match the number of clusters.
#'    When there is only one column cluster, it is not named unless
#'    `column_title` also has length 1.
#' @param row_title optional `character` string or vector
#'    to display to the left side of the heatmap row splits.
#'    The default uses `letters` to match the number of row clusters.
#'    When there is only one row cluster, it is not named unless
#'    `row_title` also has length 1.
#' @param row_title_rot `numeric` value, default 0, the rotation of
#'    `row_title` text, where `0` is not rotated, and `90` is rotated
#'    90 degrees.
#' @param colorize_by_gene `logical` default TRUE, whether to color the
#'    main heatmap body using the colors from `geneIM` to indicate
#'    each enrichment in which a given gene is involved.
#'    Colors are blended using `colorjam::blend_colors()`,
#'    using colors from `colorV` in the `geneIMcolors(Mem)`.
#' @param na_col `character` string indicating the color to use for
#'    NA or missing values. Typically this argument is only used
#'    when `colorize_by_gene=TRUE`, where entries with no color are
#'    recognized as `NA` by `ComplexHeatmap::Heatmap()`.
#' @param rotate_heatmap `logical` indicating whether the entire heatmap
#'    should be rotated so that pathway names are displayed as rows,
#'    and genes as columns.
#'    When enabled, arguments referring to columns and rows are flipped,
#'    so "column" arguments will continue to affect pathways/sets,
#'    and "row" arguments will continue to affect genes.
#'    This includes `column_method` and `row_method` as of 0.0.90.900.
#'    * Exceptions:
#'    * `row_title_rot` is only applied to rows, due to its purpose.
#'    * `column_names_rot` is only applied to columns, also due to
#'    its purpose.
#' @param colramp `character`, default "Reds", with name of color,
#'    color gradient, or a vector of colors, anything that can be converted
#'    to a color gradient by `jamba::getColorRamp()`.
#' @param column_names_max_height `grid::unit` passed to
#'    `ComplexHeatmap::Heatmap()`. When supplied as `numeric` it is
#'    converted to units in "mm". Default 180 mm.
#' @param column_names_rot `numeric` passed to
#'    `ComplexHeatmap::Heatmap()`.
#' @param show_gene_legend,show_pathway_legend `logical`,
#'    whether to show the gene IM and pathway IM legends, respectively.
#'    * The gene IM legend is FALSE by default, since it only describes the
#'    color used for each column, and is somewhat redundant with the pathway
#'    IM legend.
#'    * The pathway IM default is TRUE, it displays the color scale
#'    including the range of enrichment P-values colorized.
#' @param show_heatmap_legend `numeric` or `logical`, (default 8),
#'    the maximum number of labels to use for the heatmap color legend.
#'    When 'colorize_by_gene' is TRUE, the heatmap legend would include
#'    all possible blended colors using gene IM data, which frankly can
#'    become too much.
#'    * When `logical`, `TRUE` is converted to `8` by default.
#'    * When there are more legend items than than `show_heatmap_legend`
#'    the color legend will only display singlet colors, which means
#'    only one color per individual set defined in colorV.
#'    * This legend can be created and extracted from the
#'    output `Heatmap` object to be displayed independently for
#'    publishable figures if necessary.
#' @param use_raster `logical` passed to `ComplexHeatmap::Heatmap()`,
#'    default FALSE, whether to rasterize the heatmap body.
#'    This option is recommended FALSE when 'colorize_by_gene' is TRUE,
#'    due to the way the rasterization is handled at matrix level.
#'    For very large heatmaps you may try 'colorize_by_gene=FALSE'
#'    and 'use_raster=TRUE' to reduce the figure size with PDF and SVG
#'    output, and will improve visual fidelity in some cases.
#' @param seed `numeric` value passed to `set.seed()` to define a
#'    reproducible random seed during row and column clustering.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to `ComplexHeatmap::Heatmap()`
#'    for customization. However, if `...` causes an error, the same
#'    `ComplexHeatmap::Heatmap()` function is called without `...`,
#'    which is intended to allow overloading `...` for different
#'    functions.
#' 
#' @returns `Heatmap` object defined in `ComplexHeatmap::Heatmap()` with
#'    custom attributes with the method caption:
#'    * `"caption"`: `character` string with method details.
#'    * `"caption_legendlist"`: `ComplexHeatmap::Legends` object suitable
#'    to be included with Heatmap legends by
#'    `draw(hm, annotation_legend_list=caption_legendlist)`, or
#'    or drawn directly `grid::grid.draw(caption_legendlist)`.
#'    * `"draw_caption"` - `function` that will draw the caption in the
#'    bottom-right corner of the device by default, to be called
#'    with `attr(hm, "draw_caption")()` or `draw_caption()`.
#'
#'    The `Heatmap` row and column order can be retrieved:
#'    1. `jamba::heatmap_row_order()` - returns a `list` of vectors of
#'    rownames in the order they appear in the heatmap, with list names
#'    defined by row split.
#'    2. `jamba::heatmap_column_order()` - returns a `list` of vectors of
#'    colnames in the order they appear in the heatmap, with list names
#'    defined by row split.
#'
#'
#' @export
mem_gene_path_heatmap <- function
(mem,
 genes=NULL,
 sets=NULL,
 min_gene_ct=1,
 min_set_ct=1,
 min_set_ct_each=4,
 column_fontsize=NULL,
 column_cex=1,
 row_fontsize=NULL,
 row_cex=1,
 row_method="binary",
 column_method="binary",
 enrich_im_weight=0.3,
 gene_im_weight=0.5,
 gene_annotations=c("im",
    "direction",
    "default"),
 annotation_suffix=c(im="hit",
    direction="dir"),
 simple_anno_size=grid::unit(6, "mm"),
 cluster_columns=NULL,
 cluster_rows=NULL,
 cluster_row_slices=FALSE,
 cluster_column_slices=FALSE,
 name=NULL,
 p_cutoff=NULL,
 p_floor=1e-10,
 row_split=NULL,
 column_split=NULL,
 auto_split=TRUE,
 column_title=LETTERS,
 row_title=letters,
 row_title_rot=0,
 colorize_by_gene=TRUE,
 na_col="white",
 rotate_heatmap=FALSE,
 colramp="Reds",
 column_names_max_height=grid::unit(180, "mm"),
 column_names_rot=90,
 show_gene_legend=FALSE,
 show_pathway_legend=TRUE,
 show_heatmap_legend=8,
 use_raster=FALSE,
 seed=123,
 verbose=FALSE,
 ...)
{
   #
   if (length(seed) > 0) {
      set.seed(head(seed, 1));
   }
   
   # Recognize Mem or mem (list) input
   Mem <- NULL;
   if (inherits(mem, "Mem")) {
      Mem <- mem;
      mem <- Mem_to_list(Mem);
   } else {
      Mem <- list_to_Mem(mem);
   }
   
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
   
   if (length(colnames(mem$memIM)) == 0) {
      colnames(mem$memIM) <- rownames(mem$enrichIM);
   }
   if (length(column_names_max_height) == 0) {
      column_names_max_height <- grid::unit(180, "mm")
   } else if (!grid::is.unit(column_names_max_height)) {
      if (is.numeric(column_names_max_height)) {
         column_names_max_height <- grid::unit(column_names_max_height, "mm")
      } else {
         column_names_max_height <- grid::unit(180, "mm")
      }
   }
   memIM <- mem$memIM;
   if (length(genes) > 0) {
      memIM <- memIM[intersect(rownames(memIM), genes), , drop=FALSE];
   }
   if (length(sets) > 0) {
      memIM <- memIM[,intersect(colnames(memIM), sets),drop=FALSE];
   }
   if (any(dim(memIM) == 0)) {
      stop("No remaining data after filtering.");
   }

   if (length(name) == 0) {
      name <- "Gene Hits By\nEnrichment";
   }

   # gene_annotations: im, direction
   gene_annotations <- intersect(gene_annotations,
      c("im", "direction", "default"));
   if ("direction" %in% gene_annotations) {
      if ("default" %in% gene_annotations &&
            all(jamba::rmNA(naValue=0, geneIMdirection(Mem)) >= 0)) {
         gene_annotations <- setdiff(gene_annotations, "direction");
         if (!"im" %in% gene_annotations) {
            gene_annotations <- "im";
         }
         if (verbose) {
            jamba::printDebug("mem_gene_path_heatmap(): ",
               c("gene_annotations='direction'",
                  " was removed by 'default' ",
                  "and all non-zero geneIMdirection were positive",
                  " is not available."),
               sep="")
         }
      }
   }

   ## TODO:
   ## apply min_set_ct_each alongside p_cutoff to ensure the
   ## pathways with min_set_ct_each also meet p_cutoff
   met_p_cutoff <- (mem$enrichIM[colnames(memIM), , drop=FALSE] <= p_cutoff)
   met_min_set_ct_each <- do.call(cbind,
      lapply(jamba::nameVector(colnames(mem$geneIM)), function(icol){
         colSums(mem$geneIM[rownames(memIM), icol] * (memIM != 0)) >= min_set_ct_each;
      })
   );
   met_criteria <- rowSums(met_p_cutoff & met_min_set_ct_each) > 0;
   if (any(!met_criteria) && verbose) {
      jamba::printDebug("mem_gene_path_heatmap(): ",
         "Filtered out ",
         sum(!met_criteria),
         " of ",
         length(met_criteria),
         " sets using p_cutoff <= ",
         p_cutoff,
         " and min_set_ct_each >= ",
         min_set_ct_each);
   }
   memIM <- memIM[, met_criteria, drop=FALSE];

   ## apply min_set_ct to each enrichment test
   memIMsetct <- colSums(memIM > 0);
   if (any(memIMsetct < min_set_ct)) {
      if (verbose) {
         jamba::printDebug("mem_gene_path_heatmap(): ",
            "Filtered sets by min_set_ct:",
            min_set_ct);
      }
      sets <- colnames(memIM)[memIMsetct >= min_set_ct];
      memIM <- memIM[, sets, drop=FALSE];
   }
   memIMgenect <- rowSums(memIM > 0);
   if (any(memIMgenect < min_gene_ct)) {
      if (verbose) {
         jamba::printDebug("mem_gene_path_heatmap(): ",
            "Filtered sets by min_gene_ct:",
            min_gene_ct);
      }
      genes <- rownames(memIM)[memIMgenect >= min_gene_ct];
      memIM <- memIM[genes, , drop=FALSE];
   }

   ## Additional step to ensure columns and rows are not empty
   memIM <- memIM[, colSums(memIM > 0) > 0, drop=FALSE];
   memIM <- memIM[rowSums(memIM > 0) > 0, , drop=FALSE];
   if (any(dim(memIM) == 0)) {
      stop("No remaining data after filtering.");
   }
   if (!identical(mem$memIM, memIM)) {
      if (verbose) {
         jamba::printDebug("mem_gene_path_heatmap(): ",
            "Subsetting Mem object.");
      }
      Mem <- Mem[rownames(memIM), colnames(memIM), ];
      mem <- Mem_to_list(Mem);
      memIM <- memIM(Mem);
   }
   genes <- rownames(memIM);
   sets <- colnames(memIM);

   
   ## Handle row/column split, clustering, subsetting
   hrcs <- handle_rowcol_splits(
      Mem=Mem,
      auto_split=auto_split,
      row_split=row_split,
      row_title=row_title,
      gene_im_weight=gene_im_weight,
      cluster_rows=cluster_rows,
      row_method=row_method,
      column_split=column_split,
      column_title=column_title,
      cluster_columns=cluster_columns,
      column_method=column_method,
      enrich_im_weight=enrich_im_weight,
      p_cutoff=p_cutoff,
      p_floor=p_floor,
      ...)
   if (!identical(Mem, hrcs$Mem)) {
      if (verbose) {
         jamba::printDebug("mem_gene_path_heatmap(): ",
            "Subsetting Mem after handle_rowcol_split().");
      }
      Mem <- hrcs$Mem;
      mem <- Mem_to_list(Mem);
      genes <- genes(Mem);
      sets <- sets(Mem);
      memIM <- mem$memIM;
   }
   row_split <- hrcs$row_split;
   row_title <- hrcs$row_title;
   cluster_rows <- hrcs$cluster_rows;
   column_split <- hrcs$column_split;
   column_title <- hrcs$column_title;
   cluster_columns <- hrcs$cluster_columns;
   
   ## Automatic fontsize
   if (length(row_fontsize) == 0) {
      row_fontsize <- jamba::noiseFloor(row_cex * 60/(nrow(memIM))^(1/2),
         minimum=1,
         ceiling=18);
   }
   if (length(column_fontsize) == 0) {
      column_fontsize <- jamba::noiseFloor(column_cex * 60/(ncol(memIM))^(1/2),
         minimum=1,
         ceiling=18);
   }

   ## Apply colors to outside annotations
   ## row annotation colors (genes)
   col_iml1 <- lapply(jamba::nameVectorN(mem$colorV), function(i){
      j <- mem$colorV[[i]];
      circlize::colorRamp2(breaks=c(0,1),
         colors=jamba::fixYellow(Hrange=c(60,100), fixup=TRUE,
            jamba::getColorRamp(j, n=3)[1:2]))
   });
   ## column annotation colors (enrichments)
   col_iml4 <- lapply(jamba::nameVectorN(mem$colorV), function(i){
      j <- mem$colorV[[i]];
      p_vector <- c(-log10(p_cutoff + 1e-10),
         -log10(p_cutoff),
         (-log10(p_cutoff) - log10(p_floor)) / 2,
         -log10(p_floor));
      i_ramp <- circlize::colorRamp2(
         breaks=p_vector,
         colors=c("white",
            jamba::fixYellow(Hrange=c(60,100), fixup=TRUE,
               jamba::getColorRamp(j,
                  trimRamp=c(2,2),
                  n=3,
                  lens=3)))
      )
   });
   
   ## Optionally colorize the matrix by gene
   if (colorize_by_gene) {
      ## Convert rows to blended colors
      if (verbose) {
         jamba::printDebug("mem_gene_path_heatmap(): ",
            "colorize_by_gene");
      }
      geneIMu <- unique(mem$geneIM);
      geneIMu <- jamba::mixedSortDF(
         data.frame(check.names=FALSE,
            rowSums(geneIMu),
            -geneIMu))[,-1,drop=FALSE]*-1;
      sep <- ",\n  ";
      geneIMu_names <- jamba::cPasteS(im2list(t(geneIMu)),
         sep=sep);
      colnames(geneIMu) <- mem$colorV[colnames(mem$geneIM)];
      geneIMu_colors <- jamba::cPasteS(im2list(t(geneIMu)),
         sep=sep);
      names(geneIMu_colors) <- geneIMu_names
      if (verbose) {
         for (i in seq_along(geneIMu_colors)) {
            jamba::printDebug(jamba::nameVector(unlist(strsplit(geneIMu_colors[i], sep)),
               unlist(strsplit(geneIMu_names[i], sep))),
               indent="    ");
         }
      }

      ## blend using colorjam::blend_colors()
      gene_colorsV <- colorjam::blend_colors(strsplit(geneIMu_colors, sep));
      names(gene_colorsV) <- geneIMu_names;

      ## apply these colors to gene rows
      geneVlist <- im2list(t(mem$geneIM));
      geneV <- jamba::cPasteS(geneVlist,
         sep=sep);
      geneColors <- jamba::nameVector(gene_colorsV[geneV],
         names(geneV));

      ## convert memIM to memIM_note
      memIM_note <- apply(memIM, 2, function(i){
         ifelse(i == 0, "", geneV[rownames(memIM)]);
      });
      memIM_levels <- c("", geneIMu_names);
      gene_colorsV_idx <- c("white", gene_colorsV);
      col_hm <- circlize::colorRamp2(
         breaks=seq(from=0, to=length(gene_colorsV)),
         colors=gene_colorsV_idx);
      col_hm_at <- rev(seq_along(gene_colorsV));
      col_hm_labels <- rev(tail(memIM_levels, -1));
      if (is.logical(show_heatmap_legend)) {
         if (show_heatmap_legend > 0) {
            show_heatmap_legend <- 8;
         } else {
            show_heatmap_legend <- FALSE;
         }
      }

      if (verbose) {
         jamba::printDebug("length(col_hm_at): ", length(col_hm_at),
            ", show_heatmap_legend: ", show_heatmap_legend);
      }
      if (length(col_hm_at) > show_heatmap_legend) {
         singlet_idx <- which(lengths(strsplit(rev(geneIMu_colors), sep)) == 1);
         col_hm_at <- col_hm_at[singlet_idx];
         col_hm_labels <- col_hm_labels[singlet_idx];
      }

      memIM_idx <- as.numeric(factor(memIM_note, levels=memIM_levels)) - 1;
      memIM <- matrix(0,
         ncol=ncol(memIM),
         nrow=nrow(memIM),
         dimnames=dimnames(memIM));
      memIM[] <- memIM_idx;
   } else {
      col_hm <- jamba::getColorRamp(colramp,
         lens=5,
         trimRamp=c(0,4));
      col_hm_at <- sort(unique(as.vector(memIM[genes,sets])));
      col_hm_labels <- col_hm_at;
   }

   ## Create the heatmap
   if (length(seed) > 0) {
      set.seed(head(seed, 1));
   }

   # pathway annotation legend param
   path_annotation_legend_param <- lapply(jamba::nameVectorN(col_iml4), function(i){
      list(color_bar="discrete",
         #at=c(0, 1),
         title=paste0(i, "\n-log10P"),
         border="black")
   });
   # gene annotation legend param
   gene_annotation_legend_param <- lapply(col_iml1, function(i){
      list(color_bar="discrete",
         at=c(0, 1),
         border="black")
   });
   # heatmap legend param
   heatmap_legend_param <- list(
      color_bar="discrete",
      border="black",
      at=col_hm_at,
      labels=col_hm_labels
   );

   # generate usable caption describing the relevant parameters used
   gene_type <- "gene";
   set_type <- "set";
   caption <- c(
      paste0(jamba::formatInt(length(genes)),
         " ", gene_type, ifelse(length(genes) != 1, "s")),
      paste0(jamba::formatInt(length(sets)),
         " ", set_type, ifelse(length(genes) != 1, "s")));
   if (TRUE %in% rotate_heatmap) {
      caption <- paste(caption, c("(columns)", "(rows)"));
   } else {
      caption <- paste(caption, c("(rows)", "(columns)"));
   }
   ## Generate caption as Legend
   caption_clustering <- c(
      paste0(gene_type, "s = '", row_method, "'"),
      paste0(set_type, "s = '", column_method, "'"))
   caption_im_weights <- c(
      paste0(gene_type, "IM weight: ", format(gene_im_weight, digits=3)),
      paste0(set_type, " IM weight: ", format(enrich_im_weight, digits=3)))
   if (TRUE %in% rotate_heatmap) {
      caption_clustering <- rev(caption_clustering);
      caption_im_weights <- rev(caption_im_weights);
   }
   caption_list <- list(
      Summary=caption,
      Clustering=caption_clustering,
      Filtering=jamba::rmNA(c(
         paste0("enrichment P <= ", p_cutoff),
         ifelse(min_set_ct_each > 1,
            paste0(gene_type, "s per ",
               set_type, " >= ", min_set_ct_each,
               "\n(within enrichment)"),
            NA),
         ifelse(min_gene_ct > 1,
            paste0(set_type, "s per ",
               gene_type, " >= ", min_gene_ct),
            NA),
         ifelse(min_set_ct > min_set_ct_each,
            paste0(gene_type, "s per ",
               set_type, " >= ", min_set_ct),
            NA))),
      `IM weights`=caption_im_weights)
   # make convenient text summary
   caption <- paste0(collapse="\n",
      paste0(names(caption_list), ":\n"),
      jamba::cPaste(caption_list, sep="\n"));
   # Convert to list of Legend objects
   caption_legends <- lapply(names(caption_list), function(lname){
      caption_grob <- grid::textGrob(
         label=jamba::cPaste(sep="\n",
            caption_list[[lname]]),
         gp=grid::gpar(fontsize=10,
            lineheight=1),
         hjust=0, vjust=0)
      text_legend2 <- ComplexHeatmap::Legend(
         title=paste0(lname),
         # title_gp=grid::gpar(font=2),
         grob=caption_grob,
         title_position="topleft",
         type="grid")
   })
   caption_legendlist <- ComplexHeatmap::packLegend(
      list=caption_legends,
      direction="vertical")
   # jamba::nullPlot(plotAreaTitle="");ComplexHeatmap::draw(caption_legendlist);
   ## Unclear if we will use draw_caption() in future,
   ## since we prefer caption_legendlist as ComplexHeatmap legend
   draw_caption <- function
   (use_caption=caption_legendlist,
    x=grid::unit(1, "npc"),
    y=grid::unit(0, "npc"),
    just=NULL,
    ...)
   {
      if (length(x) == 0) {
         x <- 0.5;
      }
      if (length(y) == 0) {
         y <- 0.5;
      }
      if (!grid::is.unit(x)) {
         x <- grid::unit(x, "npc");
      }
      if (!grid::is.unit(y)) {
         y <- grid::unit(y, "npc");
      }
      xv <- (as.numeric(x));
      yv <- (as.numeric(y));
      caption_legendlist <- use_caption;
      caption_legendlist@grob$vp$justification <- c(xv, yv)
      caption_legendlist@grob$vp$valid.just <- c(xv, yv)
      if (length(just) == 0) {
         just <- c("right", "bottom");
         if (xv == 0) {
            just[1] <- "left";
         } else if (xv > 0 && xv < 1) {
            just[1] <- "center";
         }
         if (yv == 1) {
            just[2] <- "top";
         } else if (yv < 1 && yv > 0) {
            just[2] <- "center";
         }
      }
      just <- rep(just, length.out=2);
      brc <- grid::viewport(x=x,
         y=y,
         width=0.03, height=0.03,
         default.units="npc",
         just=just)
      grid::pushViewport(brc);
      grid::grid.draw(caption_legendlist@grob);
      grid::popViewport();
   }

   # raster_device workaround
   if (jamba::check_pkg_installed("ragg")) {
      raster_device <- "agg_png"
   } else {
      raster_device <- "png"
   }

   # default orientation, gene rows, pathway columns
   if (!rotate_heatmap) {
      top_annotation <- ComplexHeatmap::HeatmapAnnotation(
         which="column",
         border=TRUE,
         show_legend=show_pathway_legend,
         annotation_legend_param=path_annotation_legend_param,
         col=col_iml4,
         df=data.frame(check.names=FALSE,
            -log10(mem$enrichIM[sets,,drop=FALSE])),
         gap=grid::unit(0, "mm")
      );
      # scalable approach to annotations
      gene_anno_list <- list();
      gene_color_list <- list();
      gene_param_list <- list();
      gene_anno_gap <- numeric(0);
      for (gene_ann in gene_annotations) {
         if ("direction" %in% gene_ann) {
            if (length(gene_anno_list) > 0) {
               gene_anno_gap <- c(2, gene_anno_gap)
            }
            dir_matrix <- mem$geneIMdirection[genes, , drop=FALSE];
            if ("direction" %in% names(annotation_suffix) &&
                  any(nchar(annotation_suffix[["direction"]])) > 0) {
               colnames(dir_matrix) <- paste(colnames(dir_matrix),
                  rep(annotation_suffix[["direction"]],
                     length.out=ncol(dir_matrix)));
            }
            gene_anno_list <- c(
               list(direction=dir_matrix),
               gene_anno_list)
            gene_color_list <- c(list(
               direction=colorjam::col_div_xf(1.2)),
               gene_color_list);
            gene_param_list <- c(
               list(direction=list(color_bar="discrete",
                  at=c(-1, 0, 1),
                  labels=c("down", "no change", "up"),
                  border="black")),
               gene_param_list)
         }
         if ("im" %in% gene_ann) {
            if (length(gene_anno_list) > 0) {
               gene_anno_gap <- c(2, gene_anno_gap);
            }
            gene_anno_gap <- c(
               rep(0, ncol(mem$geneIM) - 1),
               gene_anno_gap)
            im_matrix <- mem$geneIM[genes, , drop=FALSE];
            if ("im" %in% names(annotation_suffix) &&
                  any(nchar(annotation_suffix[["im"]])) > 0) {
               colnames(im_matrix) <- paste(colnames(im_matrix),
                  rep(annotation_suffix[["im"]],
                     length.out=ncol(im_matrix)));
               names(col_iml1) <- paste(names(col_iml1),
                  rep(annotation_suffix[["im"]],
                     length.out=length(col_iml1)));
               names(gene_annotation_legend_param) <- paste(
                  names(gene_annotation_legend_param),
                  rep(annotation_suffix[["im"]],
                     length.out=length(gene_annotation_legend_param)));
            }
            gene_anno_list <- c(
               as.list(data.frame(check.names=FALSE, im_matrix)),
               gene_anno_list)
            gene_color_list <- c(
               col_iml1,
               gene_color_list);
            gene_param_list <- c(
               gene_annotation_legend_param,
               gene_param_list)
         }
      }
      # put it together
      if (length(gene_anno_list) == 0) {
         left_annotation <- NULL;
      } else {
         if (length(gene_anno_gap) == 0) {
            gene_anno_gap <- 0;
         }
         gene_alist <- alist(
            simple_anno_size=simple_anno_size,
            # gap=ComplexHeatmap::ht_opt("ROW_ANNO_PADDING"),
            col=gene_color_list,
            gap=grid::unit(gene_anno_gap, "mm"),
            annotation_legend_param=gene_param_list,
            show_legend=show_gene_legend,
            # show_annotation_name=show_gene_annotation_name,
            annotation_name_rot=column_names_rot,
            annotation_name_gp=grid::gpar(fontsize=column_fontsize),
            border=TRUE
         );
         gene_arglist <- c(
            gene_alist,
            gene_anno_list);
         left_annotation <- do.call(ComplexHeatmap::rowAnnotation,
            gene_arglist);
      }
      hm <- tryCatch({
         jamba::call_fn_ellipsis(ComplexHeatmap::Heatmap,
            matrix=memIM[genes, sets, drop=FALSE],
            border=TRUE,
            name=name,
            na_col=na_col,
            cluster_columns=cluster_columns,
            cluster_column_slices=cluster_column_slices,
            cluster_rows=cluster_rows,
            cluster_row_slices=cluster_row_slices,
            clustering_distance_columns=column_method,
            clustering_distance_rows=row_method,
            top_annotation=top_annotation,
            col=col_hm,
            show_heatmap_legend=show_heatmap_legend,
            left_annotation=left_annotation,
            row_names_gp=grid::gpar(fontsize=row_fontsize),
            column_names_gp=grid::gpar(fontsize=column_fontsize),
            column_names_rot=column_names_rot,
            column_names_max_height=column_names_max_height,
            column_title=column_title,
            heatmap_legend_param=heatmap_legend_param,
            row_split=row_split,
            row_title=row_title,
            row_title_rot=row_title_rot,
            column_split=column_split,
            use_raster=use_raster,
            raster_device=raster_device,
            ...);
      }, error=function(e){
         # same as above but without ...
         if (TRUE || verbose) {
            jamba::printDebug("mem_gene_path_heatmap(): ",
               "Error in Heatmap(), calling without '...', error is shown below:");
            print(e);
         }
         ComplexHeatmap::Heatmap(
            matrix=memIM[genes, sets, drop=FALSE],
            border=TRUE,
            name=name,
            na_col=na_col,
            cluster_columns=cluster_columns,
            cluster_column_slices=cluster_column_slices,
            cluster_rows=cluster_rows,
            cluster_row_slices=cluster_row_slices,
            clustering_distance_columns=column_method,
            clustering_distance_rows=row_method,
            top_annotation=top_annotation,
            col=col_hm,
            show_heatmap_legend=show_heatmap_legend,
            left_annotation=left_annotation,
            row_names_gp=grid::gpar(fontsize=row_fontsize),
            column_names_gp=grid::gpar(fontsize=column_fontsize),
            column_names_rot=column_names_rot,
            column_names_max_height=column_names_max_height,
            column_title=column_title,
            heatmap_legend_param=heatmap_legend_param,
            row_split=row_split,
            row_title=row_title,
            row_title_rot=row_title_rot,
            column_split=column_split,
            use_raster=use_raster,
            raster_device=raster_device)
      })
      # add caption in attributes
      attr(hm, "caption") <- caption;
      attr(hm, "caption_legendlist") <- caption_legendlist;
      attr(hm, "draw_caption") <- draw_caption;
   } else {
      ########################################################################
      # rotate heatmap 90 degrees so pathway names are rows, genes are columns
      left_annotation <- ComplexHeatmap::rowAnnotation(
         border=TRUE,
         annotation_legend_param=path_annotation_legend_param,
         annotation_name_rot=column_names_rot,
         col=col_iml4,
         df=-log10(mem$enrichIM[sets,,drop=FALSE]),
         gap=grid::unit(0, "mm")
      );
      top_annotation <- ComplexHeatmap::HeatmapAnnotation(
         col=col_iml1,
         border=TRUE,
         annotation_legend_param=gene_annotation_legend_param,
         #gp=grid::gpar(col="#00000011"), # per-cell border
         df=mem$geneIM[genes,,drop=FALSE],
         gap=grid::unit(0, "mm")
      );
      hm <- tryCatch({
         jamba::call_fn_ellipsis(ComplexHeatmap::Heatmap,
            matrix=t(memIM[genes,sets,drop=FALSE]),
            border=TRUE,
            name=name,
            na_col=na_col,
            cluster_rows=cluster_columns,  # flipped
            cluster_column_slices=cluster_row_slices, # flipped
            cluster_columns=cluster_rows,  # flipped
            cluster_row_slices=cluster_column_slices, # flipped
            clustering_distance_columns=row_method, # flipped
            clustering_distance_rows=column_method, # flipped
            top_annotation=top_annotation,
            col=col_hm,
            show_heatmap_legend=show_heatmap_legend,
            left_annotation=left_annotation,
            row_names_gp=grid::gpar(fontsize=column_fontsize), # flipped
            column_names_gp=grid::gpar(fontsize=row_fontsize), # flipped
            column_names_rot=column_names_rot, # NOT flipped
            row_names_max_width=column_names_max_height, # flipped
            row_title=column_title, # flipped
            column_title=row_title, # flipped
            heatmap_legend_param=heatmap_legend_param,
            row_title_rot=row_title_rot, # NOT flipped
            column_split=row_split, # flipped
            row_split=column_split, # flipped
            use_raster=use_raster,
            raster_device=raster_device,
            ...);
      }, error=function(e){
         # same as above but without ...
         if (verbose) {
            jamba::printDebug("mem_gene_path_heatmap(): ",
               "Error in Heatmap(), calling without '...', error is shown below:");
            print(e);
         }
         ComplexHeatmap::Heatmap(
            matrix=t(memIM[genes,sets,drop=FALSE]),
            border=TRUE,
            name=name,
            na_col=na_col,
            cluster_rows=cluster_columns,
            cluster_column_slices=cluster_row_slices,
            cluster_columns=cluster_rows,
            cluster_row_slices=cluster_column_slices,
            clustering_distance_columns=row_method,
            clustering_distance_rows=column_method,
            top_annotation=top_annotation,
            col=col_hm,
            show_heatmap_legend=show_heatmap_legend,
            left_annotation=left_annotation,
            row_names_gp=grid::gpar(fontsize=column_fontsize),
            column_names_gp=grid::gpar(fontsize=row_fontsize),
            column_names_rot=column_names_rot,
            row_names_max_width=column_names_max_height, # flipped
            row_title=column_title,
            column_title=row_title,
            heatmap_legend_param=heatmap_legend_param,
            row_title_rot=row_title_rot,
            column_split=row_split,
            row_split=column_split,
            use_raster=use_raster,
            raster_device=raster_device)
      })
      # add caption in attributes
      attr(hm, "caption") <- caption;
      attr(hm, "caption_legendlist") <- caption_legendlist;
      attr(hm, "draw_caption") <- draw_caption;
   }

   return(hm);
}

#' MultiEnrichment Heatmap of enrichment P-values
#'
#' MultiEnrichment Heatmap of enrichment P-values
#'
#' Note: It is recommended to call `mem_plot_folio()` with `do_which=1`
#' in order to utilize the gene-pathway content during clustering,
#' which is more effective at clustering similar pathways by gene
#' content. Otherwise pathways are clustered using only the
#' `-log10(p)` enrichment P-value.
#'
#' This function is a lightweight wrapper to `ComplexHeatmap::Heatmap()`
#' intended to visualize the enrichment P-values from multiple
#' enrichment results. The P-value threshold is used to colorize
#' every cell whose P-value meets the threshold, while all other
#' cells are therefore white.
#'
#' ## Drawing Dotplot with Point Legend
#' The `style` argument controls whether a heatmap or dotplot is
#' created.
#' * `style="dotplot"`: each heatmap cell is not filled, and the color
#' is drawn as a circle with size proportional to the number of
#' genes involved in enrichment. A separate point legend is returned
#' as an attribute of the heatmap object.
#' * `style="dotplot_inverted"`: each heatmap cell is filled, and
#' a circle is drawn with size proportional to the number of
#' genes involved in enrichment. A separate point legend is returned
#' as an attribute of the heatmap object.
#'
#' To draw the dotplot heatmap including the point legend,
#' use this command:
#'
#' ```R
#' ComplexHeatmap::draw(hm,
#'    annotation_legend_list=attr(hm, "annotation_legend_list"))
#' ```
#'
#' Generally, the clustering
#' using the gene-pathway incidence matrix is more effective at
#' representing biologically-driven pathway clusters.
#'
#' @family jam plot functions
#'
#' @param mem `Mem` or legacy `list` mem created by `multiEnrichMap()`.
#' @param style `character` string indicating the style of heatmap:
#'    `"heatmap"` produces a regular heatmap, shaded by `log10(Pvalue)`;
#'    `"dotplot"` produces a dotplot, where the dot size is proportional
#'    to the number of genes. See function description for details on
#'    how to include the point size legend beside the heatmap.
#'    The main benefit of using "dotplot" style is that it also indicates
#'    the relative number of genes involved in the enrichment.
#' @param apply_direction `logical`, default FALSE, whether to define
#'    a bivariate color scheme which uses `mem$enrichIMdirection`
#'    when defined. The color scheme is intended to indicate both the
#'    directional strength (usually with some type of z-score) and
#'    the statistical enrichment (usually with the enrichment P-value
#'    or adjusted P-value).
#' @param p_cutoff `numeric` value of the enrichment P-value cutoff,
#'    by default this value is obtained from `mem$p_cutoff` to be
#'    consistent with the original `multiEnrichMap()` analysis.
#'    P-values above `p_cutoff` are not colored, and are therefore white.
#'    This behavior is intended to indicate pathways with P-value above
#'    this threshold did not meet the threshold, instead of
#'    pathways with similar P-values displaying with similar color.
#' @param min_count `numeric` number of genes required for a pathway
#'    to be considered dysregulated.
#' @param p_floor `numeric` minimum P-value used for the color gradient.
#'    P-values below this floor are colored with the maximum color gradient.
#'    This value is intended to be used in cases where one enrichment
#'    P-value is very low (e.g. 1e-36) to prevent all other P-values from
#'    being colored pale red-white and not be noticeable.
#' @param point_size_factor `numeric` used to adjust the legend point size,
#'    since the heatmap point size is dependent upon the number of rows,
#'    the legend may require some manual adjustment to make sure the
#'    legend matches the heatmap.
#' @param point_size_min,point_size_max `numeric` values which define
#'    the minimum and maximum point sizes, effectively defining the range
#'    permitted when `style="dotplot"`.
#' @param row_method `character` string of the distance method
#'    to use for row and column clustering. The clustering is performed
#'    by `amap::hcluster()`.
#' @param name `character` value passed to `ComplexHeatmap::Heatmap()`,
#'    used as a label above the heatmap color legend.
#' @param row_dend_reorder `logical` indicating whether to reorder the
#'    row dendrogram using the method described in
#'    `ComplexHeatmap::Heatmap()`. The end result is minor reshuffling of
#'    leaf nodes on the dendrogram based upon mean signal in each row,
#'    which can sometimes look cleaner.
#' @param row_fontsize,column_fontsize optional `numeric` arguments passed to
#'    `ComplexHeatmap::Heatmap()` to size row and column labels.
#' @param cluster_columns `logical` indicating whether to cluster heatmap
#'    columns, by default columns are not clustered.
#' @param sets `character` vector of sets (pathways) to include in the heatmap,
#'    all other sets will be excluded.
#' @param color_by_column `logical` indicating whether to colorize the
#'    heatmap using `mem$colorV` as defined for each comparison. This
#'    option is currently experimental, and produces a base R heatmap
#'    using `jamba::imageByColors()`.
#' @param cex.axis `numeric` adjustment for axis labels, passed to
#'    `jamba::imageByColors()` only when `color_by_column=TRUE`.
#' @param lens `numeric` value used in color gradients, defining the extent
#'    the color gradient is enhanced in the mid-ranges (positive `lens`),
#'    or diminished in the mid-ranges (negative `lens`). There is no
#'    quantitative standard measure for color gradient changes, so this
#'    option is intended to help adjust and improve the visual perception
#'    of the current data.
#' @param cexCellnote `numeric` character expansion value used only
#'    when `color_by_column=TRUE`, used to adjust the P-value label size
#'    inside each heatmap cell.
#' @param column_title optional `character` string with title to display
#'    above the heatmap.
#' @param row_names_max_width,column_names_max_height,heatmap_legend_param
#'    arguments passed to `ComplexHeatmap::Heatmap()` and provided here
#'    for convenience.
#' @param hm_cell_size `grid::unit` or `numeric`, default `NULL`, to define
#'    an optional fixed heatmap cell size, useful to define absolute
#'    square heatmap cells. When `numeric` it is interpreted as
#'    "mm" units. Note that the heatmap total height is determined by
#'    the number of cells, and the total row gaps defined by the
#'    number of row gaps with `row_split` multiplied by `row_gap`.
#' @param legend_height `grid::unit`, default 6 cm (60 mm), to define
#'    the absolute height of the color gradient in the color key.
#'    This value is only used when `heatmap_legend_param` is not defined.
#' @param legend_cex `numeric` default 1, used to scale the legend
#'    fontsize relative to the default fontsize `10`.
#'    This value is only used when `heatmap_legend_param` is not defined.
#' @param top_annotation `HeatmapAnnotation` as produced by
#'    `ComplexHeatmap::HeatmapAnnotation()` or `NULL`, used to display
#'    customized annotation at the top of the heatmap. The order of
#'    columns must match the order of columns in the data displayed
#'    in the heatmap.
#' @param outline `logical` default TRUE, whether to draw an outline
#'    for each heatmap cell. Note: The outline is not drawn for
#'    `style="dotplot"` which already adds lines through the middle
#'    of each cell, not the border of each cell.
#' @param show_enrich `numeric` default NULL, indicating which of the
#'    enrichment metrics to show as a label in each cell.
#'    When only one type is shown, there is no prefix, but for multiple
#'    types, a prefix is shown for each. The metrics in order include:
#'    1. "-log10P"
#'    2. "z-score"
#'    3. "number of genes"
#' @param use_raster `logical` passed to `ComplexHeatmap::Heatmap()`
#'    indicating whether to rasterize the heatmap output, used when
#'    `style="heatmap"`. Rasterization is not relevant to dotplot output
#'    since dotplot is drawn using an individual circle in each heatmap cell.
#' @param do_plot `logical` indicating whether to display the plot with
#'    `ComplexHeatmap::draw()` or `jamba::imageByColors()` as relevant.
#'    The underlying data is returned invisibly.
#' @param ... additional arguments are passed to `ComplexHeatmap::Heatmap()`
#'    for customization.
#'
#' @export
mem_enrichment_heatmap <- function
(mem,
 style=c("dotplot_inverted",
    "dotplot",
    "heatmap"),
 apply_direction=FALSE,
 p_cutoff=NULL,
 min_count=1,
 p_floor=1e-10,
 point_size_factor=1,
 point_size_max=8,
 point_size_min=2,
 row_method="euclidean",
 column_method="euclidean",
 name="-log10P",
 row_dend_reorder=TRUE,
 row_dend_width=grid::unit(30, "mm"),
 row_fontsize=NULL,
 row_cex=1,
 row_split=NULL,
 row_gap=grid::unit(2, "mm"),
 cluster_rows=TRUE,
 column_fontsize=NULL,
 column_cex=1,
 cluster_columns=FALSE,
 sets=NULL,
 color_by_column=FALSE,
 cex.axis=1,
 lens=3,
 cexCellnote=1,
 column_title=NULL,
 row_names_max_width=grid::unit(300, "mm"),
 column_names_max_height=grid::unit(300, "mm"),
 heatmap_legend_param=NULL,
 hm_cell_size=NULL,
 legend_height=grid::unit(6, "cm"),
 legend_cex=1,
 direction_cutoff=0,
 gene_count_max=NULL,
 top_annotation=NULL,
 outline=TRUE,
 show_enrich=NULL,
 use_raster=FALSE,
 do_plot=TRUE,
 ...)
{
   #
   style <- match.arg(style);
   Mem <- NULL;
   if (inherits(mem, "Mem")) {
      Mem <- mem;
      mem <- Mem_to_list(Mem);
   } else {
      Mem <- list_to_Mem(mem);
   }
   
   if (length(p_cutoff) == 0) {
      if ("p_cutoff" %in% names(thresholds(Mem))) {
         p_cutoff <- thresholds(Mem)$p_cutoff;
      } else if ("cutoffRowMinP" %in% names(thresholds(Mem))) {
         p_cutoff <- thresholds(Mem)[["cutoffRowMinP"]];
      } else {
         p_cutoff <- 1;
      }
   }
   
   col_logp <- circlize::colorRamp2(
      breaks=c(-log10(p_cutoff + 1e-10),
         seq(from=-log10(p_cutoff),
            to=-log10(p_floor),
            length.out=25)),
      colors=c("white",
         jamba::getColorRamp("Reds",
            lens=lens,
            n=25,
            trimRamp=c(2, 2)))
   )
   if (apply_direction) {
      col_logp <- colorjam::col_div_xf(-log10(p_floor),
         floor=-log10(p_cutoff),
         colramp="RdBu_r",
         trimRamp=c(1, 1),
         lens=lens);
   }

   if (length(sets) > 0) {
      # version 0.0.76.900 change order to retain sets ordering by default
      # sets <- intersect(rownames(enrichIM(Mem)), sets);
      sets <- intersect(sets, sets(Mem));
      Mem <- Mem[, sets, ]
      mem <- Mem_to_list(Mem)
   } else {
      sets <- sets(Mem);
   }
   if (any(dim(Mem) == 0)) {
      stop("No remaining data after filtering.");
   }
   if (ncol(enrichIM(Mem)) > 1) {
      if (is.logical(cluster_rows) && TRUE %in% cluster_rows) {
         cluster_rows <- amap::hcluster(
            link="ward",
            jamba::noiseFloor(
               -log10(enrichIM(Mem)[sets,,drop=FALSE]),
               minimum=-log10(p_cutoff+1e-5),
               newValue=0,
               ceiling=-log10(p_floor)),
            method=row_method);
         cluster_rows <- as.dendrogram(cluster_rows);
         if (length(row_dend_width) == 0) {
            row_dend_width <- grid::unit(30, "mm");
         }
      }
      if (is.logical(cluster_columns) && TRUE %in% cluster_columns) {
         cluster_columns <- amap::hcluster(
            link="ward",
            jamba::noiseFloor(
               t(-log10(enrichIM(Mem)[sets,,drop=FALSE])),
               minimum=-log10(p_cutoff+1e-5),
               newValue=0,
               ceiling=-log10(p_floor)),
            #ceiling=3),
            method=column_method);
      }
   } else {
      if (is.logical(cluster_rows)) {
         cluster_rows <- FALSE;
      }
      cluster_columns <- FALSE;
      if (length(row_dend_width) == 0) {
         row_dend_width <- grid::unit(10, "mm");
      }
   }

   ## Automatic fontsize
   if (length(column_fontsize) == 0) {
      row_fontsize <- jamba::noiseFloor(
         row_cex * 60/(nrow(enrichIM(Mem)))^(1/2),
         minimum=1,
         ceiling=18);
   }
   if (length(column_fontsize) == 0) {
      column_fontsize <- jamba::noiseFloor(
         column_cex * 60/(ncol(enrichIM(Mem)))^(1/2),
         minimum=1,
         ceiling=20);
   }

   if (length(heatmap_legend_param) == 0) {
      heatmap_legend_param <- list(
         border="black",
         labels_gp=grid::gpar(fontsize=10 * legend_cex),
         title_gp=grid::gpar(fontsize=10 * legend_cex),
         legend_height=legend_height);
   }

   # optionally apply direction
   has_negative <- any(
      jamba::rmNA(naValue=0, enrichIMdirection(Mem)) < 0);
   if (TRUE %in% apply_direction && has_negative) {
      use_matrix <- -log10(enrichIM(Mem));
      # use_direction contains z-score values at or above direction_cutoff
      # otherwise it is set to zero
      use_direction <- (
         (abs(enrichIMdirection(Mem)) >= direction_cutoff) *
         enrichIMdirection(Mem));
   } else {
      use_matrix <- -log10(enrichIM(Mem));
      use_direction <- NULL;
      apply_direction <- FALSE;
   }

   # raster_device workaround
   # disabled with version 0.0.78.900
   # if (jamba::check_pkg_installed("ragg")) {
   #    raster_device <- "agg_png"
   # } else {
   raster_device <- "png"
   # }

   ## Experimental: set heatmap size with fixed cell dimensions
   hm_width <- NULL;
   hm_height <- NULL;
   if (length(hm_cell_size) > 0) {
      if (length(hm_cell_size) == 1) {
         hm_cell_size <- rep(hm_cell_size, 2);
      } else {
         hm_cell_size <- head(hm_cell_size, 2);
      }
      if (!grid::is.unit(hm_cell_size)) {
         hm_width <- grid::unit(hm_cell_size[1] * ncol(use_matrix), "mm");
         hm_height <- grid::unit(hm_cell_size[2] * nrow(use_matrix), "mm");
      } else {
         hm_width <- ncol(use_matrix) * hm_cell_size[1];
         hm_height <- nrow(use_matrix) * grid::unit(hm_cell_size[2], "mm");
      }
   }

   if ("heatmap" %in% style) {
      pch <- NULL;
   } else {
      pch <- 21;
   }
   if ("heatmap1" %in% style) {
      hm <- jamba::call_fn_ellipsis(ComplexHeatmap::Heatmap,
         matrix=use_matrix,
         name=name,
         col=col_logp,
         cluster_rows=cluster_rows,
         row_dend_reorder=row_dend_reorder,
         border=TRUE,
         row_names_gp=grid::gpar(fontsize=row_fontsize),
         row_names_max_width=row_names_max_width,
         column_names_gp=grid::gpar(fontsize=column_fontsize),
         column_names_max_height=column_names_max_height,
         cluster_columns=cluster_columns,
         row_dend_width=row_dend_width,
         column_title=column_title,
         heatmap_legend_param=heatmap_legend_param,
         use_raster=use_raster,
         raster_device=raster_device,
         ...);
   } else {
      if (length(gene_count_max) == 0) {
         ctmax <- ceiling(max(Mem@enrichIMgeneCount, na.rm=TRUE));
      } else {
         ctmax <- gene_count_max;
      }
      #jamba::printDebug("ctmax: ", ctmax);

      if (ctmax <= 1) {
         ct_ticks <- c(0, 1);
      } else {
         n <- 8;
         ct_ticks <- setdiff(unique(c(
            #1,
            min_count,
            round(pretty(c(0, ctmax), n=n)))), 0);
         ct_step <- median(diff(ct_ticks));
         if (max(ct_ticks) > ctmax) {
            ct_ticks[which.max(ct_ticks)] <- ctmax;
            if (tail(diff(ct_ticks), 1) < ceiling(ct_step / 4)) {
               ct_ticks <- head(ct_ticks, -2);
            } else if (tail(diff(ct_ticks), 1) < ceiling(ct_step / 2)) {
               ct_ticks <- c(head(ct_ticks, -2), ctmax);
            }
         }
      }
      ct_approxfun <- function(x, ...){
         approxfun(
            x=sqrt(c(min_count, ctmax)),
            yleft=0,
            ties="ordered",
            yright=point_size_max,
            y=c(point_size_min,
               point_size_max * point_size_factor))(sqrt(x), ...);
      }
      ct_tick_sizes <- ct_approxfun(ct_ticks);

      # define point size legend
      #ctbreaks <- ct_to_breaks(ctmax, n=10, maxsize=point_size_max)
      #ctbreaksize <- ct_to_size(ctbreaks, ctmax=ctmax, n=10, maxsize=point_size_max) * point_size_factor;
      pt_legend_ncol <- 1;
      if (length(ct_ticks) >= 8) {
         pt_legend_ncol <- 2;
      }
      if (any(grepl("dotplot", style))) {
         pt_legend <- ComplexHeatmap::Legend(
            labels=ct_ticks,
            title="Gene Count",
            type="points",
            pch=pch,
            ncol=pt_legend_ncol,
            labels_gp=grid::gpar(fontsize=10 * legend_cex),
            title_gp=grid::gpar(fontsize=10 * legend_cex),
            size=grid::unit(ct_tick_sizes, "mm"),
            grid_height=grid::unit(max(ct_tick_sizes) * 0.95, "mm"),
            grid_width=grid::unit(max(ct_tick_sizes) * 0.95, "mm"),
            background="transparent",
            legend_gp=grid::gpar(col="black",
               fill="grey85"));
         anno_legends <- list(pt_legend);
      } else {
         anno_legends <- list();
      }

      # custom cell label, hide 2 when directional data are not available
      if (2 %in% show_enrich && length(use_direction) == 0) {
         show_enrich <- setdiff(show_enrich, 2);
      }
      use_prefix <- NULL;
      if (length(show_enrich) > 1) {
         use_prefix <- c(
            "-log10P: ",
            "z-score: ",
            "genes: ")[show_enrich]
      }
      # improved cell_fun
      if (apply_direction) {
         tcount <- jamba::tcount;
         dir_colors <- c("royalblue4", "gold3", "firebrick3");
         dir_colors2 <- c("skyblue", "gold", "indianred1");
         dir_colors3 <- c("white", "white", "white");
         mcolor <- jamba::rbindList(list(dir_colors3, dir_colors2, dir_colors))
         # jamba::imageByColors(mcolor)
         # white_num controls the intensity of the first non-white color
         # in the color gradient
         white_num <- 2;
         mcolor2 <- matrix(ncol=3,
            c("white", "white", "white",
               colorjam::blend_colors(c(dir_colors[1], rep("white", white_num))),
               colorjam::blend_colors(c(dir_colors[2], rep("white", white_num))),
               colorjam::blend_colors(c(dir_colors[3], rep("white", white_num))),
               dir_colors),
            byrow=TRUE);
         row_breaks <- c(-log10(p_cutoff) - 1e-10,
            seq(from=-log10(p_cutoff),
               to=-log10(p_floor),
               length.out=2));
         if (p_cutoff == 1) {
            row_breaks <- tail(row_breaks, -1);
            mcolor <- mcolor[-2, , drop=FALSE]
         }
         col_bivariate <- colorRamp2D(
            column_breaks=seq(from=-2, to=2, length.out=3),
            row_breaks=row_breaks,
            mcolor);
         size_by <- match("geneCount",
            c("-log10Pvalue",
               "z-score",
               "geneCount"));
         legend_bivariate <- make_legend_bivariate(col_bivariate,
            ylab="-log10pvalue",
            xlab="z-score");
         use_col_fn <- col_bivariate;
         # if ("dotplot_inverted" %in% style) {
         #    use_col_fn <- function(x, y){
         #       rep("#FFFFFF", length.out=length(x))
         #    };
         # }
         cell_fun_custom <- cell_fun_bivariate(
            list(
               use_direction,
               use_matrix,
               Mem@enrichIMgeneCount),
            invert=grepl("invert", style),
            pch=pch,
            size_fun=ct_approxfun,
            size_by=size_by,
            outline_style="darker",
            col_hm=use_col_fn,
            show=show_enrich,
            outline=outline,
            cex=cexCellnote,
            prefix=use_prefix,
            ...
         );
         anno_legends <- c(anno_legends,
            list(legend_bivariate));
         show_heatmap_legend <- FALSE;
      } else {
         use_col_fn <- col_logp;
         # if ("dotplot_inverted" %in% style) {
         #    use_col_fn <- function(x){
         #       rep("#FFFFFF", length.out=length(x))
         #    };
         # }
         show_heatmap_legend <- TRUE;
         # remove show_enrich=2 if no supporting directional data is present
         cell_fun_custom <- cell_fun_bivariate(
            list(
               use_matrix,
               use_direction,
               Mem@enrichIMgeneCount),
            invert=grepl("invert", style),
            pch=pch,
            size_fun=ct_approxfun,
            size_by=3,
            outline_style="darker",
            col_hm=use_col_fn,
            show=show_enrich,
            outline=outline,
            cex=cexCellnote,
            type="univariate",
            prefix=use_prefix,
            ...
         );
      }
      # validate row_split
      if (length(row_split) > 0 && length(row_split) >= nrow(use_matrix)) {
         if (length(names(row_split)) > 0 &&
               all(rownames(use_matrix) %in% names(row_split))) {
            row_split <- row_split[rownames(use_matrix)];
         }
      }
      if (is.numeric(row_split) && row_split == 1) {
         row_split <- NULL
      }
      if (is.logical(cluster_rows) && FALSE %in% cluster_rows) {
         row_split <- NULL
      }

      if ("dotplot" %in% style) {
         use_raster <- FALSE
      }
      ## Add row_split to hm_height
      if (length(row_split) > 0 &&
            length(row_gap) > 0 &&
            any(as.numeric(row_gap) > 0) &&
            length(hm_height) == 1) {
         if (is.numeric(row_split) &&
               length(row_split) == 1 &&
               row_split > 1) {
            hm_height <- hm_height + (row_split - 1) * row_gap;
         }
      }

      # dot plot or heatmap style
      hm <- jamba::call_fn_ellipsis(ComplexHeatmap::Heatmap,
         matrix=use_matrix,
         name=name,
         col=col_logp,
         width=hm_width,
         height=hm_height,
         cluster_rows=cluster_rows,
         row_dend_reorder=row_dend_reorder,
         border=TRUE,
         row_names_gp=grid::gpar(fontsize=row_fontsize),
         row_names_max_width=row_names_max_width,
         row_split=row_split,
         row_gap=row_gap,
         column_names_gp=grid::gpar(fontsize=column_fontsize),
         column_names_max_height=column_names_max_height,
         cluster_columns=cluster_columns,
         row_dend_width=row_dend_width,
         column_title=column_title,
         heatmap_legend_param=heatmap_legend_param,
         rect_gp=grid::gpar(type="none"),
         cell_fun=cell_fun_custom,
         show_heatmap_legend=show_heatmap_legend,
         top_annotation=top_annotation,
         use_raster=use_raster,
         raster_device=raster_device,
         ...);
      attr(hm,
         "annotation_legend_list") <- anno_legends;
      if (do_plot) {
         ComplexHeatmap::draw(hm,
            merge_legends=TRUE,
            annotation_legend_list=anno_legends);
         # message to use draw command
         # draw(hm, annotation_legend_list=anno_legends);
      }
   }

   if ("heatmap" %in% style && color_by_column) {
      hm_sets <- rownames(enrichIM(Mem))[ComplexHeatmap::row_order(hm)];
      ## Prepare fresh image colors using p_cutoff and p_floor
      enrichIMcolors <- do.call(cbind,
         lapply(jamba::nameVector(colnames(enrichIM(Mem))), function(i){
            x <- -log10(enrichIM(Mem)[,i]);
            cr1 <- circlize::colorRamp2(
               breaks=c(-log10(p_cutoff + 1e-10),
                  seq(from=-log10(p_cutoff),
                     to=-log10(p_floor),
                     length.out=24)),
               colors=c("white",
                  getColorRamp(mem$colorV[i],
                     n=24,
                     trimRamp=c(1, 0),
                     lens=lens)));
            cr1(x);
         }));
      #enrichIMcolors <- colorjam::matrix2heatColors(
      #   x=-log10(enrichIM(Mem)),
      #   colorV=mem$colorV,
      #   baseline=-log10(p_cutoff),
      #   numLimit=-log10(p_floor),
      #   lens=lens);
      if (do_plot) {
         jamba::imageByColors(enrichIMcolors[hm_sets, , drop=FALSE],
            cellnote=sapply(enrichIM(Mem)[hm_sets, , drop=FALSE],
               base::format.pval,
               eps=1e-50,
               digits=2),
            adjustMargins=TRUE,
            flip="y",
            cexCellnote=cexCellnote,
            cex.axis=cex.axis,
            main=column_title,
            groupCellnotes=FALSE,
            ...);
      }
      retlist <- list(
         matrix=enrichIMcolors[hm_sets, , drop=FALSE],
         cellnote=sapply(enrichIM(Mem)[hm_sets, , drop=FALSE],
            base::format.pval,
            eps=1e-50,
            digits=2),
         adjustMargins=TRUE,
         flip="y",
         cexCellnote=cexCellnote,
         cex.axis=cex.axis,
         main=column_title,
         groupCellnotes=FALSE);
      return(invisible(retlist));
   } else {
      return(invisible(hm));
   }
}

#' MultiEnrichMap plot
#'
#' MultiEnrichMap plot, deprecated
#'
#' This function is replaced by `mem2emap()` and is deprecated.
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
#' @param mem legacy `list` mem from `multiEnrichMap()`, specifically
#'    containing elements `"multiEnrichMap","multiEnrichMap2"` which
#'    are expected to be `igraph` objects.
#' @param overlap numeric value between 0 and 1, indicating the Jaccard
#'    overlap filter to use for edge nodes. The value is used to
#'    delete edges whose values `E(g)$overlap` are below this threshold.
#' @param overlap_count numeric value indicating the minimum overlap count
#'    below which edges are removed. The value `E(g)$overlap_count` is
#'    used for this filter.
#' @param do_plot logical indicating whether to plot the final result.
#' @param do_legend logical indicating whether to print a color legend,
#'    valid when `do_plot=TRUE`. Arguments `...` are also passed to
#'    `mem_legend()`, for example `x="bottomleft"` can be overriden.
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
 do_legend=TRUE,
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
 edge_bundling="connections",
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
      keep_nodes <- which(igraph::V(g)$name %in% sets);
      if (length(keep_nodes) == 0) {
         stop("No sets remain after filtering V(g)$name for sets.");
      }
      g <- igraph::subgraph(g, v=keep_nodes);
   }

   ## Filter to remove nodes as needed
   if (length(overlap) > 0 && overlap > 0) {
      delete_edges <- which(igraph::E(g)$overlap < overlap);
      if (length(delete_edges) > 0) {
         g <- igraph::delete_edges(g, delete_edges);
      }
   }
   if (length(overlap_count) > 0 && "overlap_count" %in% igraph::edge_attr_names(g)) {
      delete_edges <- which(igraph::E(g)$overlap_count < overlap_count);
      if (length(delete_edges) > 0) {
         g <- igraph::delete_edges(g, delete_edges);
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
   label_edges <- head(intersect(label_edges, igraph::edge_attr_names(g)), 1);
   if (length(label_edges) > 0) {
      if (!"label" %in% label_edges) {
         edge_text <- igraph::edge_attr(g, label_edges);
         if (is.numeric(edge_text)) {
            edge_text <- format(edge_text,
               trim=TRUE,
               digits=2);
         }
         igraph::E(g)$label <- edge_text;
      }
   }
   if (length(edge_cex) > 0) {
      igraph::E(g)$label.cex <- edge_cex;
   }
   if (length(edge_color) > 0) {
      igraph::E(g)$color <- edge_color;
   }
   if (length(shape) > 0) {
      igraph::V(g)$shape <- shape;
   }
   if (length(node_size) > 0) {
      igraph::V(g)$size <- node_size;
   }
   if (length(node_cex) > 0) {
      igraph::V(g)$label.cex <- node_cex;
   }
   if (length(frame_color) > 0) {
      igraph::V(g)$frame.color <- frame_color;
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
            edge_bundling=edge_bundling,
            ...);
      }
      if (do_legend) {
         mem_legend(mem,
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
#' @param x,y,bg,box.col,title,cex,ncol,pch,pt.cex,pt.lwd,inset arguments
#'    are passed to `graphics::legend()`.
#'    Note `pt.lwd` is mostly relevant with `do_direction=TRUE`, which
#'    adds open circles to the legend, whose line width has default
#'    `pt.lwd=2`.
#' @param do_directional `logical` indicating whether to include
#'    directional colors defined in `directional_colors`, indicated only
#'    as the border color of nodes.
#' @param directional_column `character` indicating how to add the
#'    directional colors to columns of legend, with two options:
#'    * "same": appends `directional_colors` to legend colors using
#'    the defined `ncol` number of columns.
#'    * "added-bottom": appends `directional_colors` as a new column
#'    so the resulting legend with have `ncol+1` columns.
#'    In this case, intervening empty rows are filled with blank space,
#'    and the `directional_colors` are shown in the bottom-most rows in the
#'    far right column of the legend.
#'    * "added-top": appends `directional_colors` as a new column
#'    so the resulting legend with have `ncol+1` columns.
#'    In this case, intervening empty rows are filled with blank space,
#'    and the `directional_colors` are shown in the top-most rows in the
#'    far right column of the legend.
#' @param directional_colors `character` vector of R colors, named by
#'    the label to be shown in the legend, displayed in order (top to bottom,
#'    left to right) they appear in this vector.
#'    To remove the entry `"no change"`, supply a new vector:
#'    `directional_colors=c(up.regulated="firebrick3", down.regulated="blue")`
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
 pt.lwd=2,
 inset=0,
 do_directional=FALSE,
 directional_column=c("same",
    "added-bottom",
    "added-top"),
 directional_colors=c(
    `up-regulated`="firebrick3",
    `no change`="grey80",
    `down-regulated`="blue"),
 ...)
{
   ##
   Mem <- NULL;
   if (inherits(mem, "Mem")) {
      Mem <- mem;
      mem <- Mem_to_list(Mem);
   }
   directional_column <- match.arg(directional_column);
   if (!is.list(mem) || !"colorV" %in% names(mem)) {
      stop("Input mem must be a list with element 'colorV'");
   }
   colorV <- mem[["colorV"]];
   colorVb <- jamba::makeColorDarker(colorV, darkFactor=1.5);

   # directional circles for some plots
   if (do_directional) {
      # exemplar legend()
      if ("same" %in% directional_column) {
         # simply include direction at the end with same ncol
         legend_n <- length(legend) + length(directional_colors);
         legend_nrow <- ceiling(legend_n / ncol);
         pt.bg <- c(colorV,
            rep(NA, length(directional_colors)))
         col <- c(rep(NA, length(colorV)),
            directional_colors);
         pch <- c(rep(pch, length.out=length(legend)),
            rep(21, length(directional_colors)));
         legend_names <- c(names(colorV),
            names(directional_colors));
      } else {
         # add direction in its own column
         legend_nrow <- max(c(
            length(directional_colors),
            ceiling(length(legend) / ncol)));
         legend_n_diff <- (ncol * legend_nrow) - length(legend);
         legend_n_buff <- 0;
         if ("added-bottom" %in% directional_column) {
            legend_n_buff <- legend_nrow - length(legend);
         }
         pt.bg <- c(colorV,
            rep(NA, legend_n_diff + legend_n_buff),
            rep(NA, length(directional_colors)));
         col <- c(rep(NA, length(colorV)),
            rep(NA, legend_n_diff + legend_n_buff),
            directional_colors);
         pch <- c(rep(pch, length.out=length(legend)),
            rep(NA, legend_n_diff + legend_n_buff),
            rep(21, length(directional_colors)));
         legend_names <- c(names(colorV),
            rep("", legend_n_diff + legend_n_buff),
            names(directional_colors));
      }
   } else {
      legend_names <- names(colorV);
      pt.bg <- colorV;
      col <- colorVb;
   }

   tryCatch({
      legend(x=x,
         y=y,
         title=title,
         ncol=ncol,
         legend=legend_names,
         pch=pch,
         pt.cex=pt.cex,
         pt.lwd=pt.lwd,
         pt.bg=pt.bg,
         col=col,
         bg=bg,
         box.col=box.col,
         cex=cex,
         inset=inset,
         ...);
   }, error=function(e){
      legend(x=x,
         y=y,
         title=title,
         ncol=ncol,
         legend=legend_names,
         pch=pch,
         pt.cex=pt.cex,
         pt.lwd=pt.lwd,
         pt.bg=colorV,
         col=col,
         bg=bg,
         box.col=box.col,
         cex=cex,
         inset=inset);
   })
}


#' Draw Heatmap with title and subtitle using grid viewports
#'
#' Draw Heatmap with title and subtitle using grid viewports
#'
#' This function is likely to work only with `grid` or `ComplexHeatmap`
#' input. Other objects that have `draw()` defined may not work
#' properly.
#'
#' Note: This function may be superceded by the ability to include
#' an overall title when using `ComplexHeatmap::draw()`, which was not
#' originally possible.
#'
#' This function is intended to help make it easier to wrap a
#' heatmap from `ComplexHeatmap::Heatmap()` inside the relevant
#' grid viewports in order to display one large title at the
#' top of the resulting visualization, and optionally one
#' subtitle at the bottom of the visualization.
#'
#' The input can be one heatmap (`Heatmap` object),
#' or object of class `"gTree"`, or any object with a method
#' `draw()` associated with it, detected by
#' `methods::hasMethod("draw", class(object))`.
#'
#' A good example of `"gTree"` object is the result of
#' calling `draw()` inside `grid::grid.grabExpr()`, for example:
#'
#' `grid::grid.grabExpr(ComplexHeatmap::draw(...))`
#'
#' The `ComplexHeatmap::draw()` function has extended capability
#' for arranging one or more heatmaps and associated annotations.
#'
#' @family jam plot functions
#'
#' @return The byproduct of this function is to draw a grid visualization
#'    that includes a title and subtitle, then the `object` in the center.
#'
#' @param object an object with class `"Heatmap"`, `"gTree"`, or
#'    any object where `methods::hasMethod("draw", class(object))`
#'    is `TRUE`.
#' @param title character string used as title. When `NULL` or
#'    `nchar(title)==0` then no title is displayed.
#' @param title_fontsize numeric value indicating the font size,
#'    with units `"points"`.
#' @param title_just character string indicating the justification
#'    (alignment) of the title.
#' @param caption,caption_fontsize,caption_just arguments equivalent
#'    to the title_* arguments.
#' @param caption_x numeric value in grid units specifying where to
#'    position the caption text. By default when `caption_just` is `"left"`
#'    the `caption_x` is defined by `grid::unit(0.2, "npc")`, which
#'    positions the caption at the left side (20% from the left edge)
#'    of the plot, with text proceeding to the right of that point.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#'
#' @export
grid_with_title <- function
(object,
 title=NULL,
 title_fontsize=18,
 title_just="centre",
 caption=NULL,
 caption_fontsize=12,
 caption_just="left",
 caption_x=NULL,
 verbose=FALSE,
 ...)
{
   ## This function requires grid, need to verify function prefixing
   ##
   if (!any(c("gTree", "Heatmap") %in% class(object))) {
      if (!methods::hasMethod("draw", class(object))) {
         stop(paste0("Input object must have class 'Heatmap' or 'gTree', ",
            "or have method draw() detected by hasMethod('draw', class(object))"));
      }
   }
   title_nlines <- NULL;
   title_units <- NULL;
   caption_row <- 2;
   hm_row <- 1;
   if (length(title) > 0 && nchar(title) > 0) {
      title_nlines <- length(unlist(
         strsplit(title, "\n"))) * (title_fontsize / 10);
      title_units <- "lines";
      hm_row <- 2;
      caption_row <- 3;
   }
   caption_nlines <- NULL;
   caption_units <- NULL;
   if (length(caption) > 0 && nchar(caption) > 0) {
      caption_nlines <- length(unlist(
         strsplit(caption, "\n"))) * (caption_fontsize / 10);
      caption_units <- "lines";
   }
   panel_units <- c(title_units,
      "null",
      caption_units);
   panel_heights <- c(title_nlines,
      1,
      caption_nlines);
   if (verbose) {
      jamba::printDebug("grid_with_title(): ",
         "panel_heights:",
         panel_heights);
      jamba::printDebug("grid_with_title(): ",
         "panel_units:",
         panel_units);
   }

   ## Create grid layout object
   l <- grid::grid.layout(nrow=length(panel_heights),
      ncol=1,
      heights=grid::unit(panel_heights,
         panel_units));

   ## Create the combined plot
   vp <- grid::viewport(width=1,
      height=1,
      layout=l);

   grid::grid.newpage();
   grid::pushViewport(vp);

   if (verbose) {
      jamba::printDebug("grid_with_title(): ",
         "pushing viewport row:",
         hm_row);
   }
   grid::pushViewport(
      grid::viewport(
         layout.pos.row=hm_row));
   if ("gTree" %in% class(object)) {
      grid::grid.draw(object);
   } else {
      if (verbose) {
         jamba::printDebug("grid_with_title(): ",
            "creating grid object grTree.");
      }
      grobject <- grid::grid.grabExpr(
         ComplexHeatmap::draw(object));
      if (verbose) {
         jamba::printDebug("grid_with_title(): ",
            "calling grid::grid.draw()");
      }
      grid::grid.draw(
         grobject);
   }
   grid::popViewport();

   if (length(title_nlines) > 0) {
      if (verbose) {
         jamba::printDebug("grid_with_title(): ",
            "drawing title in row:",
            1);
      }
      grid::grid.text(title,
         just=title_just,
         gp=grid::gpar(fontsize=title_fontsize),
         vp=grid::viewport(layout.pos.row=1,
            layout.pos.col=1));
   }
   if (length(caption_nlines) > 0) {
      if (length(caption_x) == 0) {
         if ("left" %in% caption_just) {
            caption_x <- grid::unit(0.2, "npc");
         } else {
            caption_x <- grid::unit(0.5, "npc");
         }
      }
      if (verbose) {
         jamba::printDebug("grid_with_title(): ",
            "drawing caption in row:",
            caption_row);
         jamba::printDebug("grid_with_title(): ",
            "caption_x:",
            caption_x);
      }
      grid::grid.text(caption,
         just=caption_just,
         x=caption_x,
         gp=grid::gpar(fontsize=caption_fontsize),
         vp=grid::viewport(
            layout.pos.row=caption_row,
            layout.pos.col=1));
   }
   grid::popViewport();
}

