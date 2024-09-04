
#' MultiEnrichment Heatmap of Genes and Pathways
#'
#' MultiEnrichment Heatmap of Genes and Pathways
#'
#' This function takes the `mem` list output from
#' `multiEnrichMap()` and creates a gene-by-pathway incidence
#' matrix heatmap, using `ComplexHeatmap::Heatmap()`.
#' It uses three basic sources of data to annotate the heatmap:
#'
#' 1. `mem$memIM` the gene-set incidence matrix
#' 2. `mem$geneIM` the gene incidence matrix by dataset
#' 3. `mem$enrichIM` the pathway enrichment P-value matrix by dataset
#'
#' It will try to estimate a reasonable number of column and row
#' splits in the dendrogram, based solely upon the number of
#' columns and rows. These guesses can be controlled with argument
#' `column_split` and `row_split`, respectively.
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
#' @family jam plot functions
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
#' @param min_set_ct_each minimum number of genes required for each set,
#'    required for at least one enrichment test.
#' @param column_fontsize,row_fontsize `numeric`
#'    passed as `fontsize` to `ComplexHeatmap::Heatmap()`
#'    to define a specific fontsize for column and row labels. When
#'    `NULL` the nrow/ncol of the heatmap are used to infer a reasonable
#'    starting point fontsize, which can be adjusted with `column_cex`
#'    and `row_cex`.
#' @param row_method,column_method character string of the distance method
#'    to use for row and column clustering. The clustering is performed
#'    by `amap::hcluster()`.
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
#'    By default it uses `"im", "direction"`, and `"direction"` is removed
#'    when `mem$geneIMdirection` is not available.
#'    * `"im"` displays the gene incidence matrix `mem$geneIM` using
#'    categorical colors defined by `mem$colorV`.
#'    * `"direction"` displays the gene directionality `mem$geneIMdirection`
#'    using colors defined by `colorjam::col_div_xf(1.2)`.
#'    * When no values are given, the gene annotation is not displayed.
#'    * When two values are given, the annotations are displayed in the
#'    order they are provided.
#' @param annotation_suffix `character` vector named by values permitted
#'    by `gene_annotations`, with optional suffix to add to the annotation
#'    labels. For example it may be helpful to add "hit" or "dir" to
#'    distinguish the enrichment labels.
#' @param name character value passed to `ComplexHeatmap::Heatmap()`,
#'    used as a label above the heatmap color legend.
#' @param p_cutoff numeric value of the enrichment P-value cutoff,
#'    above which P-values are not colored, and are therefore white.
#'    The enrichment P-values are displayed as an annotated heatmap
#'    at the top of the main heatmap. Any cell that has a color meets
#'    at least the minimum P-value threshold. This value by default
#'    is taken from input `mem`, using `mem$p_cutoff`, for
#'    consistency with the input multienrichment analysis.
#' @param column_split,row_split optional arguments passed to
#'    `ComplexHeatmap::Heatmap()` to split the heatmap by columns
#'    or rows, respectively.
#'    * when `row_split` is `NULL`
#'    and `auto_split=TRUE`, it will determine an appropriate number
#'    of clusters based upon the number of rows. To turn off row `split`,
#'    use `row_split=NULL` or `row_split=0` or `row_split=1`;
#'    likewise for `column_split`.
#'    * when `row_split` or `column_split` are supplied as a named
#'    vector, the names are aligned with `sets` to be displayed
#'    in the heatmap, and will use the `intersect()` of the two.
#'    When data is clustered, `cluster_row_slices=FALSE` and
#'    `cluster_column_slices=FALSE` such that the dendrogram will
#'    be broken into separate pieces.
#' @param column_title optional character string with title to display
#'    above the heatmap.
#' @param row_title optional character string with title to display
#'    beside the heatmap. Note when `row_split` is defined, the
#'    `row_title` is applied to each heatmap section.
#' @param row_title_rot `numeric` value indicating the rotation of
#'    `row_title` text, where `0` is not rotated, and `90` is rotated
#'    90 degrees.
#' @param colorize_by_gene `logical` indicating whether to color the
#'    main heatmap body using the colors from `geneIM` which represents
#'    each enrichment in which a given gene is involved. Colors are
#'    blended using `colorjam::blend_colors()`, using colors from
#'    `mem$colorV`, applied to `mem$geneIM`.
#' @param na_col `character` string indicating the color to use for
#'    NA or missing values. Typically this argument is only used
#'    when `colorize_by_gene=TRUE`, where entries with no color are
#'    recognized as `NA` by `ComplexHeatmap::Heatmap()`.
#' @param rotate_heatmap `logical` indicating whether the entire heatmap
#'    should be rotated so that pathway names are displayed as rows,
#'    and genes as columns. Notes on how arguments are applied to rows
#'    and columns:
#'    * Column arguments applied to rows:
#'    `column_split`, `column_title`, `cluster_columns`,
#'    `column_fontsize`, `column_cex`
#'    are applied to rows since they refer to pathway data;
#'    * Row arguments applied to columns:
#'    `row_split`, `row_title`, `cluster_rows`, `row_fontsize`, `row_cex`
#'    are applied to columns since they refer to gene data;
#'    * Arguments applied directly to columns:
#'    `column_method`, `column_title_rot`
#'    are applied directly to heatmap columns since they
#'    refer to the output heatmap options.
#'    * Arguments applied directly to rows:
#'    `row_method`, `row_title_rot`
#'    are applied directly to heatmap rows since they
#'    refer to the output heatmap options.
#' @param seed `numeric` value passed to `set.seed()` to allow
#'    reproducible results, typically with clustering operations.
#' @param colramp `character` name of color, color gradient, or a
#'    vector of colors, anything compatible with input to
#'    `jamba::getColorRamp()`.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to `ComplexHeatmap::Heatmap()`
#'    for customization. However, if `...` causes an error, the same
#'    `ComplexHeatmap::Heatmap()` function is called without `...`,
#'    which is intended to allow overloading `...` for different
#'    functions.
#'
#' @returns `Heatmap` object defined in `ComplexHeatmap::Heatmap()`, with
#'    additional attributes:
#'    * `"caption"` - a `character` string with caption that described
#'    the data dimensions and clustering parameters.
#'    * `"caption_legendlist"` - a `ComplexHeatmap::Legends` object suitable
#'    to be included with Heatmap legends using
#'    `draw(hm, annotation_legend_list=caption_legendlist)`, or
#'    drawn with `grid::grid.draw(caption_legendlist)`.
#'    * `"draw_caption"` - a `function` that will draw the caption in the
#'    bottom-left corner of the heatmap by default, to be called
#'    with `attr(hm, "draw_caption")()` or `draw_caption()`.
#'
#'    In addition, the returned object can be interrogated with two helper
#'    functions that help define the row and column clusters, and the
#'    exact order of labels as they appear in the heatmap.
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
    "direction"),
 annotation_suffix=c(im="hit",
    direction="dir"),
 simple_anno_size=grid::unit(6, "mm"),
 cluster_columns=NULL,
 cluster_rows=NULL,
 cluster_row_slices=TRUE,
 cluster_column_slices=TRUE,
 name=NULL,
 p_cutoff=mem$p_cutoff,
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
 column_names_max_height=grid::unit(18, "cm"),
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

   if (length(name) == 0) {
      name <- "enrichments\nper gene";
   }

   # gene_annotations: im, direction
   gene_annotations <- intersect(gene_annotations,
      c("im", "direction"));
   if ("direction" %in% gene_annotations && !"geneIMdirection" %in% names(mem)) {
      gene_annotations <- setdiff(gene_annotations, "direction");
      if (verbose) {
         jamba::printDebug("mem_gene_path_heatmap(): ",
            c("gene_annotations='direction'",
               " was removed because ",
               "mem$geneIMdirection",
               " is not available."),
            sep="")
      }
   }

   ## TODO:
   ## apply min_set_ct_each alongside p_cutoff to ensure the
   ## pathways with min_set_ct_each also meet p_cutoff
   met_p_cutoff <- (mem$enrichIM[colnames(memIM),,drop=FALSE] <= p_cutoff)
   met_min_set_ct_each <- do.call(cbind, lapply(jamba::nameVector(colnames(mem$geneIM)), function(icol){
      colSums(mem$geneIM[rownames(memIM),icol] * (memIM != 0)) >= min_set_ct_each;
   }));
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
   memIM <- memIM[,met_criteria,drop=FALSE];

   ## apply min_set_ct to each enrichment test
   memIMsetct <- colSums(memIM > 0);
   if (any(memIMsetct < min_set_ct)) {
      if (verbose) {
         jamba::printDebug("mem_gene_path_heatmap(): ",
            "Filtered sets by min_set_ct:",
            min_set_ct);
      }
      sets <- colnames(memIM)[memIMsetct >= min_set_ct];
      memIM <- memIM[,sets,drop=FALSE];
   }
   memIMgenect <- rowSums(memIM > 0);
   if (any(memIMgenect < min_gene_ct)) {
      if (verbose) {
         jamba::printDebug("mem_gene_path_heatmap(): ",
            "Filtered sets by min_gene_ct:",
            min_gene_ct);
      }
      genes <- rownames(memIM)[memIMgenect >= min_gene_ct];
      memIM <- memIM[genes,,drop=FALSE];
   }

   ## Additional step to ensure columns and rows are not empty
   memIM <- memIM[,colSums(memIM > 0) > 0,drop=FALSE];
   memIM <- memIM[rowSums(memIM > 0) > 0,,drop=FALSE];
   if (any(dim(memIM) == 0)) {
      stop("No remaining data after filtering.");
   }
   genes <- rownames(memIM);
   sets <- colnames(memIM);

   ## Validate row_split, row_title, column_split, column_title
   if (auto_split) {
      if (length(row_split) == 0) {
         if (nrow(memIM) < 5) {
            row_split <- NULL;
         } else {
            row_split <- jamba::noiseFloor(floor(nrow(memIM)^(1/2.5)),
               ceiling=12);
            if (verbose) {
               jamba::printDebug("mem_gene_path_heatmap(): ",
                  "auto_split row_split:", row_split);
            }
         }
      }
      if (length(column_split) == 0) {
         if (ncol(memIM) < 5) {
            column_split <- NULL;
         } else {
            ncol_x <- (5 + ncol(memIM)/1.5) ^ (1/2);
            column_split <- jamba::noiseFloor(floor(ncol_x), ceiling=10);
            if (verbose) {
               jamba::printDebug("mem_gene_path_heatmap(): ",
                  "auto_split column_split:", column_split);
            }
         }
      }
   }
   if (length(row_split) == 1 && is.numeric(row_split)) {
      if (row_split > length(genes)) {
         row_split <- min(c(
            ceiling(length(genes)^(1/3)),
            10));
      }
      if (row_split <= 1) {
         row_split <- NULL;
         if (length(row_title) > 1) {
            row_title <- NULL;
         }
      } else {
         if (length(row_title) > 1) {
            row_title <- jamba::makeNames(
               rep(row_title,
                  length.out=row_split),
               ...);
         }
      }
   }
   if (length(row_split) > 0 && !is.numeric(row_split)) {
      cluster_row_slices=FALSE;
      # align row_split with genes
      if (length(names(row_split)) > 0) {
         if (!any(genes %in% names(row_split))) {
            stop("names(row_split) do not match any genes");
         }
         genes <- intersect(genes, names(row_split));
         row_split <- row_split[genes];
         memIM <- memIM[genes, , drop=FALSE];
      }
      if (length(row_title) > 0) {
         if (grepl("frame|matrix|tbl", ignore.case=TRUE, class(row_split))) {
            row_split_v <- jamba::pasteByRow(row_split);
            row_title <- jamba::makeNames(
               rep(row_title,
                  length.out=length(unique(row_split_v))),
               ...);
         } else {
            row_title <- jamba::makeNames(
               rep(row_title,
                  length.out=length(unique(row_split))),
               ...);
         }
      }
   }

   if (length(column_split) == 1 && is.numeric(column_split)) {
      if (column_split > length(sets)) {
         column_split <- length(sets);
      }
      if (column_split <= 1) {
         column_split <- NULL;
         if (length(column_title) > 1) {
            column_title <- NULL;
         }
      } else {
         if (length(column_title) == 0) {
            column_title <- LETTERS;
         }
         if (length(column_title) > 1) {
            column_title <- jamba::makeNames(
               rep(column_title,
                  length.out=column_split),
               ...);
         }
      }
   }
   if (length(column_split) > 0 && !is.numeric(column_split)) {
      cluster_column_slices=FALSE;
      # align column_split with sets
      if (length(names(column_split)) > 0) {
         if (!any(sets %in% names(column_split))) {
            stop("names(column_split) do not match any sets");
         }
         sets <- intersect(sets, names(column_split));
         column_split <- column_split[sets];
         memIM <- memIM[, sets, drop=FALSE];
      }
      if (length(column_title) > 0) {
         if (grepl("frame|matrix|tbl", ignore.case=TRUE, class(column_split))) {
            column_split_v <- jamba::pasteByRow(column_split);
            column_title <- jamba::makeNames(
               rep(column_title,
                  length.out=length(unique(column_split_v))),
               ...);
         } else {
            column_title <- jamba::makeNames(
               rep(column_title,
                  length.out=length(unique(column_split))),
               ...);
         }
      }
   }

   ## Automatic fontsize
   if (length(row_fontsize) == 0) {
      row_fontsize <- row_cex * 60/(nrow(memIM))^(1/2);
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
   ## Cluster columns
   # if cluster_columns=NULL or cluster_columns=TRUE
   if (length(cluster_columns) == 0 ||
         (length(cluster_columns) == 1 &&
               is.logical(cluster_columns) &&
               cluster_columns %in% TRUE)) {
      # Assemble the P-value matrix with gene incidence matrix
      # and cluster altogether, which has the benefit/goal of
      # accentuating similar enrichment profiles which also have
      # similar gene content.
      if (!length(enrich_im_weight) == 1 || any(enrich_im_weight > 1)) {
         enrich_im_weight <- 0.3;
      }
      enrich_weight <- round(enrich_im_weight * 10);
      im_weight <- 10 - enrich_weight;
      min_weight <- max(c(1, min(c(enrich_weight, im_weight))));
      enrich_weight <- enrich_weight / min_weight;
      im_weight <- im_weight / min_weight;
      ## 0.0.31.900 use column_matrix with enrich_im_weight adjustment
      column_matrix <- cbind(
         jamba::noiseFloor(
            -log10(mem$enrichIM[sets,,drop=FALSE]),
            minimum=-log10(p_cutoff+1e-5),
            newValue=0,
            ceiling=-log10(p_floor)) * enrich_weight,
            # ceiling=10) * enrich_weight,
         t((mem$memIM[genes, sets, drop=FALSE]) * im_weight) # non-incidence
         # t((mem$memIM[genes, sets, drop=FALSE] != 0) * 1) * im_weight # incidence
      );
      if (is.numeric(column_split)) {
         set.seed(seed);
         cluster_columns <- amap::hcluster(
            link="ward",
            column_matrix,
            method=column_method);
      } else {
         cluster_column_slices=FALSE;
         cluster_columns <- function(x, ...){
            set.seed(seed);
            userows <- rownames(x);
            use_x <- jamba::rmNA(naValue=0,
               column_matrix[match(userows, rownames(column_matrix)), , drop=FALSE]);
            amap::hcluster(use_x,
               link="ward",
               method=column_method)
         }
      }
   }
   ## Cluster rows
   if (length(cluster_rows) == 0 ||
         (length(cluster_rows) == 1 && is.logical(cluster_rows) && cluster_rows)) {
      if (!length(gene_im_weight) == 1 || any(gene_im_weight > 1)) {
         gene_im_weight <- 0.5;
      }
      gene_weight <- round(gene_im_weight * 100)/10;
      im_weight <- 10 - gene_weight;
      min_weight <- max(c(1, min(c(gene_weight, im_weight))));
      gene_weight <- gene_weight / min_weight;
      im_weight <- im_weight / min_weight;
      use_gene_im <- mem$geneIM;
      if ("geneIMdirection" %in% names(mem)) {
         # TODO: verify what to do when geneIM has 1 and geneIMdirection is NA,
         # which implies direction is "not known". For now we use 1 to maintain
         # the geneIM value.
         use_gene_im <- use_gene_im * jamba::rmNA(mem$geneIMdirection, naValue=1);
      }
      if (im_weight == 0) {
         row_matrix <- (use_gene_im[genes, , drop=FALSE]) * gene_weight;
      } else if (gene_weight == 0) {
         row_matrix <- ((mem$memIM[genes,sets,drop=FALSE]) * im_weight); # not incidence
         # row_matrix <- ((mem$memIM[genes,sets,drop=FALSE] != 0) * 1) * im_weight; # incidence
      } else {
         row_matrix <- cbind(
            (use_gene_im[genes, , drop=FALSE]) * gene_weight,
            ((mem$memIM[genes, sets, drop=FALSE])) * im_weight)
         # note: do not convert memIM to incidence
         # below: convert to incidence matrix format
         # ((mem$memIM[genes,sets,drop=FALSE] != 0) * 1) * im_weight)
      }
      if (is.numeric(row_split)) {
         set.seed(seed);
         cluster_rows <- amap::hcluster(
            link="ward",
            row_matrix,
            method=row_method);
      } else {
         cluster_row_slices=FALSE;
         cluster_rows <- function(x, ...){
            set.seed(seed);
            userows <- rownames(x);
            usematrix <- jamba::rmNA(naValue=0,
               row_matrix[x, , drop=FALSE]);
            amap::hcluster(
               usematrix,
               link="ward",
               method=row_method);
         }
      }
   }

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
   set.seed(seed);

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
   caption <- c(
      paste0(jamba::formatInt(length(genes)), " gene rows"),
      paste0(jamba::formatInt(length(sets)), " set columns"));
   if (TRUE %in% rotate_heatmap) {
      caption <- c(
         paste0(jamba::formatInt(length(sets)), " set rows"),
         paste0(jamba::formatInt(length(genes)), " gene columns"));
   }
   ## Generate caption as Legend
   caption_list <- list(
      Summary=caption,
      Clustering=c(
         paste0("column = '", column_method, "'"),
         paste0("row = '", row_method, "'")),
      Filtering=jamba::rmNA(c(
         paste0("enrichment P <= ", p_cutoff),
         ifelse(min_gene_ct > 1,
            paste0("genes per set >= ", min_gene_ct),
            NA),
         ifelse(min_set_ct > 1,
            paste0("sets per gene >= ", min_set_ct),
            NA))),
      `IM weights`=c(
         paste0("enrich_im_weight = ", format(enrich_im_weight, digits=3)),
         paste0("gene_im_weight = ", format(gene_im_weight, digits=3))))
   # make convenient text summary
   caption <- cat(paste0(collapse="\n",
      paste0(names(caption_list), ":\n"),
      jamba::cPaste(caption_list, sep="\n")))
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
         df=-log10(mem$enrichIM[sets,,drop=FALSE]),
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
            row_title=row_title,
            heatmap_legend_param=heatmap_legend_param,
            row_title_rot=row_title_rot,
            row_split=row_split,
            column_split=column_split,
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
         ComplexHeatmap::Heatmap(memIM[genes, sets, drop=FALSE],
            border=TRUE,
            name=name,
            na_col=na_col,
            cluster_columns=cluster_columns,
            cluster_rows=cluster_rows,
            clustering_distance_columns=column_method,
            clustering_distance_rows=row_method,
            top_annotation=top_annotation,
            col=col_hm,
            left_annotation=left_annotation,
            row_names_gp=grid::gpar(fontsize=row_fontsize),
            column_names_gp=grid::gpar(fontsize=column_fontsize),
            column_names_rot=column_names_rot,
            column_title=column_title,
            row_title=row_title,
            heatmap_legend_param=heatmap_legend_param,
            use_raster=use_raster,
            raster_device=raster_device,
            row_title_rot=row_title_rot,
            row_split=row_split,
            column_split=column_split);
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
            t(memIM[genes,sets,drop=FALSE]),
            border=TRUE,
            name=name,
            na_col=na_col,
            cluster_rows=cluster_columns,
            cluster_columns=cluster_rows,
            clustering_distance_columns=column_method,
            clustering_distance_rows=row_method,
            top_annotation=top_annotation,
            col=col_hm,
            left_annotation=left_annotation,
            row_names_gp=grid::gpar(fontsize=column_fontsize),
            column_names_gp=grid::gpar(fontsize=row_fontsize),
            column_names_rot=column_names_rot,
            row_title=column_title,
            column_title=row_title,
            heatmap_legend_param=heatmap_legend_param,
            row_title_rot=row_title_rot,
            column_split=row_split,
            row_split=column_split,
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
            t(memIM[genes,sets,drop=FALSE]),
            border=TRUE,
            name=name,
            na_col=na_col,
            cluster_rows=cluster_columns,
            cluster_columns=cluster_rows,
            clustering_distance_columns=column_method,
            clustering_distance_rows=row_method,
            top_annotation=top_annotation,
            col=col_hm,
            left_annotation=left_annotation,
            row_names_gp=grid::gpar(fontsize=column_fontsize),
            column_names_gp=grid::gpar(fontsize=row_fontsize),
            column_names_rot=column_names_rot,
            row_title=column_title,
            column_title=row_title,
            heatmap_legend_param=heatmap_legend_param,
            use_raster=use_raster,
            raster_device=raster_device,
            row_title_rot=row_title_rot,
            column_split=row_split,
            row_split=column_split);
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
#' This function is a lightweight wrapper to `ComplexHeatmap::Heatmap()`
#' intended to visualize the enrichment P-values from multiple
#' enrichment results. The P-value threshold is used to colorize
#' every cell whose P-value meets the threshold, while all other
#' cells are therefore white.
#'
#' The `style` argument controls whether a heatmap or dotplot is
#' created. When `style="dotplot"` the cells are colored as usual
#' but are drawn as circles sized proportional to the number of
#' genes involved in enrichment. Because `ComplexHeatmap::Heatmap()`
#' is used for this step, a separate point legend is returned
#' as an attribute of the heatmap object.
#'
#' To draw the dotplot heatmap including the point legend,
#' use this form:
#'
#' ```R
#' ComplexHeatmap::draw(hm,
#'    annotation_legend_list=attr(hm, "annotation_legend_list"))
#' ```
#'
#' Note that this function may be more easily applied through the
#' wrapper function with the format `mem_plot_folio(mem, do_which=1, ...)`.
#' The wrapper function `mem_plot_folio()` performs hierarchical
#' clustering of the underlying gene-pathway incidence matrix, which
#' informs the clustering of enrichment results shown
#' in this function `mem_enrichment_heatmap()`. Otherwise, this
#' function will cluster pathways using only the enrichment
#' P-values transformed with `-log10(p)`. Generally, the clustering
#' using the gene-pathway incidence matrix is more effective at
#' representing biologically-driven pathway clusters.
#'
#' @family jam plot functions
#'
#' @param mem `list` object created by `multiEnrichMap()`. Specifically
#'    the object is expected to contain `enrichIM`.
#' @param style `character` string indicating the style of heatmap:
#'    `"heatmap"` produces a regular heatmap, shaded by `log10(Pvalue)`;
#'    `"dotplot"` produces a dotplot, where the dot size is proportional
#'    to the number of genes. See function description for details on
#'    how to include the point size legend beside the heatmap.
#'    The main benefit of using "dotplot" style is that it also indicates
#'    the relative number of genes involved in the enrichment.
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
#' @param top_annotation `HeatmapAnnotation` as produced by
#'    `ComplexHeatmap::HeatmapAnnotation()` or `NULL`, used to display
#'    customized annotation at the top of the heatmap. The order of
#'    columns must match the order of columns in the data displayed
#'    in the heatmap.
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
 style=c("dotplot",
    "heatmap"),
 p_cutoff=mem$p_cutoff,
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
 cluster_rows=TRUE,
 column_fontsize=NULL,
 column_cex=1,
 cluster_columns=FALSE,
 sets=NULL,
 color_by_column=FALSE,
 cex.axis=1,
 lens=3,
 cexCellnote=0.0,
 column_title=NULL,
 row_names_max_width=grid::unit(30, "cm"),
 column_names_max_height=grid::unit(30, "cm"),
 heatmap_legend_param=NULL,
 legend_height=grid::unit(6, "cm"),
 legend_cex=1,
 apply_direction=FALSE,
 direction_cutoff=0,
 gene_count_max=NULL,
 top_annotation=NULL,
 show=NULL,
 use_raster=FALSE,
 do_plot=TRUE,
 ...)
{
   #
   style <- match.arg(style);

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
      # sets <- intersect(rownames(mem$enrichIM), sets);
      sets <- intersect(sets, rownames(mem$enrichIM));

      i_changes <- jamba::vigrep("^enrichIM", names(mem));
      for (i_change in i_changes) {
         i_match <- match(sets, rownames(mem[[i_change]]));
         mem[[i_change]] <- mem[[i_change]][i_match,,drop=FALSE];
      }
   } else {
      sets <- rownames(mem$enrichIM);
   }
   if (any(dim(mem$enrichIM) == 0)) {
      stop("No remaining data after filtering.");
   }
   if (ncol(mem$enrichIM) > 1) {
      if (is.logical(cluster_rows) && TRUE %in% cluster_rows) {
         cluster_rows <- amap::hcluster(
            link="ward",
            jamba::noiseFloor(
               -log10(mem$enrichIM[sets,,drop=FALSE]),
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
               t(-log10(mem$enrichIM[sets,,drop=FALSE])),
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
      row_fontsize <- jamba::noiseFloor(row_cex * 60/(nrow(mem$enrichIM))^(1/2),
         minimum=1,
         ceiling=18);
   }
   if (length(column_fontsize) == 0) {
      column_fontsize <- jamba::noiseFloor(column_cex * 60/(ncol(mem$enrichIM))^(1/2),
         minimum=1,
         ceiling=20);
   }

   if (length(heatmap_legend_param) == 0) {
      heatmap_legend_param <- list(
         border="black",
         labels_gp=gpar(fontsize=10 * legend_cex),
         title_gp=gpar(fontsize=10 * legend_cex),
         legend_height=legend_height);
   }

   # optionally apply direction
   if (apply_direction && "enrichIMdirection" %in% names(mem)) {
      use_matrix <- -log10(mem$enrichIM);
      # use_direction contains z-score values at or above direction_cutoff
      # otherwise it is set to zero
      use_direction <- (
         (abs(mem$enrichIMdirection) >= direction_cutoff) *
         mem$enrichIMdirection);
   } else {
      use_matrix <- -log10(mem$enrichIM);
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
         ctmax <- ceiling(max(mem$enrichIMgeneCount, na.rm=TRUE));
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
      if ("dotplot" %in% style) {
         pt_legend <- ComplexHeatmap::Legend(
            labels=ct_ticks,
            title="Gene Count",
            type="points",
            pch=pch,
            ncol=pt_legend_ncol,
            labels_gp=gpar(fontsize=10 * legend_cex),
            title_gp=gpar(fontsize=10 * legend_cex),
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
         cell_fun_custom <- cell_fun_bivariate(
            list(
               use_direction,
               use_matrix,
               mem$enrichIMgeneCount),
            pch=pch,
            size_fun=ct_approxfun,
            size_by=size_by,
            outline_style="darker",
            col_hm=col_bivariate,
            show=show,
            cex=cexCellnote,
            prefix=c("z-score: ",
               "-log10P: ",
               "genes: ")[show],
            ...
         );
         legend_bivariate <- make_legend_bivariate(col_bivariate,
            ylab="-log10pvalue",
            xlab="z-score");
         anno_legends <- c(anno_legends,
            list(legend_bivariate));
         show_heatmap_legend <- FALSE;
      } else {
         show_heatmap_legend <- TRUE;
         cell_fun_custom <- cell_fun_bivariate(
            list(
               use_matrix,
               use_direction,
               mem$enrichIMgeneCount),
            pch=pch,
            size_fun=ct_approxfun,
            size_by=3,
            outline_style="darker",
            col_hm=col_logp,
            show=show,
            cex=cexCellnote,
            type="univariate",
            prefix=c("z-score: ",
               "-log10P: ",
               "genes: ")[show],
            ...
         );

         cell_fun_custom_old <- function(j, i, x, y, width, height, fill) {
            cell_value <- jamba::rmNA(naValue=0,
               use_matrix[i, j]);
            cell_color <- col_logp(cell_value);
            # draw grid through center of each cell
            grid::grid.lines(x=x + width * c(-1/2, 1/2, NA, 0, 0),
               y=y + height * c(0, 0, NA, -1/2, 1/2),
               gp=grid::gpar(col="grey80"));
            if (abs(cell_value) >= -log10(p_cutoff)) {
               cell_size <- ct_approxfun(mem$enrichIMgeneCount[i, j]);
               grid::grid.points(x=x,
                  y=y,
                  pch=21,
                  default.units="mm",
                  size=grid::unit(cell_size, "mm"),
                  gp=grid::gpar(
                     col=jamba::makeColorDarker(cell_color),
                     fill=cell_color))
               if (cexCellnote > 0.01) {
                  grid::grid.text(round(mem$enrichIMgeneCount[i, j]),
                     x=x,
                     y=y,
                     gp=grid::gpar(
                        fontsize=20 * cexCellnote * 1.05,
                        fontface=2,
                        col=jamba::setTextContrastColor(jamba::setTextContrastColor(cell_color, useGrey=15)))
                  )
                  grid::grid.text(round(mem$enrichIMgeneCount[i, j]),
                     x=x,
                     y=y,
                     gp=grid::gpar(
                        fontsize=20 * cexCellnote,
                        fontface=1,
                        col=jamba::setTextContrastColor(cell_color, useGrey=15))
                  )
               }
            }
         }
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

      # dot plot or heatmap style
      hm <- jamba::call_fn_ellipsis(ComplexHeatmap::Heatmap,
         matrix=use_matrix,
         name=name,
         col=col_logp,
         cluster_rows=cluster_rows,
         row_dend_reorder=row_dend_reorder,
         border=TRUE,
         row_names_gp=grid::gpar(fontsize=row_fontsize),
         row_names_max_width=row_names_max_width,
         row_split=row_split,
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
      hm_sets <- rownames(mem$enrichIM)[ComplexHeatmap::row_order(hm)];
      ## Prepare fresh image colors using p_cutoff and p_floor
      enrichIMcolors <- do.call(cbind, lapply(jamba::nameVector(colnames(mem$enrichIM)), function(i){
         x <- -log10(mem$enrichIM[,i]);
         cr1 <- circlize::colorRamp2(
            breaks=c(-log10(p_cutoff + 1e-10),
               seq(from=-log10(p_cutoff), to=-log10(p_floor), length.out=24)),
            colors=c("white",
               getColorRamp(mem$colorV[i],
                  n=24,
                  trimRamp=c(1, 0),
                  lens=lens)));
         cr1(x);
      }));
      #enrichIMcolors <- colorjam::matrix2heatColors(
      #   x=-log10(mem$enrichIM),
      #   colorV=mem$colorV,
      #   baseline=-log10(p_cutoff),
      #   numLimit=-log10(p_floor),
      #   lens=lens);
      if (do_plot) {
         jamba::imageByColors(enrichIMcolors[hm_sets,,drop=FALSE],
            cellnote=sapply(mem$enrichIM[hm_sets,,drop=FALSE],
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
         matrix=enrichIMcolors[hm_sets,,drop=FALSE],
         cellnote=sapply(mem$enrichIM[hm_sets,,drop=FALSE],
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
#' MultiEnrichMap plot
#'
#' This function is likely to be replaced by `mem2emap()`.
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
         g <- igraph::delete.edges(g, delete_edges);
      }
   }
   if (length(overlap_count) > 0 && "overlap_count" %in% igraph::list.edge.attributes(g)) {
      delete_edges <- which(igraph::E(g)$overlap_count < overlap_count);
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
   label_edges <- head(intersect(label_edges, igraph::list.edge.attributes(g)), 1);
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
#' 1. **Enrichment Heatmap**, using enrichment P-values via
#' `mem_enrichment_heatmap()`. Plot #1.
#' 2. **Gene-Pathway Incidence Matrix Heatmap** using `mem_gene_path_heatmap()`.
#' This step visualizes the pathway clustering to be used by all
#' other plots in the folio. Plot #2.
#' 3. **Cnet Cluster Plot** representing Gene-Pathway clusters as a network,
#' created using `collapse_mem_clusters()`, then plotted with `jam_igraph()`.
#' Plots #3, #4, and #5.
#' 4. **Cnet Exemplar Plots** using exemplar pathways from each
#' gene-pathway cluster, with increase number of exemplars included
#' from each cluster (n per cluster). Cnet `igraph` objects are created
#' using `subsetCnetIgraph()`, then plotted with `jam_graph()`.
#' Plots #6, #7, and #8.
#' 5. **Cnet Individual Cluster Plots** with one plot for each gene-pathway
#' cluster defined above, including all pathways within the cluster.
#' These plots are mostly useful when a particular cluster may
#' have multiple sub-clusters included together. The plots can be useful
#' to understand the relationship between pathways in each cluster.
#' Plots #9, #10, and so on, length equal to `pathway_column_split`.
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
#' @return `list` is returned via `invisible()`, which contains each
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
#' @param do_which integer vector of plots to produce. When `do_which`
#'    is `NULL`, then all plots are produced. This argument is intended
#'    to help produce one plot from a folio, therefore each plot is referred
#'    by the number of the plot, in order.
#' @param p_cutoff numeric value indicating the enrichment P-value threshold
#'    used for `multiEnrichMap()`, but when `NULL` this value is taken
#'    from the `mem` input, or `0.05` is used by default.
#' @param p_floor numeric value indicating the lowest enrichment P-value
#'    used in the color gradient on the Enrichment Heatmap.
#' @param main character string used as a title on Cnet plots.
#' @param use_raster logical indicating whether to use raster heatmaps,
#'    passed to `ComplexHeatmap::Heatmap()`.
#' @param min_gene_ct,min_set_ct integer values passed to
#'    `mem_gene_path_heatmap()`. The `min_gene_ct` requires each set
#'    to contain `min_gene_ct` genes, and `min_set_ct` requires each gene
#'    to be present in at least `min_set_ct` sets.
#' @param min_set_ct_each minimum number of genes required for each set,
#'    required for at least one enrichment test.
#' @param column_method,row_method arguments passed to
#'    `ComplexHeatmap::Heatmap()` which indicate the distance method used
#'    to cluster columns and rows, respectively.
#' @param exemplar_range integer vector (or `NULL`) used to create Cnet
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
#' @param cex.main,cex.sub numeric values passed to `title()` which
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
#' @param do_plot `logical` indicating whether to render each plot.
#'    When `do_plot=FALSE` the plot objects will be created and returned,
#'    but the plot itself will not be rendered. This option may be
#'    useful to generate the full set of figures in one set, then
#'    review each figure one by one in an interactive session.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to downstream functions.
#'    Notably, `sets` is passed to `mem_gene_path_heatmap()` which
#'    allows one to define a specific subset of sets to use in the
#'    gene-pathway heatmap.
#'
#' @export
mem_plot_folio <- function
(mem,
 do_which=NULL,
 p_cutoff=NULL,
 p_floor=1e-10,
 main="",
 use_raster=TRUE,
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
 style="dotplot",
 enrich_im_weight=0.3,
 gene_im_weight=0.5,
 colorize_by_gene=TRUE,
 cluster_color_min_fraction=0.4,
 byCols=c("composite_rank", "minp_rank", "gene_count_rank"),
 edge_bundling="connections",
 apply_direction=NULL,
 do_plot=TRUE,
 verbose=TRUE,
 ...)
{
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
   # First step: Define the Gene-Pathway Heatmap
   #
   # All subsequent plots depend upon mem_gene_path_heatmap()
   if (verbose) {
      jamba::printDebug("mem_plot_folio(): ",
         "Gene-pathway heatmap (pre-emptive)");
   }
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
      ...);

   # Extract hclust object for re-use in the enrichment heatmap
   gp_hm_hclust <- gp_hm@column_dend_param$obj;
   if (length(gp_hm_hclust) == 0) {
      gp_hm_hclust <- gp_hm@column_dend_param$fun(t(gp_hm@matrix));
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
      mem_hm <- jamba::call_fn_ellipsis(
         mem_enrichment_heatmap,
         mem=mem,
         sets=colnames(gp_hm@matrix),
         p_cutoff=p_cutoff,
         p_floor=p_floor,
         color_by_column=color_by_column,
         row_cex=row_cex,
         row_method=row_method,
         row_split=pathway_column_split,
         column_cex=column_cex,
         column_method=column_method,
         style=style,
         apply_direction=apply_direction,
         use_raster=use_raster,
         do_plot=do_plot,
         cluster_rows=gp_hm_hclust,
         ...);
      # mem_hm <- mem_enrichment_heatmap(mem,
      #    p_cutoff=p_cutoff,
      #    p_floor=p_floor,
      #    color_by_column=color_by_column,
      #    row_cex=row_cex,
      #    row_method=row_method,
      #    column_cex=column_cex,
      #    column_method=column_method,
      #    style=style,
      #    ...);
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
      return(invisible(ret_vals));
   }
   ## All subsequent plots depend upon mem_gene_path_heatmap()
   if (verbose) {
      jamba::printDebug("mem_plot_folio(): ",
         "Gene-pathway heatmap");
   }
   if (FALSE) {
      # This section is performed as the first step and is
      # no longer required at this point.
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
         ...);
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
   clusters_mem <- jamba::heatmap_column_order(gp_hm);
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
         return(list(mem=mem, clusters_mem=clusters_mem, ret_vals=ret_vals))
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
         cnet_collapsed %>%
            subsetCnetIgraph(remove_blanks=TRUE,
               repulse=repulse,
               verbose=verbose>1);
      }, error=function(e){
         if (verbose) {
            jamba::printDebug("mem_plot_folio(): ",
               "subsetCnetIgraph() error during ",
               "subsetCnetIgraph(..., remove_blanks=TRUE):");
            print(e);
         }
         cnet_collapsed <- tryCatch({
            cnet_collapsed %>%
               subsetCnetIgraph(remove_blanks=FALSE,
                  repulse=repulse,
                  verbose=verbose>1);
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
            ...);
         cnet_exemplar <- cnet %>%
            subsetCnetIgraph(includeSets=clusters_mem_n$set,
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
            cnet_cluster <- cnet %>%
               subsetCnetIgraph(includeSets=cluster_sets,
                  repulse=repulse,
                  ...);
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

   invisible(ret_vals);
}

#' Draw Heatmap with title and subtitle using grid viewports
#'
#' Draw Heatmap with title and subtitle using grid viewports
#'
#' This function is likely to work only with `grid` or `ComplexHeatmap`
#' input. Other objects that have `draw()` defined may not work
#' properly.
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

