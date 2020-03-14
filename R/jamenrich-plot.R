
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
 min_gene_ct=2,
 min_set_ct=3,
 column_fontsize=NULL,
 row_fontsize=NULL,
 row_method="binary",
 column_method="binary",
 cluster_columns=NULL,
 cluster_rows=NULL,
 name="gene_ct",
 p_cutoff=0.05,
 row_split=NULL,
 column_split=NULL,
 auto_split=TRUE,
 column_title=NULL,
 row_title=NULL,
 row_title_rot=90,
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
   memIMgenect <- rowSums(memIM > 0);
   genes <- rownames(memIM)[memIMgenect >= min_gene_ct];
   memIMsetct <- colSums(memIM > 0);
   sets <- colnames(memIM)[memIMsetct >= min_set_ct];
   memIM <- memIM[genes,sets,drop=FALSE];
   ## Additional step to ensure columns and rows are not empty
   memIM <- memIM[,colSums(memIM > 0) > 0,drop=FALSE];
   memIM <- memIM[rowSums(memIM > 0) > 0,,drop=FALSE];
   if (any(dim(memIM) == 0)) {
      stop("No remaining data after filtering.");
   }
   genes <- rownames(memIM);
   sets <- colnames(memIM);
   ## Optional automatic row and column split
   if (auto_split) {
      if (nrow(memIM) < 5) {
         row_split <- NULL;
      } else if (length(row_split) == 0) {
         row_split <- jamba::noiseFloor(floor(nrow(memIM)^(1/2.5)),
            ceiling=12);
         if (length(row_title) == 0) {
            row_title <- letters[seq_len(row_split)];
            row_title_rot <- 0;
         }
      }
      if (ncol(memIM) < 5) {
         column_split <- NULL;
      } else if (length(column_split) == 0) {
         ncol_x <- (5 + ncol(memIM)/1.5) ^ (1/2);
         column_split <- jamba::noiseFloor(floor(ncol_x), ceiling=10);
         if (length(column_title) == 0) {
            column_title <- LETTERS[seq_len(column_split)];
         }
      }
      if (verbose) {
         printDebug("mem_gene_path_heatmap(): ",
            "auto_split row_split:", row_split,
            ", column_split:", column_split);
      }
   }
   if (length(column_split) > 0 && is.numeric(column_split) && length(column_title) == 0) {
      column_title <- LETTERS[seq_len(column_split)];
   }
   if (length(row_split) > 0 && is.numeric(row_split) && length(row_title) == 0) {
      row_title <- LETTERS[seq_len(row_split)];
   }

   ## Automatic fontsize
   if (length(column_fontsize) == 0) {
      row_fontsize <- 60/(nrow(memIM))^(1/2);
   }
   if (length(column_fontsize) == 0) {
      column_fontsize <- 60/(ncol(memIM))^(1/2);
   }

   ## Apply colors to outside annotations
   col_iml1 <- lapply(nameVectorN(mem$colorV), function(i){
      j <- mem$colorV[[i]];
      circlize::colorRamp2(breaks=c(0,1),
         colors=jamba::fixYellow(Hrange=c(60,100), fixup=TRUE,
            jamba::getColorRamp(j, n=3)[1:2]))
   });
   col_iml4 <- lapply(nameVectorN(mem$colorV), function(i){
      j <- mem$colorV[[i]];
      circlize::colorRamp2(
         #breaks=c(0,4),
         breaks=c(-log10(p_cutoff+1e-5), -log10(p_cutoff), 4, 6),
         colors=c("white",
            jamba::fixYellow(Hrange=c(60,100), fixup=TRUE,
               jamba::getColorRamp(j,
                  trimRamp=c(4,2),
                  n=3,
                  lens=3)
            )
         )
      )
   });
   ## Cluster columns and rows
   if (length(cluster_columns) == 0) {
      cluster_columns <- amap::hcluster(
         link="ward",
         cbind(
            jamba::noiseFloor(
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
   hm <- ComplexHeatmap::Heatmap(memIM[genes,sets,drop=FALSE],
      border=TRUE,
      name=name,
      cluster_columns=cluster_columns,
      cluster_rows=cluster_rows,
      clustering_distance_columns=column_method,
      clustering_distance_rows=row_method,
      top_annotation=ComplexHeatmap::HeatmapAnnotation(
         which="column",
         border=TRUE,
         #gp=grid::gpar(col="black"),
         col=col_iml4,
         df=-log10(mem$enrichIM[sets,,drop=FALSE])),
      col=getColorRamp("Reds", lens=2),
      left_annotation=ComplexHeatmap::HeatmapAnnotation(which="row",
         border=TRUE,
         col=col_iml1,
         #gp=grid::gpar(col="black"),
         df=mem$geneIM[genes,,drop=FALSE]),
      row_names_gp=grid::gpar(fontsize=row_fontsize),
      column_names_gp=grid::gpar(fontsize=column_fontsize),
      column_names_rot=90,
      column_title=column_title,
      row_title=row_title,
      row_title_rot=row_title_rot,
      row_split=row_split,
      column_split=column_split,
      ...);
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
      row_names_gp=grid::gpar(fontsize=row_fontsize),
      column_names_gp=grid::gpar(fontsize=column_fontsize),
      cluster_columns=cluster_columns,
      row_dend_width=grid::unit(30, "mm"),
      row_names_max_width=grid::unit(8, "cm"),
      column_title=column_title,
      ...);
   if (color_by_column) {
      hm_sets <- rownames(mem$enrichIM)[ComplexHeatmap::row_order(hm)];
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
#' You can use argument `label_factor_l=list(nodeType=c(Gene=0.01, Set=1))`
#' to hide labels for Gene nodes, and display labels for Set nodes.
#' Note that due to a quirk in `igraph`, setting `label.cex=0`
#' will revert the font to default size, and will not hide the
#' label.
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
#' @param node_factor numeric value multiplied by `V(x)$size` to adjust
#'    the relative size of all nodes by a common numeric scalar value.
#' @param edge_factor numeric value multiplied by `E(x)$width` to adjust
#'    the relative width of all edges by a common numeric scalar value.
#' @param label_factor numeric value multiplied by `V(x)$label.cex`
#'    and `E(x)$label.cex` to adjust the relative size of all labels on
#'    nodes and edges by a common numeric scalar value.
#' @param label_dist_factor numeric value multiplied by `V(x)$label.dist`
#'    to adjust the relative distance of all nodes labels from the node center
#'    by a common numeric scalar value.
#' @param node_factor_l,label_factor_l,label_dist_factor_l `list`
#'    of vectors, where the names of the `list` are attribute
#'    names, and the names of each vector are attributes values.
#'    The vector values are used as scalar multipliers, analogous to
#'    `node_factor`. The purpose is to apply scalar values to different
#'    subsets of nodes. For example, consider:
#'    `node_factor_l=list(nodeType=c(Gene=1, Set=2)`. The list name
#'    `"nodeType"` says to look at `vertex_attr(x, "nodeType")`. Nodes
#'    where `nodeType="Gene"` will use `1`, and where `nodeType="Set"`
#'    will use `2` as the scalar value.
#'
#' @examples
#' ## example showing how to use the list form
#' ## This form resizes nodes where V(g)$nodeType %in% "Gene" by 2x,
#' ## and resizes nodes where V(g)$nodeType %in% "Set" by 3x.
#' node_factor_l <- list(nodeType=c(Gene=2, Set=3));
#'
#' ## This form multiplies label.dist for nodeType="Gene" nodes by 1,
#' ## and multiplies label.dist for nodeType="Set" nodes by 0.5
#' label_dist_factor_l <- list(nodeType=c(Gene=1, Set=0.5))
#'
#' # jam_igraph(g, node_factor_l=node_factor_l, label_dist_factor_l=label_dist_factor_l);
#'
#' @export
jam_igraph <- function
(x,
 ...,
 xlim=c(-1,1),
 ylim=c(-1,1),
 expand=0.03,
 rescale=FALSE,
 node_factor=1,
 node_factor_l=NULL,
 edge_factor=1,
 label_factor=1,
 label_factor_l=NULL,
 label_dist_factor=1,
 label_dist_factor_l=1,
 use_shadowText=FALSE,
 plot_function=jam_plot_igraph,
 verbose=FALSE,
 debug=NULL)
{
   ##
   params <- igraph:::i.parse.plot.params(x, list(...));
   layout <- params("plot", "layout");
   vertex.size <- params("vertex", "size");
   vertex.size2 <- params("vertex", "size2");
   vertex.label.dist <- params("vertex", "label.dist") * label_dist_factor;
   edge.width <- params("edge", "width") * edge_factor;

   if (is.function(label_factor)) {
      if (verbose) {
         jamba::printDebug("jam_igraph(): ",
            "Applying ", "label_factor(label.cex)");
      }
      vertex.label.cex <- label_factor(params("vertex", "label.cex"));
      edge.label.cex <- label_factor(params("edge", "label.cex"));
   } else {
      if (verbose) {
         jamba::printDebug("jam_igraph(): ",
            "Applying ", "label.cex * label_factor");
      }
      vertex.label.cex <- params("vertex", "label.cex") * label_factor;
      edge.label.cex <- params("edge", "label.cex") * label_factor;
   }
   if (is.function(node_factor)) {
      vertex.size <- node_factor(vertex.size);
      vertex.size2 <- node_factor(vertex.size2);
   } else {
      vertex.size <- vertex.size * node_factor;
      vertex.size2 <- vertex.size2 * node_factor;
   }

   ## Handle input as list type
   handle_factor_list <- function(x, attr, factor_l, i_values) {
      if (length(factor_l) > 0 && is.numeric(factor_l)) {
         if (verbose) {
            jamba::printDebug("jam_igraph(): ",
               "Applying ['", attr, "'] * factor");
         }
         i_values <- factor_l * i_values;
      } else if (length(factor_l) > 0 && is.function(factor_l)) {
         if (verbose) {
            jamba::printDebug("jam_igraph(): ",
               "Applying factor(", attr, ")");
         }
         i_values <- factor_l(i_values);
      } else if (length(factor_l) > 0 && is.list(factor_l)) {
         for (i in names(factor_l)) {
            j <- factor_l[[i]];
            for (k in names(j)) {
               x_factor <- factor_l[[i]][[k]];
               if (i %in% list.vertex.attributes(x)) {
                  i_nodes <- (vertex_attr(x, i) %in% k);
                  if (any(i_nodes)) {
                     if (verbose) {
                        jamba::printDebug(sep="",
                           c("jam_igraph(): ",
                              "Applying ",
                              " factor_l[['",i,"']][['",k,"']] to ",
                              jamba::formatInt(sum(i_nodes)),
                              " nodes: ['", attr, "'] * (", x_factor, ")"));
                     }
                     #vertex.label.cex[i_nodes] <- vertex.label.cex[i_nodes] * x_factor;
                     i_values[i_nodes] <- i_values[i_nodes] * x_factor;
                  }
               }
            }
         }
      }
      return(i_values);
   }
   ## label_factor_l=list(nodeType=c(Gene=1, Set=2))
   vertex.label.cex <- handle_factor_list(x, attr="label.cex", label_factor_l, i_values=vertex.label.cex);
   vertex.label.dist <- handle_factor_list(x, attr="label.dist", label_dist_factor_l, i_values=vertex.label.dist);
   vertex.size <- handle_factor_list(x, attr="size", node_factor_l, i_values=vertex.size);

   if (!rescale) {
      dist_factor <- 4;
      if (min(xlim) <= min(layout[,1]) && max(xlim) >= max(layout[,1])) {
         xlim_asis <- TRUE;
      } else {
         xlim <- range(layout[,1]);
         xlim_asis <- FALSE;
      }
      if (length(debug) > 0 && any(c("vertex.label.dist","label.dist","labels") %in% debug)) {
         jamba::printDebug("jam_igraph(): ",
            "xlim before:",
            xlim);
         jamba::printDebug("jam_igraph(): ",
            "head(vertex.size, 20) before:",
            head(vertex.size, 20));
      }
      vertex.size <- vertex.size * diff(xlim) / 2;
      vertex.size2 <- vertex.size2 * diff(xlim) / 2;
      vertex.label.dist <- vertex.label.dist * diff(xlim) / dist_factor;
      if (!xlim_asis) {
         xlim <- xlim + diff(xlim) * c(-1,1) * expand;
      }
      if (min(ylim) <= min(layout[,2]) && max(ylim) >= max(layout[,2])) {
         ylim_asis <- TRUE;
      } else {
         ylim <- range(layout[,2]);
         ylim_asis <- FALSE;
         ylim <- ylim + diff(ylim) * c(-1,1) * expand;
      }
   }
   if (length(debug) > 0 && any(c("vertex.label.dist","label.dist","labels") %in% debug)) {
      jamba::printDebug("jam_igraph(): ",
         "xlim after:",
         xlim);
      jamba::printDebug("jam_igraph(): ",
         "head(vertex.label.dist, 20):",
         head(vertex.label.dist, 20));
      jamba::printDebug("jam_igraph(): ",
         "head(vertex.size, 20) after:",
         head(vertex.size, 20));
   }
   plot_function(x=x,
      ...,
      rescale=rescale,
      vertex.size=vertex.size,
      vertex.size2=vertex.size2,
      vertex.label.dist=vertex.label.dist,
      vertex.label.cex=vertex.label.cex,
      edge.label.cex=edge.label.cex,
      edge.width=edge.width,
      use_shadowText=use_shadowText,
      xlim=xlim,
      ylim=ylim,
      debug=debug);
}

#' Multienrichment folio of summary plots
#'
#' Multienrichment folio of summary plots
#'
#' This function is intended to create multiple summary plots
#' from a Multienrichment result, the `list` output from
#' `multiEnrichMap()`.
#'
#' The plots created:
#'
#' 1. Heatmap using enrichment P-values using `mem_enrichment_heatmap()`
#' 2. Gene-Pathway Heatmap using `mem_gene_pathway_heatmap()`
#' 3. Cnet plot using Gene-Pathway Heatmap collapsed clusters
#' 4. Cnet plot using Gene-Pathway Heatmap cluster examplars (n per cluster)
#' 5. Cnet plots one per cluster
#'
#' @family jam plot functions
#'
#' @return `list` is returned via invisible, which contains each
#'    relevant object.
#'
#' @param mem `list` object created by `multiEnrichMap()`. Specifically
#'    the object is expected to contain `colorV`, `enrichIM`,
#'    `memIM`, `geneIM`.
#' @param p_cutoff numeric value indicating the enrichment P-value threshold
#'    used for `multiEnrichMap()`, but when `NULL` this value is taken
#'    from the `mem` input, or `0.05` is used by default.
#' @param main character string used as a title on Cnet plots.
#' @param p_floor numeric value indicating the lowest enrichment P-value
#'    used in the color gradient on the Enrichment Heatmap.
#' @param use_raster logical indicating whether to use raster heatmaps,
#'    passed to `ComplexHeatmap::Heatmap()`.
#' @param min_gene_ct,min_set_ct integer values passed to
#'    `mem_gene_pathway_heatmap()`. The `min_gene_ct` requires each set
#'    to contain `min_gene_ct` genes, and `min_set_ct` requires each gene
#'    to be present in at least `min_set_ct` sets.
#' @param column_method,row_method arguments passed to
#'    `ComplexHeatmap::Heatmap()` which indicate the distance method used
#'    to cluster columns and rows, respectively.
#' @param exemplar_range integer vector (or `NULL`) used to create Cnet
#'    exemplar plots, using this many exemplars per cluster.
#' @param pathway_column_split integer value passed as `column_split`
#'    to `mem_gene_path_heatmap()`, indicating the number of pathway
#'    clusters to create in the gene-pathway heatmap.
#' @param cex.main,cex.sub numeric values passed to `title()` which
#'    size the default title and sub-title in Cnet plots.
#' @param do_which integer vector of plots to produce. When `do_which`
#'    is `NULL`, then all plots are produced. This argument is intended
#'    to help produce one plot from a folio, therefore each plot is referred
#'    by the number of the plot, in order.
#' @param ... additional arguments are passed to downstream functions.
#'
#' @export
mem_plot_folio <- function
(mem,
 p_cutoff=NULL,
 main="",
 p_floor=1e-6,
 use_raster=TRUE,
 min_gene_ct=2,
 min_set_ct=2,
 column_method="euclidean",
 row_method="euclidean",
 exemplar_range=c(1, 2, 3),
 pathway_column_split=NULL,
 pathway_row_split=NULL,
 edge_color=NULL,
 cex.main=2,
 cex.sub=1.5,
 do_which=NULL,
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


   #############################################################
   ## Enrichment P-value heatmap
   jamba::printDebug("mem_plot_folio(): ",
      "Enrichment P-value heatmap");
   plot_num <- plot_num + 1;
   if (length(do_which) == 0 || plot_num %in% do_which) {
      if (verbose) {
         jamba::printDebug("mem_plot_folio(): ",
            c("plot_num ", plot_num, ": "),
            c("Enrichment P-value Heatmap"),
            sep="");
      }
      mem_hm <- mem_enrichment_heatmap(mem,
         p_cutoff=p_cutoff,
         ...);
      grid_with_title(mem_hm,
         title=main);
      ret_vals$enrichment_hm <- mem_hm;
   }


   #############################################################
   ## Gene-Pathway Heatmap
   if (length(do_which) > 0 && !any(do_which > plot_num)) {
      return(invisible(ret_vals));
   }
   ## All subsequent plots depend upon mem_gene_path_heatmap()
   jamba::printDebug("mem_plot_folio(): ",
      "Gene-pathway heatmap");
   gp_hm <- mem_gene_path_heatmap(mem,
      p_cutoff=p_cutoff,
      row_method=row_method,
      column_split=pathway_column_split,
      row_split=pathway_row_split,
      column_method=column_method,
      use_raster=use_raster,
      min_gene_ct=min_gene_ct,
      min_set_ct=min_set_ct,
      ...);
   ## Obtain heatmap pathway clusters
   clusters_mem <- heatmap_column_order(gp_hm);
   ## Get number of pathway clusters
   pathway_clusters_n <- length(clusters_mem);
   if (verbose) {
      jamba::printDebug("mem_plot_folio(): ",
         c("Defined ", pathway_clusters_n, " pathway clusters."),
         sep="");
   }

   plot_num <- plot_num + 1;
   if (length(do_which) == 0 || plot_num %in% do_which) {
      if (verbose) {
         jamba::printDebug("mem_plot_folio(): ",
            c("plot_num ", plot_num, ": "),
            c("Gene-Pathway Heatmap"),
            sep="");
      }
      caption <- paste0("Hierarchical clustering: distance metric '",
         column_method, "'\n",
         "Data filtering: enrichment P-value <= ", p_cutoff,
         "; genes per set >= ", min_gene_ct,
         "; sets per gene >= ", min_set_ct);
      grid_with_title(gp_hm,
         title=main,
         caption=caption);
      ret_vals$gp_hm <- gp_hm;
   }




   #############################################################
   ## Cnet collapsed
   if (any(c(plot_num + c(1, 2, 3)) %in% do_which)) {
      jamba::printDebug("mem_plot_folio(): ",
         "Preparing Cnet collapsed");
      cnet_collapsed <- tryCatch({
         collapse_mem_clusters(mem,
            clusters_mem,
            verbose=verbose>1,
            return_type="cnet");
      }, error=function(e){
         NULL
      });
      if (length(cnet_collapsed) == 0) {
         return(list(mem=mem, clusters_mem=clusters_mem, ret_vals=ret_vals))
      }
      V(cnet_collapsed)$pie.color <- lapply(V(cnet_collapsed)$pie.color, function(i){
         j <- ifelse(names(i) %in% names(mem$colorV) & !isColorBlank(i),
            mem$colorV[names(i)],
            i);
      });
      V(cnet_collapsed)$coloredrect.color <- lapply(V(cnet_collapsed)$coloredrect.color, function(i){
         j <- ifelse(names(i) %in% names(mem$colorV) & !isColorBlank(i),
            mem$colorV[names(i)],
            i);
      });

      if (verbose) {
         jamba::printDebug("mem_plot_folio(): ",
            "subsetCnetIgraph()");
      }
      cnet_collapsed <- cnet_collapsed %>%
         subsetCnetIgraph(remove_blanks=TRUE,
            verbose=verbose>1);
      if (length(edge_color) > 0) {
         E(cnet_collapsed)$color <- edge_color;
      }
      plot_num <- plot_num + 1;
      if (length(do_which) == 0 || plot_num %in% do_which) {
         ret_vals$cnet_collapsed <- cnet_collapsed;
         if (verbose) {
            jamba::printDebug("mem_plot_folio(): ",
               c("plot_num ", plot_num, ": "),
               c("Cnet collapsed ", "with gene and cluster labels"),
               sep="");
         }
         ## Draw Cnet collapsed
         jam_igraph(cnet_collapsed);
         mem_legend(mem);
         title(sub="Cnet plot using collapsed clusters",
            main=main,
            cex.main=cex.main,
            cex.sub=cex.sub);
      }
      ## Draw Cnet collapsed with top n labels
      #isset <- (V(cnet_collapsed)$nodeType %in% "Set");
      if ("set_labels" %in% igraph::list.vertex.attributes(cnet_collapsed)) {
         V(cnet_collapsed)$label <- ifelse(
            nchar(rmNA(naValue="", igraph:::V(cnet_collapsed)$set_labels)) > 0,
            V(cnet_collapsed)$set_labels,
            V(cnet_collapsed)$name);
      }
      plot_num <- plot_num + 1;
      if (length(do_which) == 0 || plot_num %in% do_which) {
         if (verbose) {
            jamba::printDebug("mem_plot_folio(): ",
               c("plot_num ", plot_num, ": "),
               c("Cnet collapsed ", "with gene and set labels"),
               sep="");
         }
         ret_vals$cnet_collapsed_set <- cnet_collapsed;
         jam_igraph(cnet_collapsed);
         mem_legend(mem);
         title(sub="Cnet plot using collapsed clusters\nlabeled by set",
            main=main,
            cex.main=cex.main,
            cex.sub=cex.sub);
      }
      plot_num <- plot_num + 1;
      if (length(do_which) == 0 || plot_num %in% do_which) {
         if (verbose) {
            jamba::printDebug("mem_plot_folio(): ",
               c("plot_num ", plot_num, ": "),
               c("Cnet collapsed ", "with set labels, without gene labels"),
               sep="");
         }
         V(cnet_collapsed)$label <- ifelse(V(cnet_collapsed)$nodeType %in% "Gene",
            "",
            V(cnet_collapsed)$label);
         ret_vals$cnet_collapsed_set2 <- cnet_collapsed;
         jam_igraph(cnet_collapsed);
         mem_legend(mem);
         title(sub="Cnet plot using collapsed clusters\nlabeled by set\ngene labels hidden",
            main=main,
            cex.main=cex.main,
            cex.sub=cex.sub);
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
            "Preparing cnet for subsetting.");
      }
      cnet <- memIM2cnet(mem,
         ...);
      ## Freshen pie.color by using the original colorV value by name
      V(cnet)$pie.color <- lapply(V(cnet)$pie.color, function(i){
         j <- ifelse(names(i) %in% names(mem$colorV) & !isColorBlank(i),
            mem$colorV[names(i)],
            i);
      });
      ## Freshen coloredrect.color by using the original colorV value by name
      V(cnet)$coloredrect.color <- lapply(V(cnet)$coloredrect.color, function(i){
         j <- ifelse(names(i) %in% names(mem$colorV) & !isColorBlank(i),
            mem$colorV[names(i)],
            i);
      });
      cnet <- cnet %>%
         removeIgraphBlanks();
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
            ## Draw Cnet exemplar
            jam_igraph(cnet_exemplar);
            title(
               sub=paste0("Cnet plot using ",
                  exemplar_n,
                  " exemplar",
                  pluralized,
                  " per cluster"),
               main=main,
               cex.main=cex.main,
               cex.sub=cex.sub);
            mem_legend(mem);
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
                  ...);
            ## Draw Cnet cluster
            jam_igraph(cnet_cluster);
            title(sub=paste0("Cnet plot for cluster ",
               cluster_name),
               main=main,
               cex.main=cex.main,
               cex.sub=cex.sub);
            mem_legend(mem);
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
#' @param ... additional arguments are ignored.
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

   grid::pushViewport(
      grid::viewport(
         layout.pos.row=hm_row));
   if ("gTree" %in% class(object)) {
      grid::grid.draw(object);
   } else {
      grid::grid.draw(
         grid::grid.grabExpr(
            draw(object)));
   }
   grid::popViewport();

   if (length(title_nlines) > 0) {
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

