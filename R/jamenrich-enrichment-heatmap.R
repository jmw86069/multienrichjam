# jamenrich-enrichment-heatmap.R

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
#' @family custom plot functions
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
#' @param color_non_hits `logical` default FALSE, whether to colorize
#'    cells in the heatmap which do not meet P-value threshold defined
#'    by `p_cutoff.`
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
#'    to use for row and column clustering.
#'    The clustering is performed by `amap::hcluster()`.
#' @param column_method `character` string of the distance method
#'    to use for row and column clustering.
#'    The clustering is performed by `amap::hcluster()`.
#' @param name `character` value passed to `ComplexHeatmap::Heatmap()`,
#'    used as a label above the heatmap color legend.
#' @param row_dend_reorder `logical` indicating whether to reorder the
#'    row dendrogram using the method described in
#'    `ComplexHeatmap::Heatmap()`. The end result is minor reshuffling of
#'    leaf nodes on the dendrogram based upon mean signal in each row,
#'    which can sometimes look cleaner.
#' @param row_dend_width `grid::unit` width of the row dendrogram,
#'    the default is 30 mm in absolute units.
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
#' @param direction_cutoff `numeric` default 0, with directional score
#'    cutoff required for an entry to be colorized by direction.
#' @param gene_count_max `integer` number of genes in the incidence
#'    matrix pathway-enrichment with gene counts, stored as
#'    'enrichIMgeneCount'. Default NULL uses the highest observed
#'    value. This ceiling is used to define the max color.
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
#'    2. "direction"
#'    3. "z-score"
#'    4. "number of genes"
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
#' @returns `ComplexHeatmap::Heatmap` object.
#'
#' @export
mem_enrichment_heatmap <- function(
   mem,
   style=c("dotplot_inverted",
      "dotplot",
      "heatmap"),
   apply_direction=FALSE,
   p_cutoff=NULL,
   color_non_hits=FALSE,
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
   cluster_genes=TRUE,
   column_fontsize=NULL,
   column_cex=1,
   cluster_columns=FALSE,
   cluster_pathways=FALSE,
   sets=NULL,
   color_by_column=FALSE,
   cex.axis=1,
   lens=3,
   cexCellnote=1,
   placementCellnote=c("center",
      "bottomleft", "bottom", "bottomright",
      "right",
      "topleft", "top", "topright", "left"),
   fontCellnote=2,
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
   ...
) {
   #
   style <- match.arg(style);
   Mem <- NULL;
   if (inherits(mem, "Mem")) {
      Mem <- mem;
      mem <- Mem_to_list(Mem);
   } else {
      Mem <- list_to_Mem(mem);
   }
   placementCellnote <- match.arg(placementCellnote);

   if (length(p_cutoff) == 0) {
      if ("p_cutoff" %in% names(thresholds(Mem))) {
         p_cutoff <- thresholds(Mem)$p_cutoff;
      } else if ("cutoffRowMinP" %in% names(thresholds(Mem))) {
         p_cutoff <- thresholds(Mem)[["cutoffRowMinP"]];
      } else {
         p_cutoff <- 1;
      }
   }
   
   # Optional two-tone color scale
   if (p_cutoff < 1 && isTRUE(color_non_hits)) {
      col1 <- colorjam::col_linear_xf(
         -log10(p_floor),
         colramp="Reds",
         lens=lens + 0,
         floor=0
      )
      col2 <- colorjam::col_linear_xf(
         -log10(p_floor),
         colramp="Reds",
         lens=lens,
         floor=0
      )
      # showColors(c(col1(c(0, 0.5, 1)), col2(c(seq(1.5, 10, by=0.5)))))
      col1breaks <- attributes(col1)$breaks;
      col1breaks <- col1breaks[col1breaks > -log10(p_cutoff)];
      use_colors <- c(
         col1(0),
         col1(-log10(p_cutoff + 1e-10)),
         col2(-log10(p_cutoff)),
         col2(col1breaks)
      )
      use_breaks <- c(
         0,
         -log10(p_cutoff + 1e-10),
         -log10(p_cutoff),
         col1breaks
      )
      col_logp <- circlize::colorRamp2(
         breaks=use_breaks,
         colors=use_colors
      )
   } else {
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
   }
   if (apply_direction) {
      col_logp <- colorjam::col_div_xf(
         -log10(p_floor),
         open_floor=color_non_hits,
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
         ctmax <- ceiling(max(enrichIMgeneCount(Mem), na.rm=TRUE));
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
      if (any(grepl("dotplot", fixed=TRUE, style))) {
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
      if (any(c(2, 3) %in% show_enrich) && length(use_direction) == 0) {
         show_enrich <- setdiff(show_enrich, c(2, 3));
      }
      use_prefix <- NULL;
      if (length(show_enrich) > 1) {
         use_prefix <- c(
            "-log10P: ",
            "direction: ",
            "z-score: ",
            "genes: ")[show_enrich]
      }
      # improved cell_fun
      if (apply_direction) {
         # tcount <- jamba::tcount;
         dir_colors <- c("royalblue4", "gold3", "firebrick3");
         dir_colors2 <- c("skyblue", "gold", "indianred1");
         dir_colors3 <- c("white", "white", "white");
         dir_colors2.5 <- sapply(dir_colors2, function(i){
            colorjam::blend_colors(c(i, "white"))
         })
         mcolor <- jamba::rbindList(list(
            dir_colors3,
            dir_colors2.5,
            dir_colors2,
            dir_colors))
         # jamba::imageByColors(mcolor)
         # white_num controls the intensity of the first non-white color
         # in the color gradient
         # white_num <- 2;
         # mcolor2 <- matrix(ncol=3,
         #    c("white", "white", "white",
         #       colorjam::blend_colors(c(dir_colors[1], rep("white", white_num))),
         #       colorjam::blend_colors(c(dir_colors[2], rep("white", white_num))),
         #       colorjam::blend_colors(c(dir_colors[3], rep("white", white_num))),
         #       dir_colors),
         #    byrow=TRUE);
         p_cut_lvl <- 10^(-1 * (ceiling(-log10(p_cutoff)) + 5));
         row_breaks <- c(
            -log10(p_cutoff),
            seq(from=-log10(p_cutoff),
               to=-log10(p_floor),
               length.out=3));
         if (p_cutoff == 1) {
            row_breaks[2] <- row_breaks[1] + p_cut_lvl;
         } else {
            row_breaks[1] <- row_breaks[2] - p_cut_lvl;
         }
         if (p_cutoff < 1) {
            if (isTRUE(color_non_hits)) {
               row_breaks[1] <- 0;
            } else {
               row_breaks <- c(0, row_breaks);
               mcolor <- rbind(mcolor[1, , drop=FALSE], mcolor)
            }
         }
         # if (p_cutoff == 1) {
         #    row_breaks <- tail(row_breaks, -1);
         #    mcolor <- mcolor[-2, , drop=FALSE]
         # }
         col_bivariate <- colorRamp2D(
            column_breaks=seq(from=-2, to=2, length.out=3),
            row_breaks=row_breaks,
            mcolor=mcolor,
            ...);
         size_by <- match("geneCount",
            c("-log10Pvalue",
               # "direction",
               "z-score",
               "geneCount"));
         legend_bivariate <- make_legend_bivariate(
            col_bivariate,
            p_cutoff=p_cutoff,
            ylab="-log10pvalue",
            xlab="direction"
         );
         use_col_fn <- col_bivariate;
         # if ("dotplot_inverted" %in% style) {
         #    use_col_fn <- function(x, y){
         #       rep("#FFFFFF", length.out=length(x))
         #    };
         # }
         ## Switch show_enrich 1:2 to 2:1;
         show_enrich <- ifelse(
            show_enrich %in% 1, 2,
            ifelse(
               show_enrich %in% 2, 1,
               show_enrich
            )
         )
         cell_fun_custom <- cell_fun_bivariate(
            list(
               use_direction,
               use_matrix,
               enrichIMgeneCount(Mem)),
            invert=grepl("invert", fixed=TRUE, style),
            pch=pch,
            size_fun=ct_approxfun,
            size_by=size_by,
            outline_style="darker",
            col_hm=use_col_fn,
            show=show_enrich,
            outline=outline,
            cex=cexCellnote,
            placement=placementCellnote,
            font=fontCellnote,
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
         # remove show_enrich 2 or 3 if no supporting directional data is present
         cell_fun_custom <- cell_fun_bivariate(
            list(
               use_matrix,
               use_direction,
               enrichIMgeneCount(Mem)),
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
