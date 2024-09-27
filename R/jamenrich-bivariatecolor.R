
#' Color ramp for bivariate colors
#'
#' Color ramp for bivariate colors
#'
#' @family jam utility functions
#'
#' @examples
#' mcolor <- matrix(ncol=3,
#'    c("seashell", "salmon1", "firebrick3",
#'       "gray99", "lightgoldenrod1", "gold",
#'       "aliceblue", "skyblue", "dodgerblue3"));
#' row_breaks <- c(0, 0.5, 1);
#' column_breaks <- c(-1, 0, 1);
#' rownames(mcolor) <- row_breaks;
#' colnames(mcolor) <- column_breaks;
#' jamba::imageByColors(mcolor);
#' title(
#'    cex.main=1.5,
#'    cex.lab=1.5,
#'    main="Bivariate color scale",
#'    xlab="Directionality",
#'    ylab="Score");
#'
#' col_fun <- colorRamp2D(column_breaks=column_breaks,
#'    row_breaks=row_breaks,
#'    mcolor=mcolor)
#' display_colorRamp2D(col_fun, pretty.n=c(4, 5));
#' display_colorRamp2D(col_fun, pretty.n=NULL);
#'
#' mcolor1 <- matrix(ncol=3,
#'    c("white", "red",
#'       "white", "gold",
#'       "white", "blue3"));
#' row_breaks1 <- c(0, 1);
#' column_breaks1 <- c(-1, 0, 1);
#' rownames(mcolor1) <- row_breaks1;
#' colnames(mcolor1) <- column_breaks1;
#' jamba::imageByColors(mcolor1);
#'
#' col_fun1 <- colorRamp2D(column_breaks=column_breaks1,
#'    row_breaks=row_breaks1,
#'    space="LUV",
#'    mcolor=mcolor1)
#' display_colorRamp2D(col_fun1);
#'
#' @export
colorRamp2D <- function
(column_breaks,
 row_breaks,
 mcolor,
 na_color="grey15",
 return_rgb=FALSE,
 transparency=0,
 space="sRGB",
 verbose=FALSE,
 ...)
{
   #
   use_space <- space;
   x_funs <- lapply(seq_along(column_breaks), function(colnum){
      if (verbose) {
         jamba::printDebug("colnum: ", colnum,
            ", column_break: ", column_breaks[colnum],
            ", row_breaks: ", row_breaks);
      }
      circlize::colorRamp2(breaks=row_breaks,
         colors=mcolor[, colnum],
         transparency=transparency,
         space=use_space)
   });

   # attributes
   attr <- list(
      breaks=column_breaks,
      column_breaks=column_breaks,
      row_breaks=row_breaks,
      colors=mcolor,
      transparency=transparency,
      na_color=na_color,
      space=use_space)

   fun <- function
   (x=NULL,
      return_rgb=FALSE,
      max_value=1,
      y=NULL,
      ...)
   {
      if (length(x) == 0) {
         stop("x is required")
      }
      if (length(y) == 0) {
         return("transparent");
         #stop("y is required")
      }
      if (is.na(x) || is.na(y)) {
         return(na_color)
      }
      use_colors <- sapply(x_funs, function(x_fun) {
         x_fun(y,
            return_rgb=FALSE)
      });
      y_fun <- circlize::colorRamp2(breaks=column_breaks,
         colors=use_colors,
         transparency=transparency,
         space=use_space);
      y_fun(x,
         return_rgb=return_rgb);
   }
   attributes(fun) <- attr;
   return(fun);
}


#' Display colors from bivariate color function
#'
#' Display colors from bivariate color function
#'
#' This function provides a very rudimentary method for displaying
#' a bivariate 2-D color matrix out from `colorRamp2D()`.
#'
#' @family jam utility functions
#'
#' @param col_fun `function` output from `colorRamp2D()`.
#' @param pretty.n `numeric` passed to `pretty()` to determine the
#'    y-axis and x-axis label breaks, respectively, or
#'    `NULL` to use breaks encoded in `col_fun`.
#' @param xlab,ylab,main `character` strings passed to `title()`.
#' @param ... additional arguments are passed to `title()`.
#'
#' @export
display_colorRamp2D <- function
(col_fun,
 pretty.n=12,
 xlab="",
 ylab="",
 main="",
 ...)
{
   row_breaks <- attr(col_fun, "row_breaks");
   column_breaks <- attr(col_fun, "column_breaks");
   if (length(pretty.n) > 0) {
      row_breaks <- pretty(row_breaks, n=pretty.n)
      column_breaks <- pretty(column_breaks, n=pretty.n)
   }
   new_mcolor1 <- do.call(rbind, lapply(row_breaks, function(y){
      do.call(cbind, lapply(column_breaks, function(x){
         col_fun(x=x,
            y=y)
      }))
   }))
   colnames(new_mcolor1) <- column_breaks;
   rownames(new_mcolor1) <- row_breaks;
   jamba::imageByColors(new_mcolor1);
   title(
      main=main,
      xlab=xlab,
      ylab=ylab,
      ...);
}


#' ComplexHeatmap cell function with bivariant color
#'
#' ComplexHeatmap cell function with bivariant color
#'
#' This function serves as a convenient method to use a
#' bivariate color scale (biscale) to color heatmap cells.
#'
#' See:
#' * https://kwstat.github.io/pals/
#' * https://nowosad.github.io/post/cbc-bp2/
#' * https://cran.r-project.org/web/packages/biscale/vignettes/biscale.html
#'
#' This function takes two `numeric` data matrices, a color function
#' that accepts two numeric values as input and returns a color.
#'
#' This function can also optionally display a text label
#' inside each heatmap cell, use argument `show` to indicate which matrix
#' or matrices in `m` to use for the label.
#'
#' @family jam utility functions
#'
#' @param m `list` of 2 or more `matrix` objects. The first two
#'    `matrix` objects are used for the bivariate color.
#' @param prefix,suffix `character` vectors that define a prefix and
#'    suffix for each value in `m` for each cell.
#' @param cex `numeric` adjustment for the fontsize used for each label
#' @param col_hm `function` whose first two arguments accept `numeric`
#'    values, and which returns a single color. Note that when `mcolor`
#'    is provided, this argument is ignored.
#' @param outline `logical` indicating whether to draw an outline around
#'    each heatmap cell
#' @param outline_style `character` string indicating the type of outline
#'    to draw, which also requires `outline=TRUE`. Options:
#'    * none: uses no outline even when `outline=TRUE`
#'    * darker: always uses a darker color (or black)
#'    * contrast: uses a contrasting color `setTextContrastColor()`
#'    * lighter: always uses a lighter color (or white)
#'    * black: always uses `"black"`
#'    * same: use the same color as the fill color
#' @param abbrev `logical` indicating whether numeric values should
#'    be abbreviated using `jamba::asSize(..., kiloSize=1000)` which
#'    effectively reduces large numbers to `k` for thousands, `M` for
#'    millions (M for Mega), `G` for billions (G for Giga), etc.
#' @param show `numeric` indicating which list elements in `m` should
#'    be used to formulate a cell label, or `NULL` to use no label.
#' @param rot `numeric` rotation in degrees, to rotate labels inside
#'    each heatmap cell. Mainly useful for heatmaps with extremely tall
#'    cells, use `rot=90` for vertical text.
#' @param sep `character` string used as a separator between multiple
#'    labels inside each cell, used only when `show` has more than
#'    one value.
#' @param mcolor `character` matrix of R colors, with same `nrow()`
#'    and `ncol()` or each matrix in `m`. When `mcolor` is supplied,
#'    the colors are used directly, and `col_hm` is not used.
#' @param pch `numeric` point type, used only when `size_fun`
#'    is also defined. Together these arguments allow customized points.
#' @param size_fun `function` used to define point size, only when `pch`
#'    is also defined. Together these arguments allow customized points.
#' @param grid_color `character` valid R color used with `style="dotplot"`
#'    to draw lines through the center of each cell.
#' @param type `character` string indicating whether the color function
#'    uses bivariate or univariate logic. This argument is intended to
#'    allow this function to be used in both scenarios for consistency.
#' @param invert `logical` indicating whether to invert the color fill,
#'    such that each cell is filled with color, and the circle is drawn
#'    empty on top.
#' @param ... additional arguments are passed to `col_hm()` to allow
#'    custom options relevant to that function.
#'
#' @examples
#' set.seed(12);
#' m <- matrix(rnorm(36)*2.5, ncol=4)
#' colnames(m) <- LETTERS[1:4]
#' rownames(m) <- letters[1:9]
#' m2 <- m;
#' m2[] <- abs(rnorm(36)*3);
#' mcolor <- matrix(ncol=3,
#'    c("white", "white", "white",
#'    "royalblue4", "gold", "red"),
#'    byrow=TRUE);
#' col_bivariate <- colorRamp2D(
#'    column_breaks=seq(from=-2, to=2, length.out=3),
#'    row_breaks=seq(from=0, to=5, length.out=2),
#'    mcolor);
#' display_colorRamp2D(col_bivariate)
#'
#' # the heatmap can be created in one step
#' hm <- ComplexHeatmap::Heatmap(m * m2,
#'    border=TRUE,
#'    col=col_bivariate,
#'    heatmap_legend_param=list(
#'       color_bar="discrete",
#'       border=TRUE,
#'       at=-4:4),
#'    cell_fun=cell_fun_bivariate(list(m, m2),
#'       col_hm=col_bivariate,
#'       prefix=c("-log10P: ", "z-score: "),
#'       show=2:1),
#'    show_heatmap_legend=FALSE,
#' )
#'
#' lgds <- make_legend_bivariate(col_bivariate,
#'    ylab="-log10pvalue",
#'    xlab="z-score");
#' ComplexHeatmap::draw(hm, annotation_legend_list=lgds)
#'
#' lgds2 <- make_legend_bivariate(col_bivariate,
#'    row_breaks=seq(from=0, to=2, by=0.25),
#'    ylab="-log10pvalue");
#' ComplexHeatmap::draw(hm, annotation_legend_list=lgds2)
#'
#' # heatmap using point circles
#' ctmax <- 6;
#' point_size_max <- 12;
#' point_size_min <- 1;
#' size_fun_custom <- approxfun(
#'    x=c(1, ctmax),
#'    yleft=0,
#'    ties="ordered",
#'    yright=point_size_max,
#'    y=c(1,
#'       point_size_max));
#' ct_ticks <- seq(from=0, to=6);
#' ct_tick_sizes <- size_fun_custom(ct_ticks);
#'
#' hm2 <- ComplexHeatmap::Heatmap(m * m2,
#'    border=TRUE,
#'    col=col_bivariate,
#'    heatmap_legend_param=list(
#'       color_bar="discrete",
#'       border=TRUE,
#'       at=-4:4),
#'    cell_fun=cell_fun_bivariate(list(m, m2),
#'       pch=21,
#'       size_fun=size_fun_custom,
#'       size_by=2,
#'       outline_style="black",
#'       col_hm=col_bivariate,
#'       prefix=c("-log10P: ", "z-score: "),
#'       show=NULL),
#'    show_heatmap_legend=FALSE,
#' )
#' ComplexHeatmap::draw(hm2, annotation_legend_list=lgds)
#'
#' @export
cell_fun_bivariate <- function
(m,
 prefix="",
 suffix="",
 cex=1,
 col_hm,
 outline=FALSE,
 outline_style=c("none",
    "darker",
    "contrast",
    "lighter",
    "black",
    "same"),
 abbrev=FALSE,
 show=NULL,
 rot=0,
 sep="\n",
 mcolor=NULL,
 pch=NULL,
 size_fun=NULL,
 size_by=1,
 grid_color="grey80",
 type=c("bivariate",
    "univariate"),
 invert=FALSE,
 verbose=FALSE,
 ...)
{
   outline_style <- match.arg(outline_style);
   type <- match.arg(type);
   if (TRUE %in% invert) {
      grid_color <- NULL;
   }
   if (!is.list(m) || (is.list(m) && length(m) < 2)) {
      stop("Input must be a list with two or more matrix objects.")
   }
   if (length(prefix) == 0) {
      prefix <- "";
   }
   if (length(suffix) == 0) {
      suffix <- "";
   }
   if (length(mcolor) > 0) {
      if (!is.matrix(mcolor) || !"character" %in% class(mcolor)) {
         stop("mcolor must be a character matrix");
      }
      if (nrow(mcolor) == nrow(m[[1]]) && ncol(mcolor) == ncol(m[[1]])) {
         # use mcolor as direct color drop-in and ignore col_hm()
      } else {
         stop("mcolor must have identical nrow() and ncol() as m[[1]]");
      }
   }

   prefix <- rep(prefix, length.out=length(show));
   suffix <- rep(suffix, length.out=length(show));
   sep <- c("",
      head(rep(sep, length.out=length(show)), -1));
   cell_fun_return <- function(j, i, x, y, width, height, fill) {
      # convert NA to zero
      #cell_value1 <- jamba::rmNA(naValue=0, m[[1]][i, j]);
      #cell_value2 <- jamba::rmNA(naValue=0, m[[2]][i, j]);
      # do not convert NA to zero
      cell_value1 <- m[[1]][i, j];
      if (length(mcolor) > 0) {
         cell_color <- mcolor[i, j];
      } else if ("bivariate" %in% type) {
         cell_value2 <- m[[2]][i, j];
         cell_color <- col_hm(x=cell_value1,
            y=cell_value2,
            ...);
      } else {
         cell_color <- col_hm(cell_value1);
      }
      cell_label <- "";
      for (k1 in seq_along(show)) {
         k <- show[k1];
         mx <- m[[k]];
         cell_value2 <- jamba::rmNA(naValue=0, mx[i, j]);
         if (abbrev && max(cell_value2) >= 1000) {
            cell_label2 <- jamba::asSize(cell_value2,
               kiloSize=1000,
               unitType="",
               sep="");
         } else {
            cell_label2 <- format(cell_value2,
               big.mark=",",
               trim=TRUE)
         }
         cell_label <- paste(cell_label,
            paste0(prefix[k1],
               cell_label2,
               suffix[k1]),
            sep=sep[k1]);
         if (verbose && i == 1 && j == 1) {
            print("cell_label:");
            print(cell_label);
         }
         # remove leading sep
         #cell_label <- gsub(paste0("^", sep[1]),
         #   "",
         #   cell_label);
      }

      outline_col <- NA;
      if (outline && !"none" %in% outline_style) {
         if ("contrast" %in% outline_style) {
            outline_col <- jamba::setTextContrastColor(cell_color);
         } else if ("lighter" %in% outline_style) {
            outline_col <- jamba::makeColorDarker(cell_color,
               darkFactor=-1.3,
               sFactor=-1.2);
         } else if ("darker" %in% outline_style) {
            outline_col <- jamba::makeColorDarker(cell_color,
               darkFactor=1.3,
               sFactor=-1.1);
         } else if ("same" %in% outline_style) {
            outline_col <- cell_color;
         } else if ("black" %in% outline_style) {
            outline_col <- "black";
         }
      }
      if (length(pch) > 0 && "function" %in% class(size_fun)) {
         mx <- m[[size_by]];
         cell_value2 <- jamba::rmNA(naValue=0, mx[i, j]);
         cell_size <- size_fun(cell_value2);
         # draw grid through center of each cell
         if (length(grid_color) > 0) {
            grid::grid.lines(x=x + width * c(-1/2, 1/2, NA, 0, 0),
               y=y + height * c(0, 0, NA, -1/2, 1/2),
               gp=grid::gpar(col=grid_color));
         }
         if (TRUE %in% invert) {
            # fill the cell
            if (FALSE %in% outline) {
               rect_col <- NA;
            } else {
               rect_col <- outline_col;
            }
            grid::grid.rect(x=x,
               y=y,
               width=width,
               height=height,
               gp=grid::gpar(
                  col=rect_col,
                  fill=cell_color));
            text_col <- jamba::setTextContrastColor(cell_color,
               useGrey=15);
            # draw the point
            grid::grid.points(x=x,
               y=y,
               pch=pch,
               default.units="mm",
               size=grid::unit(cell_size, "mm"),
               gp=grid::gpar(
                  # col=outline_col,
                  col="black",
                  fill="#FFFFFF"));
         } else {
            grid::grid.points(x=x,
               y=y,
               pch=pch,
               default.units="mm",
               size=grid::unit(cell_size, "mm"),
               gp=grid::gpar(
                  col=outline_col,
                  fill=cell_color));
         }
         text_col <- jamba::setTextContrastColor(cell_color,
            useGrey=15);
         # text_col <- "black";
      } else {
         grid::grid.rect(x=x,
            y=y,
            width=width,
            height=height,
            gp=grid::gpar(
               col=outline_col,
               fill=cell_color));
         text_col <- jamba::setTextContrastColor(cell_color);
      }

      if (cex > 0 && !"" %in% cell_label) {
         fontsize <- (10 * cex);
         # grid::grid.text(cell_label,
         if (jamba::check_pkg_installed("shadowtext")) {
            shadowtext::grid.shadowtext(cell_label,
               x=x,
               y=y,
               rot=rot,
               bg.colour=jamba::alpha2col(alpha=0.3,
                  jamba::setTextContrastColor(text_col)),
               bg.r=0.2,
               gp=grid::gpar(
                  fontsize=fontsize,
                  fontface=1,
                  col=text_col
               ));
         } else {
            grid::grid.text(cell_label,
               x=x,
               y=y,
               rot=rot,
               gp=grid::gpar(
                  fontsize=fontsize,
                  fontface=1,
                  col=text_col
               ));
         }
      }
   }
   cell_fun_return
}


#' Display colors from bivariate color function
#'
#' Display colors from bivariate color function
#'
#' This function produces a `"Legend"` object as defined by
#' `ComplexHeatmap::Legend()`.
#'
#' @family jam utility functions
#'
#' @param col_fun `function` as defined by `colorRamp2D()`.
#' @param pretty.n `numeric` value passed to `pretty()` to help define
#'    a suitable number of labels for the x-axis and y-axis color breaks.
#'    For specific breaks, use `column_breaks`, or `row_breaks`.
#' @param name `character` string used to name the resulting
#'    `ComplexHeatmap::Legend` object, normally only useful when
#'    trying to find the `grid` object for custom modifications.
#' @param xlab,ylab `character` strings used to define x-axis and y-axis
#'    labels, effectively the units being displayed. The common
#'    values should be `xlab="z-score"` or `xlab="direction'`,
#'    and `ylab="-log10pvalue"` or `y="log10 significance"`.
#' @param title `character`, currently ignored, but may be used in future
#'    if necessary to display a title above the overall bivariate legend.
#' @param border `logical` indicating whether to draw a border around
#'    each color square in the color legend.
#'    This argument can also be a `character` R color value, which will
#'    define the color of border drawn around each color square.
#' @param digits `numeric` passed to `format()` to define the labels
#'    displayed at each position.
#' @param title_fontsize,legend_fontsize `numeric` value passed to
#'    `grid::gpar(fontsize)` to define the font sizes for
#'    legend axis labels, and numerical legend labels, respectively.
#' @param grid_height,grid_width `grid::unit()` objects to define the
#'    exact height and width of each colored square in the color legend.
#' @param row_breaks,column_breaks `numeric` optional vectors which define
#'    absolute breaks for row and column values displayed in the color
#'    legend. When not supplied, these values are defined using `pretty`
#'    and argument `pretty.n`.
#' @param row_gap,column_gap `grid::unit()` object to define optional
#'    visible gaps between color squares in the color legend. By default,
#'    there is no gap. These arguments are provided as a convenient way
#'    to impose a gap, since the method used does not otherwise provide
#'    a reasonable way to adjust the spacing.
#' @param ... additional arguments are ignored.
#'
#' @return `ComplexHeatmap::Legends-class` object as returned by
#'    `ComplexHeatmap::packLegend()`, specifically containing a group
#'    of legends, otherwise known as a legend list.
#'
#' @examples
#' mcolor <- matrix(ncol=3,
#'    c("white", "salmon1", "firebrick3",
#'       "white", "lightgoldenrod1", "gold",
#'       "white", "skyblue", "dodgerblue3"));
#' row_breaks <- c(0, 0.5, 5);
#' column_breaks <- c(-1, 0, 1);
#' rownames(mcolor) <- row_breaks;
#' colnames(mcolor) <- column_breaks;
#' jamba::imageByColors(mcolor);
#'
#' col_fun <- colorRamp2D(column_breaks=column_breaks,
#'    row_breaks=row_breaks,
#'    mcolor=mcolor)
#' lgds <- make_legend_bivariate(col_fun,
#'    ylab="-log10pvalue",
#'    xlab="z-score",
#'    pretty.n=5);
#' jamba::nullPlot(doBoxes=FALSE);
#' ComplexHeatmap::draw(lgds)
#'
#' # same as above with slightly larger grid size
#' # and slightly larger font sizes
#' lgds <- make_legend_bivariate(col_fun,
#'    ylab="-log10pvalue",
#'    xlab="z-score",
#'    title_fontsize=14,
#'    legend_fontsize=12,
#'    grid_height=grid::unit(7, "mm"),
#'    pretty.n=5);
#' jamba::nullPlot(doBoxes=FALSE);
#' ComplexHeatmap::draw(lgds)
#'
#' lgds <- make_legend_bivariate(col_fun,
#'    ylab="-log10pvalue",
#'    xlab="z-score",
#'    pretty.n=NULL);
#' jamba::nullPlot(doBoxes=FALSE);
#' ComplexHeatmap::draw(lgds)
#'
#' lgds <- make_legend_bivariate(col_fun,
#'    ylab="-log10pvalue",
#'    xlab="z-score",
#'    column_breaks=c(-1, -0.5,  0, 0.5, 1),
#'    row_breaks=c(0, 0.25, 0.5, 0.75, 1),
#'    column_gap=grid::unit(1, "mm"),
#'    row_gap=grid::unit(1, "mm"),
#'    pretty.n=5);
#' jamba::nullPlot(doBoxes=FALSE);
#' ComplexHeatmap::draw(lgds)
#'
#' @export
make_legend_bivariate <- function
(col_fun,
 pretty.n=5,
 name="bivariate",
 xlab="",
 ylab="",
 title="",
 border=TRUE,
 digits=3,
 title_fontsize=11,
 legend_fontsize=10,
 grid_height=grid::unit(5, "mm"),
 grid_width=grid_height,
 row_breaks=NULL,
 column_breaks=NULL,
 row_gap=grid::unit(0, "mm"),
 column_gap=grid::unit(0, "mm"),
 ...)
{
   if (length(pretty.n) > 0) {
      pretty.n <- rep(pretty.n, length.out=2);
   }

   if (length(column_breaks) == 0) {
      column_breaks <- attr(col_fun, "column_breaks");
      if (length(pretty.n) > 0) {
         column_breaks <- pretty(column_breaks,
            n=pretty.n[1])
      }
   }
   if (length(row_breaks) == 0) {
      row_breaks <- attr(col_fun, "row_breaks");
      if (length(pretty.n) > 0) {
         row_breaks <- pretty(row_breaks,
            n=pretty.n[2]);
      }
   }
   row_breaks <- rev(row_breaks);

   legend_seq <- c(0, seq_along(column_breaks), length(column_breaks) + 1);
   legend_columns <- lapply(legend_seq, function(i){
      labels <- rep("", length(row_breaks));
      maintitle <- "";
      title_position <- "topcenter";

      # First, last column are dedicated to blank legend for test labels only
      title_gp <- grid::gpar(fontsize=legend_fontsize,
         fontface="bold");
      labels_gp <- grid::gpar(fontsize=legend_fontsize)

      if (i %in% c(0)) {
         title_gp <- grid::gpar(fontsize=title_fontsize,
            fontface="bold");
         maintitle <- paste0(ylab,
            paste0(rep(" ", 10), collapse=""))
         title_position <- "leftcenter-rot";
         column_break <- column_breaks[1];
         igrid_width <- grid::unit(0.1, "mm");
         use_border <- "transparent";
         row_colors <- rep("transparent", length(row_breaks));
      } else if (i %in% c(length(column_breaks) + 1)) {
         maintitle <- ".\n.";
         title_gp <- grid::gpar(fontsize=title_fontsize,
            fontface="bold",
            col="transparent");
         labels <- format(row_breaks,
            digits=digits,
            trim=TRUE);
         column_break <- column_breaks[1];
         igrid_width <- grid::unit(0.1, "mm");
         use_border <- "transparent";
         row_colors <- rep("transparent", length(row_breaks));
      } else {
         column_break <- column_breaks[i];
         igrid_width <- grid_width;
         use_border <- border;
         title_gp <- grid::gpar(fontsize=legend_fontsize)
         #
         if (i == ceiling(length(column_breaks) / 2)) {
            maintitle <- ComplexHeatmap::gt_render(
               margin=grid::unit(c(0, 0, 0, 0), "pt"),
               padding=grid::unit(c(0, 0, 0, 0), "pt"),
               r=grid::unit(0, "pt"),
               paste0("<b><span style='font-size:", title_fontsize, "pt'>",
                  xlab, "</span></b><br>",
                  "<span style='font-size:", legend_fontsize, "pt'>",
                  column_break, "</span>"))
            title_position <- "topcenter";
         } else {
            maintitle <- ComplexHeatmap::gt_render(
               margin=grid::unit(c(0, 0, 0, 0), "pt"),
               padding=grid::unit(c(0, 0, 0, 0), "pt"),
               r=grid::unit(0, "pt"),
               paste0("<b><span style='font-size:", title_fontsize, "pt'>",
                  " </span></b><br>",
                  "<span style='font-size:", legend_fontsize, "pt'>",
                  column_break, "</span>"))
            title_position <- "topcenter";
         }
         row_colors <- sapply(row_breaks, function(row_break) {
            col_fun(x=column_break,
               y=row_break)
         })
      }

      # create Legend object
      lgd <- ComplexHeatmap::Legend(
         name=name,
         at=row_breaks,
         labels=labels,
         labels_gp=labels_gp,
         title=maintitle,
         title_position=title_position,
         grid_height=grid_height,
         grid_width=igrid_width,
         legend_gp=grid::gpar(fill=row_colors),
         title_gp=title_gp,
         border=use_border,
         row_gap=row_gap,
         ...)
      lgd;
   });

   legend_list <- do.call(function(...){
      ComplexHeatmap::packLegend(...,
         direction="horizontal",
         column_gap=column_gap + grid::unit(-1, "mm"))
      },
      legend_columns);

   # Note: This section is disabled, because although it works
   # for draw(legend_list) alone, it does not work with
   # draw(heatmaplist, annotation_legend_list=legend_list)
   #
   # now assemble a larger legend grob
   if (FALSE &&
         ((length(xlab) == 1 && nchar(xlab) > 0) ||
         #(length(ylab) == 1 && nchar(ylab) > 0) ||
         (length(title) == 1 && nchar(title) > 0))) {
      legend_body <- legend_list@grob;

      xlab_gp <- grid::gpar(fontsize=11, fontface=2);
      if (length(xlab) == 0 || nchar(xlab) == 0) {
         xlab_height <- grid::unit(0, "mm");
         xlab_width <- grid::unit(0, "mm");
      } else {
         xlab_grob <- grid::textGrob(xlab, gp=xlab_gp)
         xlab_height <- grid::convertHeight(grid::grobHeight(xlab_grob), "mm")
         xlab_width <- grid::convertWidth(grid::grobWidth(xlab_grob), "mm")
      }
      title_gp <- grid::gpar(fontsize=12, fontface=2);
      if (length(title) == 0 || nchar(title) == 0) {
         title_height <- grid::unit(0, "mm");
         title_width <- grid::unit(0, "mm");
      } else {
         title_grob <- grid::textGrob(title, gp=title_gp)
         title_height <- grid::convertHeight(grid::grobHeight(title_grob), "mm")
         title_width <- grid::convertWidth(grid::grobWidth(title_grob), "mm")
      }

      legend_width <- grid::convertWidth(grid::grobWidth(legend_body), "mm")
      legend_height <- grid::convertHeight(grid::grobHeight(legend_body), "mm")

      total_width <- grid::unit(max(c(xlab_width,
         title_width,
         legend_width)), "mm");

      total_height <- legend_height + xlab_height + xlab_height + title_height + title_height;

      gf <- grid::grobTree(
         grid::textGrob(xlab,
            x=grid::unit(0.5, "npc"),
            y=grid::unit(0, "npc") + xlab_height * 0.8,
            just=c("centre", "bottom"),
            gp=xlab_gp,
            rot=0),
         grid::textGrob(title,
            x=grid::unit(0, "npc"),
            y=grid::unit(1, "npc") + title_height * 1.8,
            just=c("left", "top"),
            gp=title_gp,
            rot=0),
         ComplexHeatmap:::edit_vp_in_legend_grob(legend_body,
            x=grid::unit(1, "npc"),
            y=grid::unit(1, "npc"),
            valid.just=c(1, 1)),
         # grid::grid.rect(x=grid::unit(0.5, "npc"),
         #    y=grid::unit(0.5, "npc"),
         #    width=grid::unit(1, "npc"),
         #    height=grid::unit(1, "npc"),
         #    gp=grid::gpar(col="red", fill=NA)),
         vp=grid::viewport(
            width=total_width,
            height=total_height,
            gp=grid::gpar(fontsize=11,
               lineheight=0.8)),
         cl="legend"
      );
      legend_list@grob <- gf;
   }
   if (length(name) == 1 && nchar(name) > 0) {
      legend_list@name <- name;
   }

   return(legend_list);
}


cell_fun_points <- function
(j,
   i,
   x,
   y,
   width,
   height,
   fill,
   col_hm)
{
   cell_value <- jamba::rmNA(naValue=0,
      use_matrix[i, j]);
   cell_color <- col_hm(cell_value);
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

