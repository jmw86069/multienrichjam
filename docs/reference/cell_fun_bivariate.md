# ComplexHeatmap cell function with bivariant color

ComplexHeatmap cell function with bivariant color

## Usage

``` r
cell_fun_bivariate(
  m,
  prefix = "",
  suffix = "",
  cex = 1,
  col_hm,
  outline = FALSE,
  outline_style = c("none", "darker", "contrast", "lighter", "black", "same"),
  abbrev = FALSE,
  show = NULL,
  rot = 0,
  sep = "\n",
  mcolor = NULL,
  pch = NULL,
  size_fun = NULL,
  size_by = 1,
  grid_color = "grey80",
  type = c("bivariate", "univariate"),
  invert = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- m:

  `list` of 2 or more `matrix` objects. The first two `matrix` objects
  are used for the bivariate color.

- prefix, suffix:

  `character` vectors that define a prefix and suffix for each value in
  `m` for each cell.

- cex:

  `numeric` adjustment for the fontsize used for each label

- col_hm:

  `function` whose first two arguments accept `numeric` values, and
  which returns a single color. Note that when `mcolor` is provided,
  this argument is ignored.

- outline:

  `logical` indicating whether to draw an outline around each heatmap
  cell

- outline_style:

  `character` string indicating the type of outline to draw, which also
  requires `outline=TRUE`. Options:

  - none: uses no outline even when `outline=TRUE`

  - darker: always uses a darker color (or black)

  - contrast: uses a contrasting color `setTextContrastColor()`

  - lighter: always uses a lighter color (or white)

  - black: always uses `"black"`

  - same: use the same color as the fill color

- abbrev:

  `logical` indicating whether numeric values should be abbreviated
  using `jamba::asSize(..., kiloSize=1000)` which effectively reduces
  large numbers to `k` for thousands, `M` for millions (M for Mega), `G`
  for billions (G for Giga), etc.

- show:

  `numeric` indicating which list elements in `m` should be used to
  formulate a cell label, or `NULL` to use no label.

- rot:

  `numeric` rotation in degrees, to rotate labels inside each heatmap
  cell. Mainly useful for heatmaps with extremely tall cells, use
  `rot=90` for vertical text.

- sep:

  `character` string used as a separator between multiple labels inside
  each cell, used only when `show` has more than one value.

- mcolor:

  `character` matrix of R colors, with same
  [`nrow()`](https://rdrr.io/r/base/nrow.html) and
  [`ncol()`](https://rdrr.io/r/base/nrow.html) or each matrix in `m`.
  When `mcolor` is supplied, the colors are used directly, and `col_hm`
  is not used.

- pch:

  `numeric` point type, used only when `size_fun` is also defined.
  Together these arguments allow customized points.

- size_fun:

  `function` used to define point size, only when `pch` is also defined.
  Together these arguments allow customized points.

- grid_color:

  `character` valid R color used with `style="dotplot"` to draw lines
  through the center of each cell.

- type:

  `character` string indicating whether the color function uses
  bivariate or univariate logic. This argument is intended to allow this
  function to be used in both scenarios for consistency.

- invert:

  `logical` indicating whether to invert the color fill, such that each
  cell is filled with color, and the circle is drawn empty on top.

- ...:

  additional arguments are passed to `col_hm()` to allow custom options
  relevant to that function.

## Details

This function serves as a convenient method to use a bivariate color
scale (biscale) to color heatmap cells.

See:

- https://kwstat.github.io/pals/

- https://nowosad.github.io/post/cbc-bp2/

- https://cran.r-project.org/web/packages/biscale/vignettes/biscale.html

This function takes two `numeric` data matrices, a color function that
accepts two numeric values as input and returns a color.

This function can also optionally display a text label inside each
heatmap cell, use argument `show` to indicate which matrix or matrices
in `m` to use for the label.

## See also

Other jam utility functions:
[`ashape()`](https://jmw86069.github.io/multienrichjam/reference/ashape.md),
[`avg_angles()`](https://jmw86069.github.io/multienrichjam/reference/avg_angles.md),
[`avg_colors_by_list()`](https://jmw86069.github.io/multienrichjam/reference/avg_colors_by_list.md),
[`collapse_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/collapse_mem_clusters.md),
[`colorRamp2D()`](https://jmw86069.github.io/multienrichjam/reference/colorRamp2D.md),
[`curateIPAcolnames()`](https://jmw86069.github.io/multienrichjam/reference/curateIPAcolnames.md),
[`deconcat_df2()`](https://jmw86069.github.io/multienrichjam/reference/deconcat_df2.md),
[`display_colorRamp2D()`](https://jmw86069.github.io/multienrichjam/reference/display_colorRamp2D.md),
[`enrichList2geneHitList()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2geneHitList.md),
[`find_colname()`](https://jmw86069.github.io/multienrichjam/reference/find_colname.md),
[`find_enrich_colnames()`](https://jmw86069.github.io/multienrichjam/reference/find_enrich_colnames.md),
[`get_hull_data()`](https://jmw86069.github.io/multienrichjam/reference/get_hull_data.md),
[`get_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/get_igraph_layout.md),
[`gsubs_remove()`](https://jmw86069.github.io/multienrichjam/reference/gsubs_remove.md),
[`handle_igraph_param_list()`](https://jmw86069.github.io/multienrichjam/reference/handle_igraph_param_list.md),
[`isColorBlank()`](https://jmw86069.github.io/multienrichjam/reference/isColorBlank.md),
[`make_legend_bivariate()`](https://jmw86069.github.io/multienrichjam/reference/make_legend_bivariate.md),
[`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md),
[`mem_find_overlap()`](https://jmw86069.github.io/multienrichjam/reference/mem_find_overlap.md),
[`order_colors()`](https://jmw86069.github.io/multienrichjam/reference/order_colors.md),
[`rank_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/rank_mem_clusters.md),
[`rotate_coordinates()`](https://jmw86069.github.io/multienrichjam/reference/rotate_coordinates.md),
[`summarize_node_spacing()`](https://jmw86069.github.io/multienrichjam/reference/summarize_node_spacing.md),
[`with_ht_opts()`](https://jmw86069.github.io/multienrichjam/reference/with_ht_opts.md),
[`xyAngle()`](https://jmw86069.github.io/multienrichjam/reference/xyAngle.md)

## Examples

``` r
set.seed(12);
m <- matrix(rnorm(36)*2.5, ncol=4)
colnames(m) <- LETTERS[1:4]
rownames(m) <- letters[1:9]
m2 <- m;
m2[] <- abs(rnorm(36)*3);
mcolor <- matrix(ncol=3,
   c("white", "white", "white",
   "royalblue4", "gold", "red"),
   byrow=TRUE);
col_bivariate <- colorRamp2D(
   column_breaks=seq(from=-2, to=2, length.out=3),
   row_breaks=seq(from=0, to=5, length.out=2),
   mcolor);
display_colorRamp2D(col_bivariate)


# the heatmap can be created in one step
hm <- ComplexHeatmap::Heatmap(m * m2,
   border=TRUE,
   col=col_bivariate,
   heatmap_legend_param=list(
      color_bar="discrete",
      border=TRUE,
      at=-4:4),
   cell_fun=cell_fun_bivariate(list(m, m2),
      col_hm=col_bivariate,
      prefix=c("-log10P: ", "z-score: "),
      show=2:1),
   show_heatmap_legend=FALSE,
)

lgds <- make_legend_bivariate(col_bivariate,
   ylab="-log10pvalue",
   xlab="z-score");
#> Loading required namespace: gridtext
ComplexHeatmap::draw(hm, annotation_legend_list=lgds)


lgds2 <- make_legend_bivariate(col_bivariate,
   row_breaks=seq(from=0, to=2, by=0.25),
   ylab="-log10pvalue");
ComplexHeatmap::draw(hm, annotation_legend_list=lgds2)


# heatmap using point circles
ctmax <- 6;
point_size_max <- 12;
point_size_min <- 1;
size_fun_custom <- approxfun(
   x=c(1, ctmax),
   yleft=0,
   ties="ordered",
   yright=point_size_max,
   y=c(1,
      point_size_max));
ct_ticks <- seq(from=0, to=6);
ct_tick_sizes <- size_fun_custom(ct_ticks);

hm2 <- ComplexHeatmap::Heatmap(m * m2,
   border=TRUE,
   col=col_bivariate,
   heatmap_legend_param=list(
      color_bar="discrete",
      border=TRUE,
      at=-4:4),
   cell_fun=cell_fun_bivariate(list(m, m2),
      pch=21,
      size_fun=size_fun_custom,
      size_by=2,
      outline_style="black",
      col_hm=col_bivariate,
      prefix=c("-log10P: ", "z-score: "),
      show=NULL),
   show_heatmap_legend=FALSE,
)
ComplexHeatmap::draw(hm2, annotation_legend_list=lgds)

```
