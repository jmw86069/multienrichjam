# Color ramp for bivariate colors

Color ramp for bivariate colors, using x- and y-axis input values to
apply colors.

## Usage

``` r
colorRamp2D(
  column_breaks = c(-2, -1, 0, 1, 2),
  row_breaks = c(0, 5, 10),
  mcolor = NULL,
  na_color = "grey15",
  return_rgb = FALSE,
  transparency = 0,
  space = "sRGB",
  verbose = FALSE,
  ...
)
```

## Arguments

- column_breaks:

  `numeric` vector of column breaks.

- row_breaks:

  `numeric` vector of row breaks.

- mcolor:

  `character` matrix of R colors, with dimensions equal to the number of
  `column_breaks` and `row_breaks`.

- na_color:

  `character` color used for NA values.

- return_rgb:

  `logical` default FALSE, whether the color function should return RGB
  color format, or `character` hex colors.

- transparency:

  `numeric` default 0, passed to
  [`circlize::colorRamp2()`](https://rdrr.io/pkg/circlize/man/colorRamp2.html).

- space:

  `character` string with color space, default 'sRGB'.

- verbose:

  `logical` whether to print verbose output.

- ...:

  additional arguments are ignored.

## Value

`function` that takes input with arguments:

- x: numeric x-axis value. It only accepts one value.

- return_rgb: `logical` whether to return RGB color

- max_value: maximum value, ignored

- y: numeric y-axis value. It only accepts one value.

- ...: additional arguments are ignored

## Details

This function is intentionally similar to
[`circlize::colorRamp2()`](https://rdrr.io/pkg/circlize/man/colorRamp2.html),
in that it defines a specific color to a specific numeric value. In this
case, it uses two numeric input values, intended to be x- and y-axis
values, columns and rows, respectively.

The default condition it to use divergent colors as x-axis (columns),
and a white-to-dark gradient as y-axis (rows). In this way the y-axis
represents P-value significance, and x-axis represents direction, with
-2 (down), 0 (no direction), and +2 (up).

This function is not intended to be particularly user-friendly, it
simply takes column_breaks, row_breaks, and a matrix of colors, and
returns a color function that generates one color for each x,y value
provided.

## See also

Other jam utility functions:
[`ashape()`](https://jmw86069.github.io/multienrichjam/reference/ashape.md),
[`avg_angles()`](https://jmw86069.github.io/multienrichjam/reference/avg_angles.md),
[`avg_colors_by_list()`](https://jmw86069.github.io/multienrichjam/reference/avg_colors_by_list.md),
[`cell_fun_bivariate()`](https://jmw86069.github.io/multienrichjam/reference/cell_fun_bivariate.md),
[`col_pvdir_xf()`](https://jmw86069.github.io/multienrichjam/reference/col_pvdir_xf.md),
[`collapse_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/collapse_mem_clusters.md),
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
mcolor <- cbind(
   jamba::getColorRamp("Reds", n=3, trimRamp=c(0, 2)),
   jamba::getColorRamp("Oranges", n=3, trimRamp=c(0, 3)),
   jamba::getColorRamp("gold", n=3, trimRamp=c(0, 3)),
   jamba::getColorRamp("Greens", n=3, trimRamp=c(0, 3)),
   jamba::getColorRamp("Blues", n=3, trimRamp=c(0, 2))
)
row_breaks <- c(0, 0.5, 1);
column_breaks <- c(-2, -1, 0, 1, 2);
rownames(mcolor) <- row_breaks;
colnames(mcolor) <- column_breaks;
jamba::imageByColors(mcolor);
title(
   cex.main=1.5,
   cex.lab=1.5,
   main="Bivariate color scale",
   xlab="Directionality",
   ylab="Score");


col_fun <- colorRamp2D(
   column_breaks=column_breaks,
   row_breaks=row_breaks,
   mcolor=mcolor)
display_colorRamp2D(col_fun, pretty.n=c(8, 5));

display_colorRamp2D(col_fun, pretty.n=NULL);


mcolor1 <- matrix(ncol=3,
   c("white", "red",
      "white", "gold",
      "white", "blue3"));
row_breaks1 <- c(0, 1);
column_breaks1 <- c(-1, 0, 1);
rownames(mcolor1) <- row_breaks1;
colnames(mcolor1) <- column_breaks1;
jamba::imageByColors(mcolor1);


col_fun1 <- colorRamp2D(column_breaks=column_breaks1,
   row_breaks=row_breaks1,
   space="LUV",
   mcolor=mcolor1)
display_colorRamp2D(col_fun1);

```
