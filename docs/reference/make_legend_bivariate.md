# Display colors from bivariate color function

Display colors from bivariate color function

## Usage

``` r
make_legend_bivariate(
  col_fun,
  pretty.n = 5,
  name = "bivariate",
  xlab = "",
  ylab = "",
  title = "",
  border = TRUE,
  digits = 3,
  title_fontsize = 11,
  legend_fontsize = 10,
  grid_height = grid::unit(5, "mm"),
  grid_width = grid_height,
  row_breaks = NULL,
  column_breaks = NULL,
  row_gap = grid::unit(0, "mm"),
  column_gap = grid::unit(0, "mm"),
  ...
)
```

## Arguments

- col_fun:

  `function` as defined by
  [`colorRamp2D()`](https://jmw86069.github.io/multienrichjam/reference/colorRamp2D.md).

- pretty.n:

  `numeric` value passed to
  [`pretty()`](https://rdrr.io/r/base/pretty.html) to help define a
  suitable number of labels for the x-axis and y-axis color breaks. For
  specific breaks, use `column_breaks`, or `row_breaks`.

- name:

  `character` string used to name the resulting
  [`ComplexHeatmap::Legend`](https://rdrr.io/pkg/ComplexHeatmap/man/Legend.html)
  object, normally only useful when trying to find the `grid` object for
  custom modifications.

- xlab, ylab:

  `character` strings used to define x-axis and y-axis labels,
  effectively the units being displayed. The common values should be
  `xlab="z-score"` or `xlab="direction'`, and `ylab="-log10pvalue"` or
  `y="log10 significance"`.

- title:

  `character`, currently ignored, but may be used in future if necessary
  to display a title above the overall bivariate legend.

- border:

  `logical` indicating whether to draw a border around each color square
  in the color legend. This argument can also be a `character` R color
  value, which will define the color of border drawn around each color
  square.

- digits:

  `numeric` passed to [`format()`](https://rdrr.io/r/base/format.html)
  to define the labels displayed at each position.

- title_fontsize, legend_fontsize:

  `numeric` value passed to `grid::gpar(fontsize)` to define the font
  sizes for legend axis labels, and numerical legend labels,
  respectively.

- grid_height, grid_width:

  [`grid::unit()`](https://rdrr.io/r/grid/unit.html) objects to define
  the exact height and width of each colored square in the color legend.

- row_breaks, column_breaks:

  `numeric` optional vectors which define absolute breaks for row and
  column values displayed in the color legend. When not supplied, these
  values are defined using `pretty` and argument `pretty.n`.

- row_gap, column_gap:

  [`grid::unit()`](https://rdrr.io/r/grid/unit.html) object to define
  optional visible gaps between color squares in the color legend. By
  default, there is no gap. These arguments are provided as a convenient
  way to impose a gap, since the method used does not otherwise provide
  a reasonable way to adjust the spacing.

- ...:

  additional arguments are ignored.

## Value

`ComplexHeatmap::Legends-class` object as returned by
[`ComplexHeatmap::packLegend()`](https://rdrr.io/pkg/ComplexHeatmap/man/packLegend.html),
specifically containing a group of legends, otherwise known as a legend
list.

## Details

This function produces a `"Legend"` object as defined by
[`ComplexHeatmap::Legend()`](https://rdrr.io/pkg/ComplexHeatmap/man/Legend.html).

## See also

Other jam utility functions:
[`ashape()`](https://jmw86069.github.io/multienrichjam/reference/ashape.md),
[`avg_angles()`](https://jmw86069.github.io/multienrichjam/reference/avg_angles.md),
[`avg_colors_by_list()`](https://jmw86069.github.io/multienrichjam/reference/avg_colors_by_list.md),
[`cell_fun_bivariate()`](https://jmw86069.github.io/multienrichjam/reference/cell_fun_bivariate.md),
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
mcolor <- matrix(ncol=3,
   c("white", "salmon1", "firebrick3",
      "white", "lightgoldenrod1", "gold",
      "white", "skyblue", "dodgerblue3"));
row_breaks <- c(0, 0.5, 5);
column_breaks <- c(-1, 0, 1);
rownames(mcolor) <- row_breaks;
colnames(mcolor) <- column_breaks;
jamba::imageByColors(mcolor);


col_fun <- colorRamp2D(column_breaks=column_breaks,
   row_breaks=row_breaks,
   mcolor=mcolor)
lgds <- make_legend_bivariate(col_fun,
   ylab="-log10pvalue",
   xlab="z-score",
   pretty.n=5);
jamba::nullPlot(doBoxes=FALSE);
ComplexHeatmap::draw(lgds)


# same as above with slightly larger grid size
# and slightly larger font sizes
lgds <- make_legend_bivariate(col_fun,
   ylab="-log10pvalue",
   xlab="z-score",
   title_fontsize=14,
   legend_fontsize=12,
   grid_height=grid::unit(7, "mm"),
   pretty.n=5);
jamba::nullPlot(doBoxes=FALSE);
ComplexHeatmap::draw(lgds)


lgds <- make_legend_bivariate(col_fun,
   ylab="-log10pvalue",
   xlab="z-score",
   pretty.n=NULL);
jamba::nullPlot(doBoxes=FALSE);
ComplexHeatmap::draw(lgds)


lgds <- make_legend_bivariate(col_fun,
   ylab="-log10pvalue",
   xlab="z-score",
   column_breaks=c(-1, -0.5,  0, 0.5, 1),
   row_breaks=c(0, 0.25, 0.5, 0.75, 1),
   column_gap=grid::unit(1, "mm"),
   row_gap=grid::unit(1, "mm"),
   pretty.n=5);
jamba::nullPlot(doBoxes=FALSE);
ComplexHeatmap::draw(lgds)

```
