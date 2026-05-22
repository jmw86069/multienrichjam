# Display colors from bivariate color function

Display colors from bivariate color function

## Usage

``` r
display_colorRamp2D(
  col_fun,
  pretty.n = 12,
  xlab = "",
  ylab = "",
  main = "",
  ...
)
```

## Arguments

- col_fun:

  `function` output from
  [`colorRamp2D()`](https://jmw86069.github.io/multienrichjam/reference/colorRamp2D.md).

- pretty.n:

  `numeric` passed to [`pretty()`](https://rdrr.io/r/base/pretty.html)
  to determine the y-axis and x-axis label breaks, respectively, or
  `NULL` to use breaks encoded in `col_fun`.

- xlab, ylab, main:

  `character` strings passed to
  [`title()`](https://rdrr.io/r/graphics/title.html).

- ...:

  additional arguments are passed to
  [`title()`](https://rdrr.io/r/graphics/title.html).

## Details

This function provides a very rudimentary method for displaying a
bivariate 2-D color matrix out from
[`colorRamp2D()`](https://jmw86069.github.io/multienrichjam/reference/colorRamp2D.md).

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
