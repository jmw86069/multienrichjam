# Bivariate color interpolation function with optional floors

Bivariate color interpolation function with optional floors

## Usage

``` r
col_pvdir_xf(
  p_cutoff = 0.01,
  p_floor = 1e-10,
  direction_cutoff = 0,
  direction_max = 2,
  color_non_hits = FALSE,
  pretty.n = c(10, 5),
  ...
)
```

## Arguments

- p_cutoff:

  `numeric` P-value threshold, default 0.01

- p_floor:

  `numeric` P-value lowest recognized value, default 1e-10

- direction_cutoff:

  `numeric` directional score cutoff, default 0 applies no cutoff.

- direction_max:

  `numeric` directional maximum recognized, default 2.

- color_non_hits:

  `logical` whether to apply color when the P-value does not meet the
  `p_cutoff`, default FALSE uses white for non-hits.

- pretty.n:

  `integer` number of labels used with
  [`pretty()`](https://rdrr.io/r/base/pretty.html) for rows and columns,
  respectively. Default is `10` and `5`.

- ...:

  additional arguments are ignored.

## Value

`function`

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
col_fun <- col_pvdir_xf()
attr(col_fun, "row_labels")
#> [1] "0"  "<2" "2"  "6"  "10"
display_colorRamp2D(col_fun, pretty.n=c(10, 5))


col_fun1 <- col_pvdir_xf(direction_cutoff=1)
attr(col_fun1, "row_labels")
#> [1] "0"  "<2" "2"  "6"  "10"
display_colorRamp2D(col_fun1, pretty.n=c(10, 5))


col_fun <- col_pvdir_xf(
   p_cutoff=0.05,
   direction_cutoff=1
)
display_colorRamp2D(col_fun)

```
