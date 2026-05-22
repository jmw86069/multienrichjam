# Withr mimic of with_options() for ComplexHeatmap options

Withr mimic of with_options() for ComplexHeatmap options

## Usage

``` r
with_ht_opts(new, code)

local_ht_opts(.new = list(), ..., .local_envir = parent.frame())
```

## Arguments

- new:

  `list` named by heatmap option, with corresponding values as described
  in
  [`ComplexHeatmap::ht_opt()`](https://rdrr.io/pkg/ComplexHeatmap/man/ht_opt.html).

- code:

  `any` valid R code to execute in the temporary environment.

## Value

`any` The results of the evaluation of the `code` argument.

## Details

This function uses the same mechanism as used by
[`withr::with_options()`](https://withr.r-lib.org/reference/with_options.html)
except focuses on options defined by
[`ComplexHeatmap::ht_opt()`](https://rdrr.io/pkg/ComplexHeatmap/man/ht_opt.html).

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
[`xyAngle()`](https://jmw86069.github.io/multienrichjam/reference/xyAngle.md)

## Examples

``` r
ComplexHeatmap::ht_opt("COLUMN_ANNO_PADDING")
#> [1] 1mm
#> 1mm
with_ht_opts(list(COLUMN_ANNO_PADDING=grid::unit(3, "mm")),
   ComplexHeatmap::ht_opt("COLUMN_ANNO_PADDING"))
#> [1] 3mm
#> 3mm
ComplexHeatmap::ht_opt("COLUMN_ANNO_PADDING")
#> [1] 1mm
#> 1mm

test_local <- function() {
   local_ht_opts(list(COLUMN_ANNO_PADDING=grid::unit(3, "mm")))
   print(ComplexHeatmap::ht_opt("COLUMN_ANNO_PADDING"))
}
test_local()
#> [1] 3mm
#> 3mm
ComplexHeatmap::ht_opt("COLUMN_ANNO_PADDING")
#> [1] 1mm
#> 1mm
```
