# Filter mem multienrichment object by Set names

Filter mem multienrichment object by Set names

## Usage

``` r
filter_mem_sets(mem, includeSets = NULL, ...)
```

## Arguments

- mem:

  `list` object returned by
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

- includeSets:

  `character` vector of sets, matched with `colnames(mem$memIM)`

- ...:

  additional arguments are ignored

## Value

`list` object of `mem` data content

## Details

This is intended to be an internal function. It simply takes a
`character` vector of set names, and subsets incidence matrix data in
the `mem` object. It performs no other filtering. This function is
called by
[`subset_mem()`](https://jmw86069.github.io/multienrichjam/reference/subset_mem.md),
the recommended method to subset a `mem` object.

## See also

Other jam utility functions:
[`ashape`](https://jmw86069.github.io/multienrichjam/reference/ashape.md)`()`,
[`avg_angles`](https://jmw86069.github.io/multienrichjam/reference/avg_angles.md)`()`,
[`avg_colors_by_list`](https://jmw86069.github.io/multienrichjam/reference/avg_colors_by_list.md)`()`,
[`bulk_cnet_adjustments`](https://jmw86069.github.io/multienrichjam/reference/bulk_cnet_adjustments.md)`()`,
[`call_fn_ellipsis_deprecated`](https://jmw86069.github.io/multienrichjam/reference/call_fn_ellipsis_deprecated.md)`()`,
[`cell_fun_bivariate`](https://jmw86069.github.io/multienrichjam/reference/cell_fun_bivariate.md)`()`,
[`collapse_mem_clusters`](https://jmw86069.github.io/multienrichjam/reference/collapse_mem_clusters.md)`()`,
[`colorRamp2D`](https://jmw86069.github.io/multienrichjam/reference/colorRamp2D.md)`()`,
[`deconcat_df2`](https://jmw86069.github.io/multienrichjam/reference/deconcat_df2.md)`()`,
[`display_colorRamp2D`](https://jmw86069.github.io/multienrichjam/reference/display_colorRamp2D.md)`()`,
[`enrichList2geneHitList`](https://jmw86069.github.io/multienrichjam/reference/enrichList2geneHitList.md)`()`,
[`filter_mem_genes`](https://jmw86069.github.io/multienrichjam/reference/filter_mem_genes.md)`()`,
[`find_colname`](https://jmw86069.github.io/multienrichjam/reference/find_colname.md)`()`,
[`find_enrich_colnames`](https://jmw86069.github.io/multienrichjam/reference/find_enrich_colnames.md)`()`,
[`get_hull_data`](https://jmw86069.github.io/multienrichjam/reference/get_hull_data.md)`()`,
[`get_igraph_layout`](https://jmw86069.github.io/multienrichjam/reference/get_igraph_layout.md)`()`,
[`gsubs_remove`](https://jmw86069.github.io/multienrichjam/reference/gsubs_remove.md)`()`,
[`handle_igraph_param_list`](https://jmw86069.github.io/multienrichjam/reference/handle_igraph_param_list.md)`()`,
[`isColorBlank`](https://jmw86069.github.io/multienrichjam/reference/isColorBlank.md)`()`,
[`make_legend_bivariate`](https://jmw86069.github.io/multienrichjam/reference/make_legend_bivariate.md)`()`,
[`make_point_hull`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md)`()`,
[`mem_find_overlap`](https://jmw86069.github.io/multienrichjam/reference/mem_find_overlap.md)`()`,
[`order_colors`](https://jmw86069.github.io/multienrichjam/reference/order_colors.md)`()`,
[`rank_mem_clusters`](https://jmw86069.github.io/multienrichjam/reference/rank_mem_clusters.md)`()`,
[`rotate_coordinates`](https://jmw86069.github.io/multienrichjam/reference/rotate_coordinates.md)`()`,
[`subgraph_jam`](https://jmw86069.github.io/multienrichjam/reference/subgraph_jam.md)`()`,
[`subset_mem`](https://jmw86069.github.io/multienrichjam/reference/subset_mem.md)`()`,
[`summarize_node_spacing`](https://jmw86069.github.io/multienrichjam/reference/summarize_node_spacing.md)`()`,
[`xyAngle`](https://jmw86069.github.io/multienrichjam/reference/xyAngle.md)`()`
