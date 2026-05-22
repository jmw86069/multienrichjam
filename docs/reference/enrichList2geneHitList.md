# Extract gene hit list from list of enrichResult

Extract gene hit list from list of enrichResult

## Usage

``` r
enrichList2geneHitList(
  enrichList,
  geneColname,
  geneDelim = "[,/ ]",
  make_unique = TRUE,
  verbose = FALSE,
  ...
)
```

## Arguments

- enrichList:

  `list` of `enrichResult` objects

- geneColname:

  `character` string with the column name containing delimited gene
  identifiers.

- geneDelim:

  `character` regular expression used to split delimited gene values. By
  default, `enrichResult` uses '/' forward slash as delimiter, however
  the default here will split any space or comma as well.

- make_unique:

  `logical` default TRUE, whether to return only unique genes per set,
  or potentially multiple genes per set. Typically there should only
  ever be one instance of a gene per set, but through a variety of other
  mechanisms they may exist, for example if two gene identifiers are
  resolved into the same gene symbol.

- verbose:

  `logical` whether to print verbose output.

- ...:

  additional arguments are ignored.

## Value

`list` of character vectors, containing the unique set of genes involved
in each enrichment.

## Details

This function is mainly for internal use in multienrichjam, it takes a
list of `enrichResult` objects, and determines the full set of genes
involved in each `enrichResult`.

Note that genes are sorted using
[`jamba::mixedSort()`](https://jmw86069.github.io/jamba/reference/mixedSort.html)
for alpha-numeric sorting, based upon version sorting.

This function also works with
[`ComplexHeatmap::HeatmapList`](https://rdrr.io/pkg/ComplexHeatmap/man/HeatmapList.html)
objects.

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
