# Deconcatenate delimited column values in a data.frame

Deconcatenate delimited column values in a data.frame

## Usage

``` r
deconcat_df2(x, column, split = "[,; |/]+", blank = "", ...)
```

## Arguments

- x:

  `data.frame` or compatible object

- column:

  `character` vector with one or more `colnames(x)` that should be
  de-concatenated.

- split:

  `character` pattern used by
  [`strsplit()`](https://rdrr.io/r/base/strsplit.html) to split multiple
  values in each column.

- blank:

  `character` string used to replace entries that would otherwise be
  zero-length as returned by
  [`strsplit()`](https://rdrr.io/r/base/strsplit.html).

- ...:

  additional arguments are ignored.

## Details

This function deconcatenates delimited values in a column of a
`data.frame` by calling
[`strsplit()`](https://rdrr.io/r/base/strsplit.html) on column values,
and repeating values in all other columns to match
[`lengths()`](https://rdrr.io/r/base/lengths.html) following
[`strsplit()`](https://rdrr.io/r/base/strsplit.html).

This function includes a correction for cases where
[`strsplit()`](https://rdrr.io/r/base/strsplit.html) would otherwise
return zero-length entries, and which would otherwise be dropped from
the output. From this function, zero-length entries are replaced with
`blank=""` so these rows are not dropped from the output.

## See also

Other jam utility functions:
[`ashape()`](https://jmw86069.github.io/multienrichjam/reference/ashape.md),
[`avg_angles()`](https://jmw86069.github.io/multienrichjam/reference/avg_angles.md),
[`avg_colors_by_list()`](https://jmw86069.github.io/multienrichjam/reference/avg_colors_by_list.md),
[`cell_fun_bivariate()`](https://jmw86069.github.io/multienrichjam/reference/cell_fun_bivariate.md),
[`collapse_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/collapse_mem_clusters.md),
[`colorRamp2D()`](https://jmw86069.github.io/multienrichjam/reference/colorRamp2D.md),
[`curateIPAcolnames()`](https://jmw86069.github.io/multienrichjam/reference/curateIPAcolnames.md),
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
df <- data.frame(one=c("AB", "BC", "AC"),
   two=c("a,b", "b,c", "a,c"));
deconcat_df2(df, column="two")
#>      one two
#> 1_v1  AB   a
#> 1_v2  AB   b
#> 2_v1  BC   b
#> 2_v2  BC   c
#> 3_v1  AC   a
#> 3_v2  AC   c
```
