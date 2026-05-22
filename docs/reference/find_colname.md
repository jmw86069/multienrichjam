# Find colname by character string or pattern matching

Find colname by character string or pattern matching

## Usage

``` r
find_colname(
  pattern,
  x,
  max = 1,
  index = FALSE,
  require_non_na = TRUE,
  verbose = FALSE,
  ...
)
```

## Arguments

- pattern:

  character vector containing text strings or regular expression
  patterns.

- x:

  input `data.frame` or other R object that contains colnames.

- max:

  integer maximum number of results to return.

- index:

  logical indicating whether to return the column index as an integer
  vector. When `index=FALSE` it returns the matching `colnames(x)`; when
  `index=TRUE` it returns the matching column numbers as an integer
  vector.

- require_non_na:

  logical indicating whether to require the column to contain non-NA
  values, default is TRUE. The intent of this function is to find
  colnames whose data will match expectations, and when require_non_na
  is TRUE, this function will continue until it finds a column with
  non-NA values.

- ...:

  additional arguments are passed to
  [`jamba::provigrep()`](https://jmw86069.github.io/jamba/reference/provigrep.html).

## Value

character vector with length `max`, or if no pattern match is found it
returns `NULL`. Also if there are no `colnames(x)` then it returns
`NULL`.

## Details

This function is intended to help find a colname given a character
vector of expected values or regular expression patterns. By default it
returns the first matching value, but can return multiple if `max=Inf`.

If there are no `colnames(x)` then `NULL` is returned.

The order of operations:

1.  Match exact string.

2.  Match exact string in case-insensitive manner.

3.  Match the start of each string using
    [`jamba::provigrep()`](https://jmw86069.github.io/jamba/reference/provigrep.html).

4.  Match the end of each string using
    [`jamba::provigrep()`](https://jmw86069.github.io/jamba/reference/provigrep.html).

5.  Match each string using
    [`jamba::provigrep()`](https://jmw86069.github.io/jamba/reference/provigrep.html).

The results from the first successful operation is returned.

When there are duplicate `colnames(x)` only the first unique name is
returned.

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
