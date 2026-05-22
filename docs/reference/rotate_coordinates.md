# Rotate numeric coordinates

Rotate numeric coordinates, optionally after reflecting coordinates
along one or more coordinate axes.

## Usage

``` r
rotate_coordinates(
  x,
  degrees = 0,
  reflect = c("none", "x", "y", "z"),
  center = NULL,
  center_rule = c("median", "origin", "mean", "min", "max"),
  rotation_axes = c(1, 2),
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  `matrix` with 2 or more columns.

- degrees:

  numeric value indicating the degrees to rotate layout coordinates,
  where 360 degrees is one complete rotation.

- reflect:

  `character` string indicating one or more axes to reflect coordinates,
  which flips the position of coordinates along that axis. It is usually
  called to flip x-axis or y-axis coordinates, for example with
  `reflect="x"` or `reflect=1`. Input is handled as follows:

  - if `reflect` contains `"none"`, then reflect is applied to none of
    the coordinate axes, therefore the default
    `reflect=c("none", "x", "y", "z")` will apply no reflection.

  - `character` input: `reflect` values are matched to `colnames(x)`.
    When there are no `colnames(x)`, then `reflect` values of
    `c("x", "y", "z")` are automatically recognized as columns
    `c(1, 2, 3)` respectively.

  - `integer` input is treated as a vector of column index positions,
    for example `reflect=c(2)` will reflect values on the second
    coordinate column.

- center:

  `numeric` coordinates to use as the center, or `center=NULL` to
  calculate the center using `center_rule`.

- center_rule:

  `character` string indicating which rule to apply to determine the
  center coordinates when `center=NULL`. Note that it has little effect
  on most downstream plotting assuming the plot function adjusts x- and
  y-axis ranges to the data range, but may modify the axis ranges as a
  result.

  - `"origin"` uses c(0, 0);

  - `"mean"` uses the mean of each axis;

  - `"median"` uses the median of each axis;

  - `"min"` uses the minimum of each axis;

  - `"max"` uses the max of each axis.

- rotation_axes:

  `integer` vector indicating which axis coordinates to rotate, by
  default `c(1, 2)` uses the first two axes in `x`. Note that
  `rotation_axes` must represent columns present in x.

- ...:

  additional arguments are ignored.

## Value

`numeric matrix` with the same number of columns as the input `x`.

## Details

This function rotates coordinates in two axes, by the angle defined in
`degrees`. It optionally reflects coordinates in one or more axes, which
occurs before rotation.

Note that the `reflect` is applied before `degrees`.

Rotation code kindly contributed by Don MacQueen to the `maptools`
package, and is reproduced here to avoid a dependency on `maptools` and
therefore the `sp` package.

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
[`summarize_node_spacing()`](https://jmw86069.github.io/multienrichjam/reference/summarize_node_spacing.md),
[`with_ht_opts()`](https://jmw86069.github.io/multienrichjam/reference/with_ht_opts.md),
[`xyAngle()`](https://jmw86069.github.io/multienrichjam/reference/xyAngle.md)

## Examples

``` r
layout <- cbind(0:10, 0:10);
layout_rot50 <- rotate_coordinates(x=layout, degrees=50);
layout_rot40_ctrmean <- rotate_coordinates(x=layout, degrees=40, center_rule="mean");
layout_reflectx_ctrmean <- rotate_coordinates(x=layout, reflect="x", center_rule="mean");
plot(rbind(layout, layout_rot50, layout_rot40_ctrmean, layout_reflectx_ctrmean),
   col=rep(c("darkorchid", "darkorange1", "dodgerblue", "red4"), each=11),
   pch=rep(c(17, 20, 18, 17), each=11),
   cex=2);

```
