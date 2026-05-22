# Alpha shape calculation

Alpha shape calculation for a set of points, and alpha threshold

## Usage

``` r
ashape(x, y = NULL, alpha, ...)
```

## Arguments

- x, y:

  `numeric` vector with coordinate points.

- alpha:

  `numeric` with alpha threshold to use.

- ...:

  additional arguments are ignored.

## Value

`ashape` object, which is a `list` containing:

- edges: x,y coordinates of Delauney triangulation of the alpha-shape.

- length: length of the alpha-shape.

- alpha: value of alpha used.

- alpha.extremes: `integer` index of points which were alpha-extremes.

- delvor.obj: `delvor` object with Delauney/Voronoi supporting data.

- x: x,y coordinates of input data

## Details

This function is primarily intended to be called by
[`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md),
since that function also iterates `alpha` values until it finds a
suitable, and successful, threshold.

## See also

Other jam utility functions:
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
n <- 300
theta <- runif(n, 0, 2*pi)
r <- sqrt(runif(n, 0.25^2, 0.5^2))
x <- cbind(0.5+r*cos(theta), 0.5+r*sin(theta))
alpha <- 0.1
ashape.obj <- ashape(x, alpha=alpha)

plot(ashape.obj$x, asp=1)
segments(x0=ashape.obj$edges[, "x1"], x1=ashape.obj$edges[, "x2"],
   y0=ashape.obj$edges[, "y1"], y1=ashape.obj$edges[, "y2"])

```
