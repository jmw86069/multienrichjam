# Make alpha hull from points

Make alpha hull from points

## Usage

``` r
make_point_hull(
  x,
  expand = 0.05,
  buffer = NULL,
  alpha = NULL,
  seed = 124,
  col = "#FF000033",
  border = "#FF0000FF",
  lwd = 2,
  lty = 1,
  max_iterations = 100,
  do_plot = FALSE,
  add = FALSE,
  hull_method = c("default", "ahull", "alphahull", "igraph", "sf", "chull"),
  smooth = TRUE,
  shape = 1/2,
  label = NULL,
  label.cex = 1,
  label.x.nudge = 0,
  label.y.nudge = 0,
  label_preset = NULL,
  label_adj_preset = label_preset,
  min_points = 1,
  xy_range = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  `numeric` matrix with 2 columns that contains the coordinate of each
  point.

- expand:

  `numeric` value, default 0.05, the buffer width around each point,
  scaled based upon the total range of coordinates, used only when
  `buffer` is not supplied.

- buffer:

  `numeric` value indicating the absolute buffer width around each
  point. This value is used if provided, otherwise `expand` is used to
  derive a value for `buffer`.

- alpha:

  `numeric` value passed to
  [`ashape()`](https://jmw86069.github.io/multienrichjam/reference/ashape.md)
  when hull_method is `"alphahull"`. This value determines the level of
  detail of the resulting hull.

- seed:

  `numeric` seed used with
  [`set.seed()`](https://rdrr.io/r/base/Random.html) to define
  reproducible behavior.

- lwd, lty:

  line width and line type parameters, respectively.

- max_iterations:

  `integer` number of attempts to call
  [`ashape()`](https://jmw86069.github.io/multienrichjam/reference/ashape.md)
  with varying values of `alpha`. Each iteration checks to confirm the
  resulting polygon includes all input points.

- do_plot:

  `logical` indicating whether to plot the polygon output.

- add:

  `logical` used when `do_plot=TRUE` to indicate whether the hull should
  be drawn onto an existing plot device, or whether to open a new plot
  prior to drawing the hull.

- hull_method:

  `character` string indicating the hull method to use:

  - `"default"` - will use `"alphahull"` if the `alphahull` R package is
    available.

  - `"alphahull"` - use
    [`ashape()`](https://jmw86069.github.io/multienrichjam/reference/ashape.md)
    which is the preferred method, in fact the only available option
    that will allow a concave shape in the output.

  - `"igraph"` - calls hidden function `igraph:::convex_hull()` as used
    when drawing `mark.groups` around grouped nodes.

  - `"sf"` - calls
    [`sf::st_convex_hull()`](https://r-spatial.github.io/sf/reference/geos_unary.html),
    with same effective output as `"igraph"`.

  - `"chull"` - calls
    [`grDevices::chull()`](https://rdrr.io/r/grDevices/chull.html),
    again with same effective output as `"igraph"`, but with benefit of
    not incurring additional R package dependencies.

- smooth:

  `logical` indicating whether to smooth the final polygon using
  [`graphics::xspline()`](https://rdrr.io/r/graphics/xspline.html).

- label_preset:

  `character` (default `NULL`) indicating the side to place a label,
  when `label` is provided. Recognized values:
  `"bottom", "top", "left", "right"`. When `NULL` it detects the offset
  from the plot center.

- label_adj_preset:

  `character` (default label_preset) indicating the label adjustment
  relative to the position of the label. In most cases it should equal
  `label_preset`.

- min_points:

  `integer` minimum points to use, default 1 will create a point hull
  even around only 1 point. To require at least 3 points, use
  `min_points=3`.

- xy_range:

  `numeric` range, default NULL, to define the plot range, when
  coordinates in `x` are not sufficient to describe the span of the
  plot. It assumes 1:1 aspect ratio.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

- color, border:

  `character` colors used when `do_plot=TRUE` to draw the resulting hull
  polygon.

## Value

`numeric` matrix with polygon coordinates, where each polygon is
separated by one row that contains `NA` values. This output is
sufficient for vectorized plotting in base R graphics using
[`graphics::polygon()`](https://rdrr.io/r/graphics/polygon.html).

## Details

This function makes an alpha hull around points, calling
[`ashape()`](https://jmw86069.github.io/multienrichjam/reference/ashape.md)
then piecing together the somewhat random set of outer edges into a
coherent polygon.

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
[`mem_find_overlap()`](https://jmw86069.github.io/multienrichjam/reference/mem_find_overlap.md),
[`order_colors()`](https://jmw86069.github.io/multienrichjam/reference/order_colors.md),
[`rank_mem_clusters()`](https://jmw86069.github.io/multienrichjam/reference/rank_mem_clusters.md),
[`rotate_coordinates()`](https://jmw86069.github.io/multienrichjam/reference/rotate_coordinates.md),
[`summarize_node_spacing()`](https://jmw86069.github.io/multienrichjam/reference/summarize_node_spacing.md),
[`with_ht_opts()`](https://jmw86069.github.io/multienrichjam/reference/with_ht_opts.md),
[`xyAngle()`](https://jmw86069.github.io/multienrichjam/reference/xyAngle.md)

## Examples

``` r
set.seed(12)
n <- 22
xy <- cbind(x=sample(seq_len(n), size=n, replace=TRUE),
   y=sample(seq_len(n), size=n, replace=TRUE));
xy <- rbind(xy, xy[1, , drop=FALSE])
x4 <- sf::st_multipoint(xy)

plot(x4, col="red", pch=20, cex=3,
   main="hull_method='ahull'")
phxy <- make_point_hull(x=xy,
   do_plot=TRUE,
   label="ahull",
   min_points=24,
   add=TRUE, xpd=TRUE)
   
# test single-point hull
phxy1 <- make_point_hull(x=head(xy, 1),
   do_plot=TRUE,
   xy_range=par("usr"),
   label="ahull,\nsinglet",
   col="#4169E144",
   border="royalblue",
   add=TRUE, xpd=TRUE)

# test doublet hull
phxy2 <- make_point_hull(x=xy[c(5, 22), ],
   do_plot=TRUE,
   xy_range=par("usr"),
   label="ahull,\nsinglet",
   col="#A020F055",
   border="#A020F0",
   add=TRUE, xpd=TRUE)

# test triplet hull
phxy3 <- make_point_hull(x=xy[c(11, 18), ],
   do_plot=TRUE,
   xy_range=par("usr"),
   label="ahull,\nsinglet",
   col="#A020F055",
   border="#A020F0",
   add=TRUE, xpd=TRUE)


plot(x4, col="red", pch=20, cex=3,
   main="hull_method='chull'")
phxy2 <- make_point_hull(x=xy, expand=0.05, do_plot=TRUE,
   add=TRUE, verbose=TRUE, xpd=TRUE, hull_method="chull")

#> ##  (11:03:18) 22May2026:   make_point_hull(): label.y.nudge:0 
#> ##  (11:03:18) 22May2026:   make_point_hull(): calculated alpha: 10 
#> ##  (11:03:18) 22May2026:   make_point_hull(): Iterating alpha values up to 100 times. 
#> ##  (11:03:18) 22May2026:   make_point_hull(): Iterated 1 times, 126 mxys points. 
```
