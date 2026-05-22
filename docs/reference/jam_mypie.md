# Vectorized mypie() function for igraph vertex pie polygons

Vectorized mypie() function for igraph vertex pie polygons

## Usage

``` r
jam_mypie(
  x,
  y,
  values,
  radius,
  edges = 200,
  col = NULL,
  angle = 45,
  density = NULL,
  border = NULL,
  frame.color = NULL,
  frame.lwd = NULL,
  lty = NULL,
  lwd = NULL,
  init.angle = 90,
  inner_pie_border = getOption("inner_pie_border", TRUE),
  inner_cex = 1.01,
  xy_max = NULL,
  overclose_polygons = FALSE,
  ...
)
```

## Arguments

- x, y:

  `numeric` coordinate for the center of each `igraph` node.

- values:

  `numeric` vector of relative pie wedge sizes.

- radius:

  `numeric` radius of pie

- edges:

  `integer` number of edges to make a circle, default 200.

- col:

  `character` vector of R colors to fill each pie wedge, in order,
  default NULL.

- angle, density:

  used to draw lines to fill each pie node, passed to
  [`polygon()`](https://rdrr.io/r/graphics/polygon.html). Default
  angle=45, and density=NULL so it draws none.

- border:

  `character` vector of R colors for each pie wedge border, default is
  NULL.

- frame.color:

  `character` R color used around the entire pie circle, default is
  NULL.

- frame.lwd:

  `numeric` line width, default NULL uses `par('lwd')` if a graphics
  device is open, otherwise uses 1.

- lty:

  `numeric` or `character` line type, default NULL.

- lwd:

  `numeric` line width, default NULL will use `par('lwd')` if a graphics
  device is open, otherwise uses 1.

- init.angle:

  `numeric` angle in degrees (0 to 360) where `0` is the top of the
  circle, proceeding clockwide.

- inner_pie_border:

  `logical` whether to apply `pie.border` colors only along the inside
  of each pie wedge polygon, so that adjacent colors will be seen beside
  each other without overlapping adjacent borders. This method is
  currently in development.

- overclose_polygons:

  `logical` indicating whether to close polygon coordinates with an
  extra couple points, to ensure the line join `"ljoin"` is properly
  called when drawing each polygon.

- ...:

  additional arguments are passed to
  [`polygon()`](https://rdrr.io/r/graphics/polygon.html)

## Value

`data.frame` with columns suitable for use in
[`polygon()`](https://rdrr.io/r/graphics/polygon.html) after each column
is expanded to vector form. Columns `x` and `y` are stored `AsIs()` in
`list` format, and should be combined with one `NA` value between each
`numeric` vector. The `NA` values cause
[`polygon()`](https://rdrr.io/r/graphics/polygon.html) to draw a series
of separate polygons in one vectorized step, much faster than drawing a
series of polygons in an iterative loop.

## Details

This is a function called internally by
[`shape.jampie.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.jampie.plot.md),
not intended for direct use.

See
[`reorder_igraph_nodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)
for visual examples.

This function is a light rewrite of `igraph:::mypie()`, except that this
function determines polygon coordinates without drawing them, instead
returns the polygon coordinates to the calling function
[`shape.jampie.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.jampie.plot.md)
which in turn draws all polygons once using the vectorized approach
described for
[`graphics::polygon()`](https://rdrr.io/r/graphics/polygon.html).

See
[`shape.jampie.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.jampie.plot.md)
for more detail.

To disable `pie.border` set to `NA` with `vertex.pie.border=NA` or
`V(g)[[2]]$pie.border <- NA`.

To disable `frame.color` set to `NA` with `vertex.frame.color=NA` or
`V(g)[2]$frame.color <- NA`.

## See also

Other jam igraph shapes:
[`shape.coloredrectangle.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.coloredrectangle.plot.md),
[`shape.ellipse.clip()`](https://jmw86069.github.io/multienrichjam/reference/shape.ellipse.clip.md),
[`shape.ellipse.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.ellipse.plot.md),
[`shape.jampie.plot()`](https://jmw86069.github.io/multienrichjam/reference/shape.jampie.plot.md)

## Examples

``` r
# See reorder_igraph_nodes() for visual examples.
```
