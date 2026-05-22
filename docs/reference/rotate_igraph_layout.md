# Rotate igraph layout coordinates

Rotate igraph layout coordinates, optionally after reflecting
coordinates along one or more coordinate axes.

## Usage

``` r
rotate_igraph_layout(
  g = NULL,
  degrees = 0,
  reflect = c("none", "x", "y", "z"),
  center = NULL,
  center_rule = c("median", "origin", "mean", "min", "max"),
  rotation_axes = c(1, 2),
  spread_labels = TRUE,
  do_reorder = FALSE,
  layout = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- g:

  `igraph` object that contains layout coordinates in graph attributes,
  for example `igraph::graph_attr(g, "layout")`.

- degrees:

  numeric value indicating the degrees to rotate layout coordinates,
  where 360 degrees is one complete rotation.

- reflect:

  `character` string indicating one or more axes to reflect coordinates,
  or `"none"` to reflect no axis.

- center:

  `numeric` coordinates to use as the center, or `center=NULL` to
  calculate the center using `center_rule`.

- center_rule:

  `character` string indicating which rule to apply to determine the
  center coordinates, when `center=NULL`: `"origin"` uses c(0, 0);
  `"mean"` uses the mean of each axis; `"median"` uses the median of
  each axis; `"min"` uses the minimum of each axis; `"max"` uses the max
  of each axis.

- rotation_axes:

  `integer` vector indicating which axis coordinates to rotate, by
  default `c(1, 2)` uses the first two axes.

- spread_labels, do_reorder:

  `logical` indicating whether to call
  [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md),
  and subsequently whether to call
  [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md).

- layout:

  `matrix` with 2 or more columns, when defined this layout is used and
  not the layout from the `g` `igraph` object.

- ...:

  additional arguments are passed to
  [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)
  which calls
  [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)
  when `spread_labels=TRUE` and `do_reorder=TRUE`, or to
  [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)
  when `spread_labels=FALSE` and `do_reorder=TRUE`. Notably, the
  optional argument `sortAttributes` can be passed through those
  functions to affect the node sort priority.

## Value

`igraph` when input `g` is supplied, otherwise `numeric matrix`.

## Details

This function rotates igraph layout coordinates by calling the function
[`rotate_coordinates()`](https://jmw86069.github.io/multienrichjam/reference/rotate_coordinates.md).
The input can either be `g` as `igraph` object, or `layout` as a numeric
`matrix`.

Note that the `reflect` is applied before `degrees`. To change the
order, call this function multiple times.

When both `g` and `layout` are supplied, the coordinates are used from
`layout`, rotated, then stored in the `g` `igraph` object as a graph
attribute, using `igraph::graph_attr(g, "layout")`.

When only `g` is supplied, it is expected to contain layout coordinates
in graph attributes, obtained with `igraph::graph_attr(g, "layout")`.

When only `layout` is supplied, and no `g` `igraph` object is supplied,
this function serves only as a wrapper to
[`rotate_coordinates()`](https://jmw86069.github.io/multienrichjam/reference/rotate_coordinates.md).

Rotation code kindly contributed by Don MacQueen to the `maptools`
package, and is reproduced here to avoid a dependency on `maptools` and
therefore the `sp` package.

This function also calls other useful helper functions, when
`spread_labels=TRUE` it calls
[`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)
to position labels around each node based upon the angle of incoming
edges, which has the effect of reducing label overlaps. When
`do_reorder=TRUE` it calls
[`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)
which sorts nodes within a nodeset by color then by label, to help
visually group similar nodes together in equivalent positions in the
layout.

## See also

Other jam cnet igraph functions:
[`apply_cnet_direction()`](https://jmw86069.github.io/multienrichjam/reference/apply_cnet_direction.md)

Other jam igraph layouts:
[`layout_communities()`](https://jmw86069.github.io/multienrichjam/reference/layout_communities.md),
[`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md),
[`layout_with_qfrf()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfrf.md),
[`relayout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/relayout_with_qfr.md),
[`removeIgraphSinglets()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphSinglets.md),
[`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md),
[`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)

## Examples

``` r
layout <- cbind(0:10, 0:10);
layout_rot50 <- rotate_igraph_layout(layout=layout, degrees=50);
layout_rot40_ctrmean <- rotate_igraph_layout(layout=layout, degrees=40, center_rule="mean");
plot(rbind(layout, layout_rot50, layout_rot40_ctrmean),
   col=rep(c("darkorchid", "darkorange1", "dodgerblue"), each=11),
   pch=rep(c(17, 20, 18), each=11),
   cex=2);


if (require(igraph)) {
   g <- igraph::make_graph( ~ A-B-C-D-A, E-A:B:C:D,
      F-G-H-I-F, J-F:G:H:I,
      K-L-M-N-K, O-K:L:M:N,
      P-Q-R-S-P, T-P:Q:R:S,
      B-F, E-J, C-I, L-T, O-T, M-S,
      C-P, C-L, I-L, I-P);
   g <- relayout_with_qfr(g, repulse=8);
   g2 <- rotate_igraph_layout(g, degrees=45);
   opar <- par("mfrow"=c(1,2));
   on.exit(par(opar));
   jam_igraph(g,
      main="original layout",
      node_factor=0.6,
      label_dist_factor=7);
   jam_igraph(g2,
      main="rotated 45 degrees",
      node_factor=0.6,
      label_dist_factor=7);
}

```
