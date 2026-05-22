# Draw ellipse

Draw ellipse

## Usage

``` r
drawEllipse(
  x,
  y,
  a = 1,
  b = 1,
  angle = 0,
  segment = NULL,
  arc.only = TRUE,
  nv = 100,
  deg = TRUE,
  border = NULL,
  col = NA,
  lty = 1,
  lwd = 1,
  draw = TRUE,
  ...
)
```

## Arguments

- x, y:

  `numeric` coordinates, where x can be a two-column numeric matrix of
  x,y coordinates.

- a, b:

  `numeric` values indicating x- and y-axis radius, before rotation if
  `angle` is non-zero.

- angle:

  `numeric` value indicating the rotation of ellipse.

- segment:

  NULL or `numeric` vector of two values indicating the start and end
  angles for the ellipse, prior to rotation.

- arc.only:

  `logical` indicating whether to draw the ellipse arc without
  connecting to the center of the ellipse. Set `arc.only=FALSE` when
  segment does not include the full circle, to draw only the wedge.

- nv:

  `numeric` the number of vertices around the center to draw.

- deg:

  `logical` indicating whether input `angle` and `segment` values are in
  degrees, or `deg=FALSE` for radians.

- border, col, lty, lwd:

  arguments passed to
  [`graphics::polygon()`](https://rdrr.io/r/graphics/polygon.html).

- draw:

  `logical` indicating whether to draw the ellipse.

- ...:

  additional arguments are passed to
  [`graphics::polygon()`](https://rdrr.io/r/graphics/polygon.html) when
  `draw=TRUE`.

## Value

invisible list of x,y coordinates

## Details

This function draws an ellipse centered on the given coordinates,
rotated the given degrees relative to the center point, with give x- and
y-axis radius values.

## See also

Other jam igraph functions:
[`communities2nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/communities2nodegroups.md),
[`edge_bundle_bipartite()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_bipartite.md),
[`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md),
[`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md),
[`flip_edges()`](https://jmw86069.github.io/multienrichjam/reference/flip_edges.md),
[`get_bipartite_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_bipartite_nodeset.md),
[`highlight_edges_by_node()`](https://jmw86069.github.io/multienrichjam/reference/highlight_edges_by_node.md),
[`igraph2pieGraph()`](https://jmw86069.github.io/multienrichjam/reference/igraph2pieGraph.md),
[`label_communities()`](https://jmw86069.github.io/multienrichjam/reference/label_communities.md),
[`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md),
[`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md),
[`nodegroups2communities()`](https://jmw86069.github.io/multienrichjam/reference/nodegroups2communities.md),
[`rectifyPiegraph()`](https://jmw86069.github.io/multienrichjam/reference/rectifyPiegraph.md),
[`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md),
[`subsetCnetIgraph()`](https://jmw86069.github.io/multienrichjam/reference/subsetCnetIgraph.md),
[`subset_igraph_components()`](https://jmw86069.github.io/multienrichjam/reference/subset_igraph_components.md),
[`sync_igraph_communities()`](https://jmw86069.github.io/multienrichjam/reference/sync_igraph_communities.md)

## Examples

``` r
par("mar"=c(2, 2, 2, 2));
plot(NULL,
   type="n",
   xlim=c(-5, 20),
   ylim=c(-5, 18),
   ylab="", xlab="", bty="L",
   asp=1);
xy <- drawEllipse(
   x=c(1, 11, 11, 11),
   y=c(1, 11, 11, 11),
   a=c(5, 5, 5*1.5, 5),
   b=c(2, 2, 2*1.5, 2),
   angle=c(20, -15, -15, -15),
   segment=c(0, 360, 0, 120, 120, 240, 240, 360),
   arc.only=c(TRUE, FALSE, FALSE, TRUE),
   col=jamba::alpha2col(c("red", "gold", "dodgerblue", "darkorchid"), alpha=0.5),
   border=c("red", "gold", "dodgerblue", "darkorchid"),
   lwd=1,
   nv=99)
points(x=c(1, 11), y=c(1, 11), pch=20, cex=2)
jamba::drawLabels(x=c(12, 3, 13, 5),
   y=c(14, 10, 9, 2),
   labelCex=0.7,
   drawBox=FALSE,
   adjPreset=c("topright", "left", "bottomright", "top"),
   txt=c("0-120 degrees,\nangle=-15,\narc.only=TRUE",
      "120-240 degrees,\nangle=-15,\narc.only=TRUE,\nlarger radius",
      "240-360 degrees,\nangle=-15,\narc.only=FALSE",
      "angle=20"))

```
