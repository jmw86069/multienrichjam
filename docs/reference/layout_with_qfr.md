# igraph layout using qgraph Fruchterman-Reingold

igraph layout using qgraph Fruchterman-Reingold

## Usage

``` r
layout_with_qfr(
  g,
  repulse = 3.5,
  area = 8 * (igraph::vcount(g)^2),
  repulse.rad = (igraph::vcount(g)^repulse),
  constraints = NULL,
  constrain = NULL,
  seed = 123,
  weights = NULL,
  niter = NULL,
  max.delta = NULL,
  cool.exp = NULL,
  init = NULL,
  groups = NULL,
  rotation = NULL,
  layout.control = 0.5,
  round = TRUE,
  digits = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- g:

  igraph object

- repulse:

  exponent power used to scale the radius effect around each vertex. The
  default is slightly higher than the cube of the number of vertices,
  but as the number of vertices increases, values from 3.5 to 4 and
  higher are more effective for layout.

- area:

  The area of the plot, default is the square of the number of vertices
  times 8. Changes to plot area will also affect values of `repulse` and
  `repulse.rad`.

- repulse.rad:

  Repulse radius, defaults to the the number of vertices raised to the
  `repulse` power.

- constraints:

  optional two-column matrix with the coordinates of nodes which should
  not be modified, and NA values for nodes where the position can be
  modified.

- constrain:

  `character` optional vector of node names that should be constrained.
  This argument is a convenient shortcut for defining `constraints`,
  which is a layout coordinate matrix with `NA` values on each row where
  the coordinate is free to move, and `numeric` values where the
  coordinate is fixed. For graph `g` that contains layout in
  `igraph::graph_attr(g, "layout")`, the `init` can be defined with this
  layout, then `constraints` is defined using `constrain`.

- ...:

  other arguments are sent to
  [`qgraph::qgraph.layout.fruchtermanreingold()`](https://rdrr.io/pkg/qgraph/man/qgraph.layout.fruchtermanreingold.html)

## Value

two-column numeric matrix with coordinates for each vertex.

## Details

This function provides Fruchterman-Reingold layout for an igraph object
using the implementation from the qgraph package, which provides
important configuration options deprecated in the igraph implementation.
Notably, the repulse.rad parameter is helpful in adjusting the relative
spacing of vertices, where higher values cause tighter packing of
vertices, and lower values allows greater spacing between vertices.

## See also

[`qgraph::qgraph.layout.fruchtermanreingold()`](https://rdrr.io/pkg/qgraph/man/qgraph.layout.fruchtermanreingold.html)

Other jam igraph layouts:
[`layout_communities()`](https://jmw86069.github.io/multienrichjam/reference/layout_communities.md),
[`layout_with_qfrf()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfrf.md),
[`relayout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/relayout_with_qfr.md),
[`removeIgraphSinglets()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphSinglets.md),
[`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md),
[`rotate_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/rotate_igraph_layout.md),
[`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)

## Examples

``` r
if (suppressPackageStartupMessages(require(igraph))) {
g  <- make_graph( ~ A-B-C-D-A, E-A:B:C:D,
   F-G-H-I-F, J-F:G:H:I,
   K-L-M-N-K, O-K:L:M:N,
   P-Q-R-S-P, T-P:Q:R:S,
   B-F, E-J, C-I, L-T, O-T, M-S,
   C-P, C-L, I-L, I-P)

par("mfrow"=c(2,2));
plot(g, main="default layout\n(igraph)");

plot(g, main="layout_with_fr\n(igraph)", layout=layout_with_fr);

plot(g, main="layout_with_qfr\n(qgraph default)", layout=layout_with_qfr);

plot(g, main="layout_with_qfr, repulse=6\n(qgraph custom)",
   layout=function(g)layout_with_qfr(g, repulse=6));
}

```
