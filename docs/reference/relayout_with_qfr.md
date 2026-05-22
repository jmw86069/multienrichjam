# igraph re-layout using qgraph Fruchterman-Reingold

igraph re-layout using qgraph Fruchterman-Reingold

## Usage

``` r
relayout_with_qfr(
  g,
  repulse = 3.5,
  spread_labels = TRUE,
  seed = 123,
  init = NULL,
  constrain = NULL,
  constraints = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- g:

  `igraph` object

- repulse:

  exponent power used to scale the radius effect around each vertex. The
  default is slightly higher than the cube of the number of vertices,
  but as the number of vertices increases, values from 3.5 to 4 and
  higher are more effective for layout.

- spread_labels:

  logical indicating whether to call
  [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md),
  which places node labels at an angle offset from the node, in order to
  improve default label positions.

- ...:

  additional arguments are passed to
  [`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md)
  and
  [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)
  as needed.

## Value

`igraph` object, with layout coordinates stored in graph attribute
`"layout"`, accessible for example with `graph$layout` or
`graph_attr(graph, "layout")`. When `spread_labels=TRUE`,
`V(g)$label.degree` and `V(g)$label.dist` are updated by calling
[`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md).

## Details

This function extends
[`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md)
by applying the layout to the `igraph` object itself, while also calling
[`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)
to adjust label positions accordingly.

The main benefit to using this function is to update the layout and node
label positions in one step, while also returning the `igraph` object
ready to be plotted as-is.

## See also

Other jam igraph layouts:
[`layout_communities()`](https://jmw86069.github.io/multienrichjam/reference/layout_communities.md),
[`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md),
[`layout_with_qfrf()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfrf.md),
[`removeIgraphSinglets()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphSinglets.md),
[`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md),
[`rotate_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/rotate_igraph_layout.md),
[`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)
