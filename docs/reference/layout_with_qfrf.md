# igraph layout function using qgraph Fruchterman-Reingold

igraph layout function using qgraph Fruchterman-Reingold

## Usage

``` r
layout_with_qfrf(repulse = 3.5, seed = 123, ...)
```

## Arguments

- repulse:

  numeric value typically between 3 and 5, passed to
  [`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md),
  which in turn is passed to
  [`qgraph::qgraph.layout.fruchtermanreingold()`](https://rdrr.io/pkg/qgraph/man/qgraph.layout.fruchtermanreingold.html).

- seed:

  numeric value used to set the R random seed, in order to make layouts
  consistent.

- ...:

  additional arguments are passed to
  [`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md).

## Value

function used to calculate layout coordinates of an `igraph` object.

## Details

This function returns a layout function, which can be convenient when
calling
[`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html),
in order to set layout parameters in the same call.

## See also

Other jam igraph layouts:
[`layout_communities()`](https://jmw86069.github.io/multienrichjam/reference/layout_communities.md),
[`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md),
[`relayout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/relayout_with_qfr.md),
[`removeIgraphSinglets()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphSinglets.md),
[`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md),
[`rotate_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/rotate_igraph_layout.md),
[`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)
