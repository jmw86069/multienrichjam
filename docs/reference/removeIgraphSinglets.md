# Remove igraph singlet nodes

Remove igraph singlet nodes

## Usage

``` r
removeIgraphSinglets(g, min_degree = 1, ...)
```

## Arguments

- g:

  igraph object

- min_degree:

  numeric threshold with the minimum number of connections, also known
  as the "degree", required for each node.

- ...:

  additional arguments are ignored.

## Details

This function is a lightweight method to remove igraph nodes with no
connections. In fact, the `min_degree` can be used to require a minimum
number of connections, but the intended use is to remove the singlet
nodes that have no connections.

## See also

Other jam igraph layouts:
[`layout_communities()`](https://jmw86069.github.io/multienrichjam/reference/layout_communities.md),
[`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md),
[`layout_with_qfrf()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfrf.md),
[`relayout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/relayout_with_qfr.md),
[`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md),
[`rotate_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/rotate_igraph_layout.md),
[`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)
