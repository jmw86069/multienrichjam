# Subset igraph by connected components

Subset igraph by connected components

## Usage

``` r
subset_igraph_components(
  g,
  keep = NULL,
  min_size = 1,
  order_by_size = TRUE,
  ...
)
```

## Arguments

- g:

  igraph object

- keep:

  numeric vector indicating which component or components to keep in the
  final output. When `order_by_size=TRUE`, components are ordered by
  size, from largest to smallest, in that case `keep=1` will return only
  the one largest connected subgraph.

- min_size:

  numeric value indicating the number of nodes required in all connected
  components returned. This filter is applied after the `keep` argument.

- order_by_size:

  logical indicating whether the connected components are sorted by
  size, largest to smallest, and therefore re-numbered. Otherwise, the
  components are somewhat randomly labeled based upon the output of
  [`igraph::components()`](https://r.igraph.org/reference/components.html).

- ...:

  additional arguments are passed to
  [`igraph::components()`](https://r.igraph.org/reference/components.html).

## Details

This function is intended to help drill down into an igraph object that
contains multiple connected components.

By default, it sorts the components from largest number of nodes, to
smallest, which helps choose the largest connected component, or
subsequent components in size order.

The components can also be filtered to require a minimum number of
connected nodes.

At its core, this function is a wrapper to
[`igraph::components()`](https://r.igraph.org/reference/components.html)
and
[`igraph::subgraph()`](https://r.igraph.org/reference/subgraph.html).

## See also

Other jam igraph functions:
[`communities2nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/communities2nodegroups.md),
[`drawEllipse()`](https://jmw86069.github.io/multienrichjam/reference/drawEllipse.md),
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
[`sync_igraph_communities()`](https://jmw86069.github.io/multienrichjam/reference/sync_igraph_communities.md)
