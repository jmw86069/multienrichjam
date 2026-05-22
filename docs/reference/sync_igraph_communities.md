# Sync igraph nodes and communities

Sync igraph nodes and communities

## Usage

``` r
sync_igraph_communities(g, wc, verbose = TRUE, ...)
```

## Arguments

- g:

  `igraph` object

- wc:

  `communities` object, or `list` in form of nodegroups, which is a
  `list` of `character` vectors that contain `igraph` node names.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are passed to
  [`nodegroups2communities()`](https://jmw86069.github.io/multienrichjam/reference/nodegroups2communities.md)
  only when input `wc` is supplied in `list` nodegroups format.

## Value

`list` with two elements:

- `"g"` - the `igraph` object after subsetting to match node names
  shared with `wc`, as necessary.

- `"wc'` - the `communities` object after subsetting to match node names
  shared with `g`, as necessary. When input `wc` is in `list` nodegroups
  format, that same format is returned.

## Details

This function ensures that `igraph` nodes and corresponding community
clusters are synchronized for proper downstream use. In particular, when
using a subgraph, or when communities only assign a subset of nodes to
clusters, this function ensures the two objects are in sync, the same
order, and with the same nodes.

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
[`subset_igraph_components()`](https://jmw86069.github.io/multienrichjam/reference/subset_igraph_components.md)
