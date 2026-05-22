# Subgraph using Jam extended logic

Subgraph using Jam extended logic

## Usage

``` r
subgraph_jam(graph, v)
```

## Arguments

- graph:

  `igraph` object

- v:

  `integer` or `logical` vector indicating the nodes to retain in the
  final `igraph` object.

## Details

This function extends the
[`igraph::subgraph()`](https://r.igraph.org/reference/subgraph.html)
function to include proper subset of the graph attribute `"layout"`,
which for some unknown reason does not subset the layout matrix
consistent with the subset of `igraph` nodes.

## See also

Other jam igraph utilities:
[`color_edges_by_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodegroups.md),
[`color_edges_by_nodes()`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodes.md),
[`color_nodes_by_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/color_nodes_by_nodegroups.md),
[`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md),
[`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md),
[`launch_shinycat()`](https://jmw86069.github.io/multienrichjam/reference/launch_shinycat.md)
