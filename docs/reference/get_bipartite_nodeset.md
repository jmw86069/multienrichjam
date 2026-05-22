# Get partite/connected graph nodesets

Get partite/connected graph nodesets defined by shared connections

## Usage

``` r
get_bipartite_nodeset(
  g,
  type = "nodeType",
  set_nodes = NULL,
  sep = ",",
  return_type = c("list", "df"),
  verbose = FALSE,
  ...
)
```

## Arguments

- g:

  `igraph` object that contains one attribute column with node type.

- type:

  `character` string of the node/vertex attribute that represents the
  node type.

- set_nodes:

  `character` or `NULL`, which contains the set of node neighbors for
  the requested nodeset. For example, one might want all nodes that
  connect with `c("A", "B", "C")`. When `set_nodes=NULL` then all
  nodesets are returned.

- sep:

  `character` string used as a delimiter between node names when
  defining a nodeset name

- return_type:

  `character` string indicating the type of data to return:

  - `"list"` returns a list of nodesets, each element in the `list` is a
    `character` vector with node names.

  - `"df"` returns a `data.frame` with more detailed annotation for each
    node, including nodesets, neighbor nodes, etc.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

## Details

This method is under development, the intent is to bundle edges where a
large subset of nodes are all connected to the same node neighbors. A
typical graph may not have any two nodes with the same neighbors, but
this situation tends to happen much more often with bipartite graphs,
where nodes of one type are only permitted to have node neighbors of the
other type. It is not required for this method to work, however.

The driving scenario is with Cnet (concept network) plots, which is a
bipartite network with `"Gene"` and `"Set"` nodes. It is fairly common
to have multiple genes present in the same one or few pathways. As a
result, these nodes are most often positioned near each other as a
natural by-product of having the same connected neighbor nodes.

Identifying a nodeset with identical node neighbors enables some other
useful operations:

- re-positioning, rotating, expanding, compressing the whole nodeset
  layout to improve network graph aesthetics, node label readability,
  reducing overlaps

- edge bundling to improve visual distinctiveness between multiple
  nodesets

## See also

Other jam igraph functions:
[`communities2nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/communities2nodegroups.md),
[`drawEllipse()`](https://jmw86069.github.io/multienrichjam/reference/drawEllipse.md),
[`edge_bundle_bipartite()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_bipartite.md),
[`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md),
[`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md),
[`flip_edges()`](https://jmw86069.github.io/multienrichjam/reference/flip_edges.md),
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
