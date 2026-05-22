# Subset Cnet igraph

Subset Cnet igraph

## Usage

``` r
subsetCnetIgraph(
  gCnet,
  includeSets = NULL,
  includeGenes = NULL,
  remove_singlets = TRUE,
  minSetDegree = 1,
  minGeneDegree = 1,
  remove_blanks = FALSE,
  spread_labels = TRUE,
  force_relayout = TRUE,
  do_reorder = TRUE,
  repulse = 4,
  layout = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- gCnet:

  igraph object representing Cnet concept network data

- includeSets:

  character vector, or NULL, containing the set names or labels to
  retain.

- includeGenes:

  character vector, or NULL, containing the gene names or labels to
  retain.

- remove_singlets:

  logical whether to remove singlet graph nodes, which are nodes that
  have no remaining edges. Note that argument `"removeSinglets"` is
  deprecated but will be recognized with preference, support will be
  removed in future versions.

- minSetDegree:

  integer value indicating the minimum number of edges each Set node
  must have to be retained in the resulting igraph. Use `minSetDegree=2`
  to retain only Set nodes that have multiple connected Gene nodes.

- minGeneDegree:

  integer value indicating the minimum number of edges each Gene node
  must have to be retained in the resulting igraph. Use
  `minGeneDegree=2` to retain only Gene nodes that connect to multiple
  Set nodes.

- remove_blanks:

  logical indicating whether to call
  [`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md),
  which removes blank colors from `"pie.color"` and
  `"coloredrect.color"` attributes.

- spread_labels:

  logical indicating whether to call
  [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md),
  which re-orients the node labels after the `igraph` object has been
  subsetted. When `TRUE`, the arguments `force_relayout`, and
  `do_reorder` are passed to that function.

- force_relayout:

  logical indicating whether to re-apply the `igraph` layout function.

- do_reorder:

  logical indicating whether to call
  [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md),
  which sorts equivalent nodes by color then label, intended when there
  are large numbers of nodes with the same edges, typically most common
  in a cnet plot where many gene nodes may be connected to the same
  pathway set nodes.

- layout:

  function that takes `igraph` object and returns a numeric matrix of
  node coordinates. This function is only called when
  `force_relayout=TRUE`, and must be supplied as a function in order to
  be applied properly to the subset Cnet `igraph`. To apply layout
  before the subset operation, do so with
  `igraph::set_graph_attr(g, "layout", layout)`.

- verbose:

  logical indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

## Details

This function produces a subset of a Cnet igraph based upon supplied set
names or gene names. This function is intended to be a convenient method
of filtering a Cnet igraph to a pre-defined set of "Set" names.

The function assumes graph nodes have an attribute `"nodeType"` with
values either `"Set"` or `"Gene"` to indicate the type of node.

When `includeSets` is supplied, the graph is subsetted to include only
nodes with `nodeType="Set"` with matching `V(gCnet)$name` or
`V(gCnet)$label`. Then only neighboring nodes are retained, thus
removing any nodes with `nodeType="Gene"` that do not connect to any of
the given Set nodes. The result is a proper Cnet igraph that only
contains Gene nodes connected to the subset of Set nodes.

If `includeGenes` is supplied, the graph is subsetted to include only
nodes with `nodeType="Gene"` with matching `V(gCnet)$name` or
`V(gCnet)$label`.

When `removeSinglets=TRUE` then any nodes that have no remaining edges
are removed. Especially when supplying `includeGenes`, this option is
useful to hide any Set nodes that have no connected Gene nodes.

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
[`subset_igraph_components()`](https://jmw86069.github.io/multienrichjam/reference/subset_igraph_components.md),
[`sync_igraph_communities()`](https://jmw86069.github.io/multienrichjam/reference/sync_igraph_communities.md)
