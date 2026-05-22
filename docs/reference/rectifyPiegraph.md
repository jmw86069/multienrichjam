# Convert pie igraph node shapes to coloredrectangle

Convert pie igraph node shapes to coloredrectangle

## Usage

``` r
rectifyPiegraph(
  g,
  nrow = 2,
  ncol = 5,
  byrow = TRUE,
  whichNodes = seq_len(igraph::vcount(g)),
  ...
)
```

## Arguments

- g:

  igraph object, expected to contain one or more nodes with shape
  `"pie"`.

- nrow, ncol:

  integer values indicating the default number of rows and columns to
  use when displaying the colors for each node.

- byrow:

  logical indicating whether each vector of node colors should fill the
  nrow,ncol matrix by each row, similar to how values are filled in
  [`base::matrix()`](https://rdrr.io/r/base/matrix.html) with argument
  `byrow`.

- whichNodes:

  integer vector of nodes in `g` which should be considered. Only nodes
  with shape `"pie"` will be converted which are also within the
  `whichNodes` vector. By default, all nodes are converted, but
  `whichNodes` allows converting only a subset of nodes.

- ...:

  additional arguments are ignored.

## Value

igraph object where node shapes were changed from `"pie"` to
`"coloredrectangle"`.

## Details

This function simply converts an igraph network with `"pie"` node
shapes, to use the `"coloredrectangle"` node shape provided by the
multienrichjam package.

In the process, it transfers related node attributes:

- `"pie.color"` are copied to `"coloredrect.color"`

- `"pie.names"` are copied to `"coloredrect.names"`. The
  `"coloredrect.names"` can be used to label a color key.

- `"size"` is converted to `"size2"` after applying `sqrt(size) * 1.5`.
  The `"size2"` value is used to define the size of coloredrectangle
  nodes.

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
[`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md),
[`subsetCnetIgraph()`](https://jmw86069.github.io/multienrichjam/reference/subsetCnetIgraph.md),
[`subset_igraph_components()`](https://jmw86069.github.io/multienrichjam/reference/subset_igraph_components.md),
[`sync_igraph_communities()`](https://jmw86069.github.io/multienrichjam/reference/sync_igraph_communities.md)
