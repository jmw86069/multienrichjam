# Convert igraph to use pie node shapes

Convert igraph to use pie node shapes

## Usage

``` r
igraph2pieGraph(
  g,
  valueIM = NULL,
  valueIMcolors = NULL,
  colorV = NULL,
  updateLabels = FALSE,
  maxNchar = 62,
  backgroundColor = "white",
  seed = 123,
  defineLayout = FALSE,
  repulse = 3.6,
  removeNA = FALSE,
  NAvalues = c(NA, "transparent"),
  verbose = FALSE,
  ...
)
```

## Details

This function converts an igraph to use pie node shapes, where pie
wedges are colored using values derived from a numeric matrix `valueIM`
or pre-defined in a character matrix containing colors `valueIMcolors`.

Note that pie wedge sizes are equally-sized and do not vary by score,
instead the color intensity is applied to each pie wedge.

Node names using `V(g)$name` matching `rownames(valueIMcolors)` are
colorized and the node shape is converted to pie. All other nodes are
not modified.

When `valueIMcolors` is not defined, it is derived from `valueIM` using
[`colorjam::matrix2heatColors()`](https://jmw86069.github.io/colorjam/reference/matrix2heatColors.html).
In that case, `colorV` defines the color used for numeric values in each
column, and other options are passed to
[`colorjam::matrix2heatColors()`](https://jmw86069.github.io/colorjam/reference/matrix2heatColors.html)
via `...` arguments.

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
[`label_communities()`](https://jmw86069.github.io/multienrichjam/reference/label_communities.md),
[`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md),
[`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md),
[`nodegroups2communities()`](https://jmw86069.github.io/multienrichjam/reference/nodegroups2communities.md),
[`rectifyPiegraph()`](https://jmw86069.github.io/multienrichjam/reference/rectifyPiegraph.md),
[`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md),
[`subsetCnetIgraph()`](https://jmw86069.github.io/multienrichjam/reference/subsetCnetIgraph.md),
[`subset_igraph_components()`](https://jmw86069.github.io/multienrichjam/reference/subset_igraph_components.md),
[`sync_igraph_communities()`](https://jmw86069.github.io/multienrichjam/reference/sync_igraph_communities.md)
