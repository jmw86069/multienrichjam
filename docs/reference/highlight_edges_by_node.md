# Highlight edges connected to a node or nodes

Highlight edges connected to a node or nodes

## Usage

``` r
highlight_edges_by_node(
  g,
  node = NULL,
  highlight_color = "royalblue3",
  highlight_width = NULL,
  highlight_cex = 2,
  nonhighlight_alpha = NULL,
  ...
)
```

## Arguments

- g:

  `igraph` object

- node:

  `character` node name, or `integer` node index.

- highlight_color:

  `character` color used for highlighting, default 'royalblue'.

- highlight_width:

  `numeric` default NULL, set a fixed edge width.

- highlight_cex:

  `numeric` default 2, multiplies the edge width by this expansion
  factor.

- nonhighlight_alpha:

  `numeric` default NULL, applies an optional alpha transparency to
  non-highlighted edges. Values must be between 0 (transparent) and 1
  (opaque), where a value 0.5 is recommended when used.

- ...:

  additional arguments are ignored.

## Value

`igraph` object after adjusting edge attributes.

## Details

This function highlights edges connected to one of more nodes by
applying a color, and adjusting the edge width.

Note that when edge attributes 'width' and 'color' are not yet defined
in the `igraph` object, they will be defined and stored into the object
by using
[`default_igraph_values()`](https://jmw86069.github.io/multienrichjam/reference/default_igraph_values.md).
For example: `default_igraph_values()$edge$width` and
`default_igraph_values()$edge$color`.

## See also

Other jam igraph functions:
[`communities2nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/communities2nodegroups.md),
[`drawEllipse()`](https://jmw86069.github.io/multienrichjam/reference/drawEllipse.md),
[`edge_bundle_bipartite()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_bipartite.md),
[`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md),
[`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md),
[`flip_edges()`](https://jmw86069.github.io/multienrichjam/reference/flip_edges.md),
[`get_bipartite_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_bipartite_nodeset.md),
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

## Examples

``` r
cnet1 <- make_cnet_test(seed=1);
# highlight edges by node
cnet1h <- highlight_edges_by_node(cnet1, "CA", nonhighlight_alpha=0.5)
# plot the original
jam_igraph(cnet1)

# plot the highlighted variant
jam_igraph(cnet1h, label_dist_factor=3)


# highlight edges using nodeset
nodesets <- get_cnet_nodeset(cnet1)
cnet1hns <- highlight_edges_by_node(cnet1,
   nodesets[[grep("CA", nodesets)]],
   nonhighlight_alpha=0.5)
jam_igraph(cnet1hns, node_factor=2, label_factor=2, label_dist_factor=3,
   mark.groups=unname(nodesets[grep("CA", nodesets)]))

```
