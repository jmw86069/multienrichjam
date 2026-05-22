# Color igraph edges using node colors (deprecated)

Color igraph edges using node colors (deprecated)

## Usage

``` r
color_edges_by_nodes_deprecated(g, alpha = NULL, ...)
```

## Arguments

- g:

  `igraph` object

- alpha:

  `NULL` or numeric vector with value between 0 and 1, where 0 is
  transparent and 1 is non-transparent. When supplied, this value is
  passed to
  [`jamba::alpha2col()`](https://jmw86069.github.io/jamba/reference/alpha2col.html)
  to apply alpha transparency to each edge color.

- ...:

  additional arguments are ignored.

## Details

Note: This function is deprecated in favor of
[`color_edges_by_nodes()`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodes.md).

This function uses the average color for the two nodes involved in each
edge, and applies that as the new edge color.

The color for each node depends upon the node shape, where shape `"pie"`
uses the average color from `"pie.color"`, and shape
`"coloredrectangle"` uses the avereage color from `"coloredrect.color"`.
Everything else uses `"color"`.

This function relies upon
[`avg_colors_by_list()`](https://jmw86069.github.io/multienrichjam/reference/avg_colors_by_list.md)
to blend multiple colors together.

## See also

Other jam igraph functions:
[`cnet2df`](https://jmw86069.github.io/multienrichjam/reference/cnet2df.md)`()`,
[`cnet2im`](https://jmw86069.github.io/multienrichjam/reference/cnet2im.md)`()`,
[`cnetplotJam`](https://jmw86069.github.io/multienrichjam/reference/cnetplotJam.md)`()`,
[`cnetplot_internalJam`](https://jmw86069.github.io/multienrichjam/reference/cnetplot_internalJam.md)`()`,
[`color_edges_by_nodegroups`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodegroups.md)`()`,
[`color_edges_by_nodes`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodes.md)`()`,
[`color_nodes_by_nodegroups`](https://jmw86069.github.io/multienrichjam/reference/color_nodes_by_nodegroups.md)`()`,
[`communities2nodegroups`](https://jmw86069.github.io/multienrichjam/reference/communities2nodegroups.md)`()`,
[`drawEllipse`](https://jmw86069.github.io/multienrichjam/reference/drawEllipse.md)`()`,
[`edge_bundle_bipartite`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_bipartite.md)`()`,
[`edge_bundle_nodegroups`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md)`()`,
[`enrichMapJam`](https://jmw86069.github.io/multienrichjam/reference/enrichMapJam.md)`()`,
[`fixSetLabels`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md)`()`,
[`flip_edges`](https://jmw86069.github.io/multienrichjam/reference/flip_edges.md)`()`,
[`get_bipartite_nodeset`](https://jmw86069.github.io/multienrichjam/reference/get_bipartite_nodeset.md)`()`,
[`highlight_edges_by_node`](https://jmw86069.github.io/multienrichjam/reference/highlight_edges_by_node.md)`()`,
[`igraph2pieGraph`](https://jmw86069.github.io/multienrichjam/reference/igraph2pieGraph.md)`()`,
[`jam_igraph`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)`()`,
[`jam_plot_igraph`](https://jmw86069.github.io/multienrichjam/reference/jam_plot_igraph.md)`()`,
[`label_communities`](https://jmw86069.github.io/multienrichjam/reference/label_communities.md)`()`,
[`layout_with_qfr`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md)`()`,
[`layout_with_qfrf`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfrf.md)`()`,
[`mem2emap`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md)`()`,
[`memIM2cnet`](https://jmw86069.github.io/multienrichjam/reference/memIM2cnet.md)`()`,
[`mem_multienrichplot`](https://jmw86069.github.io/multienrichjam/reference/mem_multienrichplot.md)`()`,
[`nodegroups2communities`](https://jmw86069.github.io/multienrichjam/reference/nodegroups2communities.md)`()`,
[`rectifyPiegraph`](https://jmw86069.github.io/multienrichjam/reference/rectifyPiegraph.md)`()`,
[`relayout_with_qfr`](https://jmw86069.github.io/multienrichjam/reference/relayout_with_qfr.md)`()`,
[`removeIgraphBlanks`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md)`()`,
[`removeIgraphSinglets`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphSinglets.md)`()`,
[`reorderIgraphNodes`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)`()`,
[`rotate_igraph_layout`](https://jmw86069.github.io/multienrichjam/reference/rotate_igraph_layout.md)`()`,
[`spread_igraph_labels`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)`()`,
[`subgraph_jam`](https://jmw86069.github.io/multienrichjam/reference/subgraph_jam.md)`()`,
[`subsetCnetIgraph`](https://jmw86069.github.io/multienrichjam/reference/subsetCnetIgraph.md)`()`,
[`subset_igraph_components`](https://jmw86069.github.io/multienrichjam/reference/subset_igraph_components.md)`()`,
[`sync_igraph_communities`](https://jmw86069.github.io/multienrichjam/reference/sync_igraph_communities.md)`()`,
[`with_qfr`](https://jmw86069.github.io/multienrichjam/reference/with_qfr.md)`()`
