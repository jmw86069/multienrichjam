# cnetplot internal function

cnetplot internal function

## Usage

``` r
cnetplot_internalJam(
  inputList,
  categorySize = "geneNum",
  showCategory = 5,
  pvalue = NULL,
  foldChange = NULL,
  fixed = TRUE,
  DE.foldChange = NULL,
  categoryColor = "#E5C494",
  geneColor = "#B3B3B3",
  colorRamp = "RdBu_r",
  ...
)
```

## Details

This function is intended to mimic the `DOSE:::cnetplot_internal()`
function to support
[`cnetplotJam()`](https://jmw86069.github.io/multienrichjam/reference/cnetplotJam.md)
customizations, including not plotting the output, and including
additional custom igraph attributes.

## See also

Other jam igraph functions:
[`cnet2df`](https://jmw86069.github.io/multienrichjam/reference/cnet2df.md)`()`,
[`cnet2im`](https://jmw86069.github.io/multienrichjam/reference/cnet2im.md)`()`,
[`cnetplotJam`](https://jmw86069.github.io/multienrichjam/reference/cnetplotJam.md)`()`,
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
