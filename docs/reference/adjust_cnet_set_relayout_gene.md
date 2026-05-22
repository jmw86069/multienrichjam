# Adjust Set nodes then relayout Gene nodes

Adjust Cnet Set nodes then relayout Gene nodes

## Usage

``` r
adjust_cnet_set_relayout_gene(
  g,
  nodes = NULL,
  x = 0,
  y = 0,
  use_grep = TRUE,
  do_reorder = TRUE,
  spread_labels = TRUE,
  repulse = 4,
  verbose = FALSE,
  ...
)
```

## Arguments

- g:

  `igraph` Cnet object

- nodes:

  `character` vector of one or more Set nodes.

- x, y:

  `numeric` values passed to
  [`nudge_igraph_node()`](https://jmw86069.github.io/multienrichjam/reference/nudge_igraph_node.md).

- use_grep:

  `logical` passed to
  [`nudge_igraph_node()`](https://jmw86069.github.io/multienrichjam/reference/nudge_igraph_node.md).

- do_reorder:

  `logical` indicating whether nodes should be re-positioned within each
  subcluster by calling
  [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md).

- spread_labels:

  `logical` indicating whether labels should be re-oriented around each
  node by calling
  [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md).

- repulse:

  `numeric` value passed to
  [`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md)
  to indicate the repulsion force. Typical values range between `3` for
  loosely-packed nodes, and `5` or higher for more closely-packed nodes.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are passed to functions
  [`nudge_igraph_node()`](https://jmw86069.github.io/multienrichjam/reference/nudge_igraph_node.md),
  [`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md),
  [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md),
  [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md).

## Details

This function operates on a Cnet `igraph` object, distinguished by node
attribute `"nodeType"` with value `"Gene"` for Gene nodes, and `"Set"`
for Set nodes.

This function is intended to help move a `Set` node to improve visual
spacing between nodes, then it re-positions only the `Gene` nodes using
[`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md),
keeping the `Set` nodes in fixed positions.

It calls
[`nudge_igraph_node()`](https://jmw86069.github.io/multienrichjam/reference/nudge_igraph_node.md)
using arguments `nodes,x,y`, then defines the layout coordinates of Set
nodes as `constraints` that therefore are not allowed to change when
calling
[`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md).
Note that Gene node coordinates are allowed to change, even if Gene
nodes were included in `nodes`.

## See also

Other jam cnet utilities:
[`adjust_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_nodeset.md),
[`apply_nodeset_spacing()`](https://jmw86069.github.io/multienrichjam/reference/apply_nodeset_spacing.md),
[`bulk_cnet_adjustments()`](https://jmw86069.github.io/multienrichjam/reference/bulk_cnet_adjustments.md),
[`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md),
[`get_cnet_nodeset_vector()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset_vector.md),
[`launch_shinycat()`](https://jmw86069.github.io/multienrichjam/reference/launch_shinycat.md),
[`make_cnet_test()`](https://jmw86069.github.io/multienrichjam/reference/make_cnet_test.md),
[`relayout_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/relayout_nodegroups.md)
