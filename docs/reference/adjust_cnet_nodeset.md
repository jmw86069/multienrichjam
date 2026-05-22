# Adjust Cnet node set

Adjust Cnet node set

## Usage

``` r
adjust_cnet_nodeset(
  g,
  set_nodes = NULL,
  x = 0,
  y = 0,
  percent_spacing = NULL,
  expand = 0,
  rotate_degrees = 0,
  cnet_nodesets = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- g:

  `igraph` Cnet plot, with vertex attribute `nodeType` that contains
  values `"Gene"` and `"Set"`, as produced by
  [`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md).

- set_nodes:

  `character` vector of one or more set names, which are defined here as
  the names of the nodes with `nodeType="Set"` to which each node with
  `nodeType="Gene"` is connected. For example `set_nodes=c("A", "B")`
  will return all Gene nodes that are connected only to Set nodes `"A"`
  and `"B"`. Alternatively, one can refer to a cluster by the node name
  of one member of the cluster, as a convenience option.

- x, y:

  `numeric` indicating the relative amount to move nodes in the node
  set, based upon the overall range of x and y values in the layout
  coordinates.

- percent_spacing:

  `numeric` within node spacing as percentage of the max x/y layout
  ranges. This option is preferred to `expand` when trying to keep
  nodeset spacing consistent for the overall layout. As a general rule,
  use `summarize_node_spacing(g)$nearest_within` with column 'Median' to
  determine node spacing for each nodeset, then adjust each nodeset to
  relatively consistent median spacing.

- expand:

  `numeric` the relative amount to expand node coordinates, where
  `expand=0` does no expansion, `expand=1` will expand coordinates 100%
  of the distance to the center of the cluster, twice the original
  distance; and `expand=-1` will expand coordinates to half the distance
  to the center.

- rotate_degrees:

  `numeric` value indicating the rotation in degrees to adjust the node
  coordinates, relative to the center of the coordinates.

- cnet_nodesets:

  `list` optional to provide the known nodesets upfront, in order to
  avoid re-determining the nodesets. Otherwise `get_cnet_nodesets()` is
  called, which for large networks may introduce a small time lag.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

## Details

This function is intended to help adjust a node set, defined as a set of
`"Gene"` nodes that all connect to the same `"Set"` nodes. These nodes
are usually clustered together in a Cnet plot, and can be manipulated as
a group using this function. Use `get_cnet_nodeset(g)` to get all node
sets from a Cnet `igraph` object.

Specifically, this function can be used to expand, rotate, and move the
cluster of nodes in a Cnet network layout. The purpose is usually to
help the visual clarity of the plot, and to reduce node label overlaps.

## See also

Other jam cnet utilities:
[`adjust_cnet_set_relayout_gene()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_set_relayout_gene.md),
[`apply_nodeset_spacing()`](https://jmw86069.github.io/multienrichjam/reference/apply_nodeset_spacing.md),
[`bulk_cnet_adjustments()`](https://jmw86069.github.io/multienrichjam/reference/bulk_cnet_adjustments.md),
[`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md),
[`get_cnet_nodeset_vector()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset_vector.md),
[`launch_shinycat()`](https://jmw86069.github.io/multienrichjam/reference/launch_shinycat.md),
[`make_cnet_test()`](https://jmw86069.github.io/multienrichjam/reference/make_cnet_test.md),
[`relayout_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/relayout_nodegroups.md)

## Examples

``` r
# the examples below are for internal data
# and will be rewritten to work with example data
if (1 == 2) {
g <- mem_dmjdm_gp_20_plots_jikl$cnet_collapsed_set;
set_nodes <- c("J", "L")
g2 <- adjust_cnet_nodeset(g, set_nodes=c("J", "L"), expand=0.5)

jam_igraph(g)
g2 <- adjust_cnet_nodeset(g, set_nodes=c("I", "K", "L"), x=0.01, y=-0.01, expand=0.2)
g2 <- adjust_cnet_nodeset(g2, set_nodes=c("I", "J", "K", "L"), x=-0.08, y=-0.05, expand=0.2)
g2 <- adjust_cnet_nodeset(g2, set_nodes=c("I", "L"), x=-0.05, y=0.05, expand=0.3)
g2 <- adjust_cnet_nodeset(g2, set_nodes=c("I", "J", "L"), x=-0.02, y=0.03, expand=0.3)
g2 <- adjust_cnet_nodeset(g2, set_nodes=c("I", "K"), x=0.03, y=0.09, expand=0.)
g2 <- adjust_cnet_nodeset(g2, set_nodes=c("I", "J", "K"), x=0.02, y=-0.01, expand=0.)
g2 <- adjust_cnet_nodeset(g2, set_nodes=c("I", "J"), x=0.01, y=-0.00, expand=0.1)
g2 <- adjust_cnet_nodeset(g2, set_nodes=c("J"), x=-0.03, y=-0.02, expand=0.)
jam_igraph(g2)
}
```
