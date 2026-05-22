# Color edges by nodegroups

Color edges by nodegroups

## Usage

``` r
color_edges_by_nodegroups(
  g,
  nodegroups = NULL,
  nodegroup_colors = NULL,
  edge_alpha = NULL,
  Crange = NULL,
  Lrange = NULL,
  darkFactor = 1,
  sFactor = 1,
  naColor = "grey",
  verbose = FALSE,
  ...
)
```

## Arguments

- g:

  `igraph` object that contains vertex node attribute `"color"` as seen
  with `igraph::vertex_attr(g, "color")`.

- nodegroups:

  `list` or `communities` object that references nodes in `g` and
  assigns one or more to nodegroups. Default NULL calls
  [`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md),
  passing '...' for customization.

- nodegroup_colors:

  `character` or default NULL which assigns categorical colors by
  [`colorjam::rainbowJam()`](https://jmw86069.github.io/colorjam/reference/rainbowJam.html).

- edge_alpha:

  `numeric` of default NULL uses 1, with opacity where 0 is transparent,
  1 is completely solid.

- Crange, Lrange:

  `numeric` ranges passed to
  [`colorjam::rainbowJam()`](https://jmw86069.github.io/colorjam/reference/rainbowJam.html)
  to restrict the chroma and luminance of categorical colors. Default
  NULL uses the defaults in
  [`colorjam::rainbowJam()`](https://jmw86069.github.io/colorjam/reference/rainbowJam.html).

- darkFactor, sFactor:

  `numeric` default 1 each, used to darken and saturate colors,
  respectively. When either contain values other than 1, it calls
  [`jamba::makeColorDarker()`](https://jmw86069.github.io/jamba/reference/makeColorDarker.html).

- naColor:

  `character` color, default 'grey', used for edges that for any reason
  were not assigned to a nodegroup, to prevent them from becoming NA.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are passed to
  [`colorjam::blend_colors()`](https://jmw86069.github.io/colorjam/reference/blend_colors.html).

## Value

`igraph` object with edge color attribute updated to represent the
result of blending node colors, seen by `igraph::edge_attr(g)$color`.

## Details

This function is analogous to
[`color_nodes_by_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/color_nodes_by_nodegroups.md),
and in fact uses the same nodegroup color strategy. However, by default
'Set' nodes are not colored, and therefore edges are only colored using
the nodegroup.

Using `filter_set_only=FALSE`, which is passed to
[`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md),
will assign 'Set' nodes to their own nodegroups, giving them colors. At
that point, edges will blend colors for the 'Set' and 'Gene' nodes.

## See also

Other jam igraph utilities:
[`color_edges_by_nodes()`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodes.md),
[`color_nodes_by_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/color_nodes_by_nodegroups.md),
[`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md),
[`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md),
[`launch_shinycat()`](https://jmw86069.github.io/multienrichjam/reference/launch_shinycat.md),
[`subgraph_jam()`](https://jmw86069.github.io/multienrichjam/reference/subgraph_jam.md)

## Examples

``` r
cnet1 <- make_cnet_test()
igraph::vertex_attr(cnet1, "shape") <- "circle";
cnet2 <- color_nodes_by_nodegroups(cnet1)
cnet3 <- color_edges_by_nodegroups(cnet2)
jam_igraph(cnet3)


cnet2b <- color_nodes_by_nodegroups(cnet1, filter_set_only=FALSE)
cnet3b <- color_edges_by_nodegroups(cnet2b, filter_set_only=FALSE)
jam_igraph(cnet3b)

```
