# Color edges by nodegroups

Color edges by nodegroups

## Usage

``` r
color_nodes_by_nodegroups(
  g,
  nodegroups = NULL,
  nodegroup_colors = NULL,
  node_alpha = NULL,
  Crange = NULL,
  Lrange = NULL,
  color_attributes = c("color"),
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

- node_alpha:

  `numeric` of default NULL uses 1, with opacity where 0 is transparent,
  1 is completely solid.

- Crange, Lrange:

  `numeric` ranges passed to
  [`colorjam::rainbowJam()`](https://jmw86069.github.io/colorjam/reference/rainbowJam.html)
  to restrict the chroma and luminance of categorical colors. Default
  NULL uses the defaults in
  [`colorjam::rainbowJam()`](https://jmw86069.github.io/colorjam/reference/rainbowJam.html).

- color_attributes:

  `character` vector of node attributes to assign new colors, default
  'color' assigns color to the node itself. Note that this color is not
  seen when shape is 'pie' or 'jampie', instead use 'pie.color'.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are passed
  to[`colorjam::rainbowJam()`](https://jmw86069.github.io/colorjam/reference/rainbowJam.html),
  notably 'preset' which defines the color wheel used.

## Value

`igraph` object with edge color attribute updated to represent the
result of blending node colors, seen by `igraph::edge_attr(g)$color`.

## See also

Other jam igraph utilities:
[`color_edges_by_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodegroups.md),
[`color_edges_by_nodes()`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodes.md),
[`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md),
[`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md),
[`launch_shinycat()`](https://jmw86069.github.io/multienrichjam/reference/launch_shinycat.md),
[`subgraph_jam()`](https://jmw86069.github.io/multienrichjam/reference/subgraph_jam.md)

## Examples

``` r
cnet1 <- make_cnet_test()
igraph::vertex_attr(cnet1, "shape") <- "circle";
cnet2 <- color_nodes_by_nodegroups(cnet1);
# behold color
jam_igraph(cnet2)


# you can enable mark.groups which uses the same colorjam::rainbowJam()
jam_igraph(cnet2, mark.groups=TRUE)

```
