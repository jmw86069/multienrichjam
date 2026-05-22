# Colorize igraph edges by nodes

Colorize igraph edges using node colors

## Usage

``` r
color_edges_by_nodes(
  g,
  edge_alpha = NULL,
  Crange = c(0, 100),
  Lrange = c(0, 65),
  ...
)
```

## Arguments

- g:

  `igraph` object that contains vertex node attribute `"color"` as seen
  with `igraph::vertex_attr(g, "color")`.

- edge_alpha:

  `numeric` or `NULL`, where numeric value sets the edge alpha
  transparency, where `edge_alpha=0` is completely transparent,
  `edge_alpha=0.5` is 50% transparent, and `edge_alpha=1` is completely
  not transparent, and is opaque. When `edge_alpha=NULL` the alpha
  values are supplied by
  [`colorjam::blend_colors()`](https://jmw86069.github.io/colorjam/reference/blend_colors.html)
  which blends the two values.

- ...:

  additional arguments are passed to
  [`colorjam::blend_colors()`](https://jmw86069.github.io/colorjam/reference/blend_colors.html).

## Value

`igraph` object with edge color attribute updated to represent the
result of blending node colors, seen by `igraph::edge_attr(g)$color`.

## Details

This function colorizes edges by blending colors for the nodes involved,
by calling
[`colorjam::blend_colors()`](https://jmw86069.github.io/colorjam/reference/blend_colors.html).

The color for each node depends upon the node shape, so the color or
colors used to render each node shape will be used for the edge. For
example:

- `shape="pie"` uses the average color from `V(g)$pie.color`

- `shape="coloredrectangle"` uses the avereage color from
  `V(g)$coloredrect.color`

- everything else uses `V(g)$color`

## See also

Other jam igraph utilities:
[`color_edges_by_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodegroups.md),
[`color_nodes_by_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/color_nodes_by_nodegroups.md),
[`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md),
[`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md),
[`launch_shinycat()`](https://jmw86069.github.io/multienrichjam/reference/launch_shinycat.md),
[`subgraph_jam()`](https://jmw86069.github.io/multienrichjam/reference/subgraph_jam.md)
