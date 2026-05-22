# Spread igraph node labels by angle from node center

Spread igraph node labels by angle from node center

## Usage

``` r
spread_igraph_labels(
  g,
  layout = NULL,
  y_bias = 1,
  update_g_coords = TRUE,
  do_reorder = TRUE,
  sortAttributes = NULL,
  nodeSortBy = c("x", "-y"),
  repulse = 3.5,
  force_relayout = FALSE,
  label_min_dist = 0.5,
  label_dist_l = list(nodeType = c(Gene = 2.5, Set = 0)),
  verbose = FALSE,
  ...
)
```

## Arguments

- g:

  igraph object

- layout:

  numeric matrix representing the x and y coordinates of each node in
  `g`, in the same order as `V(g)`. When `layout` is not supplied, nodes
  are checked for attributes `c("x", "y")` which define a fixed internal
  layout. When `force_layout=TRUE` these coordinates are ignored. If
  that is not supplied, then
  [`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md)
  is called along with the `repulse` argument. Subsequent coordinates
  are stored in `V(g)$x` and `V(g)$y` when argument
  `update_g_coords=TRUE`.

- y_bias:

  numeric value indicating the tendency to spread labels on the y-axis
  rather than symmetrically around each node. This argument elongates
  the circle surrounding a node into an ellipse with this ratio.

- update_g_coords:

  logical indicating whether the layout coordinates will be stored in
  `graph_attr(g, "layout")`.

- do_reorder:

  logical indicating whether to call
  [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)
  which re-distributes equivalent nodes based upon the node color(s). A
  node is "equivalent" to another node if both nodes have identical
  edges.

- sortAttributes, nodeSortBy:

  arguments passed to
  [`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md)
  when `do_reorder=TRUE`.

- repulse:

  argument passed to
  [`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md)
  only when `layout` is not supplied, and the layout is not stored in
  `c(V(g)$x, V(g)$y)`.

- force_relayout:

  logical indicating whether the `igraph` layout should be recalculated,
  in order to override coordinates that may be previously stored in the
  `igraph` object itself. Note that when `layout` is supplied, it is
  always used.

- label_min_dist:

  `numeric` minimum distance, default 0.5, used to ensure all labels are
  at least some distance from the center. Note that this argument is not
  used when `label_dist_l` is non-NULL. The units are roughly in units
  of one line height of text.

- label_dist_l:

  `list` with default to apply distance by 'nodeType' node vertex
  attribute: Gene distance 5, Set distance 0. When supplied as NULL, or
  when there is no node matched by this argument, the label distance
  uses `label_min_dist`.

- ...:

  additional arguments are passed to
  [`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md)
  when needed.

## Details

This function uses the igraph vertex attribute `"label.degree"`, which
describes the angular offset for each vertex label. (The
`"label.degree"` values are represented as radians, not degrees,
starting at 0 for right, and proceeding clockwise starting from the
right, down, left, top, right.)

This function requires a network layout, which must be fixed in order
for the vertex labels to be properly oriented. Labels are oriented
opposite the most dominant angular mean of edges from each network node.
Typically the side of a node with the fewest edges has the most space to
place a label. No further checks are performed for overlapping labels.

Note that this function only modifies the other important attribute
`"label.dist"` when \`label_min_dist“ is defined, in order to enforce a
minimum label distance from the center of each node. There is no other
logic to position small or large labels to avoid overlapping labels.

## See also

Other jam igraph layouts:
[`layout_communities()`](https://jmw86069.github.io/multienrichjam/reference/layout_communities.md),
[`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md),
[`layout_with_qfrf()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfrf.md),
[`relayout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/relayout_with_qfr.md),
[`removeIgraphSinglets()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphSinglets.md),
[`reorderIgraphNodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md),
[`rotate_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/rotate_igraph_layout.md)
