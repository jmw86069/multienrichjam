# Jam igraph vectorized plot function (internal)

Jam igraph vectorized plot internal function called by
[`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)

## Usage

``` r
jam_plot_igraph(
  x,
  ...,
  axes = FALSE,
  add = FALSE,
  xlim = NULL,
  ylim = NULL,
  mark.groups = list(),
  mark.shape = 1/2,
  mark.col = NULL,
  mark.alpha = 0.2,
  mark.border = NULL,
  mark.expand = NULL,
  mark.lwd = 2,
  mark.lty = 1,
  mark.smooth = TRUE,
  mark.cex = 1,
  mark.x.nudge = 0,
  mark.y.nudge = 0,
  pie_to_jampie = TRUE,
  use_shadowText = NULL,
  vectorized_node_shapes = TRUE,
  edge_bundling = c("default", "connections", "none", "mark.groups", "nodegroups"),
  bundle_self = FALSE,
  nodegroups = NULL,
  render_nodes = TRUE,
  render_edges = TRUE,
  render_groups = TRUE,
  render_nodelabels = TRUE,
  params = NULL,
  plot_grid = FALSE,
  verbose = FALSE,
  debug = NULL
)
```

## Arguments

- x:

  The graph to plot.

- ...:

  Additional plotting parameters. See
  [igraph.plotting](https://r.igraph.org/reference/plot.common.html) for
  the complete list.

- axes:

  Logical, whether to plot axes, defaults to FALSE.

- add:

  Logical scalar, whether to add the plot to the current device, or
  delete the device's current contents first.

- xlim, ylim:

  default x and y axis limits. When either value is `NULL` the range is
  defined by the layout coordinate ranges, respectively, then expanded
  by adding `expand` to each side of the range.

- mark.groups:

  A list of vertex id vectors. It is interpreted as a set of vertex
  groups. Each vertex group is highlighted, by plotting a colored
  smoothed polygon around and “under” it. See the arguments below to
  control the look of the polygons.

- mark.shape:

  A numeric scalar or vector. Controls the smoothness of the vertex
  group marking polygons. This is basically the ‘shape’ parameter of the
  [`graphics::xspline()`](https://rdrr.io/r/graphics/xspline.html)
  function, its possible values are between -1 and 1. If it is a vector,
  then a different value is used for the different vertex groups.

- mark.col:

  A scalar or vector giving the colors of marking the polygons, in any
  format accepted by
  [`graphics::xspline()`](https://rdrr.io/r/graphics/xspline.html); e.g.
  numeric color ids, symbolic color names, or colors in RGB.

- mark.alpha:

  `numeric` value between 0 (transparent) and 1 (opaque) indicating the
  transparency of `mark.col` color fill values, used only when
  `mark.groups` is defined, and `mark.col` is not defined.

- mark.border:

  A scalar or vector giving the colors of the borders of the vertex
  group marking polygons. If it is `NA`, then no border is drawn.

- mark.expand:

  A numeric scalar or vector, the size of the border around the marked
  vertex groups. It is in the same units as the vertex sizes. If a
  vector is given, then different values are used for the different
  vertex groups.

- mark.lwd, mark.lty:

  line with and line type parameters for each `mark.groups` polygon.

- mark.cex:

  `numeric` adjustment for mark label font size, used when `mark.groups`
  is supplied and has `names(mark.groups)`.

- mark.x.nudge, mark.y.nudge:

  `numeric` values in units of the maximum x-axis or y-axis range for
  the layout coordinates, used to adjust each label displayed when
  `names(mark.groups)` is defined. These arguments are passed to
  [`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md)
  as `label.x.nudge`, `label.y.nudge`, respectively.

- pie_to_jampie:

  `logical` indicating whether to convert vertex shape `"pie"` to
  `"jampie"` in order to use vectorized plotting.

- use_shadowText:

  `logical` to enable
  [`jamba::shadowText()`](https://jmw86069.github.io/jamba/reference/shadowText.html)
  outlines.

  - Default NULL uses graph attribute 'use_shadowText' if it exists,
    otherwise FALSE.

  - When TRUE it uses
    [`jamba::shadowText()`](https://jmw86069.github.io/jamba/reference/shadowText.html)
    instead of
    [`graphics::text()`](https://rdrr.io/r/graphics/text.html) to render
    labels, applying an outline to contrast with the text color, either
    white or black. Some global options are recognized by
    [`jamba::shadowText()`](https://jmw86069.github.io/jamba/reference/shadowText.html)
    to customize the visual details. Consider using
    [`withr::with_options()`](https://withr.r-lib.org/reference/with_options.html)
    to set these options without persistence.

- vectorized_node_shapes:

  `logical` indicating whether to plot vertex node shapes using
  vectorized operations. It is substantially faster, however the one
  drawback is that nodes are plotted in order of their shape, which
  affects the positioning of nodes when there are node overlaps. This
  tradeoff is relatively minor, and it is recommended either to
  reposition nodes to reduce or prevent overlaps, or adjust node sizes
  to reduce overlaps.

- edge_bundling:

  `character` string or `function`, where:

  - `"default"` will try to detect an appropriate method: when
    `nodegroups` or `mark.groups` are defined, it chooses the matching
    option (see below); otherwise it chooses `"connections"`.

  - `"none"` will perform no edge bundling. This method is best when
    rendering straight edges, or for rendering multiple identical edges
    with curvature as defined by
    [`igraph::igraph.plotting()`](https://r.igraph.org/reference/plot.common.html).

  - `"connections"` will perform graph edge bundling by shared
    connections by calling
    [`edge_bundle_bipartite()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_bipartite.md)
    then
    [`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md).
    This option is particularly good for bipartite graphs such as
    concept networks (cnet plots).

  - `"mark.groups"` will perform graph edge bundling using `mark.groups`
    by calling
    [`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md).
    This option is equivalent to `"nodegroups"` except that it uses
    `mark.groups` to define node groupings.

  - `"nodegroups"` will perform graph edge bundling using `nodegroups`
    by calling
    [`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md).
    This option is equivalent to `"mark.groups"` except that it uses
    `nodegroups` to define node groupings.

  - `function` will call a custom edge bundling function using the
    `igraph` object `x` and the igraph parameters `param` as input. This
    output is currently untested, and is intended to enable alternative
    edge bundling functions which may exist outside this package. The
    custom function should be able to use the node layout coordinates in
    `graph_attr(x, "layout")`, and render edges between nodes.

- nodegroups:

  `list` object as output by
  [`edge_bundle_bipartite()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_bipartite.md)
  where each list element is a `character` vector of vertex node names
  present in `igraph::V(x)$name`. If no `"name"` vertex node attribute
  exists, then integer index values are used as names. Note that all
  vertex nodes must be represented in `nodegroup` in order for the
  corresponding edges to be plotted.

- render_nodes, render_edges:

  `logical` indicating whether to render vertex nodes, or edges,
  respectively. Sometimes it can be useful to call this function for
  other byproduct outputs, for example,
  `jam_plot_igraph(graph, add=FALSE, render_nodes=FALSE, render_edges=FALSE)`
  will create a new plot device with appropriate axis ranges, and can be
  used to render edge bundling results for example.

- render_groups:

  `logical` indicating whether to render groups when `mark.groups` is
  supplied. Groups are rendered with a shaded polygon and border.

- render_nodelabels:

  `logical` indicating whether to draw node labels, which is typically
  the last operation in the plot sequence. Note that node labels can be
  rendered without also rendering the nodes or edges.

- plot_grid:

  `logical` indicating whether to plot a background grid indicating
  units of 2% across the layout of the network graph. The units are
  calculated consistent with `nudge_igraph_nodes()`,
  [`adjust_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/adjust_cnet_nodeset.md)
  and other functions, scaled relative to the maximum x- or y-coordinate
  range of layout of the graph. Layout is obtained by
  [`get_igraph_layout()`](https://jmw86069.github.io/multienrichjam/reference/get_igraph_layout.md)
  which by default uses supplied `layout`, or graph attribute
  `igraph::graph_attr(x, "layout")`. Note that by default,
  [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  represents the layout with aspect ratio = 1, so x-coordinates and
  y-coordiantes are represented with the same spacing per unit. This
  function calls
  [`plot_layout_scale()`](https://jmw86069.github.io/multienrichjam/reference/plot_layout_scale.md)
  to render the grid lines.

- verbose:

  `logical` indicating whether to print verbose output.

- debug:

  `logical` or `character` vector of attributes for which debug output
  will be plotted onscreen.

## Details

Note that this function is intended to be called by
[`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md),
and is an internal function not intended to be called directly.

The
[`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
handles the overall plot equivalent of
[`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html),
however it calculates layout coordinates, and defines more useful x- and
y-axis ranges, and then adjusts node and label sizes relevant to the
layout data range. Specifically `vertex.size=15` is only useful when the
layout range is rescaled between -1 and 1; however when using
[`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
the vertex is scaled relative to the actual layout ranges.

The steps here are a reproduction of `igraph:::plot.igraph()` with four
changes:

1.  Default `rescale=FALSE`, and `asp=1` which means igraph layout is
    drawn true to the layout coordinates without distortion. To use
    default
    [`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html)
    behavior, use `rescale=TRUE`. The new default may not be appropriate
    for bipartite layout algorithms that generate two columns, and seems
    most useful with organic layouts where aspect ratio 1 helps convey
    important meaning in the graph structure, namely by enforcing
    consistent x- and y-axis visual distance between nodes.

    - Related: the `xlim` and `ylim` values are automatically adjusted
      to include the layout coordinate range. The default
      `igraph::plot.igraph(..., rescale=FALSE)` does not adjust the
      `xlim` and `ylim` ranges, which can be problematic when supplying
      layout as a function, and therefore the output node coordinates
      are not known until the plot rendering step.

When `vectorized_node_shapes=TRUE` by default:

1.  When there are multiple different vertex `"shape"` attributes, the
    nodes are rendered vectorized one shape at a time. The original
    [`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html)
    draws each individual vertex one by one, which is substantially
    slower (minutes compared to 1-2 seconds) for large `igraph` objects.

2.  When there are multiple font families, the default plot function
    draws each label one by one. The `jam_plot_igraph()` draws labels in
    groups of font family, in order to comply with limitations in
    [`graphics::text()`](https://rdrr.io/r/graphics/text.html). This
    situation is fairly rare, however the speed improvement is
    substantial, again roughly minutes down to 1-2 seconds.

The fourth difference involves edge bundling:

1.  When `edge_bundling` is used, it renders edges differently than the
    approach in
    [`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html),
    by drawing curved splines for each bundle of edges.

Some other distinctive features include:

When `use_shadowText=TRUE` node labels call
[`jamba::shadowText()`](https://jmw86069.github.io/jamba/reference/shadowText.html)
which draws a small partly transparent outline around labels, making
them more legible when they overlap colored nodes. This step effectively
draws each label `n` times, which can slightly slow the rendering of the
overall figure.

When `pie_to_jampie=TRUE`, any nodes with `shape="pie"` are changed to
`shape="jampie"` for the purpose of rendering pie shapes in vectorized
fashion, instead of being drawn for each node separately. This change is
a substantial improvement in rendering time.

Default colors for marked node groups `mark.col` and `mark.border` when
not defined upfront, will call
[`colorjam::rainbowJam()`](https://jmw86069.github.io/colorjam/reference/rainbowJam.html)
and not
`grDevices::rainbow(). The `colorjam::rainbowJam()`produces more visually distinct categorical colors. This behavior can be controlled by supplying a`character`vector with specific colors for`mark.col`and`mark.border`. Note that the border should match the colors, or it can be set to `"grey45"\`
for a generally visible border.

Optional argument `nodegroups` can be supplied, which is a `list` of
vectors, where each vector represents a group of nodes. The `nodegroups`
can be used with `edge_bundling="nodegroups"` to define custom edge
bundling.

Finally, individual plot components can be individually disabled:

- `render_nodes=FALSE`

- `render_edges=FALSE`

- `render_groups=FALSE`

- `render_nodelabels=FALSE`

Note that this function is not called by default, and is only called by
[`multienrichjam::jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md).

All other arguments are documented in
[`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html).

## See also

Other jam igraph internal functions:
[`default_igraph_values()`](https://jmw86069.github.io/multienrichjam/reference/default_igraph_values.md),
[`get_igraph_arrow_mode()`](https://jmw86069.github.io/multienrichjam/reference/get_igraph_arrow_mode.md),
[`jam_igraph_arrows()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph_arrows.md)
