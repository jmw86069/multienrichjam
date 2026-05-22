# Jam custom function to plot an igraph network

Jam custom function to plot an igraph network

## Usage

``` r
jam_igraph(
  x,
  ...,
  xlim = NULL,
  ylim = NULL,
  expand = 0.03,
  rescale = FALSE,
  node_factor = 1,
  node_factor_l = NULL,
  edge_factor = 1,
  edge_factor_l = NULL,
  label_factor = 1,
  label_factor_l = NULL,
  border_factor = 1,
  label_fontsize_l = NULL,
  label_dist_factor = 1,
  label_dist_factor_l = 1,
  use_shadowText = NULL,
  edge_bundling = c("default", "connections", "none", "mark.groups", "nodegroups"),
  bundle_self = FALSE,
  nodegroups = NULL,
  render_nodes = TRUE,
  render_edges = TRUE,
  render_nodelabels = TRUE,
  render_groups = TRUE,
  vectorized_node_shapes = TRUE,
  plot_grid = FALSE,
  plot_function = jam_plot_igraph,
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
  verbose = FALSE,
  debug = NULL
)
```

## Arguments

- x:

  `igraph` object to be plotted

- ...:

  additional arguments are passed to
  [`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html)

- xlim, ylim:

  default x and y axis limits. When either value is `NULL` the range is
  defined by the layout coordinate ranges, respectively, then expanded
  by adding `expand` to each side of the range.

- expand:

  `numeric` value used to expand the x and y axis ranges, where `0.03`
  expands each size `3%`.

- rescale:

  `logical` indicating whether to rescale the layout coordinates to
  `c(-1, 1)`. When `rescale=FALSE` the original layout coordinates are
  used as-is without change.

- node_factor:

  `numeric` value multiplied by `igraph::V(x)$size` to adjust the
  relative size of all nodes by a common numeric scalar value.

- node_factor_l, label_factor_l, label_dist_factor_l:

  `list` of vectors, where the names of the `list` are attribute names,
  and the names of each vector are attributes values. These values are
  applied in addition to `node_factor`, `label_factor`,
  `label_dist_factor`, respectively. The vector values are used as
  scalar multipliers, analogous to `node_factor`. The purpose is to
  apply scalar values to different subsets of nodes. For example,
  consider: `node_factor_l=list(nodeType=c(Gene=1, Set=2)`. The list
  name `"nodeType"` says to look at
  `igraph::vertex_attr(x, "nodeType")`. Nodes with `nodeType="Gene"`
  will use `1`, and `nodeType="Set"` will use `2` as the scalar value.

- edge_factor:

  `numeric` value multiplied by `igraph::E(x)$width` to adjust the
  relative width of all edges by a common numeric scalar value.

- label_factor:

  `numeric` value multiplied by `igraph::V(x)$label.cex` and
  `igraph::E(x)$label.cex` to adjust the relative size of all labels on
  nodes and edges by a common numeric scalar value.

- border_factor:

  `numeric` value multiplied by frame.width and pie.lwd to adjust the
  relative width of node borders. Note that border line width is also
  affected by `node_factor` and `node_factor_l` when provided. To
  prevent borders being adjusted by `node_factor`, set `border_factor`
  to the reciprocal, for example `border_factor=1/node_factor`.

- label_dist_factor:

  `numeric` value multiplied by `igraph::V(x)$label.dist` to adjust the
  relative distance of all nodes labels from the node center by a common
  numeric scalar value.

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

- render_nodelabels:

  `logical` indicating whether to draw node labels, which is typically
  the last operation in the plot sequence. Note that node labels can be
  rendered without also rendering the nodes or edges.

- render_groups:

  `logical` indicating whether to render groups when `mark.groups` is
  supplied. Groups are rendered with a shaded polygon and border.

- vectorized_node_shapes:

  `logical` indicating whether to plot vertex node shapes using
  vectorized operations. It is substantially faster, however the one
  drawback is that nodes are plotted in order of their shape, which
  affects the positioning of nodes when there are node overlaps. This
  tradeoff is relatively minor, and it is recommended either to
  reposition nodes to reduce or prevent overlaps, or adjust node sizes
  to reduce overlaps.

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
  `jam_igraph()` represents the layout with aspect ratio = 1, so
  x-coordinates and y-coordiantes are represented with the same spacing
  per unit. This function calls
  [`plot_layout_scale()`](https://jmw86069.github.io/multienrichjam/reference/plot_layout_scale.md)
  to render the grid lines.

- plot_function:

  `function` that renders the graph, not intended to be changed except
  for very customized uses. By default `plot_function=jam_plot_igraph()`
  which calls a modified variant of `igraph:::plot.igraph()`.

- mark.groups:

  one of the following:

  - `list` with vertex (node) names

  - `list` output from
    [`igraph::components`](https://r.igraph.org/reference/components.html)
    which contains 'membership', 'csize', and 'no' elements.

  - `TRUE` will attempt
    [`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md)
    when the graph is bipartite and contains vertex attribute name
    'nodeType' as with Cnet data. **To turn off labels, set
    `mark.cex=0`.**

  - `FALSE` will not apply any 'mark.groups' even when provided in graph
    attributes.

  - `communities` object from one of the igraph cluster functions such
    as
    [`igraph::cluster_walktrap()`](https://r.igraph.org/reference/cluster_walktrap.html).

  - `function` that takes `igraph` input and returns either `list` or
    `communities` output. For example, using
    `mark.groups=igraph::cluster_walktrap` will apply the function on
    the data dynamically.

  - When not defined, if graph attributes include 'mark.groups' it will
    be used. Graph attributes are defined by
    [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md)
    for example, `igraph::graph_attr(x, "mark.groups")`. Note
    'mark.colors' is also used when defined together with 'mark.groups'.

- mark.shape:

  `numeric` default 1/2 ultimately passed to
  [`graphics::xspline()`](https://rdrr.io/r/graphics/xspline.html) when
  smoothing the alpha hull produced by
  [`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md).

- mark.col:

  `character` color assigned to each group in `mark.groups`. When NULL,
  default, it calls
  [`colorjam::rainbowJam()`](https://jmw86069.github.io/colorjam/reference/rainbowJam.html)
  for categorical colors, which is applied together with `mark.alpha`
  for transparency.

- mark.alpha:

  `numeric` value between 0 (transparent) and 1 (opaque) indicating the
  transparency of `mark.col` color fill values, used only when
  `mark.groups` is defined, and `mark.col` is not defined.

- mark.border:

  `character` color assigned to each group in `mark.groups`. When NULL,
  default, it calls
  [`colorjam::rainbowJam()`](https://jmw86069.github.io/colorjam/reference/rainbowJam.html)
  for categorical colors.

- mark.expand:

  `numeric` default NULL will choose a sensible default. The expand is a
  fraction of point range, typically `0.05` is 5% which is reasonable.
  The expansion is used by
  [`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md)
  in a non-exact process by which the region is adjusted until a
  suitable alpha hull solution is found for the points provided.

- mark.lwd, mark.lty:

  line with and line type parameters for each `mark.groups` polygon.

- mark.smooth:

  `logical` default TRUE, whether to smooth the point hull boundary.

- mark.cex:

  `numeric` default 1, used to resize the mark.group labels. **To turn
  off labels, set `mark.cex=0`.**

- mark.x.nudge, mark.y.nudge:

  `numeric` values in units of the maximum x-axis or y-axis range for
  the layout coordinates, used to adjust each label displayed when
  `names(mark.groups)` is defined. These arguments are passed to
  [`make_point_hull()`](https://jmw86069.github.io/multienrichjam/reference/make_point_hull.md)
  as `label.x.nudge`, `label.y.nudge`, respectively.

- verbose:

  `logical` indicating whether to print verbose output.

- debug:

  `logical` or `character` vector of attributes for which debug output
  will be plotted onscreen.

## Details

This function is a drop-in replacement of
[`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html),
intended to provide substantially faster vectorized plotting, to render
bundled edges when requested, and to handle `rescale=FALSE` without
requiring further adjustments. Note that this function focuses on
recognizing graph options and settings, then passes the work off to
[`jam_plot_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_plot_igraph.md)
which performs the heavy work of rendering the graph.

It also provides some convenient methods to adjust node size, label font
size, and label distance from node center, based upon node attributes.

### vectorized plots

This function calls
[`jam_plot_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_plot_igraph.md)
as a replacement for
[`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html),
and that function implements vectorized plot features when
`vectorized_node_shapes=TRUE` by default:

1.  When there are multiple different vertex `"shape"` attributes, this
    function renders nodes vectorized, one shape at a time. In this
    scenario, the original
    [`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html)
    draws each individual vertex, which is substantially slower (minutes
    compared to 1-2 seconds) for large `igraph` objects.

2.  When there are multiple font families, labels are rendered in groups
    by font family, in order to comply with limitations in
    [`graphics::text()`](https://rdrr.io/r/graphics/text.html). This
    situation is fairly rare, however the speed improvement is
    substantial, again roughly minutes down to 1-2 seconds. The
    [`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html)
    renders each node label individually when there are multiple font
    families.

### rescale=FALSE

The default
[`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html)
uses `rescale=TRUE`, which can distort the layout coordinates to fit
within a fixed x- and y-axis range `c(-1, 1)`. When using
`rescale=FALSE` the `xlim` and `ylim` values are not adjusted to the
actual coordinate range. The desired effect from this function
`jam_igraph()` is to apply `aspect=1` (`asp=1`) which fixes the aspect
ratio so the coordinates represent visual Euclidean distance between
nodes, and to define `xlim` and `ylim` to accomodate the full layout.
This function also adjusts node `vertex.size` and `vertex.label.dist`
proportionally.

### edge bundling

When `edge_bundling` is something other than `edge_bundling="none"`,
edge connections between nodes are rendered using a specific function by
drawing curved splines for each bundle of edges. The approach in
[`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html)
only draws straight edges between nodes. The recommended method is
`edge_bundling="default"` which will try to detect an appropriate method
to bundle edges. When `mark.groups` and `nodegroups` are not defined,
the default method is `"connections"` which bundles edges only among
nodes that share the same connections. The assumption is that nodes that
share the same connections usually have very similar layout coordinates,
so edge bundling is usually intuitive. In fact, for a very large set of
nodes, they are often in a round cluster

### Adjust node size, label size, label distance

The following arguments apply scaling to all nodes or edges:

- `node_factor` - `numeric` multiplied by node size

- `edge_factor` - `numeric` multiplied by edge size

- `label_factor` - `numeric` multiplied by label font size

- `label_dist_factor` - `numeric` multiplied by label distance from node
  center

New attributes `vertex.label.fontsize` and `edge.label.fontsize` which
define fixed fontsize in points for nodes and edges, respectively. These
values are not modified by `vertex.label.cex` nor `edge.label.cex` and
are intended to allow control over specific fonts used in the final
figure. Note their calculations are based upon `par("ps")` which should
represent device-dependent point size. If this value is inappropriate,
it should be adjusted to control the font sizing.

The following arguments apply scale factor based upon node attribute:

- `node_factor_l` - `list` of named vectors applied to node size

- `label_factor_l` - `list` of named vectors applied as `label_factor` .
  as a multiplier to label font size.

- `label_fontsize_l` - `list` of named vectors applied to define a
  specific, fixed label fontsize in points, which is not modified . by
  `vertex.label.cex` nor `label_factor`.

- `label_dist_factor_l` - `list` of named vectors applied to label
  distance from node center

The `factor_l` technique is as follows:

    node_factor_l = list(node_attr_name = c(
       attr_value1 = factor1,
       attr_value2 = factor2))

A specific example:

    node_factor_l = list(nodeType = c(
       Gene=1.5,
       Set=2))

In this case, nodes with attribute `igraph::V(g)$nodeType == "Gene"`
will use factor `1.5` Nodes with attribute
`igraph::V(g)$nodeType == "Set"` will use factor `2` All other nodes
will not be adjusted.

### Other features

The plot layout by default is **not** rescaled to `c(-1, 1)`, therefore
allowing direct control over plot dimensions and node sizes. The plot
aspect ratio is fixed at 1, which renders many network layouts in their
intended form, as opposed to scaling each axis to `c(-1, 1)`, which can
impose distortion of intended layout node distances.

When `use_shadowText=TRUE` node labels call
[`jamba::shadowText()`](https://jmw86069.github.io/jamba/reference/shadowText.html)
which draws a small partly transparent outline around labels, making
them more legible when they overlap colored nodes. This step effectively
draws each label `n` times, which can slightly slow the rendering of the
overall figure.

When `pie_to_jampie=TRUE`, any nodes with `shape="pie"` are changed to
`shape="jampie"` for the purpose of rendering pie shapes in vectorized
fashion, instead of being drawn for each node separately. This change is
a substantial improvement in rendering time. In addition, optional node
attributes are available:

- `pie.border` to control individual pie wedge borders, which are drawn
  as inner borders so each pie wedge border is visible without

- `pie.lwd` to control line width of pie wedge borders.

- `pie.lty` to control the line type of pie wedge borders

- `frame.color` to control the frame border color drawn around the full
  circular pie node. This border is drawn as an outer border, so it will
  not overlap any internal pie wedge border colors.

Default colors for marked node groups `mark.col` and `mark.border` when
not defined upfront, will call
[`colorjam::rainbowJam()`](https://jmw86069.github.io/colorjam/reference/rainbowJam.html)
and not
[`grDevices::rainbow()`](https://rdrr.io/r/grDevices/palettes.html). The
[`colorjam::rainbowJam()`](https://jmw86069.github.io/colorjam/reference/rainbowJam.html)
produces more visually distinct categorical colors. This behavior can be
controlled by supplying a `character` vector with specific colors for
`mark.col` and `mark.border`. Note that the border should match the
colors, or it can be set to `"grey45"` for a generally visible border.

When `mark.groups` is `NULL` and graph attribute 'mark.groups' is
defined, it will be used.

When `names(mark.groups)` is defined, the values are used as labels,
positioned at the outer edge of each polygon. The label text size is
adjusted with `label.cex`, and the position can be adjusted with
`mark.x.nudge`, `mark.y.nudge`, in units of fraction of the maximum x-
or y-axis range (effectively fraction of the layout size).

When `mark.groups=TRUE` it will attempt to call
[`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md)
to define nodegroups, and use them as the `mark.groups`. Note that '...'
are passed to this function, meaning `min_size` can be used to filter
nodesets to minimum size, and `min_points` can be used to filter the
alpha hull to require this many points.

Optional argument `nodegroups` can be supplied, which is a `list` of
vectors, where each vector represents a group of nodes. The `nodegroups`
can be used with `edge_bundling="nodegroups"` to define custom edge
bundling. This option is useful for defining a group of nodes for edge
bundling, when those nodes should not be used to render group borders as
with `mark.groups`.

Finally, individual plot components can be individually disabled:

- `render_nodes=FALSE`

- `render_edges=FALSE`

- `render_groups=FALSE`

- `render_nodelabels=FALSE`

## See also

Other multienrichjam core functions:
[`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md),
[`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

Other jam igraph utilities:
[`color_edges_by_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodegroups.md),
[`color_edges_by_nodes()`](https://jmw86069.github.io/multienrichjam/reference/color_edges_by_nodes.md),
[`color_nodes_by_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/color_nodes_by_nodegroups.md),
[`get_cnet_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_cnet_nodeset.md),
[`launch_shinycat()`](https://jmw86069.github.io/multienrichjam/reference/launch_shinycat.md),
[`subgraph_jam()`](https://jmw86069.github.io/multienrichjam/reference/subgraph_jam.md)

## Examples

``` r
# Example with karate
karate <- igraph::make_graph("Zachary");
cl <- igraph::cluster_louvain(karate);
jam_igraph(karate,
   layout=layout_with_qfrf(repulse=3.5),
   mark.groups=cl,
  mark.lwd=c(1:4),
  mark.lty=1:4, mark.shape=1,
  edge_bundling="default")


# create example cnet data
cnet <- make_cnet_test(num_sets=3)

## example showing how to use the list form
## This form resizes nodes where igraph::V(g)$nodeType %in% "Gene" by 2x,
## and resizes nodes where igraph::V(g)$nodeType %in% "Set" by 3x.
node_factor_l <- list(nodeType=c(Gene=1.2, Set=2));

## This form multiplies label.dist for nodeType="Gene" nodes by 1,
## and multiplies label.dist for nodeType="Set" nodes by 0.5
label_dist_factor_l <- list(nodeType=c(Gene=1, Set=0.5))

par("mar"=c(0, 0, 0, 0) + 0.5);
jam_igraph(cnet,
   use_shadowText=TRUE,
   node_factor_l=node_factor_l,
   label_factor=0.6,
   label_factor_l=list(nodeType=c(Gene=1, Set=2)))

par("mar"=c(2, 2, 2, 2));

# Example using edge bundling by community detection
g <- igraph::make_graph("Zachary");
gcom <- igraph::cluster_leading_eigen(g);

jam_igraph(g,
   layout=layout_with_qfr,
   edge_bundling="nodegroups",
   mark.groups=gcom,
   nodegroups=gcom,
   vertex.color=colorjam::group2colors(gcom$membership))


cfuncs <- list(cluster_leading_eigen=igraph::cluster_leading_eigen,
   cluster_edge_betweenness=igraph::cluster_edge_betweenness,
   cluster_fast_greedy=igraph::cluster_fast_greedy,
   cluster_spinglass=igraph::cluster_spinglass)
for (i in seq_along(cfuncs)) {
   cfunc <- cfuncs[[i]];
   gcom <- cfunc(g);
   igraph::V(g)$color <- colorjam::group2colors(gcom$membership);
   g <- color_edges_by_nodes(g);
   set.seed(123);
   jam_igraph(g,
      layout=layout_with_qfr,
      edge_bundling="nodegroups",
      nodegroups=gcom,
      mark.groups=gcom)
   title(main=names(cfuncs)[i]);
}





# fancy example showing mark.groups and colorizing
# edges using node colors
gcom <- igraph::cluster_spinglass(g);
igraph::V(g)$color <- colorjam::group2colors(gcom$membership);
g <- color_edges_by_nodes(g);
jam_igraph(g,
   layout=layout_with_qfrf(repulse=3.2),
   edge_bundling="nodegroups",
   nodegroups=gcom,
   mark.groups=gcom)
title(main=paste0("cluster_spinglass()\n",
   "edge_bundling='nodegroups'"))


# same but different edge_style
jam_igraph(g,
   layout=layout_with_qfrf(repulse=3.2),
   edge_bundling="nodegroups",
   nodegroups=gcom,
   mark.groups=gcom,
   bundle_style="xspline",
   detail=14)
title(main="bundle_style='xspline'")


# same but using node connections
jam_igraph(g,
   layout=layout_with_qfrf(repulse=3.2),
   edge_bundling="connections",
   nodegroups=gcom,
   mark.groups=gcom)
title(main="edge_bundling='connections'")

```
