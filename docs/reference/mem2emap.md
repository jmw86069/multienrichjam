# Convert multiEnrichMap mem output to EnrichmentMap emap

Convert multiEnrichMap mem output to EnrichmentMap emap

## Usage

``` r
mem2emap(
  mem,
  overlap = 0.2,
  p_cutoff = NULL,
  min_count = NULL,
  colorV = NULL,
  cluster_function = igraph::cluster_walktrap,
  cluster_list = NULL,
  num_keep_terms = 3,
  keep_terms_sep = ",\n",
  repulse = 3.3,
  remove_singlets = FALSE,
  color_by_nodes = FALSE,
  size_by_genes = TRUE,
  median_size = 5,
  max_size = 25,
  apply_edge_width = TRUE,
  max_edge_width = 20,
  min_edge_width = 0.5,
  apply_edge_color = apply_edge_width,
  apply_direction = TRUE,
  direction_max = 2,
  direction_floor = 0.5,
  border_lwd = 3,
  seed = 123,
  spread_labels = TRUE,
  label_min_dist = 2.5,
  y_bias = 10,
  vertex.label.font = 2,
  use_layout_components = FALSE,
  use_shadowText = TRUE,
  do_express = FALSE,
  do_plot = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- mem:

  `Mem` output from
  [`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)

- overlap:

  `numeric`, default 0.2, value between 0 and 1 with Jaccard overlap
  coefficient required between any two pathways in order to create a
  network edge connecting these two pathways.

  - Jaccard overlap coefficient is reciprocal of Jaccard distance, using
    `(1 - Jaccard_distance)`, where Jaccard_distance is calculated with
    `dist(x, method='binary')` with incidence matrix.

  - `overlap=0.2` is default, which specifies roughly 20% overlap in
    genes shared between two pathway nodes, relative to the union of
    genes in those two pathway nodes.

  - Note these genes must be involved in enrichment, and therefore does
    not use all possible genes annotated to a pathway. Therefore,
    connections are only created with enriched genes are shared between
    pathways.

- p_cutoff:

  `numeric` threshold used for significant enrichment P-value. The
  default NULL uses the value in the `mem` object. This cutoff, together
  with `min_count`, define which nodes are included, and the color fill
  for those nodes.

- min_count:

  `integer` threshold for minimum genes involved in enrichment in order
  for a pathway to be considered significant during this analysis. When
  NULL it uses the threshold in `mem`.

- colorV:

  `character` vector of R colors used for each enrichment. Default NULL
  uses colors defined in `mem`.

- cluster_function:

  `function`, default
  [`igraph::cluster_walktrap`](https://r.igraph.org/reference/cluster_walktrap.html),
  used to cluster `igraph` nodes in the resulting network graph. It is
  used to enhance the visual summary.

  - This function is not used when `cluster_list` is supplied.

- cluster_list:

  `list` default NULL, optional list of pathway clusters, containing
  `character` vectors of `igraph` node/vertex names. It will be used
  instead of `cluster_function` when supplied.

- num_keep_terms:

  `integer`, default 3, number of text terms to keep for each pathway
  cluster, used together with `cluster_function`. Common words are
  removed, remaining terms are sorted by decreasing occurrence, then
  used to summarize each cluster.

- keep_terms_sep:

  `character` string, default is comma-newline, used to separate terms,
  used together with `num_keep_terms`.

- repulse:

  `numeric` value passed to
  [`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md),
  default 3.3.

  - Use repulse 'FALSE', 'NULL', or '0' will skip the layout, which is
    used by
    [`mem_find_overlap()`](https://jmw86069.github.io/multienrichjam/reference/mem_find_overlap.md)
    to iterate numerous overlaps without spending time on layout each
    iteration.

- remove_singlets:

  `logical`, default FALSE, whether to remove pathway singlets which
  have no connections to other pathways.

  - Using TRUE will help simplify busy figures, at the expense of being
    less complete.

- color_by_nodes:

  `logical`, default TRUE, whether to colorize pathway clusters using
  node colors in each cluster. Note that a mix of colors often turns
  brown, so this feature has unpredictable benefit.

- size_by_genes:

  `logical` default TRUE, whether to size nodes proportional to the
  number of genes involved in enrichment.

- median_size:

  `numeric` used when `size_by_genes=TRUE` to define the median size,
  default is 5.

- apply_edge_width:

  `logical` default TRUE, whether to apply edge width based upon the
  edge weight, which represents the Jaccard overlap.

- max_edge_width:

  `numeric` maximum edge width, used when `apply_edge_width=TRUE`.
  Default 25 is intended to represent roughly the full size of a node.

- min_edge_width:

  `numeric` minimum edge width, used after `max_edge_width` has been
  applied. The minimum is applied as a noise floor, so that any edge at
  or below the `min_edge_width` is set to this width. All edges above
  this width are left as-is.

- apply_edge_color:

  `logical` whether to apply edge color based upon the edge weight,
  which represents the Jaccard overlap. The default matches
  `apply_edge_width` so that they are applied together.

- apply_direction:

  `logical` default TRUE, whether to apply direction to the node border
  color, when supplied. Used with `direction_max` and `direction_floor`.

- direction_max:

  `numeric`, default 2, indicating the directional score at which the
  maximum color is applied. Typically when using some derivative of
  z-score to represent directionality, values at or above 2 are
  effectively "maximum".

- direction_floor:

  `numeric` default 0.5, the minimum directional score for any color to
  be applied, below which values are considered "non-directional".
  Typically when using some form of z-score, a value at or below 0.5 may
  be considered "not directional."

- border_lwd:

  `numeric` line width of node borders, used to convey any direction
  associated with enrichment.

- seed:

  `numeric` passed to [`set.seed()`](https://rdrr.io/r/base/Random.html)
  via the
  [`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md)
  layout algorithm, default 123.

- spread_labels:

  `logical` default TRUE, whether to call
  [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md).

- label_min_dist:

  `numeric` default 2.5 defines the label distance used in
  [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md).

- y_bias:

  `numeric` default 10, passed to
  [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md).

- vertex.label.font:

  `integer` to define the font style, default 2 sets bold font;
  alternatively use 1 for normal font.

- use_layout_components:

  `logical` default FALSE, whether to call
  [`igraph::layout_components()`](https://r.igraph.org/reference/merge_coords.html)
  to combine multiple network components, rather than calling
  [`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md)
  on the overall graph object.

- use_shadowText:

  `logical` default TRUE, whether to render node labels using
  [`jamba::shadowText()`](https://jmw86069.github.io/jamba/reference/shadowText.html).
  This setting is also stored in the graph attributes and is used by
  [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md).

- do_express:

  `logical` default FALSE, when TRUE it skips a number of aesthetic
  steps. This option is intended mainly for
  [`mem_find_overlap()`](https://jmw86069.github.io/multienrichjam/reference/mem_find_overlap.md)
  to perform more rapid iterative evaluation of overlap thresholds.

- do_plot:

  `logical`, default FALSE, whether to render the resulting plot using
  [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md).

- verbose:

  `logical` whether to print verbose output.

- ...:

  additional arguments are passed to
  [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  to customize the network plot, used when `do_plot=TRUE`.

## Value

`igraph` object with an Enrichment Map, nodes are pathways, and edges
are defined by the Jaccard overlap of genes involved in enrichment of
each pair of pathways.

- Edges are defined only when the Jaccard overlap is at least `overlap`.

- Node sizes are proportional to the number of genes involved in
  enrichment, with median node size defined by `median_size`.

- By default, edge width and color is defined by the overlap, ranging
  from 0 to 1. Graph attributes store the `max_edge_width` with edge
  width at Jaccard overlap 1; and `edge_colramp` as a function to
  convert weight to color, and can be used to create a color legend.

## Details

This function is currently In development.

This function is intended to convert `mem` output from
[`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
into an EnrichmentMap `igraph` format which represents the statistical
enrichment support from each pathway enrichment.

This function can apply P-value thresholds using the input `mem`, or
using a custom value.

A node community detection function is applied by default which may help
visualize subgroups of nodes. When defined, the communities are used to
define nodegroups for edge bundling.

As of version 0.0.116.900, node size also scales the border line width
for 'frame.width' and 'pie.border.lwd', so that the line width defined
by `border_lwd` matches `median_size`, and is scaled proportionally by
node size.

## See also

Other jam Mem utilities:
[`Mem-class`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md),
[`MemPlotFolio-class`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md),
[`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md),
[`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md)

Other jam igraph functions:
[`communities2nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/communities2nodegroups.md),
[`drawEllipse()`](https://jmw86069.github.io/multienrichjam/reference/drawEllipse.md),
[`edge_bundle_bipartite()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_bipartite.md),
[`edge_bundle_nodegroups()`](https://jmw86069.github.io/multienrichjam/reference/edge_bundle_nodegroups.md),
[`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md),
[`flip_edges()`](https://jmw86069.github.io/multienrichjam/reference/flip_edges.md),
[`get_bipartite_nodeset()`](https://jmw86069.github.io/multienrichjam/reference/get_bipartite_nodeset.md),
[`highlight_edges_by_node()`](https://jmw86069.github.io/multienrichjam/reference/highlight_edges_by_node.md),
[`igraph2pieGraph()`](https://jmw86069.github.io/multienrichjam/reference/igraph2pieGraph.md),
[`label_communities()`](https://jmw86069.github.io/multienrichjam/reference/label_communities.md),
[`mem2cnet()`](https://jmw86069.github.io/multienrichjam/reference/mem2cnet.md),
[`nodegroups2communities()`](https://jmw86069.github.io/multienrichjam/reference/nodegroups2communities.md),
[`rectifyPiegraph()`](https://jmw86069.github.io/multienrichjam/reference/rectifyPiegraph.md),
[`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md),
[`subsetCnetIgraph()`](https://jmw86069.github.io/multienrichjam/reference/subsetCnetIgraph.md),
[`subset_igraph_components()`](https://jmw86069.github.io/multienrichjam/reference/subset_igraph_components.md),
[`sync_igraph_communities()`](https://jmw86069.github.io/multienrichjam/reference/sync_igraph_communities.md)

## Examples

``` r
emap <- mem2emap(Memtest, min_count=2, max_edge_width=5, overlap=0.33)
jam_igraph(emap)


emap2 <- relayout_nodegroups(emap, final_repulse=4, y_bias=10, label_min_dist=2.5)
jam_igraph(emap2)


emapB <- mem2emap(fixSetLabels(Memtest), min_count=3,
   overlap=0.35, median_size=10,
   apply_edge_width=TRUE, max_edge_width=15)
jam_igraph(emapB)

```
