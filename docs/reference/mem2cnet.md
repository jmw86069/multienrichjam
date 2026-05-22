# Convert MultiEnrichment incidence matrix to Cnet plot

Convert MultiEnrichment incidence matrix to Cnet plot

## Usage

``` r
mem2cnet(
  memIM,
  categoryShape = c("pie", "coloredrectangle", "circle", "ellipse"),
  geneShape = c("pie", "coloredrectangle", "circle", "ellipse"),
  forceColors = FALSE,
  categoryColor = "#E5C494",
  geneColor = "#B3B3B3",
  categoryLabelColor = "darkblue",
  geneLabelColor = "grey25",
  categorySize = 12,
  geneSize = 6,
  categoryCex = 1,
  geneCex = 0.8,
  frame_darkFactor = NULL,
  geneIM = NULL,
  geneIMcolors = NULL,
  geneIMdirection = NULL,
  enrichIM = NULL,
  enrichIMcolors = NULL,
  enrichIMdirection = NULL,
  coloredrect_nrow = 1,
  coloredrect_ncol = NULL,
  coloredrect_byrow = TRUE,
  colorV = NULL,
  direction_col_fn = NULL,
  direction_cutoff = NULL,
  direction_max = 2,
  gene_direction_cutoff = 0,
  gene_direction_max = 1.2,
  gene_direction_col_fn = NULL,
  hide_solo_pie = TRUE,
  remove_blanks = TRUE,
  remove_singlet_genes = TRUE,
  do_reorder = TRUE,
  spread_labels = FALSE,
  vertex.label.font = 2,
  use_shadowText = TRUE,
  repulse = 3.5,
  verbose = FALSE,
  ...
)

memIM2cnet(
  memIM,
  categoryShape = c("pie", "coloredrectangle", "circle", "ellipse"),
  geneShape = c("pie", "coloredrectangle", "circle", "ellipse"),
  forceColors = FALSE,
  categoryColor = "#E5C494",
  geneColor = "#B3B3B3",
  categoryLabelColor = "darkblue",
  geneLabelColor = "grey25",
  categorySize = 12,
  geneSize = 6,
  categoryCex = 1,
  geneCex = 0.8,
  frame_darkFactor = NULL,
  geneIM = NULL,
  geneIMcolors = NULL,
  geneIMdirection = NULL,
  enrichIM = NULL,
  enrichIMcolors = NULL,
  enrichIMdirection = NULL,
  coloredrect_nrow = 1,
  coloredrect_ncol = NULL,
  coloredrect_byrow = TRUE,
  colorV = NULL,
  direction_col_fn = NULL,
  direction_cutoff = NULL,
  direction_max = 2,
  gene_direction_cutoff = 0,
  gene_direction_max = 1.2,
  gene_direction_col_fn = NULL,
  hide_solo_pie = TRUE,
  remove_blanks = TRUE,
  remove_singlet_genes = TRUE,
  do_reorder = TRUE,
  spread_labels = FALSE,
  vertex.label.font = 2,
  use_shadowText = TRUE,
  repulse = 3.5,
  verbose = FALSE,
  ...
)
```

## Arguments

- memIM:

  one of the following:

  - `Mem` S4 object, preferred. In this case, the arguments regarding
    'geneIM' and 'enrichIM' data are taken directly from 'mem'.

  - legacy `list` mem object, for backward compatibility,

  - `numeric` matrix in the form of a `memIM` gene-pathway incidence
    matrix. In this case, other arguments involving `geneIM` and
    `enrichIM` matrices are required for correct behavior of this
    function. When `mem` format is supplied, relevant arguments which
    are empty will use corresponding data from `mem`, for example
    `geneIM`, `geneIMcolors`, `enrichIM`, `enrichIMcolors`.

- categoryShape, geneShape:

  `character` string with node shape, default 'pie' uses pie nodes. Note
  that 'pie' shapes with only one segment are converted to 'circle' for
  convenience. When using
  [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md)
  to plot, it already treats single-segment pie nodes as circle anyway,
  and vectorizes pie node plotting to make the overall rendering time
  substantially faster.

- forceColors:

  `logical` default FALSE, whether to force node colors to use
  `categoryColor` and `geneColor` arguments, otherwise colors are
  defined using
  [`enrichIMcolors()`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md)
  and
  [`geneIMcolors()`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md),
  respectively. Enabling this option will skip all pie and directional
  styling.

- categoryColor, geneColor:

  `character` R color for default node colors, used when
  `geneIMcolors`,`enrichIMcolors` is not supplied, respectively.

- categoryLabelColor, geneLabelColor:

  `character` R color used as default node label color.

- categorySize, geneSize:

  `numeric` default node size.

- categoryCex, geneCex:

  `numeric` adjustment to default node label font size.

- frame_darkFactor:

  `numeric` passed to
  [`jamba::makeColorDarker()`](https://jmw86069.github.io/jamba/reference/makeColorDarker.html)
  so the frame color is slightly darker than the node fill color.

- geneIM, geneIMcolors, geneIMdirection, enrichIM, enrichIMcolors,
  enrichIMdirection:

  `matrix` data used only when the input 'mem' is not an `Mem` or legacy
  `list` mem object which already contains these data. When input 'mem'
  is supplied as a `matrix`, these value enable this function to operate
  on nearly any custom data.

  - `geneIMcolors` is used to define gene (row) node colors

  - `enrichIMcolors` is used to define set (column) node colors

  - `geneIMdirection`,`enrichIMdirection` is used to define optional
    border colors defined by the direction, where -1 is down, 0 is no
    change, and +1 is up. Use `direction_col_fn` to define a custom
    color function, however the default uses the reversed "RdBu" Brewer
    color ramp with blue (down), white (no change) and red (up).

- coloredrect_nrow, coloredrect_ncol, coloredrect_byrow:

  arguments used when `geneShape="coloredrectangle"`, to define layout
  and placement of colors across columns in `geneIMcolors`. By default,
  one row of colors is used.

- colorV:

  `character` optional vector of R colors, taken from 'mem' for `Mem` or
  legacy `list` mem objects. When supplied, it should be a vector named
  to match `colnames(enrichIM)`. When defined, these colors override
  `enrichIMcolors`.

- direction_col_fn:

  `function` used to colorize 'Set' node borders via enrichIMdirection.
  The default uses
  [`colorjam::col_div_xf()`](https://jmw86069.github.io/colorjam/reference/col_div_xf.html)
  which applies reverse Brewer "RdBu" for blue (down), white (no
  change), and red (up) with maximum color at `2.0`. When not supplied,
  `direction_cutoff` and `direction_max` are used with
  `colorjam::col_div_xf().`

- direction_cutoff:

  `numeric` default is taken from `memIM` when provided as `Mem`,
  typical default is `1.0` to require directional z-score at least `1.0`
  in order to apply any directional border color.

- direction_max:

  `numeric` default 2, the numeric value to apply maximum color from the
  color gradient function.

- gene_direction_cutoff:

  `numeric` default 0, the minimum value in geneIMdirection to apply
  up/down color. Using 0 will colorize any non-zero color.

- gene_direction_max:

  `numeric` default 1.2, the numeric value to apply maximum color from
  the color gradient function. Typical values are `-1` and `1`, so using
  `1.2` will apply nearly the maximum color.

- gene_direction_col_fn:

  `function` used to colorize 'Gene' node borders via enrichIMdirection.
  The default uses
  [`colorjam::col_div_xf()`](https://jmw86069.github.io/colorjam/reference/col_div_xf.html)
  which applies reverse Brewer "RdBu" for blue (down), white (no
  change), and red (up) with maximum color at `1.2`. \*When not
  supplied, `gene_direction_cutoff` and `gene_direction_max` are used
  with `colorjam::col_div_xf().`

- hide_solo_pie:

  `logical` default TRUE, passed to
  [`apply_cnet_direction()`](https://jmw86069.github.io/multienrichjam/reference/apply_cnet_direction.md)
  to determine whether to display border only as one outer frame color
  when all colors are identical. When FALSE, all pie wedges are
  individually colored.

- remove_blanks:

  `logical` default TRUE, whether to remove blank color subsections from
  each node, using
  [`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md).
  This argument is useful for 'pie', 'jampie', or 'coloredrectangle'
  node shapes, so they will only indicate the relevant color.

- remove_singlet_genes:

  `logical` default TRUE, whether to remove singlet genes, which are
  genes (rows) not represented in any pathway gene sets (columns).

- do_reorder:

  `logical` default TRUE used to re-order nodes in "equivalent
  positions" by border color, color, etc. by calling
  [`reorder_igraph_nodes()`](https://jmw86069.github.io/multienrichjam/reference/reorderIgraphNodes.md).
  Equivalent nodes are common in Cnet (bipartite) networks, for example
  numerous 'Gene' nodes may be connected to the same 2 pathway 'Set'
  nodes, and therefore the node positions are interchangeable. Sorting
  by node attributes (notably color, then alphabetically by label) helps
  visual review.

- spread_labels:

  `logical` default FALSE, whether to spread node labels away from
  incoming edges, also adding label distance via vertex attribute
  'label.dist'. This step calls
  [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md),
  and can be customized further by passing arguments through '...'
  ellipses.

- vertex.label.font:

  `integer` default 2, used to set default font to use bold face
  (font=2) for labels. Use 1 for normal font.

- use_shadowText:

  `logical` default TRUE, applied to graph attributes to apply
  shadowText by default when plotting with
  [`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md).

- repulse:

  `numeric` value passed to
  [`layout_with_qfr()`](https://jmw86069.github.io/multienrichjam/reference/layout_with_qfr.md)
  when either `do_reorder` or `spread_labels` is TRUE. Otherwise there
  is no layout applied. The 'repulse' value effectively defines node
  spacing, usually with values between 3 and 4, which 3 giving broad
  spacing, and 4 giving very close spacing of nodes in equivalent
  network positions, and broad spacing otherwise. Higher values tend to
  "clump" nodes closer together.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are passed to downstream functions:

  - when `remove_blanks=TRUE` it is passed to
    [`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md)

  - when `spread_labels=TRUE` it is passed to
    [`spread_igraph_labels()`](https://jmw86069.github.io/multienrichjam/reference/spread_igraph_labels.md)

## Value

`igraph` object with Concept network data, containing pathways connected
to genes. Each node has attribute `"nodeType"` of either `"Set"` or
`"Gene"`.

## Details

This function takes `Mem` S4 object, or legacy `list` mem data, or even
bare `matrix` data, and produces a Concept network (Cnet) in the form of
an `igraph` object.

Cnet network data are a form of bipartite graph, where every node is
either 'Gene' or 'Set', and where nodes are only connected across
`'Gene' <--> 'Set'`. See
[`igraph::graph_from_biadjacency_matrix()`](https://r.igraph.org/reference/graph_from_biadjacency_matrix.html).
A Cnet plot has particularly useful characteristic that Set (pathway)
nodes tend to be connected to very many Gene nodes, creating an
imbalance which is useful when defining a visual layout.

A vertex attribute `'nodeType'` stores the type of node: 'Gene' or
'Set'.

The `memIM` defines the network, and non-zero value constitutes an edge
(connection) between row (Gene) and column (Set).

Other data are associated with nodes when available.

1.  `geneIM`, `geneIMcolors`, `geneIMdirection`: define 'Gene' node fill
    and border colors.

2.  `enrichIM`, `enrichIMcolors`, `enrichIMdirection`: define 'Set' node
    fill and border colors.

Everything else helps customize the output, and can be used as-is.

### Generic Cnet

The `memIM` can be supplied as an `matrix` of `integer` values, with no
other data supplied, and this is sufficient to create a Cnet plot.
Therefore, any incidence matrix can be used to create a network.

If no other data are provided, nodes are colored by `categoryColor`
(Set) and `geneColor` (Gene), and shapes defined by `categoryShape` and
`geneShape`.

Be aware that `remove_singlet_genes=TRUE` default will hide rows (Genes)
which have no connection to any columns (Sets). This default ensures
that an incidence matrix can be subset for columns (Sets) of interest
without also filtering rows.

However, columns (Sets) are not removed, because they typically are not
empty in the most common workflows. Singlet columns, Set nodes with no
connected Gene nodes, can be filtered by calling
[`removeIgraphSinglets()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphSinglets.md).

## See also

Other jam Mem utilities:
[`Mem-class`](https://jmw86069.github.io/multienrichjam/reference/Mem-class.md),
[`MemPlotFolio-class`](https://jmw86069.github.io/multienrichjam/reference/MemPlotFolio-class.md),
[`fixSetLabels()`](https://jmw86069.github.io/multienrichjam/reference/fixSetLabels.md),
[`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md)

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
[`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md),
[`nodegroups2communities()`](https://jmw86069.github.io/multienrichjam/reference/nodegroups2communities.md),
[`rectifyPiegraph()`](https://jmw86069.github.io/multienrichjam/reference/rectifyPiegraph.md),
[`removeIgraphBlanks()`](https://jmw86069.github.io/multienrichjam/reference/removeIgraphBlanks.md),
[`subsetCnetIgraph()`](https://jmw86069.github.io/multienrichjam/reference/subsetCnetIgraph.md),
[`subset_igraph_components()`](https://jmw86069.github.io/multienrichjam/reference/subset_igraph_components.md),
[`sync_igraph_communities()`](https://jmw86069.github.io/multienrichjam/reference/sync_igraph_communities.md)

## Examples

``` r
use_sets <- c("eNOS Signaling",
   "Growth Hormone Signaling",
   "mTOR Signaling")
jam_igraph(mem2cnet(Memtest[, use_sets, ]),
   use_shadowText=TRUE)

```
