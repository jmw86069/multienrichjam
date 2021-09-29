# multienrichjam 0.0.50.900

## changes to existing functions

* `mem_legend()` new argument `inset` passed to `legend()`.
Also this function uses `tryCatch()` to try to pass `...`
arguments, and if they fail it tries again without `...`.
Fun.


# multienrichjam 0.0.49.900

## extended function help text

* `jam_igraph()` is a drop-in replacement for `igraph::plot.igraph()`,
and its options were described in much more detail. In brief:

* it plots nodes vectorized which is substantially faster
* it changes the default to `rescale=FALSE` and maintains
aspect ratio 1:1, so layout coordinates are rendered without
distorting the x- and y-axis ranges
* it optionally bundles edge connections, which helps particularly
with large bipartite graphs (especially gene-to-pathway graphs.)
* it allows bulk adjustment to node size, label size, and label
distance from node center.

## changes to existing functions

* `mem_plot_folio()` new argument `edge_bundling="connections"`
that by default will bundle edges appropriate for Cnet plots.
Disable with `edge_bundling="none"`.
* `mem_enrichment_heatmap()` new argument `style="dotplot"` will
create a dotplot styled heatmap, whose points are sized proportional
to the number of genes involved in enrichment. This style is
under development and may require additional customization options
such as setting a max gene count for point sizes.


## bug fixes

* `mem_gene_path_heatmap()` fixed bug where enrichment names were
mangled by `data.frame(...)`, avoided by using
`data.frame(check.names=FALSE, ...)`.

Several functions were updated to include proper package prefix
for function calls:

* `matrixStats::colMins()`
* `igraph::V()`
* `igraph::E()`
* `igraph::degree()`
* `igraph::components()`
* `igraph::vcount()`
* `igraph::neighbors()`
* `igraph::vertex_attr_names()`
* `igraph::list.graph.attributes()`
* `igraph::ego()`
* `igraph::set_graph_attr()`
* `igraph::graph_attr()`

affecting functions:

* `subsetCnetIgraph()`
* `apply_nodeset_spacing()`
* `removeIgraphBlanks()`
* `removeIgraphSinglets()`
* `mem_find_overlap()`
* `spread_igraph_labels()`
* `subgraph_jam()`



# multienrichjam 0.0.48.900

## bug fixes

* `mem_gene_path_heatmap()` was updated to handle edge cases:

   * `gene_im_weight` 0 or 1; `enrich_im_weight` 0 or 1
   * user-supplied `cluster_rows` or `cluster_columns` with `colorize_by_gene`
   * The `colorize_by_gene` logic was changed to convert to an integer
   matrix that refers to colors and labels by factor levels, which helps
   for user-defined clustering methods, also helps labeling the color
   legend.
   * new argument `colramp` for user-defined color gradient.
   * `row_title` is no longer assigned when user-defined argument is `NULL`.
   This workaround helps the edge case with `gene_im_weight=1` where
   the `row_split` is forced to have a limited number of values regardless
   what `row_split` integer value is sent, thus causing `row_title`
   mismatch in length.
   * Fixed egregious type in second part of an `if` statement, only
   called when user supplies a custom cluster function.
   

# multienrichjam 0.0.47.900

## bug fixes

* `topEnrichBySource()` was not correctly handling multiple
`sourceColnames` values, instead was only using the first value.
This bug has been corrected.

## changes to existing functions

* * `topEnrichBySource()` and `topEnrichListBySource()` now
return `enrichResult` when supplied with `enrichResult`,
instead of coercing to `data.frame` which would then need
to be converted back to `enrichResult`.
* `topEnrichBySource()` and `topEnrichListBySource()` argument
default `sourceColnames` was changed to reflect
usage with R package `msigdbr` for MSigDB gene set data,
specifically `sourceColnames=c("gs_cat", "gs_subcat")`.
Similarly, default values were removed from `curateFrom`
and `curateTo`, since these defaults imposed a specific
outcome.
* `multiEnrichMap()` default values for `topEnrichCurate*` arguments
were changed to NULL; default argument values changed to
`topEnrichSources=c("gs_cat", "gs_subat")`.


# multienrichjam 0.0.46.900

## changes to existing functions

* `mem_gene_path_heatmap()` new argument `rotate_heatmap=TRUE`
will rotate the heatmap layout so pathway names are displayed
as rows, and genes are displayed as columns.

I still need example data to use for function document examples,
to show utility of rotating the gene-pathway heatmap.


# multienrichjam 0.0.45.900

## changes to existing functions

* `color_edges_by_nodes()` was implemented twice (smh) so the
older function was renamed `color_edges_by_nodes_deprecated()`.
The older function blended colors using a simpler approach
with `avg_colors_by_list()` that took a very fast hue average;
while the new function uses `colorjam::blend_colors()` that
uses red-yellow-blue additive color blending model.
Honestly the function `color_edges_by_nodes()` is not used much yet,
but is likely to be used more with community detection,
and Cnet/community edge bundling.
The new `color_edges_by_nodes()` also applies alpha to the
intermediate blended colors, to retain the relative weight
of colors during the blending step.


# multienrichjam 0.0.44.900

The edge bundling update!

The new functions are under rapid development, and represent potentially
quite useful functions for visualizing complex `igraph` networks,
specifically Cnet plots with Gene and Set node types.

The general technique of bundling edges between two groups
of nodes is relatively stable. The details and extent that edges are
bundled, including visual display of bundled edges, is
under active evaluation and development. So far it doesn't
seem to take much to improve the output figure, compared
to using straight edges.

In the near future, Cnet plots may by default enable
some form of edge bundling, as Cnet plots were in fact
the motivating example.

See pkgdown docs for `edge_bundle_nodegroups()` for visual
examples using the Karate network.

## current issues in edge bundling

`edge_bundle_nodegroups()` is the core function, and it does not
yet implement:

* edge clipping based upon igraph vertex shape, for example
using `igraph:::.igraph.shapes[["circle"]]$clip()`.
* edge arrows, which requires edge clipping for proper usage,
otherwise edge arrows will be underneath the node shape itself.
* edge labels, which could be implemented, except that bundled
edges are by definition much more likely to be too close for
effective labels. I rarely use edge labels, so I will target
the simplest thing that works.


## new igraph edge bundling functions

* `get_bipartite_nodeset()` - is a general use function to return
nodes that are grouped by having identical node neighbors. This
situation usually happens rarely, except with bipartite graphs
where there are often clusters of nodes with the same neighbors,
especially common for Cnet plot data.
* `edge_bundle_bipartite()` is the first edge bundling function,
it is actually a light wrapper around `edge_bundle_nodegroups()`.
Bipartite bundling connects a nodeset (defined as having the
same neighbor nodes) to each neighbor node.
* `edge_bundle_nodegroups()` is a general edge bundling technique
that bundles edges between node groups. These node groups can
be defined by any relevant technique, typically a
community detection algorithm such as `igraph::cluster_walktrap()`
or any other of several `igraph::cluster_*` functions.

## new igraph functions

* `color_edges_by_nodes()` is a simple function that blends the colors
of the two nodes involved in each edge, using `colorjam::blend_colors()`
since it uses RYB color blending. It is optimized so it only
blends unique combinations of colors.


## new utility functions

* `deconcat_df2()` is a necessary utility function that "expands"
a `data.frame` that has one or more columns with multiple
delimited values. It simply expands the rows to represent one
individual value per row for those columns. This function
will very likely be moved into the `"jamba"` package for
much broader use.
* `handle_igraph_param_list()` is a helper function for `jam_igraph()`
which is intended to update node and edge attributes in bulk based
upon another attribute. For example "make all the Gene nodes small,
and all the Set nodes large." Same with labels, colors, etc.


## updates to existing functions

* `jam_igraph()` has new arguments:

* `edge_bundling` - for `"nodegroups"`, `"connections"` and `"none"`
This argument will enable edge bundling by calling `edge_bundle_nodegroups()`
and `get_bipartite_nodeset()` when needed. See examples.
* `render_nodes,render_edges,render_nodelabels,render_nodegroups` are
options to skip various aspects of `igraph` plot features, either to
save time, or to allow more detailed visual layering.
* `vectorized_node_shapes` toggle the vectorized node shape plotting
when there is more than one node shape in the `igraph` object. Note
this feature is *substantially* faster, but changes the order that
nodes are rendered, by design. As a result, nodes are drawn in order
of shape, so each render is a bulk operation. For nodes that overlap,
or partially overlap, this feature will visibly change the ordering
of nodes. In general, speed is still ideal, and reducing node layout
overlaps should be a separate step.


# multienrichjam 0.0.43.900

## changes to existing functions

* `adjust_cnet_nodeset()` logic was updated for `expand`
which was not properly applying negative values to
compress the spacing between nodes in a node set. This
use is relatively rare but feasible.
* `adjust_cnet_nodeset()` argument `set_nodes` now optionally
allows referring to a nodeset using the name of one node
contained in the nodeset, as a convenience.

## new functions

* `apply_nodeset_spacing()` attempts to automate the process
of applying `adjust_cnet_nodeset()` for each subcluster of
nodes in a Cnet plot to enforce a minimum spacing between
nodes. The Fruchterman-Reingold layout algorithm
might offer a minimum distance threshold, but I could
not find it. Incidentally, this new function can also
compress node spacing, which might be useful when
gene labels are not shown.



# multienrichjam 0.0.42.900

## changes to existing functions

* `mem_plot_folio()` has much more detailed help text,
describing more details about clustering parameters,
and describing the returned objects in detail.
* `mem_plot_folio()` new argument `do_plot=TRUE` determines
whether each plot is rendered, or just returned as a plot
object to be reviewed separately.
* `mem_plot_folio()` the returned data now
includes the gene-pathway heatmap caption, which descibes
the gene/set filtering criteria, and the row/column
distance methods used for clustering. Since these arguments
have substantial effect on the pathway clusters, it is
helpful to keep this information readily available.
Similarly, the `cnet` plots as `igraph` objects store
the title as a graph attribute, accessible using
`graph_attr(cnet, "title")`.
* `mem_plot_folio()` argument `colorize_by_gene=TRUE` now
uses `colorjam::blend_colors()`.
* `mem_gene_pathway_heatmap()` is tolerant of `'...'` entries
not valid with `ComplexHeatmap::Heatmap()`, in order to
allow overloading `'...'` for other function arguments.


## new functions

* `nudge_igraph_node()` is intended to help move individual
nodes in an `igraph` layout, useful for adjusting nodes
to reduce label overlaps.
* `get_cnet_nodeset()` is a helper function to get the nodes
in a nodeset, defined as the `"Gene"` nodes that all connect
to the same `"Set"` nodes. It is helpful when looking at
a Cnet plot, and wanting easy access to the `"Gene"` nodes
in a particular cluster of nodes in the plot.
* `adjust_cnet_nodeset()` is intended to manipulate the
layout coordinates for nodes in a nodeset, useful to
expand, shift, or rotate nodes to help with visibility.
* `adjust_cnet_set_relayout_gene()` is a very useful function
with complicated name. It is intended to help move `Set`
nodes in a Cnet `igraph` layout, then re-positions all
the `Gene` nodes while keeping the `Set` nodes in fixed
positions.
* `rotate_igraph_layout()` which allows rotating an `igraph`
layout. It can also reflect layout coordinates across one
or more axes. It also calls other helper functions
as needed, `spread_igraph_labels()` and `reorderIgraphNodes()`.
* `rotate_coordinates()` is the underlying function to
rotate coordinates, it operates on a numeric matrix
and is called by `rotate_igraph_layout()`.

The next version should include `plot_cnet_heatmaps()` which
creates nice Cnet cluster plots where each cluster has a
corresponding expression heatmap displayed at the edge
of the figure.


# multienrichjam 0.0.41.900

## bug fixes

* `mem_plot_folio()` returns `NULL` when there is an error
during the Cnet collapse step, and does not report the error.
Now it at least reports the error. The underlying cause in
this case was ComplexHeatmap not accepting horizontal
color legend orientation in R-3.6.1 (grid package version 3.6.1)
because of some new unit arithmetic only available in
grid 4.0.0+. Hiding the color legend, or using vertical
color legend fixed the issue.

# multienrichjam 0.0.40.900

## bug fixes

* `collapse_mem_clusters()` - fixed error when using `apply()`
on a matrix that had one row or one column.
* `subsetCnetIgraph()` fixed issue when the Cnet igraph object
has no internal layout, it runs `relayout_with_qfr()` by default.

## changes to existing functions

* `mem_plot_folio()` now returns `clusters_mem` in the `list`,
which includes the pathway set names represented in each
cluster shown in the gene-pathway incidence matrix heatmap.

# multienrichjam 0.0.39.900

## changes to existing functions

* `cnetplot_internalJam()` was updated to include `"nodeType"` as
a node attribute, which helps distinguish `"Gene"` and `"Set"`
nodes in downstream operations. This change helps address
#5 to hide the gene labels on Cnet plots.

# multienrichjam 0.0.38.900

## bug fixes

* Update to address issue #4 in `enrichDF2enrichResult()` with
new argument `descriptionColname` which will force the resulting
colname to be `"Description"` to fit expectations of
`enrichplot:::fortify.internal()` which requires `df$Description`.
When `descriptionColname` is not supplied, or not found in the
input `enrichDF` a warning is issued that describes the problem.

## enhancements

* `heatmap_row_order()` and `heatmap_column_order()` now also
work with `HeatmapList` objects. By default they use the first
heatmap in the list, which should be consistent with all other
heatmaps.

# multienrichjam 0.0.37.900

## bug fixed

* `multiEnrichMap()` argument `subsetSets` was never implemented;
existing arguments `descriptionGrep` and `nameGrep` were used for
similar but insufficient purpose. The `subsetSets` argument defines
specific pathway names to retain for analysis, and is implemented
through `topEnrichBySource()` and by proxy `topEnrichListBySource()`.
These functions might better be called `subsetEnrichResult()` and
`subsetEnrichList()`, but will not rename these functions.

## enhancements

* `mem_plot_folio()` includes the number of rows (genes)
and columns (pathways) displayed in the gene-pathway incidence
matrix heatmap.

# multienrichjam 0.0.36.900

## bug fixed

* `enrichDF2enrichResult()` was updated to fix #3, thanks to
@john-lee-johnson for reporting. Issue arose because
`enrichplot::cnetplot()` expected enrichResult rownames to
be equal to values in the `"ID"` column of the enrichment
`data.frame`.

## changes to existing functions

* `multiEnrichMap()` new argument `min_count` requires a pathway
to contain at least `min_count` genes in order to be considered
a "hit". This filter is mostly important when used with `topEnrichN`
to use the top pathways -- it will only sort pathways then take
the top `topEnrichN` number of pathways that also contain at
least `min_count` genes.
* `topEnrichBySource()` and `topEnrichListBySource()` new arguments
`min_count`, `p_cutoff` require pathways to contain at least
`min_count` genes, and have no higher than `p_cutoff` enrichment
P-Value. Previously these functions only took the top pathways,
regardless of these filters (which were applied later). This change
allows the filters to be applied before taking the top `topEnrichN`
pathways, which mostly helps when `min_count` is greater than 1 --
sometimes pathways with only one gene involved in enrichment
have statistically significant P-value, but are not biologically
relevant for interpretation or follow-up experiments.
For `topEnrichListBySource()` the filter is applied to each
enrichment in the list, and any pathways meeting the criteria are
taken for all enrichment lists. So a pathway must 
be present in the top `topEnrichN` entries which meet both
the `p_cutoff` and `min_count` criteria to be retained by these
functions.

# multienrichjam 0.0.35.900

## bug fixes

* `mem_gene_path_heatmap()` was updated to fix a small issue
with `min_set_ct_each` which requires at least one enrichment
to have `min_set_ct_each` genes. However, this filter was
not applied alongside the `p_cutoff` -- therefore some pathways
with enough genes which were not significantly enriched were
fulfilling these criteria, as along as another enrichment was
significant. The new behavior (as expected) is to requires
an pathway to meet both the `min_set_ct_each` and `p_cutoff`
thresholds in the same enrichment in order to be retained
in the gene-pathway incidence matrix.
* `mem_gene_path_heatmap()` fixed edge case where genes present in
multiple enrichments counted more toward `min_set_ct_each`,
but should only count once per gene.


## changes to existing functions

* `mem_gene_path_heatmap()` new arguments `column_title`, and
`row_title` allow custom cluster names, which are also
carried to Cnet cluster names as needed.
* `mem_gene_path_heatmap()` was updated to handle `row_split`
and `column_split` arguments more intuitively, and to allow
`row_split` values `FALSE`,`0`,`1` to inactivate the split
completely. Previously it was not possible to turn off split.
* `mem_plot_folio()` argument was renamed from `pathway_row_split`
to `gene_row_split` to reflect the intent of this argument
more accurately. This change is early in the function lifecycle
and not in broader use yet -- better to change it now. Otherwise,
future argument name changes will not occur without some type
of backward compatibility.
* `mem_gene_path_heatmap()`, `mem_enrichment_heatmap()` and
`mem_plot_folio()` new arguments `row_cex` and `column_cex`
used to adjust row and column heatmap labels, which is helpful
when used with auto-sized axis labels to make minor adjustments.
* `mem_enrichment_heatmap()` was updated to make the color ramp
more consistent with `mem_gene_path_heatmap()`, so the color
scale atop the gene-pathway heatmap more accurately reflects the
color scale used in the enrichment heatmap itself.

# multienrichjam 0.0.34.900

## changes to existing functions

* `mem_gene_path_heatmap()` new argument `colorize_by_gene=TRUE`
will color the heatmap body using blended colors from the
`geneIMcolors` which represents the enrichments in which the
gene is involved. The default `colorize_by_gene=FALSE`
instead colors the heatmap body by the number of enrichments,
which can be confusing if one of the enrichment colors is
also red. The goal is to make it visually apparent when
a gene is involved in one enrichment, by using the color
from that enrichment. This feature is still in development
and testing.
* `mem_gene_path_heatmap()` changed argument name from `enrich_gene_weight`
to `enrich_im_weight`, before this argument is in wider use. Added
new argument `gene_im_weight`. These arguments more accurately reflect
the relative weight between enrichment and incidence matrix for
`enrich_im_weight`.

# multienrichjam 0.0.33.900

## changes to existing functions

* `collapse_mem_clusters()` was updated with new argument
`cluster_color_min_fraction` to help filter the enrichment
colors to include in each resulting cnet cluster. The intent
is not to represent colors where the number of significant
pathways is below this threshold. For example, a cluster of
10 pathways may have only one significant pathway for a
given enrichment set -- therefore this enrichment color
would not be included in the cnet cluster colors.

# multienrichjam 0.0.32.900

## changes to existing functions

* `mem_gene_pathway_heatmap()` cleaned up the heatmap overall:

   * left and top annotations have zero gap between columns and rows
   * left and top annotation legends have discrete color bars; for
   left it is always `c(0,1)`; for the top it uses -log10 integer steps.
   * top annotation legend appends `"-log10P"` to the label, to indicate
   that the color values are based upon `log10()` transform of the
   enrichment P-value.
   * the color legend has discrete color bar steps, indicating the
   number of enrichments where each gene is involved

# multienrichjam 0.0.31.900

## changes to existing functions

* `mem_gene_pathway_heatmap()` new argument `enrich_gene_weight`
used to adjust the relative influence of the enrichment `-log10 P-value`
and the gene incidence matrix on the column clustering. The
effect is to adjust how much the enrichment P-values, or the
gene content, affects the clusters. In principle both should
have similar effects, but sometimes it helps to favor gene
incidence or pathway enrichment.

# multienrichjam 0.0.30.900

## changes to existing functions

* `mem_gene_pathway_heatmap()` was altered to allow filtering
pathways for a minimum number of genes represented by at least
one enrichment result. For example, a pathway may be required to
have at least 4 genes, despite having a statistically significant
enrichment P-value. The new argument `min_set_ct_each` tests
each enrichment to see if any one has enough genes per pathway.
There are two main effects of this filter: The gene-pathway
heatmap will display fewer pathways; and any resulting
Cnet plots will have the same pathways removed.

# multienrichjam 0.0.29.900

## changes to existing functions

* `mem_plot_folio()` new argument `byCols` which is passed along
to `rank_mem_clusters()`, and is used to sort the pathway names
within each Cnet cluster. The default uses `"composite_rank"`,
which sorts by the `floor(-log10(pvalue))` then by the highest
number of genes per pathway. The `floor()` function effectively
sorts by the order of magnitude of the enrichment P-value,
dropping the details. The obvious alternative is `"minp_rank"`
which uses the `-log10(pvalue)` directly, which therefore does
not effectively sort by gene count since enrichment P-values
rarely tie. In our experience, between two pathways with
reasonably similar enrichment P-value (within one order of magnitude
such as 2.4e-5 and 3.5e-5) the pathway with more genes was
usually the more interesting/relevant biological pathway.

# multienrichjam 0.0.28.900

## new functions

`im2list()` and `imSigned2list()` are the reciprocal functions
to `list2im()` and `list2imSigned()`. The new functions convert
incidence matrix to list, or signed incidence matrix to list,
respectively. They're notable because they're blazing fast
in our testing, thanks to efficient methods from the `arules`
R package for interconverting list to compressed logical
matrix, and vice versa.

## changes to existing functions

* `mem_enrichment_heatmap()` was modified to use `p_floor`
as the ceiling for the row hierarchical clustering, previously
it used `ceiling=3` which caused the dendrogram to have `height=0`
when all enrichment results were lower than 0.001.
* `mem_enrichment_heatmap()` was modified to calculate its own
matrix colors when `color_by_column=TRUE`, until the
`colorjam::matrix2heatColors()` can be updated to the improved
method.
* `mem_plot_folio()` now correctly honors `p_floor`.

# multienrichjam 0.0.27.900

Added a new TODO.md file to track some new feature ideas.

## new functions

* `colors_from_list()` infers the proper order of colors from
a list of color vectors. It is called by `reorderIgraphNodes()`
where attributes `"pie.color"` and `"coloredrect.color"` contain
subsets of colors, but in the proper order. This function returns
all colors in their proper order.
* `colors_from_list()` takes a list of colors, and returns the
unique colors, ordered by the overall order inferred from the
list. It is internally called by `reorderIgraphNodes()` when
argument `colorV` is not suppleid, since the attribute `"pie.color"`
can be used to infer the correct order of colors.

## Changes to existing functions

* `collapse_mem_clusters()` now uses `"; "` as delimiter for the
`"set_names"` attribute, to make it easier to distinguish individual
pathway names.
* `reorderIgraphNodes()` was refactored to handle specific color
order defined either by argument `colorV` or by inferring the correct
order of colors from attributes such as `"pie.color"`. This change
should help order node colors more consistent to the original input
colors to `multiEnrichMap()`, specifically the argument `colorV`.
* `reorderIgraphNodes()` now has an example that shows the effects
of ordering pie nodes by color.

# multienrichjam 0.0.26.900

## Changes to existing functions

* `collapse_mem_clusters()` now includes the cluster name
by default with the `set_labels` which are used to display
the top `n` pathway names for each collapsed Cnet node.
For example `"Cluster A: Aryl Hydrocarbons; Granzyme Signaling"`.
* `collapse_mem_clusters()` arguments `max_labels` and `max_char_labels`
accept multiple values, which are recycled and applied to filter each
cluster. For example `max_labels=c(5,2)` will filter the first
cluster to display up to 2 labels, and the second cluster to display
up to 5 labels.

# multienrichjam 0.0.25.900

## Changes to existing functions

* Thanks to `simang5c` the function calls to `renameColumn()` were
corrected to `jamba::renameColumn()`, to resolve issue #1. Also fixed
several other references to `jamba` functions: `nameVector()`, `nameVectorN()`,
`makeNames()`, `printDebug()`, `mixedSortDF()`, `cPaste()`,
`igrepHas()`, `provigrep()`, `vigrep()`, `unvigrep()`, `rbindList()`,
`noiseFloor()`, `rmNA()`, `rmNULL()`, `normScale()`, `getColorRamp()`,
`deg2rad()`. Omg! There were so many more cases than I thought!
(Package-building options to import functions from packages seem
to go too far, the imported functions appear to be contained
in the new package which seems misleading... Ah well.)
* Also made several small changes to handle single-enrichment
input to multiEnrichMap(). It still provides some useful graphical
benefits even with only one enrichment input.

# multienrichjam 0.0.24.900

## Changes to existing functions

* `multiEnrichMap()` was updated to enable analysis with only one
enrichment result. Testing whether the downstream capabilities are
useful in the context of one enrichment result.

# multienrichjam 0.0.23.900

## Changes to existing functions

* Several functions have additional `jamba::printDebug()` output
when `verbose=TRUE`.
* `jam_plot_igraph()` and `jam_igraph()` new argument
`use_shadowText=TRUE` will enable `jamba::shadowText()` to replace
`graphics::text()` for igraph text, which affects node and edge
labels. The goal is to make labels more widely legible when they
are placed on top of light and dark colors, common when using dark
colored nodes on a white background.
* `jam_igraph()` now tries harder to respect a pre-defined xlim and ylim,
before using the range of layout coordinates
* `mem_gene_path_heatmap()` slightly adjusted the number of pathway
clusters to use based upon the number of columns in the gene-pathway
matrix, slightly increasing the number for low number of columns.
* `mem_plot_folio()` new argument pathway_column_split allows setting
a specific number of pathway clusters.
* `mem_plot_folio()` argument `do_which` allows creating only
the requested plots from the full sequence of plots. This function
is likely to become the core part of the analysis workflow:
 
> gene-pathway heatmap -> pathway clusters -> Cnet using pathway clusters

# multienrichjam 0.0.22.900

## New functions

* `shape.jampie.plot()` provides a new `igraph` vertex shape
`"jampie"`. It is a complete clone to the default `igraph`
shape `"pie"` except that it offers vectorized plotting, which
can be substantially faster for large `igraph` objects that
use `"pie"` nodes.
* `jam_plot_igraph()` is a similar clone to `igraph::plot.igraph()`
except that it uses vectorized plotting specifically when there
are multiple shapes in the same `igraph` object. Formerly, each
node is drawn individually which is substantially slower, especially
for `igraph` objects with more than 100 nodes. It also by default
converts shape `"pie"` to `"jampie"` in order to use vectorized
plotting via `shape.jampie.plot()` above.

## Changes to existing functions

* `jam_igraph()` new argument `plot_function` allows use of
custom plot function, specifically to use `jam_plot_igraph()`.

## Other changes

* Functions for igraph vertex shapes were moved into a new .R file,
and into a new function family `"jam igraph shapes"`.

# multienrichjam 0.0.21.900

## New functions

* `grid_with_title()` draws a grid object with title and optional
subtitle, allowing for grid objects such as `ComplexHeatmap::Heatmap`
objects, or even multiple heatmaps. That said, any grid `"gTree"`
object should work.

# multienrichjam 0.0.20.900

## New functions

* `mem_plot_folio()` is a new function that creates the current
recommended folio of multienrichment plots:

1. Enrichment P-value heatmap
2. Gene-pathway heatmap -- importantly with pathway clustering
3. Cnet using the pathway clusters, collapsed by cluster
4. Cnet using the pathway clusters, exemplar pathways per cluster
5. Cnet using each individual pathway cluster

## Changes to existing functions

* `multiEnrichMap()` now includes `p_cutoff` in the output,
which represents the enrichment P-value threshold used for
the analysis. This cutoff is useful in making other color
gradients respect the same threshold required for significant
enrichment results, so P-values that do not meet this threshold
can be colored white (or the background color.)
* `mem_plot_folio()` new argument (minor release) `do_which` to
help produce selected plot pages from a folio of plots.

# multienrichjam 0.0.19.900

## New functions

* `collapse_mem_clusters()` is an intriguing new method of simplifying
the numerous pathways, by using each pathway cluster from
`mem_gene_pathway_heatmap()`. Each cluster is condensed to one
result, combining all genes in each cluster. The results are
surprisingly insightful, especially when numerous pathways
are present per cluster. This function also calls
`rank_mem_clusters()`, so it's possible to pick a handful of the
top pathways per cluster, as relevant.

# multienrichjam 0.0.18.900

## New functions

* `rank_mem_clusters()` is a convenience function which ranks
pathways/sets within clusters, useful with a list of clusters
following `mem_gene_pathway_heatmap()`. It makes it easier to
rank pathways within a cluster, potentially choosing one exemplar
pathway to represent each cluster. This function is part of more
effort to streamline the overall analysis workflow.

## Bug fixes

* Fixed small issue with `mem_gene_path_heatmap()` which was
properly applying `min_gene_ct` and `min_path_ct` however it
did not check the resulting data to remove empty columns and
rows. This change has been made.

# multienrichjam 0.0.17.900

## New functions

* `color_edges_by_nodes()` takes an igraph, determines a
color for each node based upon the shape (using `avg_colors_from_list()`
for shape `"pie"` and `"coloredrectangle"`), then creates an
average color between two nodes, and uses that as the edge color.

## Changes to existing functions

* `jam_igraph()` has more configurable arguments, in the form
of named lists. The driving example is
`label_factor_l=list(nodeType=c(Gene=0.01, Set=1))`, that
will apply `label.cex*0.01` to "Gene" nodes,
and `label.cex*1.5` to "Set" nodes. Similar arguments:
`node_factor_l` to apply to node size, and
`label_dist_factor_l` to apply to label distance from node center.
* `shape.coloredrectangle.plot()` was updated to fix a small bug
in order of rendering colored rectangles, it was somehow drawing the
frame after the fill colors, which caused the frames to overlap
each other and appear transparent.


# multienrichjam 0.0.16.900

The vignette was updated to use the new functions, for
a much cleaner overall workflow.

## Changes to existing functions

* `reorderIgraphNodes()` uses new function `avg_colors_by_list()`
to create one color to represent each node, then `sort_colors()`
to sort nodes by color hue, instead of previous behavior which sorted
by the hex color string. This change affects attributes `"pie.color"`,
`"coloredrect.color"` and `"color"`. Recall that this function
`reorderIgraphNodes()` is intended to help visually organize
groups of nodes by their colors, to make it easier to tell how
many nodes are each color.
* `avg_colors_by_list()` now uses `weighted.mean()` based upon the
`"c"` channel (chroma, color saturation) in order to down-weight
the effect of grey colors, which really should have no hue. The default
weight for grey is 0.1, while the maximum value for fully saturated
colors is 100. Blending `"green"` with `"grey30"` yields slightly less
saturated green, as expected.
* `subsetCnetIgraph()` now includes convenience calls to
`removeIgraphBlanks()`, `spread_igraph_labels()`, `reorderIgraphNodes()`,
and `relayout_with_qfr()`. In almost all cases, creating a subset of
a Cnet necessitates a new layout, which therefore requires new
node ordering (sorting groups of equivalent nodes by color), and
then label orientation (the angle used to offset the node label from
the center of each node).

## New functions

* `sort_colors()` and `order_colors()` uses the `farver::encode_colour()`
function to convert colors to HCL, then sorts by hue, and returns the
original sorted vector or corresponding order, respectively. For
sorting colors inside a `data.frame`, convert colors to a factor
whose levels are `sort_colors(unique(colors))`, in order to maintain
ties during multi-column sorting.
* `apply_color_cap()` imposes a numeric range onto a color channel
for a vector or list of colors, using `jamba::noiseFloor()`. Values
outside the allowed range are forced to the range.
* `subgraph_jam()` is a custom version of `igraph::induced_subgraph()`,
that correctly handles subsetting the graph layout when the graph itself
is subsetted.
* `jam_igraph()` is a custom `igraph::plot()` with proper sizing when
`rescale=FALSE`, allowing better aspect ratios for certain network
layout coordinates.
* `with_qfr()` is a wrapper to `layout_with_qfr()` that returns
a layout specification object, for use with `igraph::add_layout_()`
in a more confusing workflow pattern than I expected.

# multienrichjam 0.0.15.900

## Changes to existing functions

* `multiEnrichMap()` several updates to improve robustness
in the overall workflow.

## new functions

* `heatmap_row_order()` and `heatmap_column_order()` provide enhanced
output from `ComplexHeatmap::row_order()`, mainly that they return
the actual rownames and colnames of data in the heatmap. Very useful
when the data used for the heatmap has been filtered internal
to the function.

# multienrichjam 0.0.14.900

## Changes to existing functions

* `list2im()` removed argument `makeUnique` because the
underlying conversion to matrix no longer requires that step,
and it was a performance hit for extremely large lists.
The argument `keepCounts` is the only requirement to
maintain the count of each entry per list.
* `multiEnrichMap()` now subsets the `geneIM` gene incidence
matrix to match the genes after `topEnrichN` filtering is
applied.
* `layout_with_qfr()` now by default will use edge attribute
`"weight"` as the `weights` argument when calling
`qgraph::qgraph.layout.fruchtermanreingold()`.
* `importIPAenrichment()` was enhanced for more robust import
conditions, specifically for different variations of missing
IPA enrichment results. Also, empty colnames are removed, to
help recognize the proper identifier in each scenario.

## New functions

* `enrichList2geneHitList()` takes a list of `enrichResult`
and returns the list of genes represented in each `enrichResult`.
Intended mainly for internal use by `multiEnrichMap()`.

# multienrichjam 0.0.13.900

## New plotting functions

* `mem_multienrichplot()` allows customized enrichMap-style plotting
of igraph objects. Notably, you can filter by Jaccard overlap,
or by overlap count -- the number of genes involved in the overlap.
There are otherwise a lot of overlaps that involve only one gene,
which is not the best way to build this type of network.
* `relayout_with_qfr()` is a light extension of `layout_with_qfr()`,
but adds default calls to `removeIgraphBlanks()`,
and `spread_igraph_labels()` to help things look pretty.
* `mem_legend()` draws a color legend in the corner of a figure,
using the colors defined in the `multiEnrichMap()` output.

## Other new functions

* `subset_igraph_components()` subsets an igraph based upon
connected components -- i.e. distinct subclusters.
* `memIM2cnet()` takes the pathway-gene incidence matrix and
produces a Cnet plot `igraph` object. If given the output from
`multiEnrichMap()` it will also color nodes using the gene
and enrichment incidence matrix colors.
* `avg_colors_by_list()` and `avg_angles()` are intended for very
rapid color blending.


# multienrichjam 0.0.11.900

## Bug fixed

* `subsetCnetGraph()` was fixed to handle rare cases where
pathway set name is identical to one or more genes, which
happens with IPA pathway analysis, in the "Upstream Regulators"
output.

# multienrichjam 0.0.10.900

## New functions

* `mem_gene_pathway_heatmap()` is a wrapper to `ComplexHeatmap::Heatmap()`
which takes output from `multiEnrichMap()` and produces a
heatmap styled similar to a pathway-gene incidence matrix. It includes
row and column annotations to help interpret the results.
* `mem_enrichment_heatmap()` produces a heatmap of the pathways
and enrichment P-values for each comparison.
* `spread_igraph_labels()` is a general use function for igraph
objects, it takes an igraph network and associated layout, and
arranges igraph nodes at an angel opposite the majority of edges
from each node. The result arranges labels around the outside of
the network, and typically away from other nodes. It isn't perfect,
but is visually a big step in the right direction.
* `removeIgraphSinglets()` is a simple function that removes nodes that
have no edges. I found myself doing it manually enough times I wanted
something quick and easy.

## enhancements

* `multiEnrichMap()` now includes two new edge attributes on MultiEnrichMap
igraphs: overlap_count which contains the number of shared genes
between two sets, and overlap_max_pct which contains the max percent
overlap between two sets -- based upon the # overlapped divided
by the smaller of the two sets. The overlap_count is useful as
an optional edge label, to show how many genes are involved.
* `enrichDF2enrichResult()` now forces `"setSize"` to contain integer
values, useful when pathway size is inferred from something like
gene ratio.

## Bug fixes

* `multiEnrichMap()` was not properly assigning `colnames(memIM)`,
although they were inferred from the `rownames(enrichIM)`. The bug
was because `as.character()` removes names from a character vector,
and the conversion was done to enforce proper handling by `strsplit()`
which gives unexpected results when the list may contain a factor.

# multienrichjam 0.0.9.900

## bug fixes

* `importIPAenrichment()` no longer calls commandline `grep`
to remove blank rows, for now it uses `readLines()` and
`jamba::vigrep()` to select rows with at least one character.

# multienrichjam 0.0.8.900

## bug fixes

* Fixed several calls to `strsplit()` that failed when sent factors,
which only happens when R `options("stringsAsFactors"=TRUE)` which
is the default in base R. Now all calls to `strsplit()` enforce
`as.character()` unless character type has already been enforced.
* Added numerous package prefixes to functions, to avoid importing
all package dependencies.


# multienrichjam 0.0.7.900

## bug fixes

* `topEnrichListBySource()` and `topEnrichBySource()` were updated
to handle more default colnames for P-values, now including the
default colnames for `enrichResult` objects `c("pvalue","padjust")`.
When no sort columns are found, a warning message is printed.

## changes

* Added documentation for `multiEnrichMap()` and other functions.
* Added vignette describing the workflow starting with Ingenuity
IPA enrichment results.

# multienrichjam 0.0.6.900

## new functions

* `importIPAenrichment()` imports Ingenuity IPA enrichment results,
by default splitting each enrichment table into its own `data.frame`.
It curates colnames to be consistent with downstream analyses.
* `curateIPAcolnames()` will curate the colnames of IPA data,
and ensures the values in the gene column are consistently delimited.
* `gsubs()` is similar to `base::gsub()` except that it applies
a vector of pattern-replacement operations in order. This function
may be moved to the `"jamba"` package.
* `find_colname()` is a helper function to find a colname given
a vector of expected values, which is matched directly, then
case-insensitively, then as a vector of patterns to match the
start, end, then any part of the colnames. By default the first
matching value from the first successful method is returned.
Good for matching c("P-Value", "pvalue", "Pval")
* `topEnrichBySource()` and `topEnrichListBySource()` subset the
input pathway enrichment results by taking the top `n` result,
then making sure the overall selected pathways are retained for
all enrichment tables.

# multienrichjam 0.0.5.900

## changes

* `list2imSigned()` is no longer dependent upon the input list `x`
having names, nor having unique names.

## enhancements

* `isColorBlank()` now handles list input, which is helpful when
applied to igraph objects where `"pie"` colors are accessed as a list.
* `removeIgraphBlanks()` was updated to use vectorized logic for `"pie"`
vertex attributes.

# multienrichjam 0.0.4.900

## bug fixed

* Fixed small issues with `fixSetLabels()`.

## changes

* Updated `fixSetLabels()` to handle words the should be kept
uppercase, with some defaults pre-configured, e.g. "mRNA".

# multienrichjam 0.0.3.900

## changes

* removed `mergeAllXY()` and `unnestList()` and moved them to the `jamba`
package. Added corresponding version requirement on jamba.

# multienrichjam 0.0.2.900

## changes

* `fixSetLabels()` was updated to include a wider set of canonical
pathway prefixes used in MsigDB v6.1.

## new functions

* `cnet2im()` and `cnet2df()` are helper functions used to convert
a Cnet igraph to either a incidence matrix, or a table summary useful
for extracting subsets of pathways using various network descriptors.
* `drawEllipse()` and `shape.ellipse.plot()` add an igraph vertex
shape "ellipse" whose shape is controlled by vertex.ellipse.ratio,
where `1` creates a circular node.

# multienrichjam 0.0.1.900

## new functions

* `multiEnrichMap()` the workhorse core function, implementing the full
workflow to pre-process the results for later inspection.
* `enrichDF2enrichResult()` converts a data.frame into an
`enrichResult` object as used in `clusterProfiler` and `DOSE`.
* `cnetplotJam()` custom function to create a cnet plot igraph object.
* `mergeAllXY()` which merges a list of data.frames while keeping
all rows.
* `unnestList()` which un-nests a list of lists, resulting in a flattened
list. It directly supports `mergeAllXY()` in order to provide a simple
list of data.frame objects from potentially nested list of lists of
data.frame objects.
* `enrichList2IM()` converts a list of enrichResult objects (or data.frames)
into an indidence matrix of gene rows and pathway columns.
* `enrichList2df()` combines a list of enrichResult objects (or data.frames)
into a single data.frame, using the union of genes, and the best P-value.
* `fixSetLabels()` to update pathway/gene set labels using some
small set of logic.
* `isColorBlank()` checks if a color is a blank color, either by comparison
to known colors, or transparency, or saturation/brightness.
* `igraph2pieGraph()` converts an igraph into one with pie node shapes,
optionally coloredrectangle node shapes.
* `layout_with_qfrf()` is an extension to `layout_with_qfr()` that returns
a function, therefore convenient to use in `igraph::plot()` with custom
arguments.
* `removeIgraphBlanks()` removes blank colors in either pie or
coloredrectangle node shapes, thereby making non-blank colors easier to
see in plots.
* `subsetCnetIgraph()` takes an igraph with Set and Gene nodeType, and
subsets based upon a fixed list of Set or Gene nodes, removing
singlet disconnected nodes from the output.

## changes

* added `dplyr` to package dependencies.
* added `DOSE` to package dependencies.
