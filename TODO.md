# TODO

This document describes plans for enhancements to the
multienrichjam R package.

## 21jun2021

* COMPLETE: Allow rotating gene-pathway incidence matrix when using
`mem_gene_path_heatmap()`, so the pathway names are rows.
* Add test data object `mem` to be convenient for function examples,
and testthis test suite.


## Smaller usability items

* edge_bundling="connections" should also allow `render_groups=TRUE`
to work, by returning the `"nodegroups"` required for that step.

   * It should probably filter out singlet nodes by default, since
   the typical use case is "nodeset-to-node" for Cnet plots.
   * Main use case is with a Cnet plot, edge bundling by connected nodes
   should allow drawing an optional border around each group of nodes.



## COMPLETE: 10apr2021: implement edge bundling techniques

Priority: high 

Status: Implemented, early active testing of usability and functionality

* "connections" - Cnet edge bundles - special case of nodeset-to-node edge bundling
* "nodegroups": Node group bundles - general case of nodeset-to-nodeset edge bundling

Useful to improve readability/aesthetics of collapsed Cnet plots,
especially Cnet plots with many gene nodes where it is
difficult to tell which pathway clusters are connected
to each gene.


## 10apr2021: enhancements to `jam_igraph()` (must be done inside `jam_plot_igraph()`)

* Consider silently returning the `igraph` object plotted
where the vertex, edge, and graph attributes are updated
to the values used at runtime. For example, if user overrides
any attributes, those attributes will be present in the
object returned, so the next iteration someone could
just call `jam_igraph()` without any custom attributes
and it would produce the same plot again.

* See https://igraph.org/r/doc/plot.common.html
and `igraph::igraph_options()`



## 10apr2021: implement edge bundling in `jam_igraph()` and `jam_plot_igraph()`

Priority: high

Status: Implemented

* This task involves making edge bundling accessible as a simple
option during plotting, to prevent having to run 3 or 4 functions:

   * `jam_igraph()` with hidden edges
   * `bundle_node_edges()` to display edges on top of nodes
   * (optional) `jam_igraph()` with hidden edges, nodes with solid white
   color so they fully cover edges.
   * `jam_igraph()` with hidden edges, to display nodes on top
   of bundled edges while allowing nodes to have alpha transparency

* Note there is a package based on `ggraph` that implements
edge bundling - if `jam_ggraph()` is going to be our future,
maybe we should implement edge bundling using a similar approach
used by that R package.

   * https://github.com/schochastics/edgebundle


## 10apr2021: add `plot_cnet_heatmaps()` to multienrichjam

Priority: high

* `plot_cnet_heatmaps()` is in development, and arranges expression
heatmaps around a central Cnet plot, using genes in each cnet
cluster. It is intended to be used with Cnet clusters from
`collapse_mem_clusters()`, by-product of `mem_plot_folio()`.

This figure seems very useful because it integrates expression changes
alongside gene-pathway connections.

The multienrichment workflow would likely become:

* Prepare enrichment data
* `multiEnrichMap()`
* `mem_plot_folio()`
* `plot_cnet_heatmap()`



## 10apr2021: enable plotting node community data in `jam_igraph()`

Priority: low

* Analogous to `igraph:::plot.communities()` with `mark.groups`
* for edge bundling, this step could optionally display boundaries
around each node group.

Value is reinforcing node groups with a boundary. In testing,
the boundary actually made plots look more complex. And with simple
plots, it does not seem necessary since nodes are already
well-spaced.



## 10apr2021: include geneCount matrix output from `multiEnrichMap()`

Priority: medium

* Purpose is to complement the "enrichIM" matrix, with P-values
for each set. However "dysregulated pathways" may also be required
to have N number of genes, yet this data is not easily available.
* Rows are pathways, columns are enrichments, values are the number
of genes involved in enrichment.

   * Consider a matrix whose cells contain delimited gene symbols

Useful to apply filtering at matrix-level for pathways that
meet enrichment P-value and gene count thresholds.


## 2020: apply node color order based upon shape of node cluster

Priority: medium

`reorderIgraphNodes()`

* when a cluster of nodes is short-wide, the order should be left-to-right
* when a cluster of nodes is tall-skinny, the order should be top-to-bottom
* all else should be sorted left-to-right (or user-defined default order)

   * it is visually confusing when tall-skinny nodes are sorted left-to-right

Useful to automate the ordering of node colors


## 2020: apply subsetSets in `multiEnrichMap()`

Priority: low

* apparently it had never been implemented? User-defined
subset of pathways/sets to include in `multiEnrichMap()`,
as an alternative to using `topEnrichBySource()`.

Users can currently perform this subset step before running
`multiEnrichMap()`.


## COMPLETE: Apply min_count in `multiEnrichMap()`

* Currently `multiEnrichMap()` does not filter by number of
genes involved in enrichment, it only filters by enrichment P-value.
New argument `min_count` is applied only when `topEnrichN` is used,
but nothing else downstream is aware of filtering by `min_count`.
The corresponding argument in `mem_plot_folio()` is
`min_set_ct_each`, which requires a set (pathway) to contain
at least this many entries in at least one enrichment result
which also meets `p_cutoff` criteria for enrichment P-value.


## 2020: Optional highlight genes

Priority: low

* `mem_plot_folio()` and subsequent plots, optional argument
`highlight_genes` which would effectively hide all gene labels
except `highlight_genes` -- to help especially crowded plots.

   * Option 1: User-supplied genes
   * Option 2: Genes defined by one or more pathways/sets to highlight
   * Implement as wrapper to `jam_igraph()` to hide/display relevant node labels
   * `plot_cnet_heatmap()` could label the heatmap rows using `anno_mark()`

Useful in Cnet plots with too many genes to label, it allows
only a subset of genes to be highlighted.
Will become high priority if a manuscript decides to use
this technique.


## Consider a strict output format for `mem_plot_folio()`

Priority: medium

* Goal is to run `mem_plot_folio()` once then be able to plot
any component separately. Would help keep all plots in sync
with the settings used.
* Benefit: Always have ability to see the exact gene-pathway
incidence matrix heatmap, and its clustering, which is used
for a collapsed Cnet plots.
* Design idea: List format

   * Named by type of output:
   
      * enrichment P-values
      * gene-pathway heatmap
      
         * Heatmap object
         * gene clusters
         * pathway clusters
         * filter/cluster settings used
         
      * collapsed cnet igraphs
      
         * with each labeling option implemented
      
      * cnet exemplar igraphs
      
         * 1 exemplar per pathway cluster
         * 2 exemplars per pathway cluster
         * 3 exemplars per pathway cluster
         
      * each cnet cluster igraph
      
         * cluster 1 cnet igraph
         * cluster 2 cnet igraph
         * cluster 3 cnet igraph
         * ...
      
   * Each type of output contains a `list` of relevant components


## 10apr2021: make plot commands word for ggraph

Priority: low

* Goal is to use ggraph for ggplot2-type plotting instead
of using base R for `igraph` plots.
* `jam_ggraph()` ?

   * Test `scatterpie` R package which implements pie
   node shapes for use in ggraph. Unsure the data content
   requirement.


## 10apr2021: document examples for `reorderIgraphNotes()`

* Examples should show how to specify color order, or change
the order if necessary.


## COMPLETE: 2020: Improve `reorderIgraphNotes()`

* `reorderIgraphNodes()` when it encounters attributes with multiple
colors per node, such as `"pie.colors"` and `"coloredrect.colors"`,
calls `avg_colors_by_list()` to generate one blended color per node,
then it sorts those colors using `sort_colors()`. However, the
color order ultimately does not match the order in the color legend,
and other plots such as heatmaps.
* Future idea is to convert node attributes with multiple colors per
node into a `data.frame` where each color is in a separate column.
Then convert each column to a factor whose levels match the color
order from `colorV` (argument to `multiEnrichMap()`). The end result
should sort nodes consistent with the order of colors.
* Note this process assumes that all node colors are from a limited
set, which ideally should match `colorV`. Therefore if any color
gradient is applied to nodes with `nodeType="gene"`, this process
will not work.
