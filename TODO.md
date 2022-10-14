# TODO

## 13oct2022

* S4 objects?

   * `mem`:

      * Idea is to have one object type to handle output from `multiEnrichMap()`
      * benefit is mostly to have default behaviors for things like `print()`
      but also for clarity during the workflow. The `mem` object can be
      clear input to other functions.
   
   * `cnet`:
   
      * inherits `igraph` and extends assumptions, again mainly for clarity
      as input argument to other functions

* `mem_gene_path_heatmap()`

   * option to display dot plot format at the top of heatmaps, equivalent
   to calling `mem_enrichment_heatmap()`. Benefit would be to indicate
   direction (color) and number of genes (size) in addition to
   P-value (intensity).

* fix `reorderIgraphNodes()`

   * DONE: add "frame.color":
   `sortAttributes=c("pie.color", "pie.color.length", "coloredrect.color", "color", "pie.border", "frame.color", "label", "name")`
   * future: make sort more intelligent, so it uses appropriate
   color based upon node shape during the sort.

      * `shape="pie"`: use "pie.color", "pie.border", "frame.color", "label", "name"
      * `shape="coloredrect"`: use "coloredrect.color", "coloredrect.border", "frame.color", "label", "name"
      * all others assumed to be `shape="circle"` or similar:
      use "color", "frame.color" and ignore "pie.color", "pie.border",
      "coloredrect.color", "coloredrect.border", "label", "name"

* DONE: Include gene direction in the workflow:

   * DONE: Add `geneIMdirection` to mem object.
   
      * Requires `geneHitIM` input to `multiEnrichMap()`.
      * Alternate backwards compatibility is function to add
      `geneIMdirection` to `mem` object, and/or `cnet` `igraph` object.
      * DONE: `memIM2cnet()` should optionally use `enrichIMdirection` to apply
      border color.
   
   * DONE: `mem2cnet()` as its own function, minor variation and extension to
   `memIM2cnet()` which only uses the `memIM` portion of the data.
   * `jam_igraph()` should recognize `lwd` when plotting nodes. This change
   diverges from default behavior in `igraph` which only uses `par("lwd")`
   as a global adjustment for all lines while plotting nodes, therefore
   does not allow any nodes to have different line widths.
   * DONE: Outline of Cnet gene nodes by direction.
   * Gene-Path heatmap rows should optionally indicate direction,
   possibly another stripe using geneIMdirection.

* `multiEnrichMap()`

   * should call `mem2cnet()` - or just avoid this step altogether since
   the `cnet` process is better done with `mem_plot_folio()`
   * consider removing the `enrichMap` and `cnet` steps altogether.

## 26sep2022

* `mem_plot_folio()`

   * FIXED: `mem_enrichment_heatmap()` does not honor `do_plot=FALSE`.

* `mem_enrichment_heatmap()`

   * DONE: add argument `cluster_rows` to allow `cluster_rows=FALSE`.
   * DONE for manual plot calls: one goal is to draw this heatmap
   using the order from the gene-pathway heatmap clustering.
   * Future idea: allow plotting this data using the same order
   as the gene-pathway heatmap. This process would require
   running the gene-pathway clustering first, determining the
   column order, then using it to order the rows in this heatmap.

## 31aug2022

* Now that `jam_igraph()` and node shapes `"jampie"` render the
`pie.border` and `frame.color`, which can indicate the direction
of change for each gene, in context of each enrichment test, some
other changes should be made to the workflow:

   * Consider adding `pie.lwd` to recognized attributes in
   `shape.jampie.plot()`, but confirm that this argument can be used
   when rendering `"jampie"` node shapes. Currently the line width
   can only be adjusted with `par("lwd")` which is a global setting.
   * The default `frame.color` (`vertex.frame.color`) should probably
   be set to `NA` or `"transparent"` so the frame color is not visible
   for pie nodes. It shows up as a small white line now, which was not
   visible previously.
   * `reorderIgraphNodes()` should include `pie.border` and `frame.color`
   in sensible default locations, so nodes will also be sorted by these
   values when relevant, without the user having to add these columns.
   * New function to populate `frame.color` and `pie.border` for Cnet plots.
   
      * When `pie` has only one value, apply color to `frame.color`
      * When `pie` has multiple values, and different directions, apply
      colors to `pie.border`
      * When `pie` has multiple values, and all have the same direction, apply
      colors to `frame.color`
      * In absence of any directional data, set all `frame.color` to default,
      and set all `pie.border` to `"transparent"` or `NA`.

* Need to include gene direction of change in the workflow:

   * Some easy method to include direction of change in the `multiEnrichMap()`
   workflow, for example argument `geneHitList` currently uses a `character`
   vector, but could accept `integer` vector direction with genes stored as
   `character` labels, same as used with `venndir` signed input lists.
   * When gene direction is available, the `frame.color` and/or
   `pie.border` colors are defined.
   * `mem_gene_path_heatmap()` option to represent direction of change
   in the `gene_im` incidence matrix.

## aug2022

* Edge bundling makes assumptions for bipartite graphs (cnet)
that are difficult to use with normal graphs.

   * `edge_bundle_nodegroups()` probably needs to subset edges
   involved in each nodegroup before bundling all edges connected
   to these nodes. Currently, nodes are assumed to share all
   connections, but the effect should occur for edges where both
   ends of the edge are contained in the nodegroup entries.
   * consider node attribute "`nodegroup`".
   * consider edge attribute `"edgegroup"`, which would probably be the
   optimal approach for edge bundling, apart from implementing another
   technique similar to force-directed, hierarchically-defined, or
   density-directed edge bundling.
   * consider an option to specify the "midpoint" for a nodegroup,
   to allow some control on the spline curvature. Bonus points for
   allowing multiple points in order, to influence a path.

* `jam_igraph()`

   * fancy effect: allow edge colors to have multiple values, then
   interpolate color along each edge. For bundled edges, make
   a gradient equal to the number of line segments. For straight
   edges, break into `detail` number of pieces. The `detail` argument
   is also used by the `edge_bundle_*()` functions.

* subset `igraph` object by nodes, with added benefit that the
`graph_attr(g, "layout")` will also be subset, if present. It is
odd that the default `igraph::subgraph()` function would subset
the graph, leaving the layout which inevitably causes an error.

## 13jul2022

* `mem_gene_pathway_heatmap()`

   * would be useful to have an option to subset enrichments,
   this option would likely be useful to `mem_plot_folio()`, and
   `mem_enrichment_heatmap()`, so the effect could cascade to
   subsequent functions.

* new function idea? `subset_mem()`

   * Could be useful to accomplish the item above, to subset by enrichments
   prior to `mem_plot_folio()` related functions.
   * The challenge is that subsetting only `enrichIM` does not also
   update corresponding `memIM` data. Data to be updated:
   * `geneIM`, `geneIMcolors`, `geneIMdirection` - simple subset by colnames
   * `enrichIM`, `enrichIMcolors`, `enrichIMgeneCount` - simple subset by colnames
   * `memIM` - re-create after updating `geneIM`
   * `enrichList` - simple subset by name
   * Others not necessary for `mem_plot_folio()`:
   
      * `multiEnrichMap`, `multiCnetPlot` should be re-created using
      methods in `multiEnrichMap()`. To be fair, these objects are not
      that useful anymore, since `mem_plot_folio()` is generally preferred.

## 05jul2022

* `mem_gene_pathway_heatmap()` when supplied with custom `column_split`
throws an error when `cluster_columns=TRUE`.

   * internally it uses `amap::hcluster()` to generate a dendrogram/hclust
   * using split and dendrogram together is not allowed by `ComplexHeatmap::Heatmap()`
   * this clustering also uses the incidence matrix combined with the
   pathway enrichment annotation displayed along the top, and these values
   are weighted with `pathway_column_weight`.
   * A proper solution would be to provide a custom function for `cluster_columns`
   that internally combined the incidence matrix and enrichment matrix data
   together prior to clustering. If this function receives a
   numeric matrix with proper colnames, it should work.

* COMPLETE: `mem_plot_folio()` argument `gene_row_title=NULL` is being passed
to `mem_gene_pathway_heatmap()` and is therefore not using the default
`row_title=letters`.

## 23mar2022 issue #7 passing arguments through ... in `mem_plot_folio()`

* User reported an error when calling `mem_plot_folio(mem, node_factor=5, label_factor=1.5)`
* `mem_plot_folio()` is passing `...` to `ComplexHeatmap::Heatmap()` which
does not accept `...` and throws an error when receiving extra arguments.
* I need to limit `...` to arguments accepted by `Heatmap()`. I feel like
there is an R package function to handle this scenario, so I can avoid
writing `do.call(Heatmap, custom_args)`.

## 04feb2022

Workflow that starts with pathway-gene incidence matrix upfront,
bypassing the `multiEnrichMap()` workflow altogether.

* `memIM2cnet()` to convert pathway-to-gene incidence matrix to Cnet `igraph`.

   * Requires `geneIM` which is the enrichment-to-gene incidence matrix.
   In the driving example, each enrichment would represent a disease subgroup,
   with the full set of differential genes associated to each subgroup.
   * Requires `enrichIM` which contains enrichment-to-pathway whose values
   are P-values. When this data is not supplied, it should use 0.001 by default.
   I think these values are only used as an optional gradient color for the
   pathway nodes.
   * Optional `geneIMcolors` which is the same as `geneIM` except the cells
   contain colors to use for each gene. Currently the function does not
   fill these colors, the best method is to use `colorjam::matrix2heatColors()`
   
   ```R
   geneIMcolors <- colorjam::matrix2heatColors(x=geneIM,
      colorV=colorV)
   ```


## 03feb2022

* option to assign pathways to "functional themes" based upon presentation
by Adeline Chin in  Dr. Hanna Kim' group.

   * This step may "collapse" multiple pathways together, similar to
   grouping functional groups:
   union of genes,
   lowest enrichment P-value.


## 25jan2021

* `heatmap_row_order()`, `heatmap_column_order()` should return a flat
vector when there are no row_split or column_split, respectively.
Currently data is returned as a list of one-length vectors.

   * Consider moving into `jamba` package, in case this function
   needs to be re-used, it should not require loading the full
   `multienrichjam` package, which itself requires things like
   `igraph`, `clusterProfiler`, `qgraph`, `DOSE`, `matrixStats`.
   * Import `jamba::heatmap_row_order()` and `jamba::heatmap_column_order()`
   to ensure any functions or packages calling this function will
   succeed without error.


This document describes plans for enhancements to the
multienrichjam R package.

## 05nov2021

Now that directional z-score can also be associated with enrichment
P-values, heatmaps might need to use a bivariate color scale, to
indicate enrichment and directionality. See "stevens.bluered".
For example the `mem_enrichment_heatmap()` colors nodes by enrichment
P-value, more intense is more significant enrichment.
The z-score direction is used to apply red "activated" or blue "inhibited".
However, pathways with no z-score, or z-score below the threshold
are colored red by default. They should use a neutral color.



## 04nov2021

For IPA "Upstream Regulators" it sometimes offers a direction
implied by the `activation z-score`. Design idea is to implement
the directionality so it can be included in downstream analyses.

   * `mem_enrichment_heatmap()` - currently shades by the `-log10(pvalue)`
   however if there is directionality, it could be signed
   `+` for activated,
   `-` for inhibited. Then the heatmap color scale would use blue-white-red
   color gradient.
   * `mem$enrichIMdirection` contains matrix of direction, by default `1`
   means all have same direction.


During `multiEnrichMap()` it filters for `topEnrichN` entries for
each enrichment. It might be useful to retain the rank number for
each enrichment, to review when setting a different `topEnrichN`
threshold.

   * `mem$enrichIMgeneCount` contains matrix of gene counts
   * `mem$enrichIMrank` contains matrix of pathway rank (after filtering gene count)


## 01nov2021

* COMPLETE: `mem_enrichment_heatmap()` - the heatmap circles and legend circles
are not the same size - they should be fixed to the same absolute size.

   * COMPLETE: Optionally label each heatmap cell with the number of genes for
   visual reinforcement.


* COMPLETE: When there are more than 3 enrichments, the color legend on the
gene-pathway heatmap becomes unwieldy - taking over the whole figure.

   * Optionally (and by default) hide color legend for the gene-pathway heatmap.

* COMPLETE: The gene-pathway heatmap use_raster=TRUE causes artifacts in output,
it should be disabled by default. In future, debug why things go wrong.

   * I think the bug is caused by rasterization being done on underlying numbers
   before the color ramp is applied, in which case this problem cannot be
   solved when `colorize_by_gene=TRUE` since categorical values are assigned
   numbers which are not a proper color gradient.



## 20jul2021

* Replace `multiEnrichMap()` with `multienrichjam()` and simplify the arguments:

   * p_cutoff
   * min_gene_count
   * top_enrich_n
   * colnames: id, name, description, pvalue, gene
   * color_sub

* Port `mem_enrichment_heatmap()` argument `colorize_by_gene=TRUE` to
use `ComplexHeatmap::Heatmap()` instead of `jamba::imageByColors()`
for consistency, also so it can support `style="dotplot"`.

## 19jul2021

* `jam_igraph()` with `rescale=TRUE` should also scale igraph vertex
size and igraph label size according to the new axis ranges.
* Debug `jam_igraph()` when `vertex.size` is defined alongside
V(g)$size, and `node_factor`. It appears not to apply the size
properly.
* COMPLETE: `mem_enrichment_heatmap()` new option for dot plot format,
based upon `enrichplot::dotplot()` that sizes each dot by the
number of genes present.

## 24jun2021

* In `multiEnrichMap()` remove default `topEnrichSources`
and `topEnrichSourceSubset` which throw errors when not
using MSigDB data.
* `topEnrichBySource()` and `topEnrichListBySource()` should
be able to accept `enrichResult` as input, and return `enrichResult`
and not `data.frame` which is understandably lossy.
I need to understand the proper method for creating a subset
of an `enrichList` object, including its internal data.
* Streamline the `topEnrichListBySource()` workflow.


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
