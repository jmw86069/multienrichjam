# multienrichjam version 0.0.20.900

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

# multienrichjam version 0.0.19.900

## New functions

* `collapse_mem_clusters()` is an intriguing new method of simplifying
the numerous pathways, by using each pathway cluster from
`mem_gene_pathway_heatmap()`. Each cluster is condensed to one
result, combining all genes in each cluster. The results are
surprisingly insightful, especially when numerous pathways
are present per cluster. This function also calls
`rank_mem_clusters()`, so it's possible to pick a handful of the
top pathways per cluster, as relevant.

# multienrichjam version 0.0.18.900

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

# multienrichjam version 0.0.17.900

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


# multienrichjam version 0.0.16.900

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

# multienrichjam version 0.0.15.900

## Changes to existing functions

* `multiEnrichMap()` several updates to improve robustness
in the overall workflow.

## new functions

* `heatmap_row_order()` and `heatmap_column_order()` provide enhanced
output from `ComplexHeatmap::row_order()`, mainly that they return
the actual rownames and colnames of data in the heatmap. Very useful
when the data used for the heatmap has been filtered internal
to the function.

# multienrichjam version 0.0.14.900

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

# multienrichjam version 0.0.13.900

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


# multienrichjam version 0.0.11.900

## Bug fixed

* `subsetCnetGraph()` was fixed to handle rare cases where
pathway set name is identical to one or more genes, which
happens with IPA pathway analysis, in the "Upstream Regulators"
output.

# multienrichjam version 0.0.10.900

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

# multienrichjam version 0.0.9.900

## bug fixes

* `importIPAenrichment()` no longer calls commandline `grep`
to remove blank rows, for now it uses `readLines()` and
`jamba::vigrep()` to select rows with at least one character.

# multienrichjam version 0.0.8.900

## bug fixes

* Fixed several calls to `strsplit()` that failed when sent factors,
which only happens when R `options("stringsAsFactors"=TRUE)` which
is the default in base R. Now all calls to `strsplit()` enforce
`as.character()` unless character type has already been enforced.
* Added numerous package prefixes to functions, to avoid importing
all package dependencies.


# multienrichjam version 0.0.7.900

## bug fixes

* `topEnrichListBySource()` and `topEnrichBySource()` were updated
to handle more default colnames for P-values, now including the
default colnames for `enrichResult` objects `c("pvalue","padjust")`.
When no sort columns are found, a warning message is printed.

## changes

* Added documentation for `multiEnrichMap()` and other functions.
* Added vignette describing the workflow starting with Ingenuity
IPA enrichment results.

# multienrichjam version 0.0.6.900

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

# multienrichjam version 0.0.5.900

## changes

* `list2imSigned()` is no longer dependent upon the input list `x`
having names, nor having unique names.

## enhancements

* `isColorBlank()` now handles list input, which is helpful when
applied to igraph objects where `"pie"` colors are accessed as a list.
* `removeIgraphBlanks()` was updated to use vectorized logic for `"pie"`
vertex attributes.

# multienrichjam version 0.0.4.900

## bug fixed

* Fixed small issues with `fixSetLabels()`.

## changes

* Updated `fixSetLabels()` to handle words the should be kept
uppercase, with some defaults pre-configured, e.g. "mRNA".

# multienrichjam version 0.0.3.900

## changes

* removed `mergeAllXY()` and `unnestList()` and moved them to the `jamba`
package. Added corresponding version requirement on jamba.

# multienrichjam version 0.0.2.900

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

# multienrichjam version 0.0.1.900

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
