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
