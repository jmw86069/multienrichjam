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
