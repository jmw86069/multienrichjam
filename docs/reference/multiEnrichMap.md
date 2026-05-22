# Prepare MultiEnrichMap data from enrichList

Prepare MultiEnrichMap data from enrichList

## Usage

``` r
multiEnrichMap(
  enrichList,
  geneHitList = NULL,
  geneHitIM = NULL,
  colorV = NULL,
  nrow = NULL,
  ncol = NULL,
  byrow = FALSE,
  enrichLabels = NULL,
  subsetSets = NULL,
  overlapThreshold = deprecated(),
  p_cutoff = 0.05,
  cutoffRowMinP = deprecated(),
  enrichBaseline = -log10(p_cutoff),
  enrichLens = 0,
  enrichNumLimit = 4,
  nEM = 500,
  min_count = 3,
  topEnrichN = 20,
  topEnrichSources = c("gs_cat", "gs_subcat"),
  topEnrichCurateFrom = NULL,
  topEnrichCurateTo = NULL,
  topEnrichSourceSubset = NULL,
  topEnrichDescriptionGrep = NULL,
  topEnrichNameGrep = NULL,
  keyColname = c("ID", "Name", "pathway", "itemsetID", "Description"),
  nameColname = c("Name", "pathway", "Description", "itemsetID", "ID"),
  geneColname = c("geneID", "geneNames", "Genes"),
  countColname = c("gene_count", "count", "geneHits"),
  pvalueColname = c("padjust", "p.adjust", "adjp", "padj", "qvalue", "qval", "q.value",
    "pvalue", "p.value", "pval", "FDR"),
  descriptionColname = c("Description", "Name", "Pathway", "ID"),
  descriptionCurateFrom = c("^Genes annotated by the GO term "),
  descriptionCurateTo = c(""),
  directionColname = c("activation.z.{0,1}score", "NES", "direction", "z.{0,1}score"),
  direction_cutoff = 0,
  pathGenes = c("setSize", "pathGenes", "Count"),
  geneHits = c("Count", "geneHits", "gene_count"),
  geneDelim = "[,/ ]+",
  GmtTname = NULL,
  msigdbGmtT = NULL,
  returnType = c("Mem", "list"),
  verbose = FALSE,
  ...
)
```

## Arguments

- enrichList:

  `list` of `enrichResult` or `data.frame` objects.

  - The `names(enrichList)` are used in subsequent results.

  - Note that `data.frame` are converted to `enrichResult` using
    `enrichDF2enrichResult(x, ...)` where the '...' ellipses are used to
    recognize colnames.

  - Recommendation is to confirm each `data.frame` is properly converted
    to `enrichResult` upfront.

- geneHitList:

  `list` of character vectors, or `list` of `numeric` vectors whose
  names represent genes, or or `NULL`. When `NULL` the gene hit list for
  each enrichment result is inferred from the enrichment results
  themselves, however this option may incompletely represent which genes
  were statistical hits. Note that `geneHitList` and `geneHitIM` serve
  the same purpose and either can be supplied.

- geneHitIM:

  `numeric` matrix with gene rows, enrichment columns, and `numeric`
  values indicating the presence and/or direction of change for each
  gene. Note that `geneHitList` and `geneHitIM` serve the same purpose
  and either can be supplied.

- colorV:

  `character` vector of colors, length equal to `length(enrichList)`,
  used to assign specific colors to each enrichment result.

- nrow, ncol, byrow:

  optional arguments used to customize `igraph` node shape
  `"coloredrectangle"`, useful when the number of `enrichList` results
  is larger than around 4. It defines the number of columns and rows
  used for each node, to display enrichment result colors, and whether
  to fill colors by row when `byrow=TRUE`, or by column when
  `byrow=FALSE`.

- enrichLabels:

  `character` vector of enrichment labels to use, as an optional
  alternative to `names(enrichList)`.

- subsetSets:

  `character` vector of optional set names to use in the analysis,
  useful to analyze only a specific subset of known pathways.

- overlapThreshold:

  `numeric` (deprecated), value between 0 and 1, used to define the
  Enrichment Map, which is only created for legacy output enabled with
  `returnType="list"`.

  - To create a Multi-Enrichment Map `igraph` object, see
    [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md).

  - overlapThreshold is the Jaccard overlap score above which two
    pathways will be linked in the resulting network.

- p_cutoff:

  `numeric` value between 0 and 1, default 0.05, enrichment P-value
  required by at least one enrichment result to be retained in
  downstream analyses. This P-value can be reviewed in `enrichIM(Mem)`
  of the output, which is a matrix of P-values by pathway and
  enrichment. The column header assigned as P-value is stored in
  `headers(Mem)$pvalueColname`.

- cutoffRowMinP:

  (deprecated in favor of 'p_cutoff'). When it is non-NULL, it will be
  used for backward compatibility.

- enrichBaseline:

  `numeric` value, default uses `-log10(p_cutoff)`, the -log10 P-value
  threshold required for the color gradient to assign a color associated
  with enrichment.

  - In other words, P-values that do not meet this threshold are not
    colored in the color gradient for enrichment P-values.

  - To color all P-values use `enrichBaseline=0`.

- enrichLens:

  `numeric` value, default 0, indicating the "lens" to enhance the
  intensity of color gradients. Numbers above 0 make the color ramp more
  compressed, and more vivid at lower numeric values.

- enrichNumLimit:

  `numeric` value indicating the `-log10(P-value)` above which each
  color gradient is considered the maximum color, useful to apply a
  fixed threshold for each color gradient.

- nEM:

  `integer` (deprecated) maximum pathways to include in the Enrichment
  Map, which is only created for legacy output enabled with
  `returnType="list"`.

  - To create a Multi-Enrichment Map `igraph` object, see
    [`mem2emap()`](https://jmw86069.github.io/multienrichjam/reference/mem2emap.md).

- topEnrichN:

  `integer` maximum rows to retain from each `enrichResult` in
  'enrichList', for each source when supplied. Set `topEnrichN=0` or
  `topEnrichN=NULL` to retain all rows. Only rows where the pathway met
  other filtering criteria in at least one enrichment in 'enrichList'
  will be retained.

- topEnrichSources, topEnrichCurateFrom, topEnrichCurateTo,
  topEnrichSourceSubset, topEnrichDescriptionGrep, topEnrichNameGrep:

  arguments passed to
  [`topEnrichListBySource()`](https://jmw86069.github.io/multienrichjam/reference/topEnrichBySource.md)
  when `topEnrichN` is greater than `0`. The default values are used
  only when input data matches these patterns.

- keyColname, nameColname, geneColname, pvalueColname,
  descriptionColname:

  `character` vector in each case with text strings or patterns to use
  when matching or prioritizing colnames to assign to each type:

  - `key`: The primary unique key for each pathway

  - `name`: The short name for each pathway

  - `gene`: column with delimited gene symbols or identifiers tested for
    enrichment of each pathway.

  - `pvalue`: column with enrichment P-value, typically prioritizing
    either 'qvalue' or 'p.adjust'. Per `enrichResult` convention, the
    selected column is renamed to 'pvalue'.

  - `description`: longer description associated with the pathway.

  Each vector is passed to
  [`find_colname()`](https://jmw86069.github.io/multienrichjam/reference/find_colname.md)
  to find a suitable matching colname for each entry in `enrichList`.
  That function prioritizes full colname matches, then leading or
  trailing matches, then substring.

- descriptionCurateFrom, descriptionCurateTo:

  `character` vectors with patterns and replacements, passed to
  `gsubs()`, intended to help curate common descriptions to shorter,
  perhaps more user-friendly labels. One example is removing the prefix
  `"Genes annotated by the GO term "` from Gene Ontology pathways. These
  label can be manually curated later in the `Mem-class` methods,
  specifically `sets()<-` will allows assignment of custom names.

- pathGenes, geneHits:

  `character` values indicating the colnames that contain the number of
  pathway genes, and the number of gene hits, respectively. These values
  are optional, and not specifically used by multienrichjam.

- geneDelim:

  `character` pattern used with
  [`strsplit()`](https://rdrr.io/r/base/strsplit.html) to split multiple
  gene values into a list of vectors. The default (required) delimiter
  for `enrichResult` objects is `'/'`. A common alternative is `','`
  (comma-delimited).

- verbose:

  `logical` indicating whether to print verbose output. For `verbose` to
  cascade to internal functions, use `verbose=2`.

- ...:

  additional arguments are passed to various internal functions.

## Value

`list` object containing various result formats:

- colorV: named vector of colors assigned to each enrichment, where
  names match the names of each enrichment in `enrichList`.

## Details

This function performs most of the work of comparing multiple enrichment
results. This function takes a list of `enrichResult` objects, generates
an overall pathway-gene incidence matrix, assembles a pathway-to-Pvalue
matrix, creates EnrichMap `igraph` network objects, and CnetPlot
`igraph` network objects. It also applies node shapes and colors
consistent with the colors used for each enrichment result.

By default, each enrichment result table is subsetted for the top `n=20`
pathways sorted by pathway source, defined by colnames
`c("Source", "Category")`. For data without a source column, the overall
enrichment results are sorted to take the top 20. Once the top 20 from
each enrichment table are selected, the overall set of pathways are used
to retain these pathways from all enrichment tables. In this way, a
significant enrichment result from one table will still be compared to a
non-significant result from another table.

The default values for `topEnrichN` and related arguments are intended
when using enrichment results from `MSigDB`, which has colnames
`c("Source","Category")` and represents 100 or more combinations of
sources and categories. The default values will select the top 20
entries from the canonical pathways, after curating the canonical
pathway categories to one `"CP"` source value.

To disable the top pathway filtering, set `topEnrichN=0`.

Colors can be defined for each enrichment result using the argument
`colorV`, otherwise colors are assigned using
[`colorjam::rainbowJam()`](https://jmw86069.github.io/colorjam/reference/rainbowJam.html).

## See also

Other multienrichjam core functions:
[`jam_igraph()`](https://jmw86069.github.io/multienrichjam/reference/jam_igraph.md),
[`mem_plot_folio()`](https://jmw86069.github.io/multienrichjam/reference/mem_plot_folio.md)

## Examples

``` r
## See the Vignette for a full walkthrough example
```
