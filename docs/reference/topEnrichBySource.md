# Subset enrichResult for top enrichment results by source

Subset enrichResult for top enrichment results by source

Subset list of enrichResult for top enrichment results by source

## Usage

``` r
topEnrichBySource(
  enrichDF,
  n = 15,
  min_count = 1,
  p_cutoff = 1,
  sourceColnames = c("gs_cat", "gs_subcat"),
  sortColname = NULL,
  countColname = c("gene_count", "count", "geneHits"),
  pvalueColname = c("P.Value", "pvalue", "FDR", "adj.P.Val", "qvalue"),
  directionColname = c("activation.z.{0,1}score", "z.{0,1}score"),
  direction_cutoff = 0,
  newColname = "EnrichGroup",
  curateFrom = NULL,
  curateTo = NULL,
  sourceSubset = NULL,
  sourceSep = "_",
  subsetSets = NULL,
  descriptionColname = c("Description", "Name", "Pathway"),
  nameColname = c("ID", "Name"),
  descriptionGrep = NULL,
  nameGrep = NULL,
  verbose = FALSE,
  ...
)

topEnrichListBySource(
  enrichList,
  n = 15,
  min_count = 1,
  p_cutoff = 1,
  sourceColnames = c("gs_cat", "gs_subcat"),
  sortColname = c(pvalueColname, "P-value", "pvalue", "qvalue", "padjust", "-GeneRatio",
    "-Count", "-geneHits"),
  countColname = c("gene_count", "count", "geneHits"),
  pvalueColname = c("P.Value", "pvalue", "FDR", "adj.P.Val", "qvalue"),
  directionColname = c("activation.z.{0,1}score", "NES", "direction", "z.{0,1}score"),
  direction_cutoff = 1,
  newColname = "EnrichGroup",
  curateFrom = NULL,
  curateTo = NULL,
  sourceSubset = NULL,
  sourceSep = "_",
  subsetSets = NULL,
  descriptionColname = c("Description", "Name", "Pathway"),
  nameColname = c("ID", "Name"),
  descriptionGrep = NULL,
  nameGrep = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- enrichDF:

  `enrichResult` or `data.frame` with enrichment results.

- n:

  `integer` maximum number of pathways to retain, after applying
  `min_count` and `p_cutoff` thresholds if relevant.

- min_count:

  `integer` minimum number of genes involved in an enrichment result to
  be retained, based upon values in `countColname`.

- p_cutoff:

  `numeric` value indicating the enrichment P-value threshold, pathways
  with enrichment P-value at or below this threshold are retained, based
  upon values in `pvalueColname`.

- sourceColnames:

  `character` vector of colnames in `enrichDF` to consider as the
  `"Source"`. Multiple columns will be combined using delimiter argument
  `sourceSep`. When `sourceColnames` is NULL or contains no
  `colnames(enrichDF)`, then data is considered `"All"`.

- sortColname:

  `character` vector, default `NULL`, indicating the colnames to
  sort/prioritize the enrichment data rows. Please use `NULL`.

  - Default `NULL` will use `pvalueColname` and the reverse of
    `countColname`, to prioritize lowest P-value, then highest gene
    count.

  - When `FALSE` it will not perform any sorting, and will use the input
    data as-is.

  - When `character` vector is provided, its values must exactly match
    the intended colnames, with optional prefix `"-"` to indicate
    reverse sort for a particular colname. These values are passed to
    [`jamba::mixedSortDF()`](https://jmw86069.github.io/jamba/reference/mixedSortDF.html)
    argument `byCols`.

- countColname:

  `character` vector of possible colnames in `enrichDF` that should
  contain the `integer` number of genes involved in enrichment. This
  vector is passed to
  [`find_colname()`](https://jmw86069.github.io/multienrichjam/reference/find_colname.md)
  to find an appropriate matching colname in `enrichDF`.

- pvalueColname:

  `character` vector of possible colnames in `enrichDF` that should
  contain the enrichment P-value used for filtering by `p_cutoff`.

- directionColname:

  `character` vector of possible colnames in `enrichDF` which may
  contain directional z-score, or other metric used to indicate
  directionality. It is assumed to be symmetric around zero, where zero
  indicates no directional bias.

- direction_cutoff:

  `numeric` threshold (default `0`) to subset enriched sets, filtering
  by magnitude of the absolute value of the `directionColname`.

- newColname:

  `character` string with new column name when `sourceColname` matches
  multiple colnames in `enrichDF`. Values for each row are combined
  using
  [`jamba::pasteByRow()`](https://jmw86069.github.io/jamba/reference/pasteByRow.html).

- curateFrom, curateTo:

  `character` vectors with pattern,replacement values, passed to
  `gsubs()` to allow some editing of values. The default values convert
  MSigDB canonical pathways from the prefix `"CP:"` to use `"CP"` which
  has the effect of combining all canonical pathways before selecting
  the top `n` results.

- sourceSubset:

  `character` vector with a subset of sources to retain. If there are
  multiple colnames in `sourceColnames`, then column values are combined
  using
  [`jamba::pasteByRow()`](https://jmw86069.github.io/jamba/reference/pasteByRow.html)
  and delimiter `sourceSep`, prior to filtering.

- sourceSep:

  `character` string used as a delimiter when `sourceColnames` contains
  multiple colnames.

- subsetSets:

  `character` optional set names to include by exact match.

- descriptionColname, nameColname:

  character vectors indicating the colnames to consider description and
  name, as returned from
  [`find_colname()`](https://jmw86069.github.io/multienrichjam/reference/find_colname.md).
  These arguments are used only when `descriptionGrep` or `nameGrep` are
  supplied.

- descriptionGrep, nameGrep:

  `character` vector of regular expression patterns, intended to subset
  pathways to include only those matching these patterns. The
  `descriptionGrep` argument searches only `descriptionColname`. The
  `nameGrep` argument searches only `nameColname`. Note that the sets
  are combined with OR logic, such that any pathways matched by
  `descriptionGrep` OR `nameGrep` or `subsetSets` will be included in
  the output.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

- enrichList:

  `list` of `enrichDF` entries, each passed to `topEnrichBySource()`.

## Value

`data.frame` subset up to `topEnrichN` rows, after applying optional
`min_count` and `p_cutoff` filters.

## Details

This function takes one `enrichResult` object, or a `data.frame` of
enrichment results, and determines the top `n` number of pathways sorted
by P-values, within each pathway source. This function may optionally
require `min_count` genes in each pathway, and `p_cutoff` maximum
enrichment P-value, prior to taking the top `topEnrichN` entries. The
default arguments do not apply filters to `min_count` and `p_cutoff`.

When the enrichment data represents pathways from multiple sources, the
filtering and sorting is applied to each source independently. The
intent is to retain the top entries from each source, as a method of
representing each source consistently even when one source may contain
many more pathways, and importantly where the range of enrichment
P-values may be very different for each source. For example, a database
of small canonical pathways would generally provide less statistically
significant P-values than a database of dysregulated genes from gene
expression experiments, where each set contains a large number of genes.

This function can optionally apply basic curation of pathway source
names, and can optionally be applied to multiple source columns. This
feature is intended for sources like MSigDB (see
http://software.broadinstitute.org/gsea/msigdb/index.jsp) which contains
columns `"Source"` and `"Category"`, and where canonical pathways are
either represented with `"CP"` or a prefix `"CP:"`. The default
parameters recognize this case and curates all prefix `"CP:.*"` down to
just `"CP"` so that all canonical pathways are considered to be the same
source. For MSigDB there are also numerous other sources, which are each
independently filtered and sorted to the top `topEnrichN` entries.

Finally, this function is useful to subset enrichment results by name,
using `descriptionGrep`, `nameGrep`, or `subsetSets`.

`topEnrichListBySource()` extends `topEnrichBySource()` by applying
filters to each `enrichList` entry, then keeping pathways across all
`enrichList` that match the filter criteria in any one `enrichList`. It
is most useful in the context of
[`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
where a pathway must meet all criteria in at least one enrichment, and
that pathway should then be included for all enrichments for the purpose
of comparative analysis.

## See also

Other jam enrichment functions:
[`add_pathway_direction()`](https://jmw86069.github.io/multienrichjam/reference/add_pathway_direction.md)

Other jam enrichment functions:
[`add_pathway_direction()`](https://jmw86069.github.io/multienrichjam/reference/add_pathway_direction.md)
