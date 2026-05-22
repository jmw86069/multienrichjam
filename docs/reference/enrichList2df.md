# Convert enrichList to data.frame

Convert enrichList to a single data.frame, merging results

## Usage

``` r
enrichList2df(
  enrichList,
  keyColname = c("ID", "Name", "pathway", "itemsetID", "Description"),
  geneColname = c("geneID", "geneNames", "Genes"),
  geneCountColname = c("gene_count", "Count", "geneCount", "geneHits"),
  pvalueColname = c("padjust", "p.adjust", "adjp", "padj", "qvalue", "qval", "q.value",
    "pvalue", "p.value", "pval", "FDR"),
  pvalueFloor = 1e-200,
  msigdbGmtT = NULL,
  geneDelim = "[,/ ]+",
  geneSep = ",",
  verbose = FALSE,
  debug = 0,
  ...
)
```

## Arguments

- enrichList:

  `list` of `enrichResult` or `data.frame` objects.

- keyColname, geneColname, geneCountColname, pvalueColname:

  `character` string indicating each column to use, or a vector of
  column names in order of priority.

- pvalueFloor:

  `numeric` lowest permitted P-value, intended to prevent using actual
  zero `0` which may occur sometimes when upstream processing rounds the
  value.

- msigdbGmtT:

  not currently implemented.

- geneDelim:

  `character` pattern used to split genes.

- geneSep:

  `character` string used to join genes, default `","` uses
  comma-delimited values.

- verbose:

  `logical` indicating whether to print verbose output.

- debug:

  `integer` used for debugging, with values `0, 1, 2, 3` recognized.
  Value `0`, default, performs no debug.

- ...:

  additional arguments are ignored.

## Value

`data.frame` with combined results across all `enrichList`.

## Details

This function takes one or more enrichment results, and combines results
into one `data.frame`.

- Genes are combined, reporting the unique genes in sorted order.

- `numeric` columns report the lowest or highest value, as appropriate:

  - P-value columns report the lowest observed value, in `pvalueColname`

  - Count columns report the highest observed value, except the primary
    count column reports the number of unique genes.

## See also

Other jam conversion functions:
[`cnet2df()`](https://jmw86069.github.io/multienrichjam/reference/cnet2df.md),
[`cnet2im()`](https://jmw86069.github.io/multienrichjam/reference/cnet2im.md),
[`enrichDF2enrichResult()`](https://jmw86069.github.io/multienrichjam/reference/enrichDF2enrichResult.md),
[`enrichList2IM()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2IM.md)
