# Convert data.frame to enrichResult

Convert data.frame to enrichResult

## Usage

``` r
enrichDF2enrichResult(
  enrichDF = NULL,
  pvalueCutoff = 1,
  pAdjustMethod = "none",
  keyColname = c("itemsetID", "ID", "Name", "Pathway"),
  pathGenes = "pathGenes",
  geneColname = c("geneNames", "geneID", "Gene", "Genes"),
  geneHits = "geneHits",
  geneRatioColname = c("GeneRatio", "^Ratio"),
  geneDelim = "[,/ ]+",
  geneSep = ",",
  pvalueColname = c("P.Value", "Pvalue", "FDR", "adj.P.Val"),
  descriptionColname = c("Description", "Name", "Pathway", "ID"),
  readable = NULL,
  msigdbGmtT = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- enrichDF:

  `data.frame` representing gene set enrichment results.

- pvalueCutoff:

  `numeric` value range 0 to 1, to define the P-value threshold for
  enrichment results to be considered in downstream processing.

- pAdjustMethod:

  `character` string to define the P-value adjustment method, or
  `"none"` for no additional adjustment. See
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) for valid
  values.

- keyColname:

  `character` value of the `colname(enrichDF)` containing the unique row
  identifier. It can be a pathway_ID or any uniquely identifying value.

- pathGenes:

  `character` or value of the `colname(enrichDF)` containing the number
  of genes in each pathway. This value will be derived from
  `geneRatioColname` if needed.

- geneColname:

  `character` value of the `colname(enrichDF)` containing delimiited
  genes in each pathway.

- geneHits:

  `character` value of the `colname(enrichDF)` containing the integer
  count of the gene hits in each pathway. This value will be derived
  from `geneRatioColname` if needed.

- geneRatioColname:

  `character` value of the `colname(enrichDF)` containing the character
  ratio of gene hits to pathway size, in format "50/100". This value is
  used when either `"pathGenes"` or `"geneHits"` are not supplied.

- geneDelim:

  `character` regular expression pattern used to separate genes in the
  `pathGenes` column into a vector of character values.

- pvalueColname:

  `character` value of the `colname(enrichDF)` containing enrichment
  P-values to use in downstream processing.

- readable:

  `logical` default NULL, sets the 'readable' flag for the resulting
  `enrichResult` object.

- msigdbGmtT:

  optional GmtT object (not currently implemented)

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

## Value

`DOSE::enrichResult-class` object.

## Details

This function takes a `data.frame` containing gene set enrichment
results, and converts it to a proper 'enrichResult' object defined in
`DOSE::enrichResult-class`. This object is supported by other functions
in the `clusterProfiler` suite of tools.

## See also

Other jam import functions:
[`IPAlist_to_hits()`](https://jmw86069.github.io/multienrichjam/reference/IPAlist_to_hits.md),
[`importIPAenrichment()`](https://jmw86069.github.io/multienrichjam/reference/importIPAenrichment.md)

Other jam conversion functions:
[`cnet2df()`](https://jmw86069.github.io/multienrichjam/reference/cnet2df.md),
[`cnet2im()`](https://jmw86069.github.io/multienrichjam/reference/cnet2im.md),
[`enrichList2IM()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2IM.md),
[`enrichList2df()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2df.md)
