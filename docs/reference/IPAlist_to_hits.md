# Convert IPA list to a gene hit list or matrix

Convert IPA list to a gene hit list or matrix

## Usage

``` r
IPAlist_to_hits(
  IPAlist,
  signColname = c("Expr Fold Change", "fold.*change", "log.*ratio", "log.*fold",
    "log.*fc", "lfc", "ratio", "fold", "fc"),
  signThreshold = 0,
  geneColname = c("name", "symbol", "ID"),
  emptyValue = 0,
  returnType = c("list", "matrix"),
  verbose = FALSE,
  ...
)
```

## Arguments

- IPAlist:

  `list` of output from
  [`importIPAenrichment()`](https://jmw86069.github.io/multienrichjam/reference/importIPAenrichment.md)

- signColname:

  `character` vector of patterns to match the colname in the 'Analysis
  Ready Molecules' IPA xref table which should indicate the direction of
  change where present.

- signThreshold:

  `numeric` minimum absolute threshold, default 0, for a value in
  'signColname' to be considered a non-zero change.

- geneColname:

  `character` vector of patterns to match the colname in the 'Analysis
  Ready Molecules' IPA xref table which should indicate the gene
  involved. This column is validated by matching with the 'Canonical
  Pathways' geneID or geneNames column to ensure all values involved in
  pathway enrichment are also present in the xref column.

- emptyValue:

  `numeric` default 0, passed to
  [`list2imSigned()`](https://jmw86069.github.io/multienrichjam/reference/list2imSigned.md)
  to use when `returnType="matrix"` in the matrix for non-hits.

- returnType:

  `character` string

  - 'list' (default): returns a signed gene list, with `integer` vectors
    named by gene.

  - 'matrix': returns a numeric `matrix` with gene rownames.

- ...:

  additional arguments are passed to internal functions

## Value

`list` by default, or `matrix` with `returnType="matrix"`.

## Details

Given a `list` of IPA data produced by
[`importIPAenrichment()`](https://jmw86069.github.io/multienrichjam/reference/importIPAenrichment.md)
produce a single `list` or incidence `matrix` representing all gene hits
tested by IPA.

## See also

Other jam import functions:
[`enrichDF2enrichResult()`](https://jmw86069.github.io/multienrichjam/reference/enrichDF2enrichResult.md),
[`importIPAenrichment()`](https://jmw86069.github.io/multienrichjam/reference/importIPAenrichment.md)
