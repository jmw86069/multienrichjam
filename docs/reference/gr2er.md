# Convert gseaResult to enrichResult equivalent result content

Convert gseaResult to enrichResult equivalent result content, focused
mainly on providing minimal result columns for downstream work.

## Usage

``` r
gr2er(er, ...)
```

## Arguments

- er:

  `gseaResult`

- ...:

  additional arguments are ignored.

## Value

`gseaResult`

## Details

This function is intended for internal use, and is not a full conversion
of `gseaResult` to a proper `enrichResult`. It simply makes minor
modifications to colnames for convenience within
[`multiEnrichMap()`](https://jmw86069.github.io/multienrichjam/reference/multiEnrichMap.md)
and other supporting functions.

Specific changes:

- Renames `'core_enrichment'` to `'geneID'` if `'geneID'` does not
  exist.

- Uses `'geneID'` to create `'Count'` column with number of genes that
  contribute to enrichment.

- Adds new attribute `'genes_signed'` which contains a subset of the
  gene rank scores stored in attribute `'geneList'`, only including
  those genes that contribute to enrichment, which is defined by the
  presence of a gene in the column `'core_enrichment'`.
