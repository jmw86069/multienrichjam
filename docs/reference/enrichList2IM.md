# Convert enrichList to IM incidence matrix

Convert enrichList to IM incidence matrix

## Usage

``` r
enrichList2IM(
  enrichList,
  addAnnotations = TRUE,
  keyColname = c("ID", "Name", "pathway", "itemsetID"),
  valueColname = c("qvalue", "q.value", "padj", "pvalue", "p.value"),
  emptyValue = NA,
  verbose = FALSE,
  GmtT = NULL,
  ...
)
```

## Arguments

- enrichList:

  `list` of `enrichResult` objects

- addAnnotations:

  `logical` not implemented, this argument is paired with `GmtT`.

- keyColname:

  `character` used to match colnames, referring to the unique
  identifier. Values in this column will become the column headers in
  the resulting incidence matrix.

- valueColname:

  `character` used to match colnames to determine the value to place in
  each cell of the incidence matrix.

- emptyValue:

  `numeric` value used to fill empty cells in the incidence matrix,
  default NA.

  - NA is used by default, indicating absence of a value.

  - When `NULL` and `valueColname` contains "gene", "count", "num",
    "hit" then `emptyValue=0`.

  - Otherwise `emptyValue=1` is used with the assumption that
    `valueColname` refers to P-values.

- verbose:

  `logical` indicating whether to print verbose output.

- GmtT:

  (not currently implemented), alternative gene set object. In future it
  may enable passing additional annotations associated with each gene
  set, including description, source, category, etc.

- ...:

  additional arguments are ignored.

## Value

incidence `matrix` with `numeric` values using the the enrichment
P-value (or FDR or Q-value as defined). Technically it is a "value
incidence matrix" where non-NA values are derived from the source data,
and NA represents absense of a value in the source data. It is possible
to have a value '0' (pure zero) which represents presence of a value,
where that value was '0', if the P-value reported by an enrichment tool
rounded to zero.

## Details

This function takes a `list` of `enrichResult` objects and creates an
incidence matrix using the value defined by `valueColname`.

TODO: Consider using
[`venndir::list2im_value()`](https://jmw86069.github.io/venndir/reference/list2im_value.html),
which itself may be moved to its own proper R package for set and list
manipulation, without the dependencies incurred by `venndir`.

## See also

Other jam conversion functions:
[`cnet2df()`](https://jmw86069.github.io/multienrichjam/reference/cnet2df.md),
[`cnet2im()`](https://jmw86069.github.io/multienrichjam/reference/cnet2im.md),
[`enrichDF2enrichResult()`](https://jmw86069.github.io/multienrichjam/reference/enrichDF2enrichResult.md),
[`enrichList2df()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2df.md)
