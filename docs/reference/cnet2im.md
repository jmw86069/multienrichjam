# Convert Cnet igraph to incidence matrix

Convert Cnet igraph to incidence matrix

## Usage

``` r
cnet2im(g = NULL, df = NULL, ...)
```

## Arguments

- g:

  igraph object containing Cnet data, specifically vertex attribute name
  "nodeType" with values "Set" and "Gene", and where "Set" nodes are
  only connected to "Gene" nodes.

- df:

  data.frame as optional input instead of `g`, usually the result of a
  previous call to
  [`cnet2df()`](https://jmw86069.github.io/multienrichjam/reference/cnet2df.md).

- ...:

  additional arguments are ignored.

## Value

numeric matrix with colnames defined by `"Set"` node names, and rownames
defined by `"Gene"` node names.

## Details

This function takes igraph object containing "Cnet" data, including
vertex attribute `"nodeType"` with values `"Set"` and `"Gene"`, and
where `"Set"` nodes are only connected to `"Gene"` nodes. It returns an
incidence matrix whose columns are "Set" node names and whose rows are
"Gene" node names.

## See also

Other jam conversion functions:
[`cnet2df()`](https://jmw86069.github.io/multienrichjam/reference/cnet2df.md),
[`enrichDF2enrichResult()`](https://jmw86069.github.io/multienrichjam/reference/enrichDF2enrichResult.md),
[`enrichList2IM()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2IM.md),
[`enrichList2df()`](https://jmw86069.github.io/multienrichjam/reference/enrichList2df.md)
